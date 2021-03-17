"""Functions used for the hypothetical structure scoring algorithms.
"""
import pandas as pd
import numpy as np

def match_ions(spectra_df, hs_frag_df, tol):
    """Matches observed ions to a hypothetical structure ions with the provided tolerance.

    Parameters
    ----------
    spectra_df : pd.DataFrame
        Dataframe with all the ms/ms spectra provided.
    hs_frag_df : pd.DataFrame
        Dataframe with all the fragments and masses for hypothetical structures.
    tol : float
        +/- tolerance to use for matching m/z values.

    Returns
    -------
    pd.DataFrame
        Dataframe with all the observed ions that matched hypothetical structures.
    """
    # Create a couple empty lists to store indices of matched ions between hypothetical and observed
    hs_ion_idxs = []
    obs_ion_idxs = []

    # iterate through all potential ions and check if they're observed (with a tolerance), record indices of each df if they match
    for hs_ion_idx, row in hs_frag_df.iterrows():
        mw = row['hyp_mw']
        # Get all matched ions for the given hypothetical ion
        matched_ions = spectra_df[(spectra_df['m/z'] >= mw - tol) & (spectra_df['m/z'] <= mw + tol)]
        
        if matched_ions.empty is not True:
            for obs_ion_idx in matched_ions.index.values:
                hs_ion_idxs.append(hs_ion_idx)
                obs_ion_idxs.append(obs_ion_idx)
    
    # Pull out the indices of each that matched and merge into a new df
    matched_hs_ions = hs_frag_df.loc[hs_ion_idxs].reset_index(drop = True)
    matched_obs_ions = spectra_df.loc[obs_ion_idxs].reset_index(drop = True)
    matched_df = pd.concat([matched_hs_ions,matched_obs_ions ], sort = False, axis= 1)
    
    return matched_df

## Scoring functions


def _frac_scorer(matched_hs_ions_df, all_hyp_ions_df, N_spectra):
    """Fraction ion observed scorer.

    Provides a score based off of the fraction of hypothetical ions that were observed
    for a given hypothetical structure.

    Parameters
    ----------
    matched_hs_ions_df : pd.DataFrame
        Dataframe of observed ions that matched a specific hypothetical structure
    all_hyp_ions_df : pd.DataFrame
        Dataframe of all possible ions for a given hypothetical structure.
    N_spectra : int
        Number of spectra provided.

    Returns
    -------
    float
        Score for a given hypothetical structure.
    """

    # Calculate the number of matched ions observed and total possible
    N_matched_hs_ions = matched_hs_ions_df.shape[0]
    N_tot_hyp_ions = all_hyp_ions_df.shape[0]
    
    score = N_matched_hs_ions / (N_tot_hyp_ions*N_spectra)

    return score

def _magic_weights_calc(abund, bkgd):
    """Magical coefficient generating function that takes the abundance and the background value
    and generates a weighted coefficient.

    Parameters
    ----------
    abund : pd.Series
        Abundances column (should probably use the ceilinged one)
    bkgd : pd.Series
        Background column.

    Returns
    -------
    pd.Series
        Weighted coefficients corresponding to the normalized abundance values.
    """

    # Divide abundances by background to get normalized abundances
    normed_abund = abund / bkgd

    # Calculate the weight coeffsicients based on the normalized abundance values for each ion
    coeffs = normed_abund.copy()
    coeffs[coeffs < 4] = 0
    coeffs[(coeffs >=4) & (coeffs <= 5)] = 6.4*(coeffs[(coeffs >=4) & (coeffs <=5)]/10)**3
    coeffs[(coeffs >5) & (coeffs <= 6.25)] = coeffs[(coeffs >5) & (coeffs <=6.25)]/6.25
    coeffs[(coeffs >6.25) & (coeffs <= 8)] = np.sqrt((coeffs[(coeffs >6.25) & (coeffs <=8)])/6.25)
    coeffs[(coeffs >8) & (coeffs <= 9.5)] = 0.5656854*(coeffs[(coeffs >8) & (coeffs <= 9.5)])**(1/3)
    coeffs[coeffs >9.5] = 1.2
    
    return coeffs

def _magic_weights_scorer(matched_hs_ions_df, all_hyp_ions_df, N_spectra):
    """Weights scoring method

    Calculates weights for observed m/z values based off of their normalized abundances. Taken from
    the Van Der Donk paper.

    Parameters
    ----------
    matched_hs_ions_df : pd.DataFrame
        Dataframe of observed ions that matched a specific hypothetical structure
    all_hyp_ions_df : pd.DataFrame
        Dataframe of all possible ions for a given hypothetical structure.
    N_spectra : int
        Number of spectra provided.

    Returns
    -------
    float   
        Score for a given hypothetical structure.
    """

    # First calculate the magical coefficients for the observed ions
    abund = matched_hs_ions_df['abund_ceil']
    bk = matched_hs_ions_df['bkgd']

    # Sum the magical weighted coefficients for observed ions
    sum_weights = _magic_weights_calc(abund,bk).sum()
    N_tot_hyp_ions = all_hyp_ions_df.shape[0]
    
    score = sum_weights / (N_tot_hyp_ions*N_spectra)

    return score

def score_wrapper(matched_ions_df, all_ions_df, N_spectra, score_method = 'frac'):
    """Wrapper function for scoring methods

    Parameters
    ----------
    matched_ions_df : pd.DataFrame
        Dataframe of observed ions that matched hypothetical ions for a given hypothetical structure.
    all_ions_df : pd.DataFrame
        Dataframe of all hypothetical fragment ions for a given hypothetical structure.
    N_spectra : int
        Number of spectra provided.
    score_method : str, optional
        Scoring method to be used, by default 'frac'

    Returns
    -------
    pd.DataFrame
        Dataframe with scores for all the hypothetical structures

    Raises
    ------
    ValueError
        If scoring method provided is not a valid method.
    """

    # Create an empty score dataframe
    scores_df = pd.DataFrame({'hs_id': [], 'score': [], 'score_method': []})
    scores_df['hs_id'] = scores_df['hs_id'].astype(int)
    score_name = score_method # name used to refer to this scoring metric

    # Create empty lists to store hypothetical structure IDs and scores
    hs_ids = []
    scores = []

    if score_method == 'frac':
        scorer = _frac_scorer
    elif score_method == 'weights':
        scorer = _magic_weights_scorer
    else:
        raise ValueError('{method} is not a supported scoring metric'.format(method = score_method))

    # Iterate through each hypothetical structure and score it
    for hs_id in all_ions_df['hs_id'].unique():
        # Get all the matched ions and all possible ions for the hypothetical structure
        matched_hs_ions_df = matched_ions_df[matched_ions_df['hs_id'] == hs_id]
        all_hyp_ions_df = all_ions_df[all_ions_df['hs_id'] == hs_id]
        # Pass it to the scoring function
        score  = scorer(matched_hs_ions_df, all_hyp_ions_df, N_spectra)
        
        hs_ids.append(hs_id)
        scores.append(score)

    df = pd.DataFrame({'hs_id': hs_ids, 'score': scores, 'score_method': [score_name]*len(hs_ids)})
    scores_df = scores_df.append(df)
    
    return scores_df