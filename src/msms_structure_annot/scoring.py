"""Functions used for the hypothetical structure scoring algorithms.
"""
import pandas as pd
import numpy as np

def magic_coeffs(abund, bkgd):
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
    coeffs : pd.Series
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

def match_ions(spectra_df, hs_frag_df, tol):
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