"""Hypothetical structure related functions.

Generate hypothetical structures from proposed PTMs. Also fragments and
assigns hypothetical masses.

"""

import itertools
import pandas as pd

# Define molecular weights for polymerized amino acids
aa_mws = pd.DataFrame({
    'name': [
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P",
        "Q", "R", "S", "T", "V", "W", "Y","B","J"
    ],
    'mw': [
        71.0371, 103.0091, 115.0269, 129.0426, 147.0684, 57.0214, 137.0589, 
        113.0840, 128.0949, 113.0840, 131.0404, 114.0429, 97.0527, 128.0585, 
        156.1011, 87.0320, 101.0476, 99.0684, 186.0793, 163.0633, 69.03, 83.03
    ]
})

def gen_hss(ptms_df):
    """Generate hypothetical structures by marrying all combinations of all PTMs.

    Parameters
    ----------
    ptms_df : pd.DataFrame
        DataFrame with the post-translational modification information contained.
        Columns it should have: an ID "ptm_id", positions possible for modification "poss_mod_pos",
        and total number of mods observed "num_mods".

    Returns
    -------
    pd.DataFrame
        Long-form dataframe with PTMs and PTM locations for all the hypothetical structures. 
    """

    # Iterate through each PTM and get all combinations of possible modified positions for each
    all_mod_combs = list() # Store it by ptm_id number

    for i in range(0,ptms_df.shape[0]):
        ptm_id = ptms_df.loc[i]['ptm_id']
        mod_poss = ptms_df.loc[i]['poss_mod_pos']
        mod_count = ptms_df.loc[i]['num_mods']

        # Get all combinations of possible modified positions and store in dictionary
        all_mod_combs.append([comb for comb in itertools.combinations(mod_poss, mod_count)])

    # Then get all combinations of those combinations; each list item is a hypothetical structure
    mult_ptm_all_combs = list(itertools.product(*all_mod_combs))

    # Iterate through hypothetical structures list and massage into long-form dataframe
    hs_id_col = []
    ptm_id_col = []
    ptm_locs_col = []

    for hs_id, ptms in enumerate(mult_ptm_all_combs):
        for ptm_id, ptm_locs in enumerate(ptms):
            hs_id_col.append(hs_id)
            ptm_id_col.append(ptm_id)
            ptm_locs_col.append(ptm_locs)
            #print(hs_id,ptm_id,ptm_locs)

    hs_df = pd.DataFrame({'hs_id': hs_id_col, 'ptm_id': ptm_id_col, 'ptm_locs': ptm_locs_col})
    hs_df

    return hs_df


def _get_mass_seq(trunc_seq, hs_sub_df, ptms_df, trunc_type, parent_seq, N_term_mod, C_term_mod):
    """Gets the mass of a sequence given its truncation and the expected PTMs for the
    hypothetical structure.

    Parameters
    ----------
    trunc_seq : str
        Truncated sequence AA string.
    hs_sub_df : pd.DataFrame
        Hypothetical structure dataframe subsetted to only include relevant hypothetical strucuture.
    ptms_df : pd.DataFrame
        PTM information dataframe.
    trunc_type : str
        Truncation type (C-term or N-term)
    parent_seq : str
        Untruncated peptide sequence
    N_term_mod : float
        N-terminal modification mass shift.
    C_term_mod : float
        C-terminal modification mass shift.

    Returns
    ----------
    pd.DataFrame
        Dataframe row with all the info for the new hypothetical ion.
    """
    # Get the hypothetical structure ID.
    hs_id = hs_sub_df['hs_id'].values[0]
    # Function to get the mass of a sequence given its truncation and the expected ptms[]
    mass = 0
    ion_len = len(trunc_seq)
    
    for AA in trunc_seq: #iterate through each position and add the mass of each amino acid
            mass = mass + aa_mws[aa_mws['name'] == AA]['mw'].values[0]
    
    if trunc_type == 'C':
        term_mod = N_term_mod
        
        #included_ptms = [x for x in ptm_locs if x <= ion_len] # get mod sites included in the truncation
        ptm_range = (1,ion_len)
        ion_type = 'b'
        ion_name = ion_type + str(ion_len)
        #print(trunc_seq, ' bion length:',len(trunc_seq), '; ptms: ', str(included_ptms))
        
    elif trunc_type == 'N':
        #included_ptms = [x for x in ptm_locs if x > (len(parent_seq) - ion_len)]
        ptm_range = ((len(parent_seq) - ion_len),ion_len)
        
        if ion_len == len(parent_seq): # This is the parent ion
            
            term_mod = N_term_mod + C_term_mod
            ion_type = 'p'
            ion_name = ion_type
            #print(trunc_seq, ' p length:',len(trunc_seq), '; ptms: ', str(included_ptms))
            
        else: # name the gamma ion 
            term_mod = C_term_mod # + N_term_mod * added N_term_mod here for spectra that aren't deconvoluted
            ion_type = 'y'
            ion_name = ion_type + str(ion_len)
            
            #print(trunc_seq, ' yion length:',len(trunc_seq), '; ptms: ', str(included_ptms))
        
    mass = mass + term_mod # Add the terminal modification(s)
    
    # Iterate through PTMs based off PTM range and add PTM mass shifts for each PTM
    for ptm_id in hs_sub_df['ptm_id'].unique():
        # Get the shift from the ptm_id and the ptm dataframe
        ptm_shift = ptms_df[ptms_df['ptm_id'] == ptm_id]['m_shift'].values[0]
        # count the numbers of PTM in the fragment
        ptm_locs = hs_sub_df[hs_sub_df['ptm_id'] == ptm_id]['ptm_locs'].values[0]
        num_included_ptms = len([x for x in ptm_locs if (x > ptm_range[0] and x <= ptm_range[1])])
        # Shift mass by number of PTMs and the according mass shift
        mass = mass + ptm_shift*num_included_ptms

    #print(mass)

    hyp_ion_row = {'hs_id': hs_id, 'seq': trunc_seq, 'hyp_mw': mass, 'ion_name': ion_name, 'b_y_p': ion_type}
    
    return hyp_ion_row

def frag_hs(hs_df, ptms_df, parent_seq, N_term_mod, C_term_mod):
    """Fragments parental sequence for each hypothetical structure from N-term and C-term.
    Uses the PTM locations to define the masses for each fragment. 

    Parameters
    ----------
    hs_df : pd.DataFrame
        Hypothetical structure dataframe.
    ptms_df : pd.DataFrame
        PTM information dataframe.
    parent_seq : str
        Untruncated peptide sequence
    N_term_mod : float
        N-terminal modification mass shift.
    C_term_mod : float
        C-terminal modification mass shift.

    Returns
    -------
    pd.DataFrame
        Dataframe with all the fragmented hypothetical fragmented ions and masses.
    """

    columns = ['hs_id', 'seq', 'hyp_mw', 'ion_name', 'b_y_p']
    frag_df = pd.DataFrame(columns= columns)

    # Iterate through each hypothetical sequence and specify the mass
    for hs_id in hs_df['hs_id'].unique():

        # Pull PTM locations from hs_df
        hs_sub_df = hs_df[hs_df['hs_id'] == hs_id]
        
        # Truncate from the C-term
        for i in range(1,len(parent_seq)):
            trunc_seq = parent_seq[:i]        
            new_row = _get_mass_seq(trunc_seq, hs_sub_df, ptms_df, 'C', parent_seq, N_term_mod, C_term_mod)
            
            frag_df = frag_df.append(new_row, ignore_index=True)
        
        # Truncate from the N-term
        for i in range(0,len(parent_seq)):
            trunc_seq = parent_seq[i:]        
            new_row = _get_mass_seq(trunc_seq, hs_sub_df, ptms_df, 'N', parent_seq, N_term_mod, C_term_mod)
            
            frag_df = frag_df.append(new_row, ignore_index=True)
        
    return frag_df


def _add_charge(frag_df,charge_N, proton_m):
    """Returns dataframe with charged ions m/z values.

    Parameters
    ----------
    frag_df : pd.DataFrame
        Fragmented ion dataframe.
    charge_N : int
        Number of charges to add
    proton_m: float
        Mass of a proton constant
    
    Returns
    ----------
    pd.DataFrame
        Dataframe with charges added to m/z values.
    """

    charge_df = frag_df.copy()
    charge_df['charge'] = charge_N
    charge_df['hyp_mw'] = (charge_df['hyp_mw'] + charge_N*proton_m) / charge_N
    charge_df['ion_name'] = charge_df['ion_name'] + '^' + '{+'+ str(charge_N)+'}'

    return charge_df

def mk_charge_df(frag_df, charges, proton_m):
    """Makes charged versions of fragmented ions. 
    
    Charges the values provided in charges list. Also changes the ion name so it plots well
    as a label.

    Parameters
    ----------
    frag_df : pd.DataFrame
        Fragmented ion dataframe.
    charges : list
        List of integer numbers of hydrogens to add.
    proton_m : float
        Constant for mass of a proton.

    Returns
    -------
    pd.DataFrame
        Fragmented dataframe with charged versions appended.
    """
    # Add additional column for charge state; Everything that is uncharged is a deconvoluted mass
    new_frag_df = frag_df.copy()
    new_frag_df['charge'] = 0
    # Add brackets on name
    new_frag_df['ion_name'] = new_frag_df['ion_name'].apply((lambda x: x[0] + '_' + '{' + x[1:] + '}')) 

    # Create list of charged dfs and concatenate them
    charged_dfs = [_add_charge(new_frag_df, N, proton_m) for N in charges]
    frag_df_charged = pd.concat([new_frag_df] + charged_dfs).reset_index()

    # Do some formatting on the name so it displays nicely when plotted
    frag_df_charged['ion_name'] = frag_df_charged['ion_name'].apply((lambda x: '$' + x + '$'))
    # Some shenanigans to fix the parental ion names
    frag_df_charged.loc[frag_df_charged['b_y_p'] == 'p', 'ion_name'] = frag_df_charged.loc[
        frag_df_charged['b_y_p'] == 'p', 'ion_name'
        ].apply((lambda x: x[:2] + x[3:]))

    return frag_df_charged