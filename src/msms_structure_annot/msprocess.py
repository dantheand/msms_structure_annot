"""Spectra processing functions and utilities

All the mass-spec processing code is here

"""
import pandas as pd
import re
import numpy as np

def import_ms_files(ms_files_dir):
    """MS file importer

    Imports a series of msms files in csv or tab-separated format into a single dataframe.
    Files should be numbered ms1.csv, ms2.csv, ms3.csv.
    Columns of the tab-separated file should be: m/z | abundance

    Parameters
    ----------
    ms_files_dir : Path
        Path to directory with the ms/ms files in it.

    Returns
    ----------
    ms_df : pd.DataFrame
        Dataframe with all msms data concatenated
    """

    assert ms_files_dir.exists(), 'No directory found at that path!'

    # Import all msms files in the directory
    ms_files = list(ms_files_dir.glob('ms*'))
    assert len(ms_files) > 0, 'No MS files found! Did you specify the names as "ms[0-9].csv"?'

    # Iterate through ms files and import them into a dataframe
    ms_df = pd.DataFrame(columns = ['m/z', 'orig_abundance', 'spec_num'])

    for ms_file in ms_files:
        spec_num = int(re.search('ms[\d+].', str(ms_file))[0][2:-1]) # extract the number in the ms file
        # Import the csv file (tab separated works too)
        new_ms = pd.read_csv(ms_file, skiprows = [0,1], delimiter = '\t', header = 0, 
            names = ['m/z','orig_abundance'], index_col = False)
        new_ms['spec_num'] = spec_num
        ms_df = ms_df.append(new_ms, ignore_index = True)
    
    return ms_df

# Change abundances such that they cap out at a given abundance value
def abund_ceiling(abundances, upper_lim):
    """Set an abundance ceiling value

    Parameters
    ----------
    abundances : pd.Series
        Mass abundance values
    upper_lim : float
        Multiplier of the mean abundance value to set the ceiling at.

    Returns
    -------
    abundances : pd.Series
        Abundances with ceiling applied.
    """
    
    # First set all ions above a certain threshold to a single value
    max_val = upper_lim * np.mean(abundances)
    abundances.mask(abundances > max_val, max_val, inplace = True)
    
    return abundances

def bkgd_calc_ser(abundances, N):
    """Background abundance value calculator.

    Sections off the m/z axis into N different sections.
    Calculates the "background" value (threshold to be used for signal) 
    from the average of a given section and its two neighboring sections.

    Parameters
    ----------
    abundances : pd.Series
        Ceilinged abundances.
    N : int
        Number of sections to split into when calculating background.

    Returns
    -------
    bkgd
        Background value for each abundance point.
    """

    # Split indices into sections
    sections = np.array_split(abundances.index.values, N)

    bkgd = pd.Series(data = None, index = abundances.index.values, dtype = 'float64')

    # Iterate through each section and calculate background noise value based on algorithm above
    #   Edge cases are required fro the first and last sections
    for i in range(len(sections)):
        # Get indices of appropriate m/z signals for each section's calculation
        if i == 0: # If it's the first section, dont use the previous section in the calculation
            idxs = np.concatenate((sections[i],sections[i+1]))
        elif i == (len(sections)-1): # If its the last section... don't use the following section
            idxs = np.concatenate((sections[i-1],sections[i]))
        else: # Otherwise calculate background signal based on previous section, current section, and following section
            idxs = np.concatenate((sections[i-1],sections[i],sections[i+1]))
        
        bkgd_sig = np.mean(abundances.loc[idxs]) # take the mean of the abundance values in the given window
        bkgd.loc[sections[i]] = bkgd_sig

    return bkgd