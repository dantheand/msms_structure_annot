"""Spectra processing functions and utilities

All the mass-spec processing code is here

"""
import pandas as pd
import re

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
    
    return(ms_df)