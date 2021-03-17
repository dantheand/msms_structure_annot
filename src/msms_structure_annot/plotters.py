"""Plotting functions

Useful plotting functions are here.

"""

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

# Set matplotlib parameters
mpl.rcParams['figure.dpi']=150
mpl.rcParams['pdf.fonttype'] = 42 # For manipulatable fonts in pdfs
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['savefig.transparent'] = True

# Private method to label points
def _label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for _, point in a.iterrows():
        ax.text(point['x']+10, point['y'], str(point['val'])) # Add label offset here if you like


def label_spectra_plot(ms_df, matched_df, ms_file_nums, hs_id):
    """Vertical line plotting function for mass spec data. Plots
    a vertical line at each m/z value with the height according to the abundance.

    Parameters
    ----------
    ms_df : pd.DataFrame
        Datafrane with all the mass spectra m/z values and their abundances for
        all the different spectra.
    matched_df : pd.DataFrame
        Dataframe with the m/z values that matched hypothetical structures.
    ms_file_nums : list
        List of ms file numbers (e.g. [1,2,...] for ms1.txt, ms2.txt, etc.)
    hs_id : int
        Hypothetical structure ID to plot ions for (should correspond to a value in
        the matched_df 'hs_id' column).

    Returns
    -----------
    matplotlib.figure.Figure
        The whole figure that is plotted.
    np.array
        Array of axes being plotted.
    """
    # Make the plot
    fig, axs = plt.subplots(nrows = len(ms_file_nums), figsize = (8,len(ms_file_nums)*4))

    # Iterate through each spectrum and plot and label
    for ms_file, ax in zip(ms_file_nums,axs):
        # Get the observed masses / abundances for a single spectrum
        sub_df = ms_df[ms_df['spec_num'] == ms_file]

        # Get the masses that match both the spectra masses and the hypothetical structure masses
        sub_matched_df = matched_df[(matched_df['spec_num'] == ms_file) & (matched_df['hs_id'] == hs_id)]

        # Plot the original spectrum as vertical lines at each point
        mz_vals = sub_df['m/z'].values
        abunds = sub_df['orig_abundance'].values # Using original abundances here (could use ceilings too)
        ax.vlines(x = mz_vals, ymin = 0, ymax=abunds, linewidth = 1, color = 'k')

        # Plot the masses that match the hypothetical structure and color it by b/y/p ion
        
        mz_vals = sub_matched_df['m/z'].values
        abunds = sub_matched_df['orig_abundance'].values
        # Add color generation code here

        # Plot vertical lines
        ax.vlines(x = mz_vals, ymin = 0, ymax=abunds, linewidth = 1)
        # Plot some dots on top too
        ax.scatter(x = mz_vals, y=abunds, s = 5)

        # Label each point
        _label_point(x = sub_matched_df['m/z'], y = sub_matched_df['orig_abundance'], 
            val = sub_matched_df['ion_name'], ax = ax)

        # Format the plot
        ax.set_yscale('log')
        ax.set_xlim(0, 2000)
        ax.set_ylim(1e1, 1e4)
        ax.set_xlabel( "m/z", size = 10)
        ax.set_ylabel( "Abundance", size = 10)

    return(fig, axs)