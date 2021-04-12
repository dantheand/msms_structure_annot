"""Plotting functions

Useful plotting functions are here.

"""

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
from adjustText import adjust_text

# Set matplotlib parameters
mpl.rcParams['figure.dpi']=150
mpl.rcParams['pdf.fonttype'] = 42 # For manipulatable fonts in pdfs
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['savefig.transparent'] = True

# Private method to label points
def _label_point(x, y, val, ax):
    # Should probably make this take arrays instead of series
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    # Store annotations to be used by "adjust_text"
    annots = []
    for _, point in a.iterrows():
        annot = ax.text(point['x'], point['y'], str(point['val'])) # Add label offset here if you like
        annots.append(annot)
    adjust_text(annots)


def label_spectra_plot(ms_df, matched_df, ms_file_nums, hs_id, xlims = [(0,2000)], ylims= [(0,1e5)],
                        auto_yscale = True, annot_sn_lim = 0):
    """Vertical line plotting function for mass spec data. 
    
    Plots a vertical line at each m/z value with the height according to the abundance. Also
    plots the matched hypothetical ions.

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
    xlims : list
        List of x limits to use for each plot. Only specify one if you want it globally applied.
    ylims : list
        List of y limits to use for each plot. Only specify one if you want it globally applied.
    auto_yscale : bool
        Boolean to decide whether or not to autoscale the y axis.
    annot_sn_lim : float
        Signal to noise limit to exclude annotations below.

    Returns
    -----------
    matplotlib.figure.Figure
        The whole figure that is plotted.
    np.array
        Array of axes being plotted.
    """
    # Make the plot
    N_plots = len(ms_file_nums)
    
    fig, axs = plt.subplots(nrows = N_plots, figsize = (8,N_plots*4))

    # If there's only one plot, make the axes object an iterable so it can be parsed below
    if N_plots == 1:
        axs = [axs]

    # Make xlim / ylim arrays
    xlims_processed = []
    ylims_processed = []
    if len(xlims) == 1: # If only a single limit is provided, apply limits across all plots
        xlims_processed = xlims * N_plots
    else:
        assert len(xlims) == N_plots, "Number of x-axes limits specified must match number of spectra! (or be 1)"
        xlims_processed = xlims

    if len(ylims) == 1: # If only a single limit is provided, apply limits across all plots
        ylims_processed = ylims * N_plots
    else:
        assert len(ylims) == N_plots, "Number of y-axes limits specified must match number of spectra! (or be 1)"
        ylims_processed = ylims


    # Iterate through each spectrum and plot and label
    for ms_file, ax, xlim, ylim in zip(ms_file_nums, axs, xlims_processed, ylims_processed):
        # Get the observed masses / abundances for a single spectrum
        sub_df = ms_df[ms_df['spec_num'] == ms_file]

        # Get the masses that match both the spectra masses and the hypothetical structure masses
        sub_matched_df = matched_df[(matched_df['spec_num'] == ms_file) & (matched_df['hs_id'] == hs_id)]

        # Plot the original spectrum as vertical lines at each point
        mz_vals = sub_df['m/z'].values
        abunds = sub_df['orig_abundance'].values # Using original abundances here (could use ceilings too)
        ax.vlines(x = mz_vals, ymin = 0, ymax=abunds, linewidth = 1, color = 'k')

        # Define axes limits
        # Scale y value to the largest abundance value if autoscale is true
        if auto_yscale:
            upper_y_lim = 1.25*np.amax(abunds)
        else:
            upper_y_lim = ylim[1]

        # Format plot
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], upper_y_lim)
        ax.set_xlabel( "m/z", size = 10)
        ax.set_ylabel( "Abundance", size = 10)

        # Label each point

        # Truncate labels to only include labels inside the axes limits
        # Also remove labels where abundance is below s/n threshold
        trunc_labels_df = sub_matched_df[
            (sub_matched_df['m/z'] > xlim[0]) & (sub_matched_df['m/z'] < xlim[1]) & 
            (sub_matched_df['orig_abundance'] > ylim[0]) & 
            (sub_matched_df['orig_abundance'] < ylim[1]) &
            (sub_matched_df['orig_abundance'] > annot_sn_lim)
        ]
        # Plot and label hypothetical masses (if they exist in the axes ranges)
        if not trunc_labels_df.empty:
            mz_vals = trunc_labels_df['m/z'].values
            abunds = trunc_labels_df['orig_abundance'].values
            # Add color generation code here

            # Plot vertical lines
            ax.vlines(x = mz_vals, ymin = 0, ymax=abunds, linewidth = 1)
            # Plot some dots on top too
            ax.scatter(x = mz_vals, y=abunds, s = 5)
            # Label all the points
            _label_point(x = trunc_labels_df['m/z'], y = trunc_labels_df['orig_abundance'], 
                val = trunc_labels_df['ion_name'], ax = ax)

        
        

    return(fig, axs)