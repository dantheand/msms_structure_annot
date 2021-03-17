# MS/MS hypothetical structure annotation

## Overview

This package allows you to propose structures for modified peptides with unknown modification patterns using their MS/MS spectra. It also allows you to annotate the MS/MS spectra of modified peptides with these hypothetical structures. The package generates multiple possible hypothetical structures, computationally fragments them, and then aligns those hypothetical fragments to the observed ions in the MS/MS data to give a best estimate as to where the modifications are taking place on the peptide.

This package was primarily coded to support the work presented here:
*\<cite papers when they come out\>*

## Using this project

The proposed workflow is meant to balance reproducibility and tinkering required for MS/MS data.

### Proposed workflow

- Put all MS/MS data pertaining to a given experiment into a folder in `/data/<exp_name>`
- Copy a `msms_match_template.ipynb` template notebook from `/notebooks/templates` into `/notebooks` and rename it something useful like `<exp_name>.ipynb`
- Provide parameters for analysis in the new notebook
- Run the notebook and export a report into a folder in `/reports/<exp_name>/reports001/`
- If you want to try different parameters / variations on the same experiment, use the same `<exp_name>.ipynb` notebook, but increment the reports number to output to a different location
  - Take notes on what you're changing in the "Notes" section at the top of the notebook
  - At the end of the notebook, the entire notebook is exported so you preserve the modifications you made in the final report

## Quickstart / Installation

### Binder image

If you just want to take a quick peek at the code, use a Binder instance to interactively explore the example Jupyter notebooks in `/notebooks`. ([binder link](https://mybinder.org/v2/gh/dantheand/msms_structure_annot/2c3342c374ff8129687ca0300d3b27541448a8ec))

### First-time use: create and install conda environment

Download the repository (e.g., on Github use the green "Clone or download" button, then "Download ZIP").

Navigate to the project root directory and run this code in an Anaconda terminal:

```shell
conda env create -f ./environment.yml msms_structure_annot-env
pip install -e ./src
```

This will create the conda environment `msms_structure_annot-env` with all the required dependencies. It also installs the custom package in "editable mode".

### Later uses

After this initial install, you only need to activate the environment before using this package:

```shell
conda activate msms_structure_annot-env
```

### Troubleshooting

If you get the error `Cannot find module: msms_structure_annot`, make sure your Jupyter kernel is using the `msms_structure_annot-env` environment:

```shell
# assuming you have already activated your environment,
python -m ipykernel install --user --name msms_structure_annot-env
```

## More details

### How it works

- You provide:
  - Input of one or more ms/ms files
  - Linear peptide sequence
  - Knowledge of potential modification types / locations
- It does:
  - MS/MS spectra s/n filtering
  - Generates a series of hypothetical modified peptide structures
  - Generates fragmentation profiles for those hypothetical structures
  - Maps observed MS/MS peaks onto each hypothetical structure
  - Provides metrics to score which hypothetical structure is most likely
- It outputs:
  - Plots of ms/ms spectra with masses from hypothetical structure fragment masses mapped onto it
  - Tables of matched masses
  - Tables of hypothetical structures and their scores
  - The Jupyter notebook used to make a given report

### Inputs
- All currently provided in the notebook

#### File names and locations

- data files location
- export location

#### MS/MS data

- One or more MS/MS spectra as tab-separated or comma separated file for a compound:
  - Column 1: m/z
  - Column 2: Abundance
- When using multiple spectra (e.g. when trying multiple fragmentation strengths),
  - Name files as `ms[0-9].csv` and put them into the same folder in `/data/<experiment-name>`

#### Modification information

- For modification type:
  - mass shift
  - potential modified residue locations
  - total number of modifications (you should know this from the total mass shift of the selected ion)

```python
# Original AA sequence
parent_seq = 'GCMSKELEKVLESSSMAKGDGWKVMAKGDGWE' # Will be referred to as one-indexed from here on

# Define N and C-term modifications and their mass shifts
N_term_mod = 1.0078
C_term_mod = 17.0027
proton_m = 1.0078
# Define number of charges to calculate m/z values for
charges = [1,2,3]

ptm_dict = {
    'name': [
        'rSAM thioether',  # Name of the modification type
    ],
    'm_shift': [  # Mass shift for a given modification
        -1.007, 
    ], 
    'num_mods': [ # Total number of modifications observed
        2, 
    ], 
    'poss_mod_pos': [ # Potential modification positions (one-indexed)
        [18,22,15,5],
    ],
    'type': [  # Type of modification (ring or point); (ring feature not currently implemented)
        'point',
    ]
}
```

#### Spectra processing parameters

```python
tol = 0.05 # mass-deviation tolerance; Default 0.4
sn_thr = 0.1 # signal-to-noise threshold (must be sn_thr times above background for ion to count); Default 5
N = 12800 # Number of sections to split m/z datapoints into when calculating background values; Default 500
upper_lim = 10 # limit to deviance above average ms ion intensity to set ion values to a limit; Default 50
```

### Outputs

#### Ranked hypothetical structures

Based on various metrics, the likelihood of the peptide being a specific structure is scored. Higher scores are more likely structures

![Hypothetical structure ranking example.](/docs/images/hs_rank_example.png)

#### Annotated MS/MS spectra

After choosing one hypothetical structure, you can map the hypothetical ions (blue) onto the observed ions for spectra.

![MS/MS annotation.](/docs/images/msms_matched_example.png)

## Future features that would be helpful

- Expand the fragmentation code to work for ring structures.

## Acknowledgements

This work heavily leans on the methods proposed and developed by the van der Donk group here: [10.1073/pnas.1406418111](https://doi.org/10.1073/pnas.1406418111)
