# MS/MS hypothetical structure annotation

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dantheand/msms_structure_annot/HEAD)

## Overview

This package allows you to propose structures for modified peptides with unknown modification patterns using their MS/MS spectra.

Basic features:

- MS/MS spectra cleanup
- Computational generation of modified peptides "hypothetical structures" and their fragmentation patterns
- Matching of observed masses with computationally-generated spectra
- Likelihood scoring of hypothetical structures given experimental MS/MS spectra

This package was primarily coded to support the work presented here:

Glassey, E., King, A.M., Anderson, D.A., Zhang, Z., & Voigt, C.A. (2022).
**Functional expression of diverse post-translational peptide-modifying enzymes in Escherichia coli under uniform expression and purification conditions.**
PLOS ONE https://doi.org/10.1371/journal.pone.0266488

## Using this project

### Project organization

- `/data`: raw ms/ms data is stored here (never altered)
- `/notebooks`: where analysis notebooks are created and stored; most of the heavy-lifting code is factored out into the source code
- `/reports`: exported reports from analysis notebooks go here
- `/src/msms_structure_annot`: the lightweight package used to contain and organize the custom source code in modules

### Proposed workflow

The proposed workflow is meant to balance reproducibility while still allowing the tinkering required for MS/MS data:

- Put all MS/MS data pertaining to a given experiment into a folder in `/data/<exp_name>`
- Copy a `HalA2_example.ipynb` template notebook from `/notebooks` and rename it something useful like `<exp_name>.ipynb`
- Provide parameters for analysis in the new notebook
- Run the notebook and export a report into a folder in `/reports/<exp_name>/reports001/`
- If you want to try different parameters / variations on the same experiment, use the same `<exp_name>.ipynb` notebook, but increment the reports number to output to a different location
  - Take notes on what you're changing in the "Notes" section at the top of the notebook
  - At the end of the notebook, the entire notebook is exported so you preserve the modifications you made in the final report

## Quickstart / Installation

### Binder image

If you just want to take a quick peek at the functionality, use a Binder instance to interactively explore the example Jupyter notebooks in `/notebooks`. I recommend HalA2. ([binder link](https://mybinder.org/v2/gh/dantheand/msms_structure_annot/2c3342c374ff8129687ca0300d3b27541448a8ec))

### First-time use: install conda environment and package

Download the repository (e.g., on Github use the green "Clone or download" button, then "Download ZIP").

Navigate to the project root directory and run this code in an Anaconda terminal:

```shell
conda env create -f ./environment.yml msms_structure_annot-env
conda activate msms_structure_annot-env
pip install -e ./src
```

This will create the conda environment `msms_structure_annot-env` with all the required dependencies. It also installs the custom package in "editable mode".

Then open a Jupyter server in the `msms_structure_annot-env` environment and give it a go:

```shell
conda activate msms_structure_annot-env
jupyter notebook
```

### Later uses

After this initial install, you only need to activate the environment before opening a Jupyter server:

```shell
conda activate msms_structure_annot-env
jupyter notebook
```

### Troubleshooting

If you get the error `Cannot find module: msms_structure_annot`, make sure your Jupyter kernel is using the `msms_structure_annot-env` environment:

```shell
# assuming you have already activated your environment,
python -m ipykernel install --user --name msms_structure_annot-env
```

## More details

### What it's actually doing

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
- Include isotopic peaks (currently only does monoisotopic masses, which isn't the most prevalent in large fragments)

## Acknowledgements

This work heavily leans on the methods proposed and developed by the van der Donk group here: [10.1073/pnas.1406418111](https://doi.org/10.1073/pnas.1406418111)

The adjustText Python package was use to ease annotation (https://github.com/Phlya/adjustText) doi: [10.5281/zenodo.3924114](https://zenodo.org/badge/latestdoi/49349828)
