# MS/MS hypothetical structure annotation

## Overview

SIMPLE OVERVIEW HERE.

### Steps to accomplish this

- takes input of one or more ms/ms files
- takes a peptide sequence
- takes knowledge of potential modification types / locations
- generates a series of hypothetical modified peptide structures
- generates hypothetical fragmentation profiles for those hypothetical structures
- maps observed ms/ms peaks onto each hypothetical structure
- provides metrics to score which hypothetical structure is the most likely
- writes to a report

## Using this project

The methodology for applying this pipeline is meant to provide a streamlined and reproducible workflow, while still allowing the tweaking and freedom required to interface with this kind of data.

## Quickstart / Installation

## Download / install this repository / package ?

Instructions for that.

## First-time use: create and install conda environment

Navigate to the project root directory and run this code:

```shell
conda env create -f ./scripts/environment.yml msms_structure_annot-env
cd src
pip install -e
```

This will create the conda environment `msms_structure_annot-env` with all the required dependencies. It also installs the custom package in "developer mode". After this initial install, you only need to activate the environment before using this package:

```shell
conda activate msms_structure_annot-env
```

If you get the error `Cannot find module: msms_structure_annot`, make sure your Jupyter kernel is using the `msms_structure_annot-env` environment.

## How to use it

- how to use this code
- also expanded Jupyter notebook explaining what each step of the process does

### Basics

#### Inputs

##### MS/MS data

- One or more MS/MS spectra as tab-separated file for a compound (not deconvoluted):
  - Column 1: m/z; Column 2: Abundance
- When using multiple spectra (e.g. when trying multiple fragmentation strengths),
  - Name files as `ms[0-9].csv` and put them into the same folder in `/data`

##### Modification information

- For modification type:
  - mass shift
  - potential modified residue locations
  - total number of modifications (you should know this from the total mass shift of the selected ion)

### Example use

- structural annotation example

## Todo

### Near-term

- Things that need to be done in the script:
  - plots
    - fix and factor out barplots
  - factor out fragmentation code
  - outputs
    - make standard report output and save jupyter notebook outputs

### In the future

- Make it work for ring structures