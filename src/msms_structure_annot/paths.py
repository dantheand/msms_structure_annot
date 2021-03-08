# Designate universal project paths in this file
#   import them into notebooks or scripts with:
#       from msms_structure_annot.paths import {____} <- put path name here
from pyprojroot import here

root = here(project_files=[".git"])

notebooks_dir = root / "notebooks"
data_dir = root / "data"