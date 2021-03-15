"""Project paths definitions.

Universal project paths are defined in this file.

To use them in notebooks or scripts:
    from msms_structure_annot.paths import {____} <- put path name of interest here

Append new directories and use them like:
    new_dir = data_dir / "folder_name"

"""

from pyprojroot import here

# Path names
root = here(project_files=[".git"])

notebooks_dir = root / "notebooks"
script_dir = root / "scripts"
data_dir = root / "data"
reports_dir = root / "reports"  # Output directory for reports