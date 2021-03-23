"""Project paths definitions.

Universal project paths are defined in this file.

To use them in notebooks or scripts:
    from msms_structure_annot.paths import {____} <- put path name of interest here

Append new directories and use them like:
    new_dir = data_dir / "folder_name"

"""

from pyprojroot import here

# Get root directory by looking for the environment file
root = here(project_files=["environment.yml"])

notebooks_dir = root / "notebooks"
script_dir = root / "scripts"
data_dir = root / "data"
reports_dir = root / "reports"  # Output directory for reports
test_data_dir = root / "src" / "tests" / "tests_data" # Directory where example data is stored for tests