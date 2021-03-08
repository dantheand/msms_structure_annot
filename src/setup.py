import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="msms_structure_annot-env", # environment name
    version="0.0.1",         
    packages=setuptools.find_packages(),
)