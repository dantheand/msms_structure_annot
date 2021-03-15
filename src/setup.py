# This makes the project pip installable with pip install -e
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="msms_structure_annot", # package name
    version="0.0.1",         
    packages=setuptools.find_packages(),
)