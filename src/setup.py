# This makes the project pip installable with pip install -e
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="msms_structure_annot", # package name (if it were published to PyPI)
    version="0.0.1",
    long_description = long_description,
    license = 'MIT',
    packages=setuptools.find_packages(),
)