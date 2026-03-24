# Imports:
from setuptools import setup, find_packages

# Setup:
setup(
    name='lat_binned_analysis',
    version="0.1",
    url='https://github.com/ckarwin/Fermi-LAT_Binned_Analysis',
    author='Chris Karwin',
    author_email='ckarwin@clemson.edu',
    packages=find_packages(),
    description = 'perform lat binned analysis',
    entry_points = {"console_scripts":["make_analysis = lat_binned_analysis.make_new_analysis:main"]}
)
