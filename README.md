# *Fermi*-LAT Binned Analysis

### Introduction

This repository contains the main analysis script that was used for the following studies:

1. Dark Matter Interpretation of the *Fermi*-LAT Observation Toward the Galactic Center ([link](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.95.103005))

2. *Fermi*-LAT Observations of Gamma-Ray Emission Toward the Outer Halo of M31 ([link](https://iopscience.iop.org/article/10.3847/1538-4357/ab2880))


## Required Software <br />
The module requires installation of the Fermi Science Tools (available [here](https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/)).

## Purpose <br />
The code streamlines the analysis procedure that's detailed in the binned likelihood tutorial (available [here](https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/binned_likelihood_tutorial.html)). 

## Getting Help <br />
For any help/problems with running the code please contact Chris Karwin at: ckarwin@clemson.edu. 

## Quickstart Guide <br /> 
<pre>
1. Clone repository:
  $ git clone https://github.com/ckarwin/Fermi-LAT_Binned_Analysis.git

2. Install with pip:
  $ cd Fermi-LAT_Binned_Analysis
  $ pip install -e .
     
3. Start a new analysis directory, and enter the commmand-line prompt:
  $ make_analysis

4. Specify inputs in inputs.yaml </b>

5. To run the code:  </b>
  - Uncomment the functions inside run_sims.py that you want to run.
  - The code can be ran directly from the terminal or submitted to a batch system.

</pre>
