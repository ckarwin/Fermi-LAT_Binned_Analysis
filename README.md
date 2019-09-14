# Fermi-LAT Binned Analysis

### Introduction

This repositroy contains the main analysis script that was used for the following studies:

1. Dark Matter Interpretation of the *Fermi*-LAT Observation Toward the Galactic Center ([link](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.95.103005))

2. *Fermi*-LAT Observations of Gamma-Ray Emission Toward the Outer Halo of M31 ([link](https://iopscience.iop.org/article/10.3847/1538-4357/ab2880))

The module requires installation of the Fermi Science Tools (available [here](https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/)).

The code provided here essentially streamlines the analysis procedure that's detailed in the binned likelihood tutorial (available [here](https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/binned_likelihood_tutorial.html)). The code is easily generalizable, although it will require some modifications of input directories and files. 


### Module name: 
       analysis_module

### Classes: <br/>
       1) Analysis: superclass 
       2) Iteration(Analysis) 
       3) PyLikelihood(Iteration,Analysis) 
       4) Plots(Analysis) 
      
### Purpose:
      Perform a binned analysis of the Fermi-LAT data.
      
### Index of Functions:
       Analysis Class(superclass):
               -gtselect()
               -maketime()
               -cmap()
               -ccube()
               -expCube()      
               -expMap()
               -srcMaps(xml_file, outfile)
               -run_gtlike(xml_file, outfile)
               -model_map(src_file, xml_file, outfile, out_type)
               -Residual_map(cmap, model, outfile, operation)
               -Residual_Map2(cmap, model, ccube, modelcube)
               -counts_map_small
               -Make_TS(model_file_xml, xref, yref, nxpix, nypix, bin_size, outfile)
               -gtobsSim(infile, srclist, seed, name)

       Iteration(Analysis):
               -ff_source(source, ff)
               -ff_by_name(root, name_list, ff)
               -iterate(starting_xml_infile, src_file, mass) 

       PyLikelihood(Iteration, Analysis)
               -run_gtlike(src_file, xml_infile, xml_outfile, mass, pass_count)
               
       Plots(Analysis, Iteration)
               -plot_main(x_list, y_list, savefig, Title, x_label, y_label)
               -Counts_spectra_plots()
               -calculate_flux(model_file, exp_file, name, low, high)
               -pylikelihood_flux(src, xml, source_name_list)
               -pylikelihood_counts(input_list, source_plot_list)
               -pyLikelihood_fractional_residuals(src_file,xml_file)
               
### Methodology 
       For any given analysis do the following:
       
       1) Define a new directory
       2) Make a parameter card
               -This file defines all the neccessary parameters for the analysis
               -The parameter card must follow the same format as that in parameters.txt
       3) Make a TS file
               -This file contains the point sources to be iterated over
               -The TS file must follow the same format as TS.txt
       4) Make a client code
               -This is a python script that specifies the neccessary functions to be called
               -The file client.py can be used as a template

