#imports
from analysis_module import Analysis
from analysis_module import Plots
from analysis_module import pyLikelihood

#define instance of class:
iteration_instance = pyLikelihood('parameters.txt','TS.txt')
instance = Analysis('parameters.txt','TS.txt')
instance_plots = Plots('parameters.txt','TS.txt')

#uncomment to run desired function:
#instance.gtselect()
#instance.maketime()
#instance.cmap()
#instance.ccube()
#instance.expCube()
#instance.expMap()
#instance.srcMaps('M31_input.xml','M31_binned_srcmaps.fits')
#instance.run_gtlike('input_model.xml','fit.xml')
#iteration_instance.iterate('M31_input.xml', 'M31_binned_srcmaps.fits', 'NA')
#instance.model_map('M31_binned_srcmaps.fits','M31_pass10_out.xml','2d_model.fits','cmap')
#instance.model_map('M31_binned_srcmaps.fits','M31_pass10_out.xml','3d_model.fits','ccube')
#instance.Make_TS('pylikelihood_fit.xml',121.1744,-21.5729,140,140,0.2,'TS_Map_final.fits')
#instance_plots.calculate_flux('3D_New_model.fits',-1, 490000)
#instance.counts_map_small()
#instance.Residual_map(data_path + 'M31_cmap_small.fits','2d_model.fits','Residual.fits','minus')
#instance.Residual_Map2(data_path + 'M31_cmap_small.fits','2d_model.fits',data_path + 'M31_ccube.fits','3d_model.fits')
#instance_plots.pylikelihood_flux('M31_binned_srcmaps.fits','M31_pass10_out.xml',look)
#instance_plots.pyLikelihood_counts([['M31_binned_srcmaps.fits','M31_pass10_out.xml']],look)
#instance_plots.pyLikelihood_fractional_residuals('M31_binned_srcmaps.fits','M31_pass10_out.xml')
