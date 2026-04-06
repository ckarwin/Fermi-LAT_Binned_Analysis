import matplotlib.pyplot as plt
from fermipy.gtanalysis import GTAnalysis
import pyLikelihood
from BinnedAnalysis import *
import pandas as pd
import math
import numpy as np
from lat_binned_analysis import Analysis, Plots
from astropy.io import fits

# The following commands were ran from a subdirectory of the main analysis directory,
# and thus the paths have "../"

instance = Analysis('inputs.yaml')
instance_plots = Plots('inputs.yaml')

# calc total observed flux:
instance_plots.calculate_flux('../ccube_00.fits','../bexpmap_00.fits',"observed", ccube=True)

# calc total model flux:
instance_plots.calculate_flux('model_cube.fits','../bexpmap_00.fits',"model")

# calc isodiff flux:
instance.model_map('../srcmap_00.fits','../fit_model_00_isodiff.xml',\
        '../ltcube_00.fits','../bexpmap_00.fits','isodiff_3d_model.fits','ccube')
instance_plots.calculate_flux('isodiff_3d_model.fits','../bexpmap_00.fits',"isodiff")

# calc galdiff flux:
instance.model_map('../srcmap_00.fits','../fit_model_00_galdiff.xml',\
        '../ltcube_00.fits','../bexpmap_00.fits','galdiff_3d_model.fits','ccube')
instance_plots.calculate_flux('galdiff_3d_model.fits','../bexpmap_00.fits',"galdiff")

# Plot:
fig, ax = plt.subplots(figsize=(8, 6))

# observed:
hdu = fits.open('../ccube_00.fits')
energy_A = hdu[1].data
xerr_low_list = []
xerr_high_list = []
for E in range(0,len(energy_A)):
    energy = ((energy_A[E][2]/1000)*(energy_A[E][1]/1000))**0.5
    err_low = energy - energy_A[E][1]/1000
    err_high = energy_A[E][2]/1000 - energy
    xerr_low_list.append(err_low)
    xerr_high_list.append(err_high)
xerr = np.array([xerr_low_list,xerr_high_list])

df = pd.read_csv("observed_flux.dat", skiprows=2, delim_whitespace=True)
plt.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="black",ls="",marker="s", label="total observed")
plt.errorbar(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],xerr=xerr,color="black",ls="",marker="s", label="_NONE_")

# model:
df = pd.read_csv("model_flux.dat", skiprows=2, delim_whitespace=True)
plt.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="black",ls="-",label="total model")

# galdiff:
df = pd.read_csv("galdiff_flux.dat", skiprows=2, delim_whitespace=True)
plt.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="cornflowerblue",ls="--",label="galdiff")

# isodiff:
df = pd.read_csv("isodiff_flux.dat", skiprows=2, delim_whitespace=True)
plt.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="darkorange",ls="--",label="isodiff")

plt.legend(loc=1,fontsize=12,ncol=1,frameon=False)
plt.xlabel("Energy [MeV]", fontsize=16)
plt.ylabel("$\mathrm{E^2 dN/dE \  [MeV \ cm^{-2} \ s^{-1}]}$", fontsize=16)
plt.ylim(1e-9,10)
#plt.ylim(1e-3,1e-1)
plt.xlim(2.5e2,0.9e6)
plt.grid(ls="--",color="grey",alpha=0.2)
ax.tick_params(axis='both', labelsize=14)
plt.savefig("flux_spectra.png")
plt.show()
