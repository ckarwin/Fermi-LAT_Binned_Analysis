import numpy as np
import matplotlib.pyplot as plt
from fermipy.gtanalysis import GTAnalysis

# The following commands were ran from a subdirectory of the main analysis directory,
# and thus the paths have "../"

c = np.load('../fit_model.npy', allow_pickle=True).flat[0]

# Fractional count residuals:
total_model_counts = c['roi']['model_counts']
total_counts = c['roi']['counts']
energies = c['roi']['energies']

# Get energy data points.
# The energy array gives the bin edges. 
# Calculate geometric mean and xerror.
energy_mean_list = np.sqrt(energies[:-1]*energies[1:])
err_low_list = energy_mean_list - energies[:-1]
err_high_list = energies[1:] - energy_mean_list

# Get contribution from point sources:
src_list = sorted(c['sources'].keys())
for i in range(0,len(src_list)-2):
    #print(src_list[i])
    src_counts = c['sources'][src_list[i]]['model_counts']    
    src_ts = c['sources'][src_list[i]]['model_counts']
    if i == 0:
        total_src_counts = src_counts
    else:
        total_src_counts += src_counts

# Create subplots:
fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(9, 9),
    sharex=True,
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
)


# plot observed counts:
xerr = np.array([err_low_list,err_high_list])
yerr = np.sqrt(total_counts)
ax1.loglog(energy_mean_list, total_counts, marker='s', color='black', ls='', label='total observed')
ax1.errorbar(energy_mean_list, total_counts, xerr=xerr, yerr=yerr, marker='s', color='black', ls='', label='_nolabel_')

# plot model counts:
ax1.loglog(energy_mean_list,total_model_counts, ls='-', lw=2, color='black', label='total predicted')

# plot galdiff and iso:
galdiff_counts = c['sources']['galdiff']['model_counts']
isodiff_counts = c['sources']['isodiff']['model_counts']
ax1.loglog(energy_mean_list,galdiff_counts, ls='--', lw=2, color='cornflowerblue', label='galdiff')
ax1.loglog(energy_mean_list,isodiff_counts, ls='-.', lw=2, color='orange', label='isodiff')

# plot source counts:
ax1.loglog(energy_mean_list,total_src_counts, ls=':', lw=2, color='green', label='point sources')

# Sanity check:
# Plot sum:
#model_sum = galdiff_counts + isodiff_counts + total_src_counts
#ax1.loglog(energy_mean_list,model_sum, ls='--', lw=2, color='cyan', label='sum')

# fractional residuals:
resid = (total_counts - total_model_counts) / total_model_counts
resid_err = (np.sqrt(total_counts)/total_model_counts)

# Plot kwargs:
ax1.set_xlim(3e2,0.8e6)
ax1.set_ylabel("Counts", fontsize=16)
ax1.tick_params(axis='both', labelsize=14)
ax1.legend(loc=1, frameon=False, fontsize=10)
ax1.grid(ls="--",color="grey",alpha=0.2)

ax2.semilogx(energy_mean_list, resid, ls='', color='black', marker='s')
ax2.errorbar(energy_mean_list, resid, xerr=xerr, yerr=resid_err, ls='', color='black', marker='s')
ax2.axhline(y=0, ls='-',color='grey',lw=1)
ax2.set_ylim(-0.25,0.25)
ax2.set_xlabel("Energy [MeV]", fontsize=16)
ax2.set_ylabel("(Data - Model)/Model", fontsize=16)
ax2.tick_params(axis='both', labelsize=14)
ax2.grid(ls="--",color="grey",alpha=0.2)
plt.savefig("counts_spectra.png")
plt.show()
