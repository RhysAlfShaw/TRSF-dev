import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

def calculate_beam():
    '''

    Calculates the beam of the image in pixels.

    returns:
        beam sampling factor, and creates object variable BMAJp and BMINp.

    '''
    
    # get beam info from header

    BMAJ = 6
    BMIN = 6
    
    # convert to pixels
    
    arcseconds_per_pixel = 1*3600*-1
    beam_size_arcseconds = 6*3600
    BMAJ_oversampled_spacial_width = (BMAJ**2 + beam_size_arcseconds**2)**0.5
    BMAJ = BMAJ_oversampled_spacial_width/arcseconds_per_pixel
    
    beam_size_arcseconds = 6*3600
    BMIN_oversampled_spacial_width = (BMIN**2 + beam_size_arcseconds**2)**0.5
    BMIN = BMIN_oversampled_spacial_width/arcseconds_per_pixel
    
    BMAJp = BMAJ
    BMINp = BMIN
    print('BMAJp = ',BMAJp)
    print('BMINp = ',BMINp)
    return np.pi * (BMAJ)*(BMIN) / (4*np.log(2))


DRUIDcatalogue = Table.read('/data/typhon2/Rhys/data/SIMS/DRUIDsimulation_catalogue_Combined.fits')
Truth = Table.read('/data/typhon2/Rhys/data/SIMS/Simulation_Catalogue_Truth_Combined.fits')

Beam = calculate_beam()
matched = Table.read('/data/typhon2/Rhys/data/SIMS/Matched_Catalogue.fits')
matched['Flux_total'] = matched['Flux_total'] * Beam
print(len(Truth))
# plot histogram of Truth fluxes
plt.figure(figsize=(5,5))
plt.hist(Truth['flux'],bins=np.logspace(np.log10(Truth['flux'].min()),np.log10(Truth['flux'].max()),50),label='Truth',alpha=0.5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Flux')
plt.ylabel('N')
plt.savefig('data/Truth_flux_dist.png',dpi=300)
plt.show()

matched['f_diff'] = abs(matched['flux'] - matched['Flux_total'])/matched['flux']
matched = matched[matched['f_diff'] < 10]
print(len(matched))
plt.style.use('/home/rs17612/GitHub/mplrc_sotiria.mplstyle')
plt.figure(figsize=(5,5))
plt.axhline(0,color='k',linestyle='--')
plt.scatter(matched['Flux_total'], matched['f_diff'],marker='x',s=10)
#plt.axvline(matched['Noise'].mean(),color='r')
plt.xlabel('True Total Flux')
plt.ylabel('|TF - SF| / TF')
plt.xscale('log')
plt.savefig('data/flux_total_difference.png',dpi=300)
#plt.yscale('log')
plt.show()

plt.figure(figsize=(5,5))
plt.scatter(matched['flux'], matched['Flux_total'],marker='x',s=30,color='black')
plt.axline([0,0],[1,1],color='k',linestyle='--',marker='',zorder=0)
plt.xlabel('True Total Flux')
plt.ylabel('DRUID Total Flux')
plt.xscale('log')
plt.yscale('log')
plt.savefig('data/flux_total.png',dpi=300)
plt.show()

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
ax1.scatter(matched['flux'], matched['Flux_total'], marker='x',s=30,zorder=30)
ax1.axline([0,0],[1,1],color='k',linestyle='--',marker='',zorder=0)
ax1.set_ylabel('Value')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend()
ax2.scatter(matched['flux'], matched['f_diff'],s=5,marker='x',zorder=30)
ax2.axhline(0,color='k',linestyle='--',zorder=0)
ax2.set_xlabel('X-axis Label')
ax2.set_ylabel('Fractional Difference')
plt.savefig('data/Total_flux_difference.png',dpi=300)
plt.show()
matched['pf_diff'] = (matched['peak_flux'] - matched['Flux_peak'])/matched['peak_flux'] 

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
ax1.scatter(matched['peak_flux'], matched['Flux_peak'], marker='x',s=30,zorder=30)
ax1.axline([0,0],[1,1],color='k',linestyle='--',marker='',zorder=0)
ax1.set_ylabel('Value')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend()
ax2.scatter(matched['peak_flux'], matched['pf_diff'],s=5,marker='x',zorder=30)
ax2.axhline(0,color='k',linestyle='--',zorder=0)
ax2.set_xlabel('X-axis Label')
ax2.set_ylabel('Fractional Difference')
plt.savefig('data/peak_flux_difference.png',dpi=300)
plt.show()

bins = np.logspace(0,np.log10(matched['flux'].max()),50)

mathced_hist = np.histogram(matched['flux'],bins=bins)[0]
truth_hist = np.histogram(Truth['flux'],bins=bins)[0]
DRUID_hist = np.histogram(DRUIDcatalogue['Flux_total']*Beam,bins=bins)[0]

completness = mathced_hist/truth_hist
mathced_hist = np.histogram(matched['Flux_total'],bins=bins)[0]
relibility = mathced_hist/DRUID_hist

plt.figure(figsize=(5,5))
plt.plot(bins[:-1],completness,label='Completness',marker='',zorder=5)
plt.plot(bins[:-1],relibility,label='Reliability',marker='',zorder=5)
plt.axhline(1,color='k',ls='--',zorder=0)
plt.xscale('log')
plt.xlabel('Flux')
plt.ylabel('Metric')
plt.legend()
plt.savefig('data/completness_relibility.png',dpi=300)
plt.show()