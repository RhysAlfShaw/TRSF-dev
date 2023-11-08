from matplotlib.pylab import plt
import numpy as np
from astropy.table import Table

def calculate_beam():
    import numpy as np
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
    beam_size_arcseconds = 3*3600
    BMAJ_oversampled_spacial_width = (BMAJ**2 + beam_size_arcseconds**2)**0.5
    BMAJ = BMAJ_oversampled_spacial_width/arcseconds_per_pixel
    
    beam_size_arcseconds = 3*3600
    BMIN_oversampled_spacial_width = (BMIN**2 + beam_size_arcseconds**2)**0.5
    BMIN = BMIN_oversampled_spacial_width/arcseconds_per_pixel
    
    BMAJp = BMAJ
    BMINp = BMIN
    print('BMAJp = ',BMAJp)
    print('BMINp = ',BMINp)
    return np.pi * (BMAJ)*(BMIN) / (4*np.log(2))

catalogue = Table.read('DRUIDsimulation_catalogue.fits')
truth = Table.read('Simulation_Catalogue_Truth.fits')
Beam = calculate_beam()
matched = Table.read('Matched_DRUID_Truth_Sim.fits')
matched['Flux_total'] = matched['Flux_total'] * Beam
matched['f_diff'] = (matched['flux'] - matched['Flux_total'])/matched['flux']
# remove any source that arent class 0

matched = matched[matched['Class'] == 0] # maybe remove?

plt.style.use('/Users/rs17612/Documents/GitHub/style/mplrc_sotiria.mplstyle')

plt.figure(figsize=(5,5))
plt.axhline(0,color='k',linestyle='--')
plt.scatter(matched['Flux_total'], abs(matched['f_diff']),marker='x',s=10)
#plt.axvline(matched['Noise'].mean(),color='r')
plt.xlabel('True Total Flux')
plt.ylabel('|TF - SF| / TF')
plt.xscale('log')
#plt.yscale('log')
plt.show()

plt.scatter(matched['flux'], matched['Flux_total'],marker='x',s=30,color='black')
plt.plot([0,1e9],[0,1e9],color='k',linestyle='--',marker='')

plt.xlabel('True Total Flux')
plt.ylabel('DRUID Total Flux')
plt.xscale('log')
plt.yscale('log')
plt.show()


matched['pf_diff'] = (-matched['peak_flux'] + matched['Flux_peak'])/matched['peak_flux']

plt.figure(figsize=(5,5))
plt.axhline(0,color='black',linestyle='--')
plt.scatter(matched['peak_flux'], matched['pf_diff'],marker='x',s=30)
#plt.axvline(matched['Noise'].mean()*15,color='r')
plt.xlabel('True Peak Flux')
plt.ylabel('(TFp - PFp) / TFp')
plt.xscale('log')
#plt.yscale('log')
plt.show()

# completness and reliability
import numpy as np

bins = np.logspace(1.5,np.log10(matched['flux'].max()),50)

mathced_hist = np.histogram(matched['flux'],bins=bins)[0]
truth_hist = np.histogram(truth['flux'],bins=bins)[0]
DRUID_hist = np.histogram(catalogue['Flux_total']*Beam,bins=bins)[0]

completness = mathced_hist/truth_hist
mathced_hist = np.histogram(matched['Flux_total'],bins=bins)[0]
relibility = mathced_hist/DRUID_hist

# plot completness and reliability
plt.figure(figsize=(5,5))
plt.plot(bins[:-1],completness,label='Completness',marker='',zorder=5)
plt.plot(bins[:-1],relibility,label='Reliability',marker='',zorder=5)
plt.axhline(1,color='k',ls='--',zorder=0)
plt.xscale('log')
plt.xlabel('Flux')
plt.ylabel('Fraction')
plt.legend()
plt.show()