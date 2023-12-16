from DRUID import sf
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('/Users/rs17612/Documents/GitHub/style/mplrc_sotiria.mplstyle')
from astropy.io import fits

import pandas as pd

PATH_COSMOS_MIGHTEE = '/Users/rs17612/Documents/Radio_Data/MIGHTEE/MIGHTEE_Continuum_Early_Science_COSMOS_r-1p2.app.restored.circ.fits'
PATH_pb = '/Users/rs17612/Documents/Radio_Data/MIGHTEE/MIGHTEE_Continuum_Early_Science_COSMOS_r-1p2_effective_frequency_pbcut_corrected.fits'


findmysources = sf(image=None,image_PATH=PATH_COSMOS_MIGHTEE,
                   pb_PATH=PATH_pb,mode='Radio',cutup=True,
                   cutup_size=250,output=False)

cutouts = findmysources.cutouts

findmysources.set_background(detection_threshold=5,analysis_threshold=2)
findmysources.phsf()
findmysources.source_characterising()
catalogue = findmysources.catalogue
print(catalogue)
flux_corrs = catalogue['Corr_f']

plt.hist(flux_corrs,bins=100)
plt.xlim(0,20)
plt.show()

# Get Given catalogue.

PATH_Cosmo_r1p2_cat = '/Users/rs17612/Documents/Radio_Data/MIGHTEE/MIGHTEE_CONTINUUM_Early_Science_COSMOS_r-1p2.app.restored.circ.gaul.fits'

bdsf_df = pd.DataFrame(fits.open(PATH_Cosmo_r1p2_cat)[1].data)

catalogue_test = catalogue[catalogue['Class'] == 0]

DRUIDFlux = catalogue_test['Flux_total']
BDSfFlux = bdsf_df['Total_flux']

bins = np.logspace(-5,-2,30)

DRUIDFluxVals = np.histogram(DRUIDFlux,bins=bins)[0]
BDSfFluxVals = np.histogram(BDSfFlux,bins=bins)[0]

plt.step(bins[:-1],DRUIDFluxVals,marker='',label='DRUID')
plt.step(bins[:-1],BDSfFluxVals,marker='',label='BDSF')
plt.xscale('log')
plt.show()