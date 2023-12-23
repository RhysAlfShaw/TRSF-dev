from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from DRUID import sf
from matplotlib import colors

PATH_image = '/data/typhon2/Rhys/data/KiDS/ADP.2019-02-11T13:02:26.713.fits'

hdulist = fits.open(PATH_image)
image = hdulist[0].data

# source find on the above image c

xmin = 10600
xmax = 11600
ymin = 10000
ymax = 11000
image = image[xmin:xmax,ymin:ymax]
# theshold the image to the background level
findmysources = sf(image=image,image_PATH=None,mode='optical',area_limit=5,smooth_sigma=1.5,nproc=1,GPU=True)
findmysources.set_background(detection_threshold=2,analysis_threshold=2,mode='Radio')
findmysources.phsf()                                           # current bottle neck probaly in the getting mask function. (GPU?) acceleration Here?
#findmysources.source_characterising(use_gpu=True)              # GPU acceleration here?
#findmysources.plot_sources(cmap="gray",figsize=(10,10),norm=colors.LogNorm(clip=True,vmin=1E-13,vmax=1E-9))