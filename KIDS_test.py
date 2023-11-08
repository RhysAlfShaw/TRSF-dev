import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib

# Read in the data
PATH_image = '/data/typhon2/Rhys/data/KiDS/ADP.2019-02-11T13:02:26.713.fits'

# Open the fits file
hdulist = fits.open(PATH_image)
#hdulist.info()
image = hdulist[0].data

from TRSF import sf

# source find on the above image
xmin = 10600
xmax = 15600
ymin = 10000
ymax = 15000
image = image[xmin:xmax,ymin:ymax]
findmysources = sf(image=image,image_PATH=None,mode='other',output=False,area_limit=5,smooth_sigma=1.5)
findmysources.set_background(detection_threshold=1,analysis_threshold=1,set_bg=1E-12)
findmysources.phsf(lifetime_limit=0.8E-12)
findmysources.create_polygons()
ploygons_fixed = findmysources.polygons
catalogue = findmysources.catalogue
catalogue.columns
# Open the catalogue
PATH_CATALOGUE = '/data/typhon2/Rhys/data/KiDS/ADP.2019-02-11T13:02:26.716.fits'
from astropy.table import Table
REF_catalogue = Table.read(PATH_CATALOGUE, format='fits')
REF_catalogue.columns

# crop the catalogue to the x and y limits for the x_image and y_image columns
REF_catalogue = REF_catalogue[REF_catalogue['X_IMAGE'] > ymin]
REF_catalogue = REF_catalogue[REF_catalogue['X_IMAGE'] < ymax]
REF_catalogue = REF_catalogue[REF_catalogue['Y_IMAGE'] > xmin]
REF_catalogue = REF_catalogue[REF_catalogue['Y_IMAGE'] < xmax]

# offset the REF_catalogue to the image
REF_catalogue['X_IMAGE'] = REF_catalogue['X_IMAGE'] - ymin
REF_catalogue['Y_IMAGE'] = REF_catalogue['Y_IMAGE'] - xmin
print(len(REF_catalogue))

# save this part of the catalogue

REF_catalogue.write('/data/typhon2/Rhys/data/KiDS/ADP.2019-02-11T13:02:26.716_cropped.fits', format='fits',overwrite=True)

# save the DRUID catalogue at this stage its to see how many matches we get.

catalogue = catalogue.drop(columns=['enclosed_i'])
print(catalogue.columns)

# format the catalogue types 
for column in catalogue.columns:
    catalogue[column] = catalogue[column].astype(float)


catalogue = Table.from_pandas(catalogue)
catalogue.write('/data/typhon2/Rhys/data/KiDS/ADP.2019-02-11T13:02:26.716_DRUID.fits', format='fits',overwrite=True)

# transform the corrds of polygons to the correct coords based on image cutout
polygons = findmysources.polygons
for i in range(len(polygons)):
    for j in range(len(polygons[i])):
        polygons[i][j][0] = polygons[i][j][0] + xmin
        polygons[i][j][1] = polygons[i][j][1] + ymin
        
# save as ds9 region file

findmysources.polygons = polygons
findmysources.save_polygons_to_ds9('/data/typhon2/Rhys/data/KiDS/Galaxy_polys.reg')
