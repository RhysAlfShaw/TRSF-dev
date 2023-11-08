import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import pdb
from astropy.table import Table
from scipy.ndimage import gaussian_filter
from astropy.io import fits
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--N', type=int, default=1000, help='Number of sources to simulate.')
parser.add_argument('--xmax', type=int, default=3000, help='Size of the image in x.')
parser.add_argument('--ymax', type=int, default=3000, help='Size of the image in y.')
parser.add_argument('--sn', type=int, default=1, help='simulation number.')
parser.add_argument('--seed', type=int, default=0, help='seed for random number generator.')
parser.add_argument('--save_path', type=str, default='.', help='path to save the simulation to.')

args = parser.parse_args()

# set seed random but repeatable results.
np.random.seed(args.seed)

N = args.N
xmax = args.xmax
ymax = args.ymax

# source positions (x,y) these will be generated randomly

x = np.random.uniform(0, xmax, N)
y = np.random.uniform(0, ymax, N)

flux = np.random.lognormal(-3, 6, N)

Catalogue = Table([x, y, flux], names=('y', 'x', 'flux'))  # note that x and y are flipped this is intentional.

# populate the image with the sources
image = np.zeros((xmax, ymax))
for i in range(N):
    image[int(x[i]), int(y[i])] = flux[i]

# add noise to the image
gaussian_noise = np.random.normal(0, 1, (xmax, ymax))

image_radio = image + gaussian_noise

beam = 6 # FWHM in pixels assume circular beam
image_radio = gaussian_filter(image_radio, sigma=beam/2.355)

peak_flux = []
for i in range(N):
    peak_flux.append(image_radio[int(x[i]), int(y[i])])

Catalogue['peak_flux'] = peak_flux
Catalogue.write('Simulation_Catalogue_Truth'+str(args.sn)+'.fits', format='fits',overwrite=True)

# save image to fits file with astropy.


hdul = fits.PrimaryHDU(image_radio)
# create a header for the fits file
hdul.header['BMAJ'] = beam
hdul.header['BMIN'] = beam 
hdul.header['CDELT1'] = 1  # to ensure that the units are correct int source finder.
hdul.header['BPA'] = 0
hdul.writeto('Simulation_Image_radio'+str(args.sn)+'.fits', overwrite=True)