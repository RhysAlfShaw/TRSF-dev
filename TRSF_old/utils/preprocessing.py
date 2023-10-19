from astropy.io import fits
import numpy as np

'''
File: preprocessing.py
Author: Rhys Shaw
E-mail: rhys.shaw@bristol.ac.uk
Date: 26-04-2023
Description: This file contains the functions for the pre-processing of SKA images.

'''

def load_fits(filename,formatting=True):
    """
    Load a FITS file

    Parameters
    ----------
    filename: str
        The name of the FITS file

    Returns
    -------
    data: ndarray (M, N)
        The data in the FITS file

    """

    # Load the FITS file
    hdul = fits.open(filename)

    # Extract the data
    data = hdul[0].data
    header = hdul[0].header

    if formatting == True:
        # Remove the non used axis
        try: 
            data = np.squeeze(data, axis=(0,1))
        except ValueError:
            print('Error ValueError: data is not 3D!!')
            
    return data, header


class preprocess:
    '''
    Object to preprocess the image.
    - Loading the image
    - Loading the primary beam
    - Regridding the primary beam
    - Applying the primary beam correction
    - Returning Image for further analysis.
    '''
    def __init__(self,Image_Path,formatting):
       
        self.Image_Path = Image_Path
        self.Image_Path = Image_Path
        self.formatting = formatting
        # Run the main function
        self.main()

    def main(self):
        img = load_fits(self.Image_Path,formatting=self.formatting)
        #pb = load_fits(self.Pb_Path)
        self.img, self.header = load_fits(self.Image_Path,formatting=self.formatting)