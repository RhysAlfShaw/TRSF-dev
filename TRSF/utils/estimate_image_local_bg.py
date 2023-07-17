"""
author: Rhys Shaw
email: rhys.shaw@bristol.ac.uk
date: 26-06-2023
desc:
    Estimate the local background of an image using median std.
"""
from astropy.stats import mad_std
from TRSF.utils import hom_rms_estimate
import numpy as np

def estimate_bg_madstd(image):
    """
    Estimate the local background of an image using median std.

    Parameters
    ----------
    image: ndarray (M, N)
        The image to estimate the local background of.

    Returns
    -------
    local_bg: ndarray (M, N)
        The local background of the image.

    """
    # Calculate the local background
    local_bg = mad_std(image, ignore_nan=True)
    #print('local_bg', local_bg)
    return local_bg

def estimate_bg_from_homology(image):

    # evaluate madstd
    mad_std = estimate_bg_madstd(image)
    #print('mad_std', mad_std)
    if mad_std == 0 or mad_std == np.nan:
        print('mad_std is 0, evaluating homology')
        # evaluate homology
        local_bg = hom_rms_estimate.estimate_noise_level(image)[0]
        sigma = 1 
    else:
        #print('mad_std is not 0, using mad_std')
        sigma = 5 
        local_bg = mad_std
    return local_bg, sigma


