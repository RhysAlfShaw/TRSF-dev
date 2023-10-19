"""
author: Rhys Shaw
email: rhys.shaw@bristol.ac.uk
date: 26-06-2023
description: This file contains the functions for the creation of regions.

"""

from skimage.measure import label, regionprops
import numpy as np

def create_regions(img, intentsity_img):
    """
    Create a list of regions from a binary image
    
    Parameters
    ----------
    img: ndarray (M, N)
        A binary image
    min_size: int
        The minimum size of a region to be included in the list of regions
        
    Returns
    -------
    regions: list of Region objects
        A list of regions
    
    labels: ndarray (M, N)
        A labelled image
    

    """
    

    # Label regions
    labels = label(img)
    
    # Create a list of regions

    regions = regionprops(labels,intensity_image=intentsity_img)
        
    return labels, regions

def remove_low_flux_regions(img, regions, threshold=None):
    '''
    Remove regions with peak flux < threshold

    Parameters
    ----------
    img: ndarray (M, N)
        A binary image
    regions: list of Region objects
        A list of regions
    threshold: float
        The threshold for the peak flux

    Returns
    -------
    regions: list of Region objects
        A list of regions

    '''

    props = regions

    # for each region evaluate the peak flux in the original image
    peak_flux = []
    for prop in props:
        minr, minc, maxr, maxc = prop.bbox
        peak_flux.append(img[minr:maxr,minc:maxc].max())

    # create cutout from original image for each region
    cutouts = []
    for prop in props:
        minr, minc, maxr, maxc = prop.bbox
        cutouts.append(img[minr:maxr,minc:maxc])

    # evaluate the peak flux in each cutout
    peak_flux_cutout = []
    for cutout in cutouts:
        peak_flux_cutout.append(cutout.max())

    # Threshold the peak flux in the cutout
    
    # Create a suggested threshold but allow input

    if threshold is None:
        threshold = (np.sqrt(np.mean(peak_flux))**2 )*2

    print('Threshold',threshold)

 
    print(np.mean(peak_flux_cutout))
    # remove region with index of cutout with peak flux < 0.0001
    region_idexes_to_remove = []
    for i in range(len(peak_flux_cutout)):
        if peak_flux_cutout[i] < threshold:
            region_idexes_to_remove.append(i)

    # remove regions
    props = np.delete(props,region_idexes_to_remove)    

    # create new mask from regionprops
    mask = np.zeros(img.shape)
    for prop in props:
        minr, minc, maxr, maxc = prop.bbox
        mask[minr:maxr,minc:maxc] = 1


    return props, mask

