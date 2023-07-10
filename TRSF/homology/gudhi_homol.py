'''
File: persistent_digram.py
Author: Rhys Shaw
E-mail: rhys.shaw@bristol.ac.uk
Date: 26-04-2023
Description: This file contains the functions for the calculation of the persistent diagram of an image.

'''

import numpy as np
import gudhi
import pandas as pd


def gudhi_cubical_complex(img):
    '''
    Compute the persistence diagram of a single channel image using Gudhi.
    Checks if the image is positive or negative and computes the persistence diagram accordingly.

    Parameters
    ----------
    img: ndarray (M, N)
        An array of single channel image data
        Infinite entries correspond to empty pixels

    Returns
    -------
    dgm: ndarray (N, 2)
        A persistence diagram

    '''
    max_abs = np.max(np.abs(img))

    if np.max(img) == max_abs:
        # image is positive
        img = -img
    elif np.max(-img) == max_abs:
        # image is negative
        pass
    else:
        raise ValueError('Image is not positive or negative')

    return gudhi.CubicalComplex(dimensions=img.shape,top_dimensional_cells=img.flatten()).persistence()

def dgm_to_lifetime(dgm):
    '''
    Compute the lifetime of a persistence diagram

    Parameters
    ----------
    dgm: ndarray (N, 2)
        A persistence diagram

    Returns
    -------
    lifetime: ndarray (N,)
        The lifetime of each point in the persistence diagram
    '''
    return dgm[:,1] - dgm[:,0]


def format_gudhi_dgm(dgm,dim):
    '''
    function to format the output of gudhi's persistence function

    Parameters
    ----------
    dgm : list
        list of birth and death times of the persistence diagram.
    dim : int
        dimension of the persistence diagram you want.

    Returns
    -------
    dgm_alt : ndarray
        array of birth and death times of the persistence diagram.

    '''
    Birth_Death = pd.DataFrame(dgm)

    Birth_Death = Birth_Death.where(Birth_Death[0]==dim).dropna()[1].to_numpy()
    Death = []
    Birth = []

    for i in range(0, len(Birth_Death)):
        Birth.append(Birth_Death[i][0])
        Death.append(Birth_Death[i][1])

    dgm_alt = np.array([Birth, Death]).T


    # remove infinities
    dgm_alt = dgm_alt[~np.isinf(dgm_alt).any(axis=1)]

    return dgm_alt