"""
Author: Rhys Shaw
E-mail: rhys.shaw@bristol.ac.uk
Date: 26-06-2023
Description: This file contains the functions for the Gaussian fitting of images.

"""

import numpy as np
from scipy.optimize import curve_fit


def gaussian_2d(xy,amplitude,x0,y0,sigma_x,sigma_y,theta):
    '''
    2D Gaussian function
    
    Parameters
    ----------
    xy : tuple of array_like
        x and y coordinates
    amplitude : float
        Gaussian amplitude
    x0 : float
        x coordinate of the center
    y0 : float
        y coordinate of the center
    sigma_x : float
        standard deviation in x
    sigma_y : float
        standard deviation in y
    theta : float
        rotation angle in radians

    Returns
    -------
    g : array_like
        Gaussian function evaluated at (x, y)

    Examples
    --------
    >>> import numpy as np
    >>> from egauss_fit import gauss_2d
    >>> x, y = np.meshgrid(np.arange(10), np.arange(10))
    >>> g = gauss_2d((x, y), 1, 5, 5, 1, 1, 0)
    >>> g[5, 5]
    1.0
        
    '''

    x, y = xy
    
    a = np.cos(theta)**2 / (2 * sigma_x**2) + np.sin(theta)**2 / (2 * sigma_y**2)
    
    b = np.sin(2 * theta) / (4 * sigma_x**2) - np.sin(2 * theta) / (4 * sigma_y**2)
    
    c = np.sin(theta)**2 / (2 * sigma_x**2) + np.cos(theta)**2 / (2 * sigma_y**2)

    Z = amplitude * np.exp(-(a * (x - x0)**2 + 2 * b * (x - x0) * (y - y0) + c * (y - y0)**2))
    return Z




def fit_gaussian_2d(img, maxfev,inital_guess=None):
    
    '''
    Fit a 2D Gaussian function to an image using scipy.optimize.curve_fit

    Parameters
    ----------
    img : array_like
        Image to fit

    Returns
    -------
    popt : array_like
        Optimal values for the parameters so that the sum of the squared error
        of f(xdata, *popt) - ydata is minimized
    fitted_img : array_like
        Fitted image

    '''


    def gauss_2d_m(xy,amplitude,x0,y0,sigma_x,sigma_y,theta):
        # Flatten the image to 1D.
        return gaussian_2d(xy,amplitude,x0,y0,sigma_x,sigma_y,theta).ravel()

    # bring scale of image up to 1
    
    img_temp = np.copy(img/img.max() * 100)           

    #the optimzation algorithm is more stable if the image is scaled up to 0-100

    ny, nx = img_temp.shape
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))

    # initial guess of ellipse parameters
    if inital_guess is not None:
        amplitude = img_temp.max()
        x0_guess, y0_guess = img.shape[1]/2, img.shape[0]/2
        sigma_xguess = min([nx/2, ny/2])
        sigma_yguess = min([ny/2, nx/2])
        theta_guess = 0
        p0 = [amplitude, x0_guess, y0_guess, sigma_xguess,sigma_yguess, theta_guess]
    else:
        p0 = inital_guess

    # Perform the fit and return the fitted image
    popt, pcov = curve_fit(gauss_2d_m, (x, y), img_temp.ravel(), p0=p0, maxfev=maxfev)
    
    # scale amplitude back down to original image scale
    popt[0] = img.max()*popt[0] /100         
    
    fitted_img = gaussian_2d((x, y), *popt).reshape(img.shape)

    return popt, fitted_img