# author: Rhys Shaw
# date: 16/06/2023
import gudhi
# import GUDHI library this is the topology library, you will need to pip install this
# https://gudhi.inria.fr/python/latest/installation.html
import numpy as np
import scipy.stats as stats
import pandas as pd

def gudhi_cubical_complex(img):
    '''
    Compute the persistence diagram of a single channel image using Gudhi

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
    return gudhi.CubicalComplex(dimensions=img.shape,top_dimensional_cells=img.flatten()).persistence()


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


def estimate_noise_level(image,plot=False):
    '''
    Estimate the noise level of an image using the persistence diagram of the 0th homology group.

    Parameters
    ----------
    image : ndarray (M, N)
        An array of single channel image data
        Infinite entries correspond to empty pixels
    plot : bool, optional
        Whether to plot the model dist vs real. The default is False.

    Returns
    -------
    mean : float
        Mean of the fitted skew normal distribution.
    mode : float
        Mode of the fitted skew normal distribution.
    std_left : float
        Standard deviation of the left side of the fitted skew normal distribution.
    std_right : float
        Standard deviation of the right side of the fitted skew normal distribution.

    '''

    # note that here we invert the image.
    dgm = gudhi_cubical_complex(-image)
    # select only the 0th dimension homology group.
    dgm_0 = format_gudhi_dgm(dgm,0)
    # remove the infinite points
    dgm_0 = -dgm_0[dgm_0[:,0] != np.inf]

    briths = dgm_0[:,0]
    # remove any births above 0.5
    briths = briths[briths<0.5]
    # remove the nan values
    briths = briths[~np.isnan(briths)]
    # fit this to a sqewed normal distribution


    params = stats.exponnorm.fit(briths)
    mean,var,skew,kurt = stats.exponnorm.stats(*params,moments='mvsk')
    meadian = stats.exponnorm.median(*params)
    mean = mean
    var = var
    std = np.sqrt(var)

    peak_freq_of_births = np.histogram(briths, bins=100)[0].argmax()
    peak = np.histogram(briths, bins=100)[1][peak_freq_of_births]

    if plot:
        import matplotlib.pyplot as plt
        # plot the brith frequency over birth value for the 0th homology group
        plt.figure(figsize=(10,5))
        # Plot the skew-normal distribution
        x = np.linspace(np.min(briths), np.max(briths), 1000)
        pdf = stats.exponnorm.pdf(x, *params)
        plt.plot(x, pdf, 'k-', lw=2, label='skewnorm pdf')
        # plot the histogram
        plt.hist(briths, bins=100, density=True, alpha=0.6, color='g')
        # plot lin where the mean is
        plt.axvline(mean, color='k', linestyle='dashed', linewidth=1)
        plt.axvline(meadian, color='k', linestyle='dashed', linewidth=1)
        plt.xlabel('$log_{10}$($Birth Value of the 0th Homology Group)')
        plt.ylabel('Frequency')
        plt.show()

    return mean, var, std, meadian, peak
