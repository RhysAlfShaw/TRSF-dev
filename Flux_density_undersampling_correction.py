import numpy as np


def generate_2D_gaussian(A,shape, center, sigma_x, sigma_y, angle_deg=0,norm=True):
    """
        Generate a 2D elliptical Gaussian distribution on a 2D array.

        Parameters:
            shape (tuple): Shape of the output array (height, width).
            center (tuple): Center of the Gaussian distribution (x, y).
            sigma_x (float): Standard deviation along the x-axis.
            sigma_y (float): Standard deviation along the y-axis.
            angle_deg (float): Rotation angle in degrees (default is 0).

        Returns:
            ndarray: 2D array containing the Gaussian distribution.
        """
    
    x, y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
    x_c, y_c = center
    angle_rad = np.radians(angle_deg)
    
    # Rotate coordinates
    x_rot = (x - x_c) * np.cos(angle_rad) - (y - y_c) * np.sin(angle_rad)
    y_rot = (x - x_c) * np.sin(angle_rad) + (y - y_c) * np.cos(angle_rad)
    
    # Calculate Gaussian values
    gaussian = A *np.exp(-(x_rot ** 2 / (2 * sigma_x ** 2) + y_rot ** 2 / (2 * sigma_y ** 2)))
    
    if norm:
        return gaussian / (2 * np.pi * sigma_x * sigma_y) # normalize the gaussian
    else:
        return gaussian
    

def model_beam_func(peak_flux,shape,x,y,bmaj_p,bmin_p,bpa):
    """
    This function creates a model beam with the same shape as the image cutout.
    """
    model_beam = generate_2D_gaussian(peak_flux,shape,(x,y),bmaj_p/2,bmin_p/2,bpa,norm=False)
    return model_beam

def calculate_beam_correction(mask,peak_flux,shape,x,y,bmaj_p,bmin_p,bpa):
    """
    Calculates the correction factor for the model beam.
    """

    model_beam = model_beam_func(peak_flux,shape,x,y,bmaj_p,bmin_p,bpa)
    model_beam_flux = np.sum(model_beam)
    masked_beam_flux = np.sum(mask*model_beam)
    correction_factor = model_beam_flux/masked_beam_flux
    return correction_factor