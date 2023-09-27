
import numpy as np
import pandas
from skimage import measure
from TRSF.source_props.expand_region import region_expansion_downhill
from TRSF.source_props.gaussian_fitting import fit_gaussian_2d
import matplotlib.pyplot as plt
from tqdm import tqdm
from astropy.io import fits 
import time
import skimage.measure as measure
from astropy.stats import mad_std

import skimage.measure as measure

class cal_props:

    def __init__(self, pd: pandas.DataFrame, img: np.ndarray, local_bg: float, 
                 sigma: float, pbimg = None,method: str = None,
                 expsigma: float = 3,beam: float = 1,
                 bmajp: float = None, bminp: float = None, bpa: float = None):
        
        self.pd = pd
        self.img = img
        self.local_bg = local_bg
        self.sigma = sigma
        self.method = method
        self.expsigma = expsigma
        self.pbimg = pbimg
        self.beam = beam
        self.bmajp = bmajp
        self.bminp = bminp
        self.bpa = bpa

    def get_enclosing_mask(self,x, y, mask):
        '''
        Returns the mask of the enclosed area of the point (x,y) in the mask.
        '''
        
        # Ensure the point is inside the mask
        if not mask[y, x]:
            return None
        # Create a copy of the mask
        enclosed_mask = np.copy(mask)
        # Perform a flood fill starting from the point
        h, w = mask.shape
        stack = [(x, y)]
        while stack:
            x, y = stack.pop()
            if enclosed_mask[y, x]:
                enclosed_mask[y, x] = False
                if x > 0:
                    stack.append((x - 1, y))
                if x < w - 1:
                    stack.append((x + 1, y))
                if y > 0:
                    stack.append((x, y - 1))
                if y < h - 1:
                    stack.append((x, y + 1))
        
        return mask & (~enclosed_mask)

    def create_source_mask(self, idx, pd):
        '''
        
        '''
        
        point = pd.iloc[idx]
        mask = np.zeros(self.img.shape)
        mask = np.logical_or(mask,np.logical_and(self.img <= point.Birth,self.img > point.Death))
        mask = self.get_enclosing_mask(int(point.y1),int(point.x1),mask)
        mask = mask.astype(int)
        return mask

    def fit_all_single_sources(self,gaussian_fit: bool = True,expand: bool = True):

        print('Calculating properties...')
        print(self.pd)
        
        # seperate pd into classes 
        pdclass0 = self.pd[self.pd['single'] == 0] # single sources
        pdclass1 = self.pd[self.pd['single'] == 1] # child sources
        pdclass2 = self.pd[self.pd['single'] == 2] # parent sources
        pdclass3 = self.pd[self.pd['single'] == 3] # child sources of parent sources

        # process each class seperately

        # class 0

        for i in range(0,len(pdclass0)):

            mask = self.create_source_mask(i,pdclass0)
            
            region_props = measure.regionprops(mask,self.img)[0]
            plt.figure(figsize=(10,10))
            plt.imshow(mask*self.img)
            plt.scatter(region_props['weighted_centroid'][1],region_props['weighted_centroid'][0],c='r')
            plt.show()

            plt.imshow(self.img)
            plt.show()
            # calculate the flux
            flux = np.sum(self.img*mask)
            x_peak_loc, y_peak_loc = region_props['weighted_centroid'][1], region_props['weighted_centroid'][0]
            #y_peak_loc, x_peak_loc = peak_coords[0][0], peak_coords[1][0]
            print(self.bmajp,self.bminp,self.bpa)
            Model_Beam = self.model_beam_func(region_props['max_intensity'],self.img.shape,x_peak_loc,y_peak_loc,self.bmajp/2,self.bminp/2,self.bpa)
            print(flux)
            print(Model_Beam)
            plt.figure(figsize=(10,10))
            plt.imshow(Model_Beam*mask)
            plt.show()
            
            print('beam*mask sum: ',np.sum(Model_Beam*mask))
            print('sum beam: ',np.sum(Model_Beam))
            print('correction factor: ', np.sum(Model_Beam)/np.sum(Model_Beam*mask))
            print('mask*img sum:',np.sum(mask*self.img))
            print('Beam', self.beam)
            print('Flux_density: ',np.sum(mask*self.img)/self.beam)


    def generate_2d_gaussian(self,A,shape, center, sigma_x, sigma_y, angle_deg=0,norm=True):
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



    def model_beam_func(self,peak_flux,shape,x,y,bmaj,bmin,bpa):
        model_beam = np.zeros(shape)
        model_beam = self.generate_2d_gaussian(peak_flux,shape,(x,y),bmaj,bmin,bpa,norm=False)
        return model_beam