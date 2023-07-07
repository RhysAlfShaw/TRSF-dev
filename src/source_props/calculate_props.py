'''
package: TRSF-SourceProperties
author: Rhys Shaw
email: rhys.shaw@bristol.ac.uk
date: 30-06-2023
description: This file contains the functions for turning a processed persistence 
    diagram into a set of sources with properties.
'''

from typing import Any
import numpy as np
import pandas
from skimage import measure
from src.source_props.expand_region import region_expansion_downhill
from src.source_props.gaussian_fitting import fit_gaussian_2d
import matplotlib.pyplot as plt
from tqdm import tqdm
import time

class cal_props:

    def __init__(self, pd: pandas.DataFrame, img: np.ndarray, local_bg: float, sigma: float, method: str = None):
        self.pd = pd
        self.img = img
        self.local_bg = local_bg
        self.sigma = sigma
        self.method = method

    def calculate_bounding_box_of_mask(self,mask):
        '''
        Calculates the bounding box of a mask.
        Returns the bounding box.
        '''
        regionprops = measure.regionprops(mask.astype(int))
        return regionprops[0].bbox


    def fit_all_single_sources(self,gaussian_fit: bool = True,expand: bool = True):
        '''
        Fits all the single sources in the persistence diagram.
        Returns a pandas dataframe with the properties of the sources.
        '''
        # prehaps a parrellization option here. 
        #ss_pd = self.pd 
        ss_pd =self.pd[self.pd['single']!=2] # e do not what fit to extended sources.
        params = []
        for i in tqdm(range(len(ss_pd)),total=len(ss_pd),desc='Fitting single sources'):
            # create mask
            
            mask = self.create_source_mask(i,ss_pd)
            name = ss_pd.iloc[i].name
            # get region props
            regionprops = self.get_region_props(mask)
            regionprops = self.props_to_dict(regionprops[0])
        
            # expand mask to improve fitting.
            # check if mask is empty
            #
            # this takes a long time. around 0.75s for each source in a 500x500 image.
            if gaussian_fit == True:
                if expand == True:
                    mask = self.expand_mask_downhill(mask,max_iter=3)
                bbox = self.calculate_bounding_box_of_mask(mask)
                #print(mask.sum())
                
                if mask.sum() == 0:
                    print('Empty mask. {}'.format(ss_pd.iloc[i].name))
                    amp, x0, y0, sigma_x, sigma_y, theta = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
                    params.append([amp, x0, y0, sigma_x, sigma_y, theta])
                    continue
                # crop mask based on regionprops bounding box
                # bouding box need to be updated to account for the expansion.

                temp_img = mask*self.img
                temp_img = temp_img[bbox[0]:bbox[2],bbox[1]:bbox[3]]


                try:

                    amp, x0, y0, sigma_x, sigma_y, theta = self.gaussian_fit(temp_img,regionprops)
                
                except (RuntimeError, TypeError): # we will assign these types of failures to nan.
                    print('WARNING: Fitting Failure - Failed to fit source. ID {}'.format(ss_pd.iloc[i].name))
                    amp, x0, y0, sigma_x, sigma_y, theta = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
            else:
                #plt.imshow(mask)
                #plt.show()
                if expand == True:
                    mask = self.expand_mask_downhill(mask,max_iter=3)
                bbox = self.calculate_bounding_box_of_mask(mask)
                #plt.imshow(mask)
                #plt.plot([bbox[1],bbox[3],bbox[3],bbox[1],bbox[1]],[bbox[0],bbox[0],bbox[2],bbox[2],bbox[0]])
                #plt.show()
                amp, x0, y0, sigma_x, sigma_y, theta = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
            # correct x0 and y0 for the bounding box.
            #if self.fit_param_check(bbox,[amp, x0, y0, sigma_x, sigma_y, theta]):
            x0 = x0 + bbox[1]
            y0 = y0 + bbox[0]
            peak_flux = regionprops['max_intensity']
            x_c = regionprops['centroid'][0]
            y_c = regionprops['centroid'][1]
            bbox = bbox#regionprops['bbox']
            # add the pd params to the list.
            
            params.append([name, amp, x0, y0, sigma_x, sigma_y, theta, peak_flux, x_c, y_c, bbox,ss_pd.iloc[i]['single']])
            # deconstruct regionprops.
        #arams_df = self.create_params_df(params)
        return params
    
    def create_params_df(self,params):
        '''
        Creates a pandas dataframe from the parameters.
        '''
        columns = ['index','amp','x','y','sigma_x','sigma_y','theta','peak_flux','x_c','y_c','bbox','class']
        params = pandas.DataFrame(params,columns)
        return params
    
    def format_gaussian_fit(self,params):
        pass
    
    def fit_param_check(self,bbox,params):
        '''
        Checks the fitting parameters for any issues.
        '''
        # check that the x0 and y0 are within the bounding box.
        # check that the sigma_x and sigma_y are not more than 1x the size of the bounding box.
        amp, x0, y0, sigma_x, sigma_y, theta = params
        if x0 < bbox[1] or x0 > bbox[3]:
            #print('x0 out of bounds. {}'.format(x0))
            return False
        if y0 < bbox[0] or y0 > bbox[2]:
            #print('y0 out of bounds. {}'.format(y0))
            return False
        
        if sigma_x > (bbox[3]-bbox[1]):
            #print('sigma_x too big. {}'.format(sigma_x))
            return False
        if sigma_y > (bbox[2]-bbox[0]):
            #print('sigma_y too big. {}'.format(sigma_y))
            return False
        return True
        



    def gaussian_fit(self,temp_img,regionprops):
        
        sigma_x = regionprops['major_axis_length']/2
        sigma_y = regionprops['minor_axis_length']/2
        theta = regionprops['orientation']
        x0 = regionprops['centroid'][0]
        y0 = regionprops['centroid'][1]
        amp = regionprops['max_intensity']
        guess = [amp,x0,y0,sigma_x,sigma_y,theta]
        #print(guess)
        params, fitted_gauss = fit_gaussian_2d(temp_img,maxfev=12000,inital_guess=guess)
        amp, x0, y0, sigma_x, sigma_y, theta = params
        return amp, x0, y0, sigma_x, sigma_y, theta

    def get_region_props(self,mask):
        region = measure.regionprops(mask,self.img)
        return region
    
    def expand_mask_downhill(self,mask,max_iter=3):
        mask = region_expansion_downhill(mask,self.img,self.local_bg*self.sigma,method=self.method,max_iter=max_iter)
        return mask
    
    def props_to_dict(self,regionprops):
        #print(regionprops)
        # return mask here if wanted.
        dict = {
            'area': regionprops.area,
            'bbox': regionprops.bbox,
            'bbox_area': regionprops.bbox_area,
            'centroid': regionprops.centroid,
            'convex_area': regionprops.convex_area,
            'eccentricity': regionprops.eccentricity,
            'equivalent_diameter': regionprops.equivalent_diameter,
            'euler_number': regionprops.euler_number,
            'extent': regionprops.extent,
            'filled_area': regionprops.filled_area,
            'major_axis_length': regionprops.major_axis_length,
            'minor_axis_length': regionprops.minor_axis_length,
            'moments': regionprops.moments,
            'perimeter': regionprops.perimeter,
            'solidity': regionprops.solidity,
            'orientation': regionprops.orientation,
            'max_intensity':regionprops.max_intensity,
        }
        return dict
        
    def get_source(self,num):
        return self.pd.iloc[num]
    
    def create_source_mask(self, idx, pd):
        point = pd.iloc[idx]
        mask = np.zeros(self.img.shape)
        mask = np.logical_or(mask,np.logical_and(self.img <= point.Birth,self.img > point.Death))
        mask = self.get_enclosing_mask(int(point.y1),int(point.x1),mask)
        mask = mask.astype(int)
        return mask

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