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
from TRSF.source_props.expand_region import region_expansion_downhill
from TRSF.source_props.gaussian_fitting import fit_gaussian_2d
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
import skimage.measure as measure

class cal_props:

    def __init__(self, pd: pandas.DataFrame, img: np.ndarray, local_bg: float, sigma: float, method: str = None,expsigma: float = 3):
        self.pd = pd
        self.img = img
        self.local_bg = local_bg
        self.sigma = sigma
        self.method = method
        self.expsigma = expsigma

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
        ss_pd =self.pd[self.pd['single']!=2] # we do not what fit to extended components.
        params = []
        for i in tqdm(range(len(ss_pd)),total=len(ss_pd),desc='Fitting single sources',leave=False):
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
                    #print('Expanding mask')
                    mask = self.expand_mask_downhill(mask,max_iter=5)
                bbox = self.calculate_bounding_box_of_mask(mask)
                
                
                if mask.sum() == 0:
                    print('Empty mask. {}'.format(ss_pd.iloc[i].name))
                    amp, x0, y0, sigma_x, sigma_y, theta = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
                    params.append([amp, x0, y0, sigma_x, sigma_y, theta])
                    continue
                

                temp_img = mask*self.img
                temp_img = temp_img[bbox[0]:bbox[2],bbox[1]:bbox[3]]


                try:

                    amp, x0, y0, sigma_x, sigma_y, theta = self.gaussian_fit(temp_img,regionprops)
                    # apply correction to x0 and y0 and check for boundaries, or extreamly elliptical fits.
                except RuntimeError: # we will assign these types of failures to nan.
                    # print the error this is for debugging.
                    #plt.show()
                    #plt.imshow(temp_img)
                    #raise RuntimeError('Fitting Failure - Failed to fit source. ID {}'.format(ss_pd.iloc[i].name))
                    print('WARNING: Fitting Failure - Failed to fit source. ID {}'.format(ss_pd.iloc[i].name))
                    amp, x0, y0, sigma_x, sigma_y, theta = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
                except TypeError:
                    #plt.imshow(temp_img)
                    #raise TypeError('Fitting Failure - Failed to fit source. ID {}'.format(ss_pd.iloc[i].name))
                    #plt.show()
                    print('WARNING: Fitting Failure - Failed to fit source. ID {}'.format(ss_pd.iloc[i].name))
                    amp, x0, y0, sigma_x, sigma_y, theta = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
            else:
                
                if expand == True:
                    mask = self.expand_mask_downhill(mask,max_iter=3)
                bbox = self.calculate_bounding_box_of_mask(mask)
                
                amp, x0, y0, sigma_x, sigma_y, theta = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
            # correct x0 and y0 for the bounding box.
            #if self.fit_param_check(bbox,[amp, x0, y0, sigma_x, sigma_y, theta]):
            x0 = x0 + bbox[1]
            y0 = y0 + bbox[0]
            peak_flux = regionprops['max_intensity']
            x_c = regionprops['centroid'][0]
            y_c = regionprops['centroid'][1]
            bbox = bbox 
            params.append([name, amp, x0, y0, sigma_x, sigma_y, theta, peak_flux, x_c, y_c, 
                           bbox,ss_pd.iloc[i]['single'],ss_pd.iloc[i]['Birth'],ss_pd.iloc[i]['Death'],
                           ss_pd.iloc[i]['x1'],ss_pd.iloc[i]['y1'],ss_pd.iloc[i]['lifetime']])
            
        
        return self.create_params_df(params)
    


    def _polygon(self,row):
        mask = np.zeros(self.img.shape)
        mask = np.logical_or(mask,np.logical_and(self.img <= row.Birth,self.img > row.Death))
        mask = self.get_enclosing_mask(int(row.y1),int(row.x1),mask)
        mask = mask.astype(int)
        #mask = self.create_source_mask_s(row)
        # create polygon from the mask
        contour = measure.find_contours(mask, 0.5)[0]    
        # correct x,y of the contour for the image coords.
        return contour


    def create_params_df(self,params):
        '''
        Creates a pandas dataframe from the parameters.
        '''
        
        params = pandas.DataFrame(params,columns=['index','amp','x','y','sigma_x','sigma_y','theta','peak_flux','x_c','y_c','bbox','Class','Birth','Death','x1','y1','lifetime'])
        pd_alt_combine = self.pd[self.pd['single'] == 2]
        #print(pd_alt_combine)
        if len(pd_alt_combine) > 0:
            pd_alt_combine['polygon'] = pd_alt_combine.apply(self._polygon,axis=1)
        
        pd_alt_combine.rename({'single':'Class'})
        pd_alt_combine = pd_alt_combine.drop(columns=['new_row','parent_tag','len_enclosed'])
        #combined_dataframe = pandas.DataFrame(params,columns=['index','amp','x','y','sigma_x','sigma_y','theta','peak_flux','x_c','y_c','bbox','Class','single','Birth','Death','x1','y1'])
        combined_dataframe = pandas.concat([params, pd_alt_combine])
        return combined_dataframe
    
    def format_gaussian_fit(self,params):
        pass
    
    def fit_param_check(self,bbox,params):
        '''
        Checks the fitting parameters for any issues.
        '''

        amp, x0, y0, sigma_x, sigma_y, theta = params
        if x0 < bbox[1] or x0 > bbox[3]:
            return False
        if y0 < bbox[0] or y0 > bbox[2]:
            return False    
        if sigma_x > (bbox[3]-bbox[1]):
            return False
        if sigma_y > (bbox[2]-bbox[0]):
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
        mask = region_expansion_downhill(mask,self.img,self.local_bg*self.expsigma,method=self.method,max_iter=max_iter)
        return mask
    
    def props_to_dict(self,regionprops):
        '''
        Converts the regionprops to a dictionary.
        '''
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
    
    def create_source_mask_s(self, row):
        point = row
        mask = np.zeros(self.img.shape)
        mask = np.logical_or(mask,np.logical_and(self.img <= row.Birth,self.img > row.Death))
        mask = self.get_enclosing_mask(int(row.y1),int(row.x1),mask)
        mask = mask.astype(int)
        return mask
    
    def create_source_mask(self, idx, pd):
        '''
        
        '''
        
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