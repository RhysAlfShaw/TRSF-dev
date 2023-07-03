'''
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

class cal_props:

    def __init__(self, pd: pandas.DataFrame, img: np.ndarray, local_bg: float, sigma: float):
        self.pd = pd
        self.img = img
        self.local_bg = local_bg
        self.sigma = sigma

    def get_region_props(self,mask):
        region = measure.regionprops(mask,self.img)
        return region
    
    def expand_mask_downhill(self,mask,max_iter=3):
        mask = region_expansion_downhill(mask,self.img,self.local_bg,max_iter=max_iter)
        return mask
    
    def props_to_dict(self,regionprops):
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
    
    def create_source_mask(self, idx):
        point = self.get_source(idx)
        mask = np.zeros(self.img.shape)
        mask = np.logical_or(mask,np.logical_and(self.img <= point.Birth,self.img > point.Death))
        mask = self.get_enclosing_mask(int(point.y1),int(point.x1),mask)
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