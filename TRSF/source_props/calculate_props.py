'''
package: TRSF-SourceProperties
author: Rhys Shaw
email: rhys.shaw@bristol.ac.uk
date: 30-06-2023
description: This file contains the functions for turning a processed persistence 
    diagram into a set of sources with properties.
'''
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

class cal_props:
    
    '''

    Calculates source properties for each given image.

    Args:
        - Persistent Diagram (pandas.DataFrame).
        - Image (numpy array).
        - Local bg (float). 
    
    Returns:
        - 


    
    '''
    
    
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

    def calculate_bounding_box_of_mask(self,mask):
        '''
        Calculates the bounding box of a mask.
        Returns the bounding box.
        '''
        regionprops = measure.regionprops(mask.astype(int))
        return regionprops[0].bbox

    def _open_img_pbcorr(self,img_path):
        hdu = fits.open(img_path)
        img = hdu[0].data
        # measure shape 
        img = np.squeeze(img,axis=(0,1))
        
        print('pbcorr Image shape: {}'.format(img.shape))
        img = self._crop_image(img)
        print('pbcorr Image shape: {}'.format(img.shape))
        img[np.isnan(img)] = 0
        if img.shape[0] != img.shape[1]:
            img = self._make_square(img)
        return img
    
    def _make_square(self,img):
        # make image square by adding padding
        if img.shape[0] > img.shape[1]:
            diff = img.shape[0] - img.shape[1]
            pad = np.zeros((img.shape[0],diff))
            img = np.concatenate((img,pad),axis=1)
        elif img.shape[1] > img.shape[0]:
            diff = img.shape[1] - img.shape[0]
            pad = np.zeros((diff,img.shape[1]))
            img = np.concatenate((img,pad),axis=0)
        return img
    
    def _crop_image(self,arr): 
        arr = np.array(arr)  # Convert input to numpy array
        mask_rows = np.all(np.isnan(arr), axis=1)
        mask_cols = np.all(np.isnan(arr), axis=0)
        return arr[~mask_rows][:, ~mask_cols]

    def fit_all_single_sources(self,gaussian_fit: bool = True,expand: bool = True):
        '''
        Fits all the single sources in the persistence diagram.
        Returns a pandas dataframe with the properties of the sources.
        '''
        # prehaps a parrellization option here. 
        #if pbcorr == True and pbimage==True: # corrected image
        #    pbcorrimg = self._open_img_pbcorr(pbcorr_path)
            
            # open pbcorr img same as original img.
            # if image is tghe primary beam, them create the corrected image.
        #if pbcorr == True and pbimage == False:
        #    pbimg = self._open_img_pbcorr(pbcorr_path)
            # make corrected image
            # are the images the same shape?
        #    if pbimg.shape != self.img.shape:
        #        print(pbimg.shape,self.img.shape)
        #        raise ValueError('Image and primary beam image are not the same shape, make change and try again.')
            
        #    pbcorrimg = self.img/pbimg

        #ss_pd =self.pd[self.pd['single']!=2] # we do not what fit to extended components.   
        ss_pd = self.pd    
        params = []
        Flux_correction_list = []
        Flux_tot_before = []
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
                    mask = self.expand_mask_downhill(mask,max_iter=10)
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
                    #print('WARNING: Fitting Failure - Failed to fit source. ID {}'.format(ss_pd.iloc[i].name))
                    amp, x0, y0, sigma_x, sigma_y, theta = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
                except TypeError:
                    #plt.imshow(temp_img)
                    #raise TypeError('Fitting Failure - Failed to fit source. ID {}'.format(ss_pd.iloc[i].name))
                    #plt.show()
                    #print('WARNING: Fitting Failure - Failed to fit source. ID {}'.format(ss_pd.iloc[i].name))
                    amp, x0, y0, sigma_x, sigma_y, theta = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
            else:
                
                if expand == True:
                    mask = self.expand_mask_downhill(mask,max_iter=3)
                bbox = self.calculate_bounding_box_of_mask(mask)
                
                amp, x0, y0, sigma_x, sigma_y, theta = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
            
            #calculate int flux

            Beam = self.beam
            # coords where the peak is located.
            peak_coords = np.where(self.img == regionprops['max_intensity'])
            #print(peak_coords)
            y_peak_loc, x_peak_loc = peak_coords[0][0], peak_coords[1][0]
            Model_Beam = self.model_beam_func(regionprops['max_intensity'],self.img.shape,x_peak_loc,y_peak_loc,self.bmajp/2,self.bminp/2,self.bpa)
            Flux_correction_factor = self._flux_correction_factor(mask,Model_Beam)
            #print(Flux_correction_factor)
            Flux_correction_list.append(Flux_correction_factor)
            # calculate flux here # should be sum(mask*img)
            if self.pbimg is None:
                Flux_tot_corr = 0
                Flux_peak_corr = 0
            else:
                
                background = mad_std(self.pbimg)
                # if shape of mask and image are not the same, then we crop the pbimg to the mask shape.
                if mask.shape != self.pbimg.shape:
                    pbimg = np.resize(self.pbimg, mask.shape)
                else:
                    pbimg = self.pbimg
                
                Flux_tot_corr = np.nansum(mask*pbimg) / Beam 
#                if correction == True:
                #if ss_pd.iloc[i]['single'] == 0:
                Flux_before = Flux_tot_corr
                Flux_tot_corr = Flux_tot_corr * Flux_correction_factor
#                else:
                #Flux_tot_corr = Flux_tot_corr
                Flux_peak_corr = np.nanmax(mask*pbimg)      #- background
            Flux_tot_before.append(Flux_before)
                #print('Corrected Flux {}'.format(Flux_tot_corr))
            background = mad_std(self.img) 
            Flux_tot = np.sum(mask*self.img) / Beam 
            Flux_tot = Flux_tot * Flux_correction_factor
            Flux_peak = np.max(mask*self.img)  #- background
            Area = np.sum(mask)
            
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
                           ss_pd.iloc[i]['x1'],ss_pd.iloc[i]['y1'],ss_pd.iloc[i]['lifetime'],Flux_tot,Flux_peak,Flux_peak_corr,Flux_tot_corr,Area])
            
        #print(Flux_correction_list)
        return self.create_params_df(params), [Flux_correction_list,Flux_tot_before]
    

    def _flux_correction_factor(self,mask,Model_Beam):
        # calculate the correction factor
        model_beam_flux = np.sum(Model_Beam)
        masked_beam_flux = np.sum(mask*Model_Beam)
        correction_factor = model_beam_flux/masked_beam_flux
        #if correction_factor > 100:
        #    correction_factor = 100
        return correction_factor

    def model_beam_func(self,peak_flux,shape,x,y,bmaj,bmin,bpa):
        model_beam = np.zeros(shape)
        model_beam = self.generate_2d_gaussian(peak_flux,shape,(x,y),bmaj,bmin,bpa,norm=False)
        return model_beam
    
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


    def _calculate_int_flux(self,row):

        x, y = np.meshgrid(np.arange(0, 100, 1), np.arange(0, 100, 1))
        try:
            gaussian = row.amp * np.exp(-((x-50)**2/(2*row.sigma_x**2) + (y-50)**2/(2*row.sigma_y**2)))
            return gaussian.sum()
        except:
            return row.flux_tot_corr

    def _polygon(self,row):
        mask = np.zeros(self.img.shape)
        mask = np.logical_or(mask,np.logical_and(self.img <= row.Birth,self.img > row.Death))
        mask = self.get_enclosing_mask(int(row.y1),int(row.x1),mask)
        mask = mask.astype(int)
        #mask = self.create_source_mask_s(row)
        # create polygon from the mask
        contour = measure.find_contours(mask, 0.5)[0]    
        background = mad_std(self.pbimg)
        # correct x,y of the contour for the image coords.
        Flux_tot = np.sum(mask*self.img) /self.beam  
        Flux_peak = np.max(mask*self.img)
        Area = np.sum(mask)
        if mask.shape != self.pbimg.shape:
            pbimg = np.resize(self.pbimg, mask.shape)
        else:
            pbimg = self.pbimg
        Flux_tot_corr = np.nansum(mask*pbimg) /self.beam   # - background
        Flux_peak_corr = np.nanmax(mask*pbimg)            # - background 
        return contour,Flux_tot,Flux_peak,Area,Flux_tot_corr,Flux_peak_corr


    def create_params_df(self,params):

        '''
        
        Creates a pandas dataframe from the parameters.
        
        '''
        
        params = pandas.DataFrame(params,columns=['index','amp','x','y','sigma_x','sigma_y','theta','peak_flux','x_c','y_c','bbox','Class','Birth','Death','x1','y1','lifetime','flux_tot','flux_peak','flux_peak_corr','flux_tot_corr','area'])
        pd_alt_combine = self.pd[self.pd['single'] == 2]
        #print(pd_alt_combine)
        if len(pd_alt_combine) > 0:
            result = pd_alt_combine.apply(self._polygon,axis=1,result_type='expand')
            pd_alt_combine['polygon'] = result[0]
            pd_alt_combine['flux_tot'] = result[1]
            pd_alt_combine['flux_peak'] = result[2]
            
            pd_alt_combine['area'] = result[3]

            pd_alt_combine['flux_tot_corr'] = result[4]
            pd_alt_combine['flux_peak_corr'] = result[5]

        pd_alt_combine.rename({'single':'Class'})
        try:
            pd_alt_combine = pd_alt_combine.drop(columns=['new_row','parent_tag','dist','enclosed_i'])
        except KeyError:
            #print(pd_alt_combine.columns)
            if 'encloses_i' not in pd_alt_combine.columns:
                pd_alt_combine['encloses_i'] = 0
            pd_alt_combine = pd_alt_combine.drop(columns=['new_row','len_enclosed','encloses_i'])
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

        '''

        Expansion of the regions using the region_expansion_downhill 
        
        '''
        
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
        
        Creates source for each of the source finder.

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