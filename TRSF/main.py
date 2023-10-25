'''
file: main.py
author: Rhys Shaw
date: 17-10-2023
description: This file contains the main class for the new source finder.
'''

from TRSF.homology2D import compute_ph_components
import numpy as np
from skimage import measure
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas

## Add checker of path function.
## Create a similar style known flux density image and compare the corrected flux densities to know values.
## plot the difference as a funciton of correction factor.


class sf:

    


    def __init__(self,image,image_PATH,mode,pb_PATH=None,cutup=False,cutup_size=500,output=True):
        
        print('Initialising Source Finder..')
        self.cutup = cutup
        self.output = output
        self.image_PATH = image_PATH
        if image_PATH is None:
            self.image = image
        else:
            self.image,self.header = self.open_image()
        
        self.mode = mode
        self.pb_PATH = pb_PATH

        if self.pb_PATH is not None:
            self.pb_image = self.open_pb()

        if self.cutup:
            # cut the image into sizes of 500x500 pixels.
            # will need to correct the x1,y1,x2,y2 values in the catalogue. accordingly for the cutup. 
            # when calculating source_characteristics
            self.cutouts, self.coords = self.cut_image(cutup_size,self.image)
        
            # if we have a pb image, we need to cut that up too.
            if self.pb_PATH is not None:
                self.pb_cutouts, self.pb_coords = self.cut_image(cutup_size,self.pb_image)
                
        print('Done.')    
            
            
    def cut_image(self,size,image):
        # return a list of cutouts of the image. and coords of the cutout

        cutouts = []
        coords = []
        for i in range(0,image.shape[0],size):
            for j in range(0,image.shape[1],size):
                # Create a cutout
                cutout = image[i:i+size,j:j+size]
                # Append the cutout to the list
                cutouts.append(cutout)
                # Append the coordinates to the list
                coords.append((i,j))
        
        return cutouts, coords


    def phsf(self,lifetime_limit=0):

        if self.cutup == True:
            
            # loop over all cutouts
            catalogue_list = []
            for i, cutout in enumerate(self.cutouts):
                
                print('Computing for Cutout number : {}/{}'.format(i,len(self.cutouts)))
                catalogue = compute_ph_components(cutout,self.local_bg[i],lifetime_limit=lifetime_limit,output=self.output)
                # add cutout coords to catalogue
                catalogue['Y0_cutout'] = self.coords[i][0]
                catalogue['X0_cutout'] = self.coords[i][1]
                catalogue_list.append(catalogue)
            
            # combine the catalogues
            self.catalogue = pandas.concat(catalogue_list)
                
        else:   
        
            self.catalogue = compute_ph_components(self.image,self.local_bg,lifetime_limit=lifetime_limit,output=self.output)




    def source_characterising(self):

        if self.mode == 'Radio':
            # need to read beam size from fits header
        
            
            
            if self.cutup:
                for i, cutout in enumerate(self.cutouts):
                    
                    cutout_cat = self.catalogue[(self.catalogue['Y0_cutout'] == self.coords[i][0]) & (self.catalogue['X0_cutout'] == self.coords[i][1])]
                    
                    Cutout_catalogue = self.radio_characteristing(catalogue=cutout_cat,cutout=cutout,cutout_pb=self.pb_cutouts[i])  
                    # add cutout coords to catalogue
                    Cutout_catalogue['Y0_cutout'] = self.coords[i][0] # these get removed in the previous function.
                    Cutout_catalogue['X0_cutout'] = self.coords[i][1]
                    if i == 0:
                        Processed_catalogue = Cutout_catalogue
                    else:
                        Processed_catalogue = pandas.concat([Processed_catalogue,Cutout_catalogue])
                self.catalogue = Processed_catalogue
                # correct for the poistion of the cutout.
                
                self.catalogue['y1'] = self.catalogue['y1'] + self.catalogue['X0_cutout']
                self.catalogue['y2'] = self.catalogue['y2'] + self.catalogue['X0_cutout']
                self.catalogue['x1'] = self.catalogue['x1'] + self.catalogue['Y0_cutout']
                self.catalogue['x2'] = self.catalogue['x2'] + self.catalogue['Y0_cutout']
                self.catalogue['Xc'] = self.catalogue['Xc'] + self.catalogue['X0_cutout']
                self.catalogue['Yc'] = self.catalogue['Yc'] + self.catalogue['Y0_cutout']
                self.catalogue['bbox1'] = self.catalogue['bbox1'] + self.catalogue['Y0_cutout'] - 1
                self.catalogue['bbox2'] = self.catalogue['bbox2'] + self.catalogue['X0_cutout'] - 1
                self.catalogue['bbox3'] = self.catalogue['bbox3'] + self.catalogue['Y0_cutout']
                self.catalogue['bbox4'] = self.catalogue['bbox4'] + self.catalogue['X0_cutout']
                
                # add cutout coords to catalogue
            
            else:
                self.radio_characteristing()
                



    def open_pb(self):
        
        hdul = fits.open(self.pb_PATH)
        image = hdul[0].data

        image_shape = image.shape

        if len(image_shape) == 4:
        
            image = np.squeeze(image, axis=(0,1))
        
        if len(image_shape) == 3:
        
            image = np.squeeze(image, axis=(0))    

        return image



    def radio_characteristing(self,catalogue=None,cutout=None,cutout_pb=None):

        # get beam in pixels
        
        if cutout is None:
            self.Beam = self.calculate_beam()
            catalogue = self.catalogue
            image = self.image
            pb_image = self.pb_image
        
            
        else:
            self.Beam = self.calculate_beam()
            image = cutout
            pb_image = cutout_pb
        
        # for each source in the catalogue create mask and measure properties. prephorm source flux correction.
        
        self.flux_correction_list = []
        
        params = []

        for i, source in tqdm(catalogue.iterrows(),total=len(catalogue),desc='Calculating Source Properties..',disable=not self.output):
            
            mask = self.get_mask(row=source,image=image)
            #print(np.sum(mask))
            #plt.imshow(mask)
            #plt.show()
        
            source_props = self.get_region_props(mask,image=image)
            
            source_props = self.props_to_dict(source_props[0])

            peak_coords = np.where(image == source_props['max_intensity'])

            y_peak_loc = peak_coords[0][0]
            x_peak_loc = peak_coords[1][0]

            Model_Beam = self.model_beam_func(source_props['max_intensity'],image.shape,x_peak_loc,y_peak_loc,self.BMAJp,self.BMINp,self.BPA)
            Flux_correction_factor = self._flux_correction_factor(mask, Model_Beam)

            self.flux_correction_list.append(Flux_correction_factor)

            # calculate the flux of the source with option for pb correction.
            if self.pb_PATH is not None:
                
                Flux_total = np.nansum(mask*image/pb_image)/self.Beam   # may need to be altered for universality.
                Flux_peak = np.nanmax(mask*image/pb_image)              # may need to be altered for universality.
                
            else:

                Flux_total = np.nansum(mask*image)/self.Beam
                Flux_peak = np.nanmax(mask*image)
            
            Flux_total = Flux_total*Flux_correction_factor
            Area = np.sum(mask)

            Xc = source_props['centroid'][1]
            Yc = source_props['centroid'][0]

            bbox1 = source_props['bbox'][0]
            bbox2 = source_props['bbox'][1]
            bbox3 = source_props['bbox'][2]
            bbox4 = source_props['bbox'][3]

            Maj = source_props['major_axis_length']
            Min = source_props['minor_axis_length']
            Pa = source_props['orientation']

            params.append([source.name,
                           source.Birth,
                           source.Death,
                           source.x1,
                           source.y1,
                           source.x2,
                           source.y2,
                           Flux_total,
                           Flux_peak,
                           Flux_correction_factor,
                           Area,
                           Xc,
                           Yc,
                           bbox1,
                           bbox2,
                           bbox3,
                           bbox4,
                           Maj,
                           Min,
                           Pa,
                           source.parent_tag,
                           source.Class])
        if self.cutup:   
            cutup_Catalogue = self.create_params_df(params)
            return cutup_Catalogue
        else:
            self.create_params_df(params)
        

    




    def create_params_df(self,params):
       
        '''
        
        Creates a pandas dataframe from the parameters.
        
        '''
        
        params = pandas.DataFrame(params,columns=['ID',
                                                  'Birth',
                                                  'Death',
                                                  'x1',
                                                  'y1',
                                                  'x2',
                                                  'y2',
                                                  'Flux_total',
                                                  'Flux_peak',
                                                  'Corr_f',
                                                  'Area',
                                                  'Xc',
                                                  'Yc',
                                                  'bbox1',
                                                  'bbox2',
                                                  'bbox3',
                                                  'bbox4',
                                                  'Maj',
                                                  'Min',
                                                  'Pa',
                                                  'parent_tag',
                                                  'Class'])
       
        if self.cutup:
            return params
        
        else:
            self.catalogue = params



    def get_region_props(self,mask,image):
        region = measure.regionprops(mask,image)
        return region
    
    def save_polygons_to_ds9(self, filename):
        
        '''
        Saves the polygons to a ds9 region file.
        '''

        with open(filename, 'w') as f:
            f.write('# Region file format: DS9 version 4.1\n')
            f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
            for polygon in self.polygons:
                f.write('polygon(')
                for i, point in enumerate(polygon):
                    f.write('{:.2f},{:.2f}'.format(point[1]+1, point[0]+1)) # note this transformation as the index in some CARTA inmages start at -1.
                    if i < len(polygon) - 1:
                        f.write(',')
                f.write(')\n')



    def get_mask(self,row,image):
        
        mask = np.zeros(image.shape)
        mask = np.logical_or(mask,np.logical_and(image <= row.Birth, image > row.Death))
        mask_enclosed = self.get_enclosing_mask(int(row.y1),int(row.x1),mask)
        # set mask as integer
        mask_enclosed = mask_enclosed.astype(int)
        return mask_enclosed




    def open_image(self):

        hdul = fits.open(self.image_PATH)
        image = hdul[0].data
        header = hdul[0].header

        image_shape = image.shape

        if len(image_shape) == 4:
        
            image = np.squeeze(image, axis=(0,1))
        
        if len(image_shape) == 3:
        
            image = np.squeeze(image, axis=(0))    

        return image, header
    




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
    



    def plot_sources(self,vmin=0.12,vmax=1,figsize=(10,10)):

        plt.figure(figsize=figsize)
        plt.imshow(self.image,cmap='gray')
        for i, poly in enumerate(self.polygons):
            plt.plot(poly[:,1],poly[:,0])
        plt.show()





    def _get_polygons(self,x1,y1,birth,death):
        
        '''
        Returns the polygon of the enclosed area of the point (x,y) in the mask.
        '''

        mask = np.zeros(self.image.shape)
        mask = np.logical_or(mask,np.logical_and(self.image <= birth,self.image > death))
        
        mask = self.get_enclosing_mask(int(y1),int(x1),mask)
        contour = measure.find_contours(mask,0)[0]

        return contour



    def create_polygons(self):
        
        polygons = []
        for index, row in tqdm(self.catalogue.iterrows(),total=len(self.catalogue),desc='Creating polygons'):
            try:
                contour = self._get_polygons(row.x1,row.y1,row.Birth,row.Death)
                polygons.append(contour)
            except:
                continue
        self.polygons = polygons
        
        
    def _get_polygons_in_bbox(self,Xmin,Xmax,Ymin,Ymax,x1,y1,birth,death):
        image_crop = self.image[Ymin:Ymax,Xmin:Xmax]
    
        mask = np.zeros_like(image_crop)
        mask = np.logical_or(mask,np.logical_and(image_crop <= birth,image_crop > death))
        mask = self.get_enclosing_mask(int(y1-Xmin),int(x1-Ymin),mask)
        contour = measure.find_contours(mask, 0)[0]
        
        # correct the coordinates to the original image
        contour[:,0] += Ymin
        contour[:,1] += Xmin
        # plot the contour on the original image

        return contour  
  
        


    def create_polygons_fast(self):
        # since we have a bounding box, we can just create a polygon in the bounding box.
        polygons = []
        for index, row in tqdm(self.catalogue.iterrows(),total=len(self.catalogue),desc='Creating polygons'):
            contour = self._get_polygons_in_bbox(row.bbox2-10,row.bbox4+10,row.bbox1-10,row.bbox3+10,row.x1,row.y1,row.Birth,row.Death)
            polygons.append(contour)
        
        self.polygons = polygons
            

    def set_background(self,detection_threshold,set_bg=None):
        
        self.sigma = detection_threshold

        if self.mode == 'Radio':
            
            if self.cutup:
                
                # loop though each cutout and calculate the local background.
                
                local_bg_list = []
                for cutout in self.cutouts:
                    local_bg = self.radio_background(cutout)
                    local_bg_list.append(local_bg*self.sigma)
                local_bg = local_bg_list
            else:
                
                # Radio background is calculated using the median absolute deviation of the total image.

                local_bg = self.radio_background(self.image)
                local_bg = local_bg*self.sigma
            
            
        if self.mode == 'Optical':
            # Optical background is calculated using a random sample of pixels
            mean_bg, std_bg = self.optical_background(nsamples=1000)
            local_bg = mean_bg + self.sigma*std_bg
            
        if self.mode == 'X-ray':
            # no implemented X-ray specific background function.
            mean_bg,std_bg = self.xray_background()
            local_bg = mean_bg + self.sigma*std_bg

        if self.mode == 'other':
            # If the user has a custom background function, they can pass it in here.
            local_bg = set_bg
        
        self.local_bg = local_bg
        if self.cutup:
            print('Mean Background across cutouts: ', np.nanmean(self.local_bg))
        else:
            print('Background set to: ',self.local_bg)




    def radio_background(self,image):
        
        from astropy.stats import mad_std

        local_bg = mad_std(image,ignore_nan=True)
        
        return local_bg
    



    def optical_background(self,nsamples):
        
        '''
        
        Returns the mean and standard deviation of a random sample of pixels in the image.
        This function also resamples the mean by doing a 3-sigma clip.

        '''
        
        sample_box = np.array([[1,1,1],
                        [1,1,1],
                        [1,1,1]])
        image_shape = self.image.shape
        min_x = 2
        max_x = image_shape[0]-2
        min_y = 2
        max_y = image_shape[1]-2

        # list of random x and y coordinates
        x = np.random.randint(min_x, max_x, nsamples)
        y = np.random.randint(min_y, max_y, nsamples)

        background_values = []

        for i in range(len(x)):
            sample_box = np.array([[1,1,1],
                        [1,1,1],
                        [1,1,1]])
            
            sample_box = sample_box*self.image[x[i]-1:x[i]+2,y[i]-1:y[i]+2]
            background_values += sample_box.flatten().tolist()
        

        background_values = np.array(background_values)
        background_values = background_values[background_values < np.mean(background_values) + 3*np.std(background_values)]
        background_values = background_values[background_values > np.mean(background_values) - 3*np.std(background_values)]

        mean_bg = np.mean(background_values)
        std_bg = np.std(background_values)
        return mean_bg, std_bg
    
    



    def xray_background(self):
        # no implemented X-ray specific background function.
        mean_bg, std_bg = self.optical_background(nsamples=1000)
        return mean_bg, std_bg


    def _image_smoothing(self,img,smooth_param):
        
        '''
        
        Applies a Gaussian filter to the image to smooth it as desired.

        '''
        
        # import gaussian filter
        from scipy.ndimage import gaussian_filter
        # smooth image
        img = gaussian_filter(img, sigma=smooth_param) # std of gaussian kernel
        return img


    def calculate_beam(self):
        '''

        Calculates the beam of the image in pixels.

        returns:
            beam sampling factor, and creates object variable BMAJp and BMINp.
        
        '''
        
        # get beam info from header
        
        BMAJ = self.header['BMAJ']
        BMIN = self.header['BMIN']
        # convert to pixels
        arcseconds_per_pixel = self.header['CDELT1']*3600*-1
        beam_size_arcseconds = self.header['BMAJ']*3600
        BMAJ_oversampled_spacial_width = (BMAJ**2 + beam_size_arcseconds**2)**0.5
        BMAJ = BMAJ_oversampled_spacial_width/arcseconds_per_pixel
        BMIN_oversampled_spacial_width = (BMIN**2 + beam_size_arcseconds**2)**0.5
        BMIN = BMIN_oversampled_spacial_width/arcseconds_per_pixel
        self.BMAJp = BMAJ
        self.BMINp = BMIN
        try:
            self.BPA = self.header['BPA']
        except KeyError:
            self.BPA = 0
        return np.pi * (BMAJ)*(BMIN) / (4*np.log(2))