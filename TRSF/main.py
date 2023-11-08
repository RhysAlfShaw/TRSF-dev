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
from astropy.wcs import WCS
import pandas
## Add checker of path function.
## Create a similar style known flux density image and compare the corrected flux densities to know values.
## plot the difference as a funciton of correction factor.


class sf:

    
    def __init__(self,image,image_PATH,mode,pb_PATH=None,cutup=False,cutup_size=500,output=True,area_limit=1,smooth_sigma=0):
       
        print("""   
##############################################
_______   _______  _______  _______  _______  
(  __  \ (  ____ )|\     /|\__   __/(  __  \ 
| (  \  )| (    )|| )   ( |   ) (   | (  \  )
| |   ) || (____)|| |   | |   | |   | |   ) |
| |   | ||     __)| |   | |   | |   | |   | |
| |   ) || (\ (   | |   | |   | |   | |   ) |
| (__/  )| ) \ \__| (___) |___) (___| (__/  )
(______/ |/   \__/(_______)\_______/(______/ 
        
#############################################
Detector of astRonomical soUrces in optIcal and raDio images

Version: 0.1.0
For more information see:
https://github.com/RhysAlfShaw/TRSF-dev
        """)
        print('Initialising Source Finder..')
        self.cutup = cutup
        self.output = output
        self.image_PATH = image_PATH
        self.area_limit = area_limit
        self.smooth_sigma = smooth_sigma
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
            else:
                self.pb_cutouts = None
                self.pb_coords = None
        if self.smooth_sigma != 0:
            self.image = self._image_smoothing(self.image,self.smooth_sigma)
            print('Image smoothed with sigma: ',self.smooth_sigma)
      
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
                catalogue = compute_ph_components(cutout,self.local_bg[i],analysis_threshold_val=self.analysis_threshold_val[i],lifetime_limit=lifetime_limit,output=self.output,bg_map=self.bg_map,area_limit=self.area_limit)
                print('Done ')
                # add cutout coords to catalogue
                catalogue['Y0_cutout'] = self.coords[i][0]
                catalogue['X0_cutout'] = self.coords[i][1]
                catalogue_list.append(catalogue)

            # combine the catalogues
            self.catalogue = pandas.concat(catalogue_list)

        else:
            
            print(self.analysis_threshold_val)
            
            self.catalogue = compute_ph_components(self.image,self.local_bg,analysis_threshold_val=self.analysis_threshold_val,lifetime_limit=lifetime_limit,output=self.output,bg_map=self.bg_map,area_limit=self.area_limit)


    def gaussian2dkernal_convolution(self, image, sigma):
        '''
        This function takes an image and a sigma value and convolves the image with a 2D gaussian kernel.
        '''
        from scipy.ndimage import gaussian_filter
        print(sum(image.flatten()))
        image = gaussian_filter(image, sigma=sigma)
        print(sum(image.flatten()))
        return image



    def source_characterising(self):

        if self.mode == 'Radio':
            # need to read beam size from fits header



            if self.cutup:
                for i, cutout in enumerate(self.cutouts):

                    cutout_cat = self.catalogue[(self.catalogue['Y0_cutout'] == self.coords[i][0]) & (self.catalogue['X0_cutout'] == self.coords[i][1])]
                    
                    if self.pb_PATH is not None:
                        Cutout_catalogue = self.radio_characteristing(catalogue=cutout_cat,cutout=cutout,cutout_pb=self.pb_cutouts[i],background_map=self.local_bg[i])
                    else:
                        Cutout_catalogue = self.radio_characteristing(catalogue=cutout_cat,cutout=cutout,cutout_pb=None,background_map=self.local_bg[i])
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
        
        if self.mode == 'optical':
            self.optical_characteristing()
            
            
        if self.header:
            try:
                Ra, Dec = self._xy_to_RaDec(self.catalogue['Xc'],self.catalogue['Yc'])
                self.catalogue['RA'] = Ra
                self.catalogue['DEC'] = Dec
                self.catalogue['RA'] = self.catalogue['RA'].astype(float)
                self.catalogue['DEC'] = self.catalogue['DEC'].astype(float)
            except:
                pass
                
        self.set_types_of_dataframe()

    def set_types_of_dataframe(self):

        self.catalogue['ID'] = self.catalogue['ID'].astype(int)
        self.catalogue['Birth'] = self.catalogue['Birth'].astype(float)
        self.catalogue['Death'] = self.catalogue['Death'].astype(float)
        self.catalogue['x1'] = self.catalogue['x1'].astype(float)
        self.catalogue['y1'] = self.catalogue['y1'].astype(float)
        self.catalogue['x2'] = self.catalogue['x2'].astype(float)
        self.catalogue['y2'] = self.catalogue['y2'].astype(float)
        self.catalogue['Flux_total'] = self.catalogue['Flux_total'].astype(float)
        self.catalogue['Flux_peak'] = self.catalogue['Flux_peak'].astype(float)
        self.catalogue['Corr_f'] = self.catalogue['Corr_f'].astype(float)
        self.catalogue['Area'] = self.catalogue['Area'].astype(float)
        self.catalogue['Xc'] = self.catalogue['Xc'].astype(float)
        self.catalogue['Yc'] = self.catalogue['Yc'].astype(float)
        self.catalogue['bbox1'] = self.catalogue['bbox1'].astype(float)
        self.catalogue['bbox2'] = self.catalogue['bbox2'].astype(float)
        self.catalogue['bbox3'] = self.catalogue['bbox3'].astype(float)
        self.catalogue['bbox4'] = self.catalogue['bbox4'].astype(float)
        self.catalogue['Maj'] = self.catalogue['Maj'].astype(float)
        self.catalogue['Min'] = self.catalogue['Min'].astype(float)
        self.catalogue['Pa'] = self.catalogue['Pa'].astype(float)
        self.catalogue['parent_tag'] = self.catalogue['parent_tag'].astype(float)
        self.catalogue['Class'] = self.catalogue['Class'].astype(float)
        if self.cutup:   
            self.catalogue['Y0_cutout'] = self.catalogue['Y0_cutout'].astype(float)
            self.catalogue['X0_cutout'] = self.catalogue['X0_cutout'].astype(float)
        self.catalogue['SNR'] = self.catalogue['SNR'].astype(float)
        self.catalogue['Noise'] = self.catalogue['Noise'].astype(float)



    def open_pb(self):

        hdul = fits.open(self.pb_PATH)
        image = hdul[0].data

        image_shape = image.shape

        if len(image_shape) == 4:

            image = np.squeeze(image, axis=(0,1))

        if len(image_shape) == 3:

            image = np.squeeze(image, axis=(0))

        return image



    def radio_characteristing(self,catalogue=None,cutout=None,cutout_pb=None,background_map=None):
        
        # get beam in pixels

        if cutout is None:
            self.Beam = self.calculate_beam()
            catalogue = self.catalogue
            image = self.image
            if self.pb_PATH is not None:
                pb_image = self.pb_image
            else:
                pb_image = None
            background_map = self.local_bg


        else:
            self.Beam = self.calculate_beam()
            image = cutout
            if self.pb_PATH is not None:
                pb_image = cutout_pb
            background_map = background_map


        # for each source in the catalogue create mask and measure properties. prephorm source flux correction.

        self.flux_correction_list = []

        params = []

        for i, source in tqdm(catalogue.iterrows(),total=len(catalogue),desc='Calculating Source Properties..',disable=not self.output):
            #print(i)
            try:
                mask = self.get_mask(row=source,image=image)
            except:
                continue
            #print(np.sum(mask))
            #plt.imshow(mask)
            #plt.show()

            source_props = self.get_region_props(mask,image=image)

            source_props = self.props_to_dict(source_props[0])

            peak_coords = np.where(image == source_props['max_intensity'])

            y_peak_loc = peak_coords[0][0]
            x_peak_loc = peak_coords[1][0]

            Model_Beam = self.model_beam_func(source_props['max_intensity'],image.shape,x_peak_loc,y_peak_loc,self.BMAJp/2,self.BMINp/2,self.BPA)
            Flux_correction_factor = self._flux_correction_factor(mask, Model_Beam)

            self.flux_correction_list.append(Flux_correction_factor)

            # calculate the flux of the source with option for pb correction.
             
            if self.pb_PATH is not None:
                
                background_mask = mask*background_map/self.sigma                    # fixed problem with slight offset.
                Flux_total = np.nansum(mask*image/pb_image - background_mask)/self.Beam   # may need to be altered for universality.
                Flux_peak = np.nanmax(mask*image/pb_image) - background_mask[y_peak_loc,x_peak_loc]
            
                # may need to be altered for universality.
                # get location of peak in the image
                
                Flux_peak_loc = np.where(image == Flux_peak)
                Flux_peak = Flux_peak - background_mask[Flux_peak_loc[0][0],Flux_peak_loc[1][0]]
                
            else:
                
                background_mask = mask*background_map/self.sigma                 # fixed problem with slight offset.
                #print(self.Beam)
                Flux_total = np.nansum(mask*image - background_mask)/self.Beam   # may need to be altered for universality.
                Flux_peak = np.nanmax(mask*image) - background_mask[y_peak_loc,x_peak_loc]
            
            #pdb.set_trace() # for debugging
            background_mask = np.where(background_mask == 0, np.nan, background_mask) # set background mask to nan where there is no background.
            Noise = np.nanmean(background_mask)
            #print('Noise: ',Noise)
            Flux_total = Flux_total*Flux_correction_factor
            
            Area = np.sum(mask)
            SNR = Flux_total/Noise
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
                           source.Class,
                           SNR,
                           Noise])
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
                                                  'Class',
                                                  'SNR',
                                                  'Noise'])

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
                    f.write('{:.2f},{:.2f}'.format(point[1], point[0])) # note this transformation as the index in some CARTA inmages start at -1.
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
        
        #plt.imshow(mask)
        #plt.show()
        
        #plt.imshow(Model_Beam)
        #plt.show()
        
                
        #raise ValueError('This function is not working correctly. Please use the flux correction factor from the catalogue.')
                
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


    def _xy_to_RaDec(self,x,y):
        # convert x,y to ra,dec
        """
        Convert an X and Y coordinate to RA and Dec using a header file with astropy.

        Parameters:
        x (float): The X coordinate.
        y (float): The Y coordinate.
        header_file (str): The path to the FITS header file.
        stokes (int): The Stokes dimension.
        freq (int): The frequency dimension.

        Returns:
        tuple: A tuple containing the RA and Dec in degrees.
        """ 
        stokes = 0 # stokes and freq are not used in this function.
        freq = 0 # stokes and freq are not used in this function.
        wcs = WCS(self.header)
        ra, dec, _, _ = wcs.all_pix2world(x, y, stokes, freq, 0)
        return ra, dec
        



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
        image_crop = self.image[int(Ymin):int(Ymax),int(Xmin):int(Xmax)]

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
            contour = self._get_polygons_in_bbox(row.bbox2-2,row.bbox4+2,row.bbox1-2,row.bbox3+2,row.x1,row.y1,row.Birth,row.Death)
            polygons.append(contour)

        self.polygons = polygons


    def set_background(self,detection_threshold,analysis_threshold,set_bg=None,bg_map=None,box_size=10):
        self.bg_map = bg_map
        self.sigma = detection_threshold
        self.analysis_threshold = analysis_threshold
        if self.mode == 'Radio':

            if self.cutup:
                
                # loop though each cutout and calculate the local background.
                
                if bg_map is not None:
                    
                    # users wants to use background map so lets make it
                    local_bg_list = []
                    analysis_threshold_list = []
                    for i, cutout in enumerate(self.cutouts):
                        local_bg_map = self.radio_background_map(cutout, box_size)
                        
                        analysis_threshold_list.append(local_bg_map*self.analysis_threshold)
                        local_bg_list.append(local_bg_map*self.sigma)
                    
                
                else:

                    local_bg_list = []
                    analysis_threshold_list = []
                    for cutout in self.cutouts:
                        local_bg = self.radio_background(cutout)
                        analysis_threshold_list.append(local_bg*self.analysis_threshold)
                        local_bg_list.append(local_bg*self.sigma)
                local_bg = local_bg_list
                analysis_threshold = analysis_threshold_list
                    
            else:

                # Radio background is calculated using the median absolute deviation of the total image.

                local_bg_o = self.radio_background(self.image)
                local_bg = local_bg_o*self.sigma
                analysis_threshold = local_bg_o*self.analysis_threshold

        if self.mode == 'Optical':
            # Optical background is calculated using a random sample of pixels
            mean_bg, std_bg = self.optical_background(nsamples=1000)
            local_bg = mean_bg + self.sigma*std_bg
            analysis_threshold = mean_bg + std_bg*self.analysis_threshold

        if self.mode == 'X-ray':
            # no implemented X-ray specific background function.
            mean_bg,std_bg = self.xray_background()
            local_bg = mean_bg + self.sigma*std_bg

        if self.mode == 'other':
            # If the user has a custom background function, they can pass it in here.
            local_bg = set_bg*self.sigma
            analysis_threshold = local_bg*self.analysis_threshold
            print('Background set to: ',local_bg)
            print('Analysis threshold set to: ',analysis_threshold)
            
        self.analysis_threshold_val = analysis_threshold
        self.local_bg = local_bg
        
        if bg_map:
            print('Using bg_map for analysis.')
        else:
            if self.cutup:
            
                print('Mean Background across cutouts: ', np.nanmean(self.local_bg))
            
            else:
                print('Background set to: ',self.local_bg)


 


    def radio_background_map(self,cutout2, box_size):
        '''
        This function takes an image and a box size and calculates the radio_background() for each box to create a map of local background.
        '''
        from astropy.stats import mad_std
        
        # step size
        step_size = box_size//2
        
        # initialize the map
        map_shape = (cutout2.shape[0]//step_size, cutout2.shape[1]//step_size)
        bg_map = np.zeros(map_shape)
        
        # box
        box = np.ones((box_size, box_size))
        
        for i in range(0, cutout2.shape[0], step_size):
            for j in range(0, cutout2.shape[1], step_size):
                # get the box
                box_image = cutout2[i:i+box_size, j:j+box_size]
                # calculate the radio background
                local_bg = mad_std(box_image, ignore_nan=True)
                # set the value in the map
                bg_map[i//step_size, j//step_size] = local_bg
                
        # make a version of the image where each pixel is assigned a loval bg value from the map
        # set nans to 0
        bg_map[np.isnan(bg_map)] = 0
        bg_image = np.zeros(cutout2.shape)
        for i in range(bg_map.shape[0]):
            for j in range(bg_map.shape[1]):
                bg_image[i*step_size:(i+1)*step_size, j*step_size:(j+1)*step_size] = bg_map[i, j]
                
        
        return bg_image
        
    
    
    
        
        


    def radio_background(self,image,metric='mad_std'):

        if metric == 'mad_std':
                
            from astropy.stats import mad_std

            local_bg = mad_std(image,ignore_nan=True)
            
        elif metric == 'rms':
            
            local_bg = np.sqrt(np.nanmean(image**2))
            
        else:
            raise ValueError('metric not recognised. Please use mad_std or rms')
        
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
        
        beam_size_arcseconds = self.header['BMIN']*3600
        BMIN_oversampled_spacial_width = (BMIN**2 + beam_size_arcseconds**2)**0.5
        BMIN = BMIN_oversampled_spacial_width/arcseconds_per_pixel
        
        self.BMAJp = BMAJ
        self.BMINp = BMIN
        try:
            self.BPA = self.header['BPA']
        except KeyError:
            self.BPA = 0
        return np.pi * (BMAJ)*(BMIN) / (4*np.log(2))