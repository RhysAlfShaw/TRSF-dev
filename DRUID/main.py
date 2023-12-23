'''
file: main.py
author: Rhys Shaw
date: 17-10-2023
description: This file contains the main class for the new source finder.
'''

from DRUID.homology2D import compute_ph_components
import numpy as np
from skimage import measure
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import pandas
import setproctitle
from multiprocessing import Pool
import time
import os
from scipy.ndimage import label
import pdb
import pandas as pd
from scipy import ndimage


try:
    import cupy as cp
    from cupyx.scipy.ndimage import label as cupy_label

    # Get the list of available GPU devices
    available_gpus = cp.cuda.runtime.getDeviceCount()

    if available_gpus > 0:
        print(f"Number of GPUs available: {available_gpus}, DRUID with GPU Acceleration Avalible.")
        GPU_AVAILABLE = True
    else:
        
    #  from scipy.ndimage import label
        
        print('No GPU available. DRUID GPU acceleration will not be available.')
        GPU_AVAILABLE = False
    
except:
    
 #   from scipy.ndimage import label
    print('Could not import cupy. DRUID GPU acceleration will not be available.')
    GPU_AVAILABLE = False
    


# set the process title
setproctitle.setproctitle('DRUID')

## Add checker of path function.
## Create a similar style known flux density image and compare the corrected flux densities to know values.
## plot the difference as a funciton of correction factor.

class sf:

    
    def __init__(self,image,image_PATH,mode,pb_PATH=None,cutup=False,cutup_size=500,output=True,area_limit=1,smooth_sigma=1,nproc=4,GPU=False):
       
        print("""   
              
#############################################
_______   _______          _______  _______  
(  __  \ (  ____ )|\     /|\__   __/(  __  \ 
| (  \  )| (    )|| )   ( |   ) (   | (  \  )
| |   ) || (____)|| |   | |   | |   | |   ) |
| |   | ||     __)| |   | |   | |   | |   | |
| |   ) || (\ (   | |   | |   | |   | |   ) |
| (__/  )| ) \ \__| (___) |___) (___| (__/  )
(______/ |/   \__/(_______)\_______/(______/ 
        
#############################################
Detector of astRonomical soUrces in optIcal and raDio images

Version: 0.0.1
For more information see:
https://github.com/RhysAlfShaw/TRSF-dev
        """)
        print('Initialising Source Finder..')
        self.cutup = cutup
        self.output = output
        self.image_PATH = image_PATH
        self.area_limit = area_limit
        self.smooth_sigma = smooth_sigma
        self.GPU = GPU
        if self.GPU:
            print('GPU=TRUE | GPU acceleration enabled.')
        # set nproc as a env variable
        self.nproc = nproc

        if image_PATH is None:
            self.image = image
            self.header = None
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
                catalogue = compute_ph_components(cutout,self.local_bg[i],analysis_threshold_val=self.analysis_threshold_val[i],lifetime_limit=lifetime_limit,output=self.output,bg_map=self.bg_map,area_limit=self.area_limit,nproc=self.nproc,GPU=self.GPU)
                print('Done ')
                # add cutout coords to catalogue
                catalogue['Y0_cutout'] = self.coords[i][0]
                catalogue['X0_cutout'] = self.coords[i][1]
                catalogue_list.append(catalogue)

            # combine the catalogues
            self.catalogue = pandas.concat(catalogue_list)

        else:

            self.catalogue = compute_ph_components(self.image,self.local_bg,analysis_threshold_val=self.analysis_threshold_val,lifetime_limit=lifetime_limit,output=self.output,bg_map=self.bg_map,area_limit=self.area_limit,nproc=self.nproc,GPU=self.GPU)




    def gaussian2dkernal_convolution(self, image, sigma):
        '''
        This function takes an image and a sigma value and convolves the image with a 2D gaussian kernel.
        '''
        from scipy.ndimage import gaussian_filter
        print(sum(image.flatten()))
        image = gaussian_filter(image, sigma=sigma)
        print(sum(image.flatten()))
        return image





    def source_characterising(self,use_gpu=False):

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
            self.optical_characteristing(use_gpu=use_gpu)
        
            

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

    def get_bounds(self,mask):
        
        # Assuming 'mask' is your binary mask
        #labeled_mask, num_labels = ndimage.label(mask)
        bounding_boxes = ndimage.find_objects(mask)

        # Now bounding_boxes is a list of slice objects, which represent the bounding boxes of the labeled regions
        # You can access the bounding box of the first region like this:
        bbox = bounding_boxes[0]

        # bbox is a tuple of slice objects, which can be used to index numpy arrays
        # Here's how you can get the coordinates of the bounding box:
        min_row, max_row = bbox[0].start, bbox[0].stop
        min_col, max_col = bbox[1].start, bbox[1].stop
        # output as minx, maxx, miny, maxy
        return min_col, max_col, min_row, max_row
    
    
    def large_mask_red_image_procc(self,Birth,Death,x1,y1,image):
        '''
        Does all gpu processing for the large mask.
        
        return the red_image and red_mask and the bounding box.
        
        '''
        
        
        mask = cp.zeros(image.shape,dtype=cp.bool_)
        mask = cp.logical_or(mask,cp.logical_and(image <= Birth, image > Death))
       # mask_enclosed = self.get_enclosing_mask_gpu(y1,x1,mask)
        labeled_mask, num_features = cupy_label(mask)
        
        # Check if the specified pixel is within the mask
        if 0 <= y1 < mask.shape[1] and 0 <= x1 < mask.shape[0]:
            label_at_pixel = labeled_mask[x1, y1]
            
            if label_at_pixel != 0:
                # Extract the connected component containing the specified pixel
                component_mask = (labeled_mask == label_at_pixel)
                #plt.imshow(component_mask.get())
                #plt.show()
                #pdb.set_trace()
                non_zero_indices = cp.nonzero(component_mask)

                # Extract minimum and maximum coordinates
                xmin = cp.min(non_zero_indices[1])
                ymin = cp.min(non_zero_indices[0])
                xmax = cp.max(non_zero_indices[1])
                ymax = cp.max(non_zero_indices[0])
                
                # images are not being cropped?
                
                red_image = image[ymin:ymax+1, xmin:xmax+1]
                red_mask = component_mask[ymin:ymax+1, xmin:xmax+1]
                
                return red_image, red_mask, xmin, xmax, ymin, ymax
        


    def optical_characteristing(self,use_gpu,catalogue=None,cutout=None,background_map=None):
        '''
        
        Characterising the Source Assuming the input image of of the format of a optical astronomical image.
        
        # needs to work on cutup images.
        # work on the whole image.
    
        '''
        # single image only.
        if cutout is None:
            catalogue = self.catalogue
            image = self.image
            background_map = self.local_bg
            
        else:
            image = cutout
            background_map = background_map
            
        #return self.catalogue


        params = []
        
        # main loop over the catalogue.
        
        #image_gpu = cp.asarray(image, dtype=cp.float64)
        # put the image on the GPU memory
        #Deaths = catalogue['Death']
        #list_to_remove = []
        #for i in range(len(Deaths)):
        #    if type(Deaths.iloc[i]) == pd.core.series.Series:
        #        #print(Deaths.iloc[i])
        #        list_to_remove.append(i)
        #catalogue = catalogue.drop(catalogue.index[list_to_remove])
        
        self.image_gpu = cp.asarray(image, dtype=cp.float64)
        shape = self.image_gpu.shape
        ## # tuple to array
        #shape = np.asarray(shape)
        # # convert to cupy array
        #shape = cp.asarray(shape, dtype=cp.int64)
        x1 = catalogue['x1'].to_numpy()
        y1 = catalogue['y1'].to_numpy()
        x2 = catalogue['x2'].to_numpy()
        y2 = catalogue['y2'].to_numpy()
        Birth = catalogue['Birth'].to_numpy()
        Death = catalogue['Death'].to_numpy()
        parent_tag = catalogue['parent_tag'].to_numpy()
        Class = catalogue['Class'].to_numpy()
        #t0 = time.time()
       # vectorised_get_mask_gpu = cp.vectorize(self.get_mask_gpu)
        #mask = cp.zeros(shape,dtype=cp.float64)
        #num = 6
        #masks = vectorised_get_mask_gpu(Birth,Death,x1,y1,self.image_gpu,shape,num)
        #t1 = time.time()
        # print('Time to vectorise get_mask_gpu: ',t1-t0)
        # print(masks)
        #return vectorised_get_mask_gpu
    
       #print('Time to vectorise get_mask_gpu: ',t1-t0)
        print(len(catalogue))
        print(len(Birth))
        polygons = []
        for i, source in tqdm(enumerate(Birth),total=len(Birth),desc='Calculating Source Properties..',disable=not self.output):
        #     if use_gpu:
        #         #try:
        #         t0_get_mask = time.time()
        #         #print('source: ',source)
        #         mask = self.get_mask_gpu(Birth[i],Death[i],x1[i],y1[i],self.image_gpu).get()
        #         # mask to interger
        # #        mask_gpu = cp.asarray(mask, dtype=cp.int64)
        #         mask = mask.astype(int)
        #         #print('mask: ',mask)
        #         t1_get_mask = time.time()
        #         #
        #         print('Time to get mask: ',t1_get_mask-t0_get_mask)
                
        #         #except:
        #         #   continue
        #     else:
        #         try:
                    
        #             #t0_get_mask = time.time()
        #             mask = self.get_mask(row=source,image=image)
                    
        #             #t1_get_mask = time.time()
        #             #print('Time to get mask: ',t1_get_mask-t0_get_mask)
        #         except:
        #             continue
            # t0 = time.time()
            # min_col, max_col, min_row, max_row = self.get_bounds(mask)
            # t1 = time.time()
            # print('Time to get boundingbox: ',t1-t0)
            # t0_crop = time.time()
            # red_image = image[min_row:max_row,min_col:max_col]
            # red_mask = mask[min_row:max_row,min_col:max_col]
            # t1_crop = time.time()
            # print('Time to crop image and mask: ',t1_crop-t0_crop)
            #t0_large_mask = time.time()
            red_image, red_mask, xmin, xmax, ymin, ymax = self.large_mask_red_image_procc(Birth[i],Death[i],x1[i],y1[i],self.image_gpu)
            #print(xmin,xmax,ymin,ymax)
            red_mask = red_mask.astype(int)
            # make sure red_image and red_mask are on the cpu.
            red_image = red_image.get()
            red_mask = red_mask.get()
            xmin = xmin.get()
            xmax = xmax.get()
            ymin = ymin.get()
            ymax = ymax.get()
            # create the polygons here
            
            contour = self._get_polygons_in_bbox(xmin,xmax,ymin,ymax,x1[i],y1[i],Birth[i],Death[i],red_mask)
            #print(red_image.shape)
            #t1_large_mask = time.time()
            #print('Time to create large mask: ',t1_large_mask-t0_large_mask)
            #t0_region_props = time.time()
            source_props = self.get_region_props(red_mask,image=red_image)
            #t1_region_props = time.time()
            #print('Time to get region props: ',t1_region_props-t0_region_props)
            #t0_props_to_dict = time.time()
            source_props = self.props_to_dict(source_props[0])
            #t1_props_to_dict = time.time()
            #print('Time to convert props to dict: ',t1_props_to_dict-t0_props_to_dict)
            #t0_peak_coords = time.time()
            #red_image = Use the bounding box to reduce the size of the image.
            #red_image = image[int(source_props['bbox'][0]):int(source_props['bbox'][2]),int(source_props['bbox'][1]):int(source_props['bbox'][3])]
            #red_mask = mask[int(source_props['bbox'][0]):int(source_props['bbox'][2]),int(source_props['bbox'][1]):int(source_props['bbox'][3])]
            
            #t0_model_beam = time.time()
            #background_map = np.asarray(background_map, dtype=cp.float64)
            #mask = cp.asarray(mask, dtype=cp.float64)
            #background_map = cp.asarray(background_map, dtype=cp.float64)
            red_background_mask = np.where(red_mask == 0, np.nan, red_mask*background_map)
            
            #t1_model_beam = time.time()
            #print('Time to create background mask: ',t1_model_beam-t0_model_beam)
            #t0_noise = time.time()
            Noise = np.nanmean(red_background_mask)
            #t1_noise = time.time()
            #print('Time to calculate noise: ',t1_noise-t0_noise)
            
            #t0_fluxes = time.time()
            Flux_total = np.nansum(red_mask*red_image - red_background_mask)
            #t1_fluxes = time.time()
            #print('Time to calculate fluxes: ',t1_fluxes-t0_fluxes)
            #t0_define_props = time.time()
            Area = source_props['area']
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
            
            params.append([i,
                            Birth[i],
                            Death[i],
                            x1[i],
                            y1[i],
                            x2[i],
                            y2[i],
                            Flux_total,
                            Flux_total,
                            1,
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
                            parent_tag[i],
                            Class[i],
                            SNR,
                            Noise])
            polygons.append(contour)
            #t1_define_props = time.time()
            # ##
            #print('Time to define props: ',t1_define_props-t0_define_props)
        self.polygons = polygons
        if self.cutup:
            # return as we are looping over the cutouts.
            cutup_Catalogue = self.create_params_df(params)
            return cutup_Catalogue
        else:
            # return nothing as we only have one catalogue.
            #t0_create_params_df = time.time()
            self.create_params_df(params)
            #t1_create_params_df = time.time()
            #print('Time to create params df: ',t1_create_params_df-t0_create_params_df)
        
        

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


    def get_enclosing_mask_gpu(self, x, y, mask):
        '''
        
        Returns the connected components inside the mask starting from the point (x, y) using GPU.
        
        '''
        # which 
        labeled_mask, num_features = cupy_label(mask)
        
        # Check if the specified pixel is within the mask
        if 0 <= x < mask.shape[1] and 0 <= y < mask.shape[0]:
            label_at_pixel = labeled_mask[y, x]
            
            if label_at_pixel != 0:
                # Extract the connected component containing the specified pixel
                component_mask = (labeled_mask == label_at_pixel)
                return cp.asarray(component_mask)
            else:
                # The specified pixel is not part of any connected component
                return None
        else:
            # The specified pixel is outside the mask
            return None

    
    def get_mask(self,row,image):

        #t0_numpyMask = time.time()
        mask = np.zeros(image.shape)
        mask = np.logical_or(mask,np.logical_and(image <= row.Birth, image > row.Death))
        #t1_numpyMask = time.time()
        #print('Time to create numpy mask: ',t1_numpyMask-t0_numpyMask)
        #t0_enclosingMask = time.time()
        mask_enclosed = self.get_enclosing_mask(int(row.y1),int(row.x1),mask)
        #t1_enclosingMask = time.time()
        #print('Time to create enclosing mask: ',t1_enclosingMask-t0_enclosingMask)
        # set mask as integer
        mask_enclosed = mask_enclosed.astype(int)
        return mask_enclosed


    def get_mask_gpu(self,Birth,Death,x1,y1,img):
    
        # same as get_mask but uses gpu with the cupy library
        #t0 = time.time()
        
        mask = cp.zeros(img.shape,dtype=cp.float64)
        mask = cp.logical_or(mask,cp.logical_and(img <= Birth,img > Death))
        mask_enclosed = self.get_enclosing_mask_gpu(y1,x1,mask)
        
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
        # slow function when all are called.
        '''

        dict = {
            'area': regionprops.area,
            'bbox': regionprops.bbox,
            #'bbox_area': regionprops.bbox_area,
            'centroid': regionprops.centroid,
            #'convex_area': regionprops.convex_area,
            'eccentricity': regionprops.eccentricity,
            #'equivalent_diameter': regionprops.equivalent_diameter,
            #'euler_number': regionprops.euler_number,
            #'extent': regionprops.extent,
            #'filled_area': regionprops.filled_area,
            'major_axis_length': regionprops.major_axis_length,
            'minor_axis_length': regionprops.minor_axis_length,
            #'moments': regionprops.moments,
            #'perimeter': regionprops.perimeter,
            #'solidity': regionprops.solidity,
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
        
        return correction_factor



    def model_beam_func(self,peak_flux,shape,x,y,bmaj,bmin,bpa):
        model_beam = np.zeros(shape)
        model_beam = self.generate_2d_gaussian(peak_flux,shape,(x,y),bmaj,bmin,bpa,norm=False)
        return model_beam





    


    def get_enclosing_mask_new(self, x, y, mask):
        '''
        Returns the connected components inside the mask starting from the point (x, y).
        '''
        labeled_mask, num_features = label(mask)
        
        # Check if the specified pixel is within the mask
        if 0 <= x < mask.shape[1] and 0 <= y < mask.shape[0]:
            label_at_pixel = labeled_mask[y, x]
            
            if label_at_pixel != 0:
                # Extract the connected component containing the specified pixel
                component_mask = (labeled_mask == label_at_pixel)
                return component_mask
            else:
                # The specified pixel is not part of any connected component
                return None
        else:
            # The specified pixel is outside the mask
            return None



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






    def _get_polygons(self,x1,y1,birth,death):

        '''
        Returns the polygon of the enclosed area of the point (x,y) in the mask.
        '''

        mask = np.zeros(self.image.shape)
        mask = np.logical_or(mask,np.logical_and(self.image <= birth,self.image > death))

        mask = self.get_enclosing_mask(int(y1),int(x1),mask)
        contour = measure.find_contours(mask,0)[0]

        return contour
    
    
    
    
    def _get_polygons_gpu(self,x1,y1,birth,death):
        '''
        Returns the polygon of the enclosed area of the point (x,y) in the mask.
        '''
        # is the image on the GPU memory?
        #if not cp.is_cuda_array(self.image):
        #    self.image = cp.asarray(self.image, dtype=cp.float64)
            
        mask = cp.zeros(self.image.shape)
        mask = cp.logical_or(mask,cp.logical_and(self.image_gpu <= birth,self.image_gpu > death))
        mask = self.get_enclosing_mask_gpu(int(y1),int(x1),mask)
        contour = measure.find_contours(mask,0)[0]
        
        return contour 




    def _xy_to_RaDec(self,x,y):
        
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
        
        if self.mode == 'Radio': 
            
            stokes = 0  # stokes and freq are not used in this function.
            freq = 0    # stokes and freq are not used in this function.
            wcs = WCS(self.header)
            ra, dec, _, _ = wcs.all_pix2world(x, y, stokes, freq, 0)
            
        elif self.mode == 'Optical':
            # image is 2d so no stokes or freq
            ra, dec = wcs.all_pix2world(x, y, 0)
        
        return ra, dec
        



    def create_polygons(self,use_nproc=False,nproc=4):
        '''
        
        Parrallelised version of create_polygons. not recommended working.
        
        '''
        def process(i):
            
            row = self.catalogue.iloc[i]
            
            try:
            
                return self._get_polygons(row.x1, row.y1, row.Birth, row.Death)
            
            except:
            
                return None
                
        if use_nproc:
            
            self.nproc = nproc
            print('Creating polygons with {} processes'.format(self.nproc))
            t0 = time.time()
            with Pool() as pool:
                polygons = list(pool.imap(process, range(len(self.catalogue)), chunksize=len(self.catalogue)//self.nproc))
            t1 = time.time()
            print('Time to create polygons: ',t1-t0)
        
        else:
            
            polygons = []
            print(len(self.catalogue))
            for index, row in tqdm(self.catalogue.iterrows(),total=len(self.catalogue),desc='Creating polygons'):
            
                try:
            
                    contour = self._get_polygons(row.x1,row.y1,row.Birth,row.Death)
                    polygons.append(contour)
            
                except:
            
                    continue
            
        self.polygons = polygons




    def _get_polygons_in_bbox(self,Xmin,Xmax,Ymin,Ymax,x1,y1,birth,death,mask,pad=1):
        
        
        mask = np.pad(mask, pad, mode='constant', constant_values=0)    
        contour = measure.find_contours(mask, 0)[0]
        contour = measure.find_contours(mask, 0)[0]
        
        # remove the border
        contour[:,0] -= pad
        contour[:,1] -= pad
        
        # correct the coordinates to the original image
        contour[:,0] += Ymin
        contour[:,1] += Xmin
        
        return contour




    def create_polygons_fast(self):
        
        # since we have a bounding box, we can just create a polygon in the bounding box.
        
        polygons = []
        for index, row in tqdm(self.catalogue.iterrows(),total=len(self.catalogue),desc='Creating polygons'):
            contour = self._get_polygons_in_bbox(row.bbox2-2,row.bbox4+2,row.bbox1-2,row.bbox3+2,row.x1,row.y1,row.Birth,row.Death)
            polygons.append(contour)

        self.polygons = polygons
        
        
        
        
    def create_polygons_gpu(self):
        '''
        Recommended when using GPU acceleration. and not using the bounding box to simplify the polygon creation.
        '''
        polygons = []
        self.image_gpu = cp.asarray(self.image, dtype=cp.float64)
        
        for index, row in self.catalogue.iterrows():
            t0 = time.time()
            contour = self._get_polygons_gpu(row.x1,row.y1,row.Birth,row.Death)
            t1 = time.time()
            print('Time to create polygon: ',t1-t0)
            polygons.append(contour)
        self.polygons = polygons




    def set_background(self,detection_threshold,analysis_threshold,set_bg=None,bg_map=None,box_size=10,mode='Radio'):
        self.bg_map = bg_map
        self.sigma = detection_threshold
        self.analysis_threshold = analysis_threshold
        print(mode)
        if mode == 'Radio':

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
                if bg_map is not None:
                    local_bg_o = self.radio_background_map(self.image, box_size)
                    local_bg = local_bg_o*self.sigma
                    analysis_threshold = local_bg_o*self.analysis_threshold    
                else:
                    local_bg_o = self.radio_background(self.image)
                    local_bg = local_bg_o*self.sigma
                    analysis_threshold = local_bg_o*self.analysis_threshold

        if mode == 'Optical':
            # Optical background is calculated using a random sample of pixels
            mean_bg, std_bg = self.optical_background(nsamples=1000)
            local_bg = mean_bg + self.sigma*std_bg
            analysis_threshold = mean_bg + std_bg*self.analysis_threshold

        if mode == 'X-ray':
            # no implemented X-ray specific background function.
            mean_bg,std_bg = self.xray_background()
            local_bg = mean_bg + self.sigma*std_bg

        if mode == 'other':
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
    
        # now upsample the map to the original image size
        bg_map = np.repeat(bg_map, step_size, axis=0)
        bg_map = np.repeat(bg_map, step_size, axis=1)
        
        # shift the map to the correct position
            
        return bg_map
        
        
        
    
    
    
        
        


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
    
    
    def plot_sources(self,cmap,figsize=(10,10),norm='linear'):

        plt.figure(figsize=figsize)
        # scale the cmap to be logarithmic
        #cmap = plt.cm.get_cmap('viridis')
        plt.imshow(self.image,cmap=cmap,origin='lower',norm=norm)
        
        for i, poly in enumerate(self.polygons):
            if poly is not None:
                plt.plot(poly[:,1],poly[:,0])
        plt.show()