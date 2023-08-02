'''
author: Rhys Shaw
email: rhys.shaw@bristol.ac.uk
date: 14-07-2023
description: This file contains the main class for TRSF.

'''

from .utils import preprocessing, create_subimages
from .utils import estimate_image_local_bg
import numpy as np
import pandas
import matplotlib.pyplot as plt
from .homology import cripser_homol
from .source_props import calculate_props as cp
from tqdm import tqdm
import time
import warnings
from astropy.io import fits
from astropy.table import Table

class trsf:
    '''
    Main class for TRSF. This class will take an image and calculate the persistence diagrams and source properties.

    Parameters
    ----------
    img_path : str
        Path to the image to be analysed.
    cutup_img_size : int
        Size of the cutouts to be used. Default is 250.

    Attributes'
    ----------
    img_path : str
        Path to the image to be analysed.
    cutup_img_size : int
        Size of the cutouts to be used. Default is 250.
    confusion_limit : int   
        The confusion limit of the image. This is read from the header of the image.
    gaussian_fitting : bool
        If True, the sources will be fit with a 2D Gaussian. Default is False.
    region_expansion : bool
        If True, the sources will be expanded to improve the Gaussian fitting. Default is False.
    method : str
        Method to be used for calculating the persistence diagrams. Default is 'cython'.
    sigma : float
        Sigma value to be used for calculating the persistence diagrams. Default is None.

    Methods
    -------
    save_catalogue(path: str, type: str) -> None
        Saves the catalogue to the specified path. The type can be 'csv', 'hdf', 'fits' or 'json'.
        
    '''
    def __init__(self, img_path, cutup_img_size=250, method='cython', 
                 gaussian_fitting=False, region_expansion=False,cutup_img = True,
                 sigma=None,plot=True,expsigma=3,smoothing=False,smooth_param=1):
        self.img_path = img_path
        self.cutup_img = cutup_img 
        self.cutup_img_size = cutup_img_size
        self.confusion_limit = 3 # this needs to be read from the header.
        self.gaussian_fitting = gaussian_fitting
        self.region_expansion = region_expansion
        self.method = method
        self.sigma = sigma
        self.sum_plot = plot
        self.expsigma = expsigma
        self.smoothing = smoothing
        self.smooth_param = smooth_param
        #add here where it reads and prints image basic info.
        self._main()
    
    def save_catalogue(self, path: str, type: str) -> None:
        self._set_data_types()
        if type == 'csv':
            self.catalogue.to_csv(path)
            print('Catalogue saved to {}'.format(path))
        if type == 'hdf':
            self.catalogue.to_hdf(path,key='catalogue')
            print('Catalogue saved to {}'.format(path), 'with key "catalogue"')
        if type == 'fits':
            table = Table.from_pandas(self.catalogue)
            table.write(path, format='fits', overwrite=True)
            print('Catalogue saved to {}'.format(path))
        if type == 'json':
            self.catalogue.to_json(path)
            print('Catalogue saved to {}'.format(path))


    def _main(self):
        # suppress all warnings
        warnings.filterwarnings("ignore")
        t0 = time.time()
        print("""   
###########################
 _____   ___    ___    ___ 
|_   _| | _ \  / __|  | __|
  | |   |   /  \__ \  | _| 
  |_|   |_|_\  |___/  |_|  
        
###########################
Topological Radio Source Finder.
        """)
        print('-------------------')
        print('Starting TRSF')
        print('NOTICE: Image path: {}'.format(self.img_path))
        print('Attempting to open Image...')
        self._open_img()
        print('Calculating persistence diagrams and source properties..')
        #self.full_img[np.isnan(self.full_img)] = 0
        counter = 0
        for num, img in tqdm(enumerate(self.Cutouts),total=len(self.Cutouts),desc='Cutouts Completed'):
            
            # check if image is empty
            if img.sum() == 0:
                #print('Empty image. Skipping.')
                continue
            else:

                #plt.imshow(img)
                #plt.show()
                
                local_bg, sigma = self._background_estimate(img)
                if self.sigma != None:
                    sigma = self.sigma
                #print('local_bg', mad_std(img, ignore_nan=True))
                if np.isnan(local_bg):
                    continue
                img[np.isnan(img)] = 0
                
                pd = self._calculate_persistence_diagrams(img,local_bg,sigma)
                if len(pd) == 0:
                    continue
                src_cat = self._create_component_catalogue(pd,img,local_bg,sigma)
                # alter catalogue to include the cutout coordinates
                src_cat['Yc'] = src_cat['x_c'] + self.Coords[num][0]
                src_cat['Xc'] = src_cat['y_c'] + self.Coords[num][1]
                src_cat['x'] = src_cat['x'] + self.Coords[num][1]
                src_cat['y'] = src_cat['y'] + self.Coords[num][0]
                # drop the columns that are not needed
                src_cat = src_cat.drop(columns=['x_c','y_c'])
                # rename the columns
                src_cat = src_cat.rename(columns={'Xc':'x_c','Yc':'y_c'})
                # alter bbox to include the cutout coordinates
                bbox_list = src_cat['bbox'].tolist()
                
                new_bbox = []
                for i in range(len(bbox_list)):
                    if isinstance(bbox_list[i],float):
                        if np.isnan(bbox_list[i]):
                            new_bbox.append(None)
                    else:
                        bbox = list(bbox_list[i])
                        bbox[0] += self.Coords[num][0]
                        bbox[1] += self.Coords[num][1]
                        bbox[2] += self.Coords[num][0]
                        bbox[3] += self.Coords[num][1]
                        new_bbox.append(bbox)
                src_cat['bbox'] = new_bbox
                
                # see if src has column 'polygon'
                if 'polygon' not in src_cat.columns:
                    src_cat['polygon'] = None
                else:
                    polygon_list = src_cat['polygon'].tolist()
                    new_polygon = []
                    for i in range(len(polygon_list)):
                        if isinstance(polygon_list[i],float):
                            if np.isnan(polygon_list[i]):
                                new_polygon.append(None)
                        else:
                            polygon = polygon_list[i]
                            polygon[:,0] += self.Coords[num][0]
                            polygon[:,1] += self.Coords[num][1]
                            new_polygon.append(polygon)
                    src_cat['polygon'] = new_polygon

                if counter == 0 :
                    self.catalogue = src_cat
                else:
                    self.catalogue = pandas.concat([self.catalogue,src_cat])
                counter += 1
        
        if self.sum_plot == True:
            self._summary_plots()
        print('TRSF finished.')
        print('Time taken: {} seconds'.format(time.time()-t0))
        print('-------------------')

    def _set_data_types(self):
        # set any nan to 0
        self.catalogue = self.catalogue.fillna(0)
        self.catalogue['index'] = self.catalogue['index'].astype(int)
        self.catalogue['amp'] = self.catalogue['amp'].astype(float)
        self.catalogue['x'] = self.catalogue['x'].astype(float)
        self.catalogue['y'] = self.catalogue['y'].astype(float)
        self.catalogue['sigma_x'] = self.catalogue['sigma_x'].astype(float)
        self.catalogue['sigma_y'] = self.catalogue['sigma_y'].astype(float)
        self.catalogue['theta'] = self.catalogue['theta'].astype(float)
        self.catalogue['peak_flux'] = self.catalogue['peak_flux'].astype(float)
        self.catalogue['x_c'] = self.catalogue['x_c'].astype(float)
        self.catalogue['y_c'] = self.catalogue['y_c'].astype(float)
        self.catalogue['bbox'] = self.catalogue['bbox'].astype(str)
        self.catalogue['Class'] = self.catalogue['Class'].astype(int)
        self.catalogue['Birth'] = self.catalogue['Birth'].astype(float)
        self.catalogue['Death'] = self.catalogue['Death'].astype(float)
        self.catalogue['x1'] = self.catalogue['x1'].astype(int)
        self.catalogue['y1'] = self.catalogue['y1'].astype(int)
        self.catalogue['lifetime'] = self.catalogue['lifetime'].astype(float)
        self.catalogue['x2'] = self.catalogue['x2'].astype(int)
        self.catalogue['y2'] = self.catalogue['y2'].astype(int)
        self.catalogue['polygon'] = self.catalogue['polygon'].astype(str)
        self.catalogue['encloses_i'] = self.catalogue['encloses_i'].astype(str)

    def _open_img(self):
        self.full_img = preprocessing.preprocess(self.img_path,formatting=True).img 
        print('NOTICE: Input Image Size {}'.format(self.full_img.shape))
        self.full_img = self._crop_image(self.full_img)
        copy_full_img = self.full_img.copy()
        self.full_img[np.isnan(self.full_img)] = 0

        if self.smoothing == True:
            self.full_img = self._image_smoothing(img=self.full_img,smooth_param = self.smooth_param)

        print('NOTICE: Image Size with reduced padding {}'.format(self.full_img.shape))
        if len(self.full_img.shape) > 2:
            self.full_img = self.full_img[:,:,0]
        if self.cutup_img == True:
            self.Cutouts, self.Coords = create_subimages.create_cutouts(copy_full_img,size=self.cutup_img_size)
        else:
            self.Cutouts = [self.full_img]
        print('NOTICE: Image opened and cut into {} pieces.'.format(len(self.Cutouts)))
        
    def _background_estimate(self,img):
        local_bg, sigma = estimate_image_local_bg.estimate_bg_from_homology(img)
        return local_bg, sigma

    def _calculate_persistence_diagrams(self,img,local_bg,sigma):
        pd = cripser_homol.compute_ph_cripser(img,local_bg,sigma,maxdim=0)
        pd = cripser_homol.apply_confusion_limit(pd,self.confusion_limit)
        
        if len(pd) == 0:
            print('Empty persistence diagram. Skipping.')
        else:
            try:
                pd = cripser_homol.ph_precocessing(pd,img,local_bg,sigma)
            except:
                print('Error in preprocessing. Skipping.')
                print(pd)
                raise ValueError
        return pd
    
    def _create_component_catalogue(self,pd,img,local_bg,sigma):
        props = cp.cal_props(pd,img,local_bg,sigma,method=self.method,expsigma=self.expsigma)
        source_catalogue = props.fit_all_single_sources(gaussian_fit=self.gaussian_fitting,expand=self.region_expansion)
        return source_catalogue
    
    def _crop_image(self,arr):
        arr = np.array(arr)  # Convert input to numpy array
        mask_rows = np.all(np.isnan(arr), axis=1)
        mask_cols = np.all(np.isnan(arr), axis=0)
        return arr[~mask_rows][:, ~mask_cols]

    def _summary_plots(self):
        # plot image with source locations
        plt.figure(figsize=(8,8))
        plt.imshow(self.full_img,vmax=np.nanpercentile(self.full_img,95),vmin=np.nanpercentile(self.full_img,0.1))
        plt.title('Image with Source Locations')
        plt.scatter(self.catalogue['x_c'],self.catalogue['y_c'],color='red',s=1,marker='+')
        # plot polygons
        for i in range(len(self.catalogue)):
            try:
                polygon = self.catalogue.iloc[i].polygon
                plt.plot(polygon[:,1],polygon[:,0],color='red',alpha=0.5)
            except:
                pass
        plt.legend()
        plt.show()
    
    def _image_smoothing(self,img,smooth_param):
        # import gaussian filter
        from scipy.ndimage import gaussian_filter
        # smooth image
        img = gaussian_filter(img, sigma=smooth_param) # std of gaussian kernel
        return img