from .utils import preprocessing, create_subimages
from .utils import estimate_image_local_bg
import numpy as np
import pandas
import matplotlib.pyplot as plt
from .homology import cripser_homol
from .source_props import calculate_props as cp

class trsf:
    '''
    Main class for TRSF. This class will take an image and calculate the persistence diagrams and source properties.

    Parameters
    ----------
    img_path : str
        Path to the image to be analysed.
    cutup_img_size : int
        Size of the cutouts to be used. Default is 250.

    '''
    def __init__(self, img_path, cutup_img_size=250, method='cython', gaussian_fitting=False, region_expansion=False,cutup_img = True):
        self.img_path = img_path
        self.cutup_img = cutup_img 
        self.cutup_img_size = cutup_img_size
        self.confusion_limit = 3 # this needs to be read from the header.
        self.gaussian_fitting = gaussian_fitting
        self.region_expansion = region_expansion
        self.method = method
        #add here where it reads and prints image basic info.
        self._main()

    def _main(self):
        print('Starting TRSF')
        print('Image path: {}'.format(self.img_path))
        print('Attempting to open Image.')
        self._open_img()
        print('Calculating persistence diagrams and source properties.')
        print('This may take a while {}.'.format(len(self.Cutouts)))

        counter = 0
        for num, img in enumerate(self.Cutouts):
            # check if image is empty
            if img.sum() == 0:
                print('Empty image. Skipping.')
                continue
            else:

                local_bg, sigma = self._background_estimate(img)
                if np.isnan(local_bg):
                    continue
                #img[np.isnan(img)] = 0
                #img[np.isnan(img)] = 0
                pd = self._calculate_persistence_diagrams(img,local_bg,sigma)
                src_cat = self._create_component_catalogue(pd,img,local_bg,sigma)
                # alter catalogue to include the cutout coordinates
                src_cat['x'] = src_cat['x'] + self.Coords[counter][0]
                src_cat['y'] = src_cat['y'] + self.Coords[counter][1]
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

                if counter == 0 :
                    self.catalogue = src_cat
                else:
                    self.catalogue = pandas.concat([self.catalogue,src_cat])
                counter += 1
        print('TRSF finished.')

    def _open_img(self):
        self.full_img = preprocessing.preprocess(self.img_path,formatting=True).img 
        self.full_img[np.isnan(self.full_img)] = 0
        print(self.full_img.shape)
        self.full_img = self._crop_image(self.full_img)
        print(self.full_img.shape)
        if len(self.full_img.shape) > 2:
            self.full_img = self.full_img[:,:,0]
        if self.cutup_img == True:
            self.Cutouts, self.Coords = create_subimages.create_cutouts(self.full_img,size=self.cutup_img_size)
        else:
            self.Cutouts = [self.full_img]
        print('Image opened and cut into {} pieces.'.format(len(self.Cutouts)))
        
    def _background_estimate(self,img):
        local_bg, sigma = estimate_image_local_bg.estimate_bg_from_homology(img)
        return local_bg, sigma

    def _calculate_persistence_diagrams(self,img,local_bg,sigma):
        pd = cripser_homol.compute_ph_cripser(img,local_bg,sigma,maxdim=0)
        pd = cripser_homol.apply_confusion_limit(pd,self.confusion_limit)
        pd = cripser_homol.ph_precocessing(pd,img,local_bg,sigma)
        return pd
    
    def _create_component_catalogue(self,pd,img,local_bg,sigma):
        props = cp.cal_props(pd,img,local_bg,sigma,method=self.method)
        source_catalogue = props.fit_all_single_sources(gaussian_fit=self.gaussian_fitting,expand=self.region_expansion)
        return source_catalogue
    
    def _crop_image(self,image):
        # Find the indices of the non-zero elements in each row and column
        rows = np.where(image.any(axis=1))[0]
        cols = np.where(image.any(axis=0))[0]

        # Calculate the cropping boundaries
        min_row, max_row = min(rows), max(rows)
        min_col, max_col = min(cols), max(cols)

        # Crop the image
        cropped_image = image[min_row:max_row+1, min_col:max_col+1]

        return cropped_image