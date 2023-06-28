"""
Author: Rhys Shaw
E-mail: rhys.shaw@bristol.ac.uk
Date: 26-06-2023
Description: This file contains the functions for the parallelisation of the Gudhi Cubical Complex Filtration.

"""

from tqdm import *
from multiprocessing import Pool
import numpy as np
from src.TRSF import create_regions

def _iterator(args):
    i, idx, dgm, img, threshold, detection, Xx, Yy = args
    try:
            #if i == 0:
                ## Since the the first point has infinate persistence.
            #    pass
            #else:
            bidx = np.argmin(np.abs(img - dgm[idx,0]))
            didxs = np.argmin(np.abs(img - dgm[idx,1]))
            img_temp = np.ones_like(img)
            
            death = dgm[idx,1]
            if death > -threshold:
                if dgm[idx,0] < - 5*detection:
                    death = - 10*detection
                else:
                    death = -threshold
    
            elif dgm[idx,0] < - 5*detection:
                death = dgm[idx,1] - 10*detection
            
            else:
                death = dgm[idx,1] - threshold

            
            img_temp[np.logical_or(img < dgm[idx, 0], death <= img)] = 0     # setting the values that are outside of brith and death to 0.


        
            regions = create_regions.create_regions(img_temp,-img)
            for prop in regions[1]:
                props = regions[1]
                minr, minc, maxr, maxc = prop.bbox
                if Xx[i] >= minc and Xx[i] <= maxc and Yy[i] >= minr and Yy[i] <= maxr:
                    # The point (x,y) is inside the bounding box of this region, so return the region label
                    region_label_index = prop.label - 1 # labels begin at 1 rather than 0!!
            
            bounding_box = props[region_label_index].bbox
            image_mask = img_temp[bounding_box[0]:bounding_box[2],bounding_box[1]:bounding_box[3]]

            region_data = props[region_label_index]

            # get attributes from skimage regionprops object
            # speeds up moving data between processes.

            attributes = {
                'label': region_data.label,
                'area': region_data.area,
                'bbox': region_data.bbox,
                'centroid': region_data.centroid,
                'convex_area': region_data.convex_area,
                'eccentricity': region_data.eccentricity,
                'equivalent_diameter': region_data.equivalent_diameter,
                'euler_number': region_data.euler_number,
                'extent': region_data.extent,
                'filled_area': region_data.filled_area,
                'major_axis_length': region_data.major_axis_length,
                'minor_axis_length': region_data.minor_axis_length,
                'moments': region_data.moments,
                'perimeter': region_data.perimeter,
                'solidity': region_data.solidity,
                'orientation': region_data.orientation,
                'max_intensity':region_data.max_intensity,
                'mask': image_mask,
                # Add more attributes as needed
            }  
            
            return attributes
    
    except (UnboundLocalError,IndexError):
        # likely only a single point
        # NEED TO FIX THIS as its not good practice to use exceptions for flow control.
        # lets return the point as a region
        attributes = {
                    'label': 0,
                    'area': 1,
                    'bbox': (Xx[i]-1,Yy[i]-1,Xx[i]+1,Yy[i]+1),
                    'centroid': (Xx[i],Yy[i]),
                    'convex_area': 1,
                    'eccentricity': 1,
                    'equivalent_diameter': 1,
                    'euler_number': 1,
                    'extent': 1,
                    'filled_area': 1,
                    'major_axis_length': 1,
                    'minor_axis_length': 1,
                    'moments': 1,
                    'perimeter': 1,
                    'solidity': 1,
                    'orientation': 1,
                    'max_intensity':img[Xx[i],Yy[i]],
                    'mask':'single point',
                    # Add more attributes as needed
                }
        return attributes



def Persistent_point_to_mask_SNR(dgm,img,sigma,detection,num_processes=4):
    """
    Calculate Persistent Regions from a persistence diagram, with multiprocessing.

    Parameters
    ----------
    dgm : ndarray
        Persistence diagram.
    img : ndarray
        Image.
    sigma : float
        Standard deviation of the noise.
    detection : float
        Detection threshold.
    num_processes : int, optional
        Number of processes to use. The default is 4.
        Increasing this does not necessarily increase the speed of the calculation, due to the overhead of multiprocessing.
        This is especially true for small images.

    Returns
    -------
    data : list
        List of dictionaries containing the region properties for each possible source.
    
    """


    threshold = sigma * detection 
    #print('Threshold: ',threshold)
    
    # thresholds detection of sources. 
    idxs = np.where(dgm[:,0] < -threshold)[0]
    im = - np.copy(img)
    #print('Number of points: ',len(idxs))
    X, Y = np.meshgrid(np.arange(im.shape[1]), np.arange(im.shape[0]))
    X = X.flatten()
    Y = Y.flatten()
    Xx = []
    Yy = []
    for idx in idxs:
        bidx = np.argmin(np.abs(im + dgm[idx, 0]))
        Xx.append(X[bidx])
        Yy.append(Y[bidx])
    #print(Xx,Yy)

    del (im)
    
    bx = []
    by = []
    data = []
    #for i, idx in tqdm(enumerate(idxs),total=len(idxs),desc='Calculating Persistent Mask & Regions'):
    with Pool(processes=num_processes) as p:
        with tqdm(total=len(idxs),desc='Calculating Persistent Mask & Regions') as pbar:
            for _ in p.imap_unordered(_iterator,[(i,idx,dgm,img,threshold,detection,Xx,Yy) for i, idx in enumerate(idxs)],chunksize=50):
                pbar.update()
                if _ is not None:
                    data.append(_)

    
    return  data