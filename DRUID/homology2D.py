'''
file: homology2D.py
author: Rhys Shaw
date: 17-10-2023
description: This file contains the main class for the new source finder.
'''

import cripser
import numpy as np
import pandas
from tqdm import tqdm
import matplotlib.pyplot as plt
from multiprocess import Pool
import time
import setproctitle
import os
from collections import deque
from scipy.ndimage import label

try:
    import cupy as cp
    from cupyx.scipy.ndimage import label as cupy_label

    # Get the list of available GPU devices
    available_gpus = cp.cuda.runtime.getDeviceCount()

    if available_gpus > 0:
        print(f"Number of GPUs available: {available_gpus}, DRUID with GPU Acceleration Avalible.")
        GPU_AVAILABLE = True
    else:
        
#        from scipy.ndimage import label
        
        print('No GPU available. DRUID GPU acceleration will not be available.')
        GPU_AVAILABLE = False
    
except:
    
 #   from scipy.ndimage import label
    print('Could not import cupy. DRUID GPU acceleration will not be available.')
    GPU_AVAILABLE = False
    



def parent_tag_func_vectorized(df):
    # Create a dictionary to map each enclosed_i value to its corresponding index
    enclosed_i_dict = {idx: set(row['enclosed_i']) for idx, row in df.iterrows()}

    def find_parent_tag(row):
        for idx, enclosed_i_set in enclosed_i_dict.items():
            if row.name in enclosed_i_set:
                return idx
        return np.nan

    return df.apply(find_parent_tag, axis=1)




def parent_tag_func(row,pd):
    # search for yourself in anouther points enclosed_i
    for i in range(len(pd)):
        if row.name in pd.iloc[i].enclosed_i:
            return pd.iloc[i].name
    return np.nan



def classify_single(row):
    if row.new_row == 0:
        if len(row.enclosed_i) == 0: # no children
            if np.isnan(row.parent_tag): # no parent 
                return 0 # no children, no parent.
            else:
                return 1 # no child has parent.
        else:
            if np.isnan(row.parent_tag):
                return 2 # has children, no parent.
            else:
                return 3 # has children, has parent.
    else:
        return 4 # new row has children.
        



def correct_first_destruction(pd,output):

    pd['new_row'] = 0

    for i in tqdm(range(0,len(pd)),total=len(pd),desc='Correcting first destruction',disable=output):

        row = pd.iloc[i]
        
        enlosed_i = row['enclosed_i']
        
        if len(enlosed_i) >=1: # careful this lead to a bug make sure its >= 1.
            new_row = row.copy()
            new_row['Death'] = pd.loc[enlosed_i[0]]['Death']
            new_row['new_row'] = 1
            new_row.name = len(pd)+i
            pd = pandas.concat((pd,new_row.to_frame().T), ignore_index=False)
        
    return pd


def make_point_enclosure_assoc(row,pd,img):
    '''
    Returns a list of the indices of the points that are enclosed by the mask pd point.
    # Evaluate if a GPU is available and use it if it is.
    '''
    
    mask = get_mask(row,img)
    # check if any brith points from other points are inside the mask
    
    encloses = []
    for i in range(len(pd)):
        point = pd.iloc[i]
        if point.name == row.name:
            continue
        if mask[int(point.x1),int(point.y1)]:
            # add column to pd    
            encloses.append(point.name)

    return encloses


def make_point_enclosure_assoc_GPU(Birth,Death,row,pd,img,img_gpu):
    '''
    Returns a list of the indices of the points that are enclosed by the mask pd point.
    # Evaluate if a GPU is available and use it if it is.
    '''
    #img = img.get()
    #t0_gpu_mask = time.time()
    mask = get_mask_gpu(Birth,Death,row,img_gpu)
    #t1_gpu_mask = time.time()
    #print('GPU mask time: ',t1_gpu_mask-t0_gpu_mask)
    # check if any brith points from other points are inside the mask
    #t0_vectorized = time.time()
    mask_coords = np.column_stack((pd['x1'], pd['y1']))
    points_inside_mask = mask[mask_coords[:, 0].astype(int), mask_coords[:, 1].astype(int)]
    encloses_vectorized = pd.iloc[points_inside_mask].index.tolist()
    # remove self from list
    encloses_vectorized.remove(row.name)
    #t1_vectorized = time.time()

    #print('loop time: ',t1_vectorized-t0_vectorized)
    
    return encloses_vectorized



def get_enclosing_mask_gpu(x, y, mask):
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
            return cp.asnumpy(component_mask)
        else:
            # The specified pixel is not part of any connected component
            return None
    else:
        # The specified pixel is outside the mask
        return None




def get_enclosing_mask_new(x, y, mask):
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




def get_enclosing_mask(x, y, mask):
    '''
    Returns the mask of the enclosed area of the point (x,y) in the mask.
    '''
    
    if not mask[y, x]:
        return None
    enclosed_mask = np.copy(mask)
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




def get_enclosing_mask_old(x, y, mask):
    
    '''
    Returns the mask of the enclosed area of the point (x,y) in the mask.
    '''
    
    if not mask[y, x]:
        return None
    
    enclosed_mask = np.copy(mask)
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




def calculate_area(row, img):

    mask = get_mask(row,img)    
    area = np.sum(mask)
    return area


def calculate_area_GPU(Birth,Death,row, img_gpu,img):
    
    mask = get_mask_gpu(Birth,Death,row,img_gpu)
    area = np.sum(mask)
    return area


def get_mask(row, img):
    
    """
    
    get mask for a single row
    
    """
    
    mask = np.zeros(img.shape)
    mask = np.logical_or(mask,np.logical_and(img <= row.Birth,img > row.Death))
    mask_enclosed = get_enclosing_mask_new(int(row.y1),int(row.x1),mask)
    
    return mask_enclosed


def get_mask_gpu(Birth,Death,row,img):
    
    # same as get_mask but uses gpu with the cupy library
    t0 = time.time()
    mask = cp.zeros(img.shape,dtype=cp.float64)
    mask = cp.logical_or(mask,cp.logical_and(img <= Birth,img > Death))
    mask_enclosed = get_enclosing_mask_gpu(int(row.y1),int(row.x1),mask)
    t1 = time.time()
    #print('GPU mask time: ',t1-t0)
    return mask_enclosed    



def compute_ph_components(img,local_bg,analysis_threshold_val,lifetime_limit=0,output=True,bg_map=False,area_limit=3,nproc=1,GPU=False):
    
    
    global GPU_Option 
    GPU_Option = GPU
    pd = cripser.computePH(-img,maxdim=0)
    pd = pandas.DataFrame(pd,columns=['dim','Birth','Death','x1','y1','z1','x2','y2','z2'],index=range(1,len(pd)+1))
    pd.drop(columns=['dim','z1','z2'],inplace=True)
    pd['lifetime'] = pd['Death'] - pd['Birth']
    pd['Birth'] = -pd['Birth'] 
    pd['Death'] = -pd['Death'] 
    
    
    if bg_map:
        
        list_of_index_to_drop = []
    
        #print(pd)
        #print(local_bg)
        
        #plt.imshow(local_bg)
        #plt.show()
        
        for index, row in pd.iterrows():
            # check if local_bg is a map or a value
            if row['Birth'] < local_bg[int(row.x1),int(row.y1)]:
                list_of_index_to_drop.append(index)
        
        pd.drop(list_of_index_to_drop,inplace=True)
    
        
        # for each row evaluate if death is below analysis thresholdval map value at its birth point. if its below then set Death to bg map value.    
        for index, row in pd.iterrows():
            Analy_val = analysis_threshold_val[int(row.y1),int(row.x1)]
            if row['Death'] < Analy_val:
                row['Death'] = Analy_val
    
            
    else:
        
        pd = pd[pd['Birth']>local_bg] # maybe this should be at the beginning.
        pd['Death'] = np.where(pd['Death'] < analysis_threshold_val, analysis_threshold_val, pd['Death'])
       
    pd['lifetime'] = abs(pd['Death'] - pd['Birth'])
        

    print('Persis Diagram computed. Length: ',len(pd))

    if lifetime_limit > 0:
        pd = pd[pd['lifetime'] > lifetime_limit]

    
    pd.sort_values(by='Birth',ascending=False,inplace=True,ignore_index=True)
    
    if len(pd) > 0:
        
        def process_area(i):
            # handles worker function
            return calculate_area(pd.iloc[i], img)
        
        def process_assoc(i):
            # hanldes worker function
            return make_point_enclosure_assoc(pd.iloc[i], pd, img)

        setproctitle.setproctitle('DRUID')

        area_list = []
        
        # save the first 10 rows from the pd
        #rows = pd.iloc[0:10]
        # save rows to csv
        #rows.to_csv('test_pd_area.csv') 
        # save image as npy
        
        #np.save('are_test_image.npy',img)     
        # # parrallelize this loop
        
        if nproc == 1:
            
            if GPU_Option:
                
                if GPU_AVAILABLE:
                    
                    # convert img to cupy array and define type so it does not have to be converted each time.
                    img_gpu = cp.asarray(img,dtype=cp.float64)
                    # Calculate area and enforce area limit Single Process.
                    print('Calculating area with GPU...')
                    t0 = time.time()
                    
                    for i in range(0,len(pd)):
                        row = pd.iloc[i]
                        Birth = row.Birth
                        Death = row.Death
                        area = calculate_area_GPU(Birth,Death,row,img_gpu,img)
                        area_list.append(area)
                        percentage_completed = (i/len(pd))*100
                        
                        #if percentage_completed % 10 == 0:
                        #    print(percentage_completed,'%')
                       
                    print('Area calculated! t='+str(time.time()-t0)+' s')
    
                    pd['area'] = area_list
                    pd = pd[pd['area'] > area_limit]
                    
                    
                    enclosed_i_list = []
                    print('Calculating enclosed_i with GPU...')
                    t0 = time.time()
                    for i in range(0,len(pd)):
                        row = pd.iloc[i]
                        Birth = row.Birth
                        Death = row.Death
                        enclosed_i = make_point_enclosure_assoc_GPU(Birth,Death,row,pd,img,img_gpu)
                        enclosed_i_list.append(enclosed_i)
                        
                    print('enclosed_i calculated! t='+str(time.time()-t0)+' s')
            
                    pd['enclosed_i'] = enclosed_i_list
                                        
                    
            else:   
                
                
                # Calculate area and enforce area limit Single Process.
    
                t0 = time.time()   
                
                
                for i in tqdm(range(0,len(pd)),total=len(pd),desc='Calculating area',disable=not output):
                    area = calculate_area(pd.iloc[i],img)
                    #print(area)
                    area_list.append(area)
                    
                print('Area calculated! t='+str(time.time()-t0)+' s')
            
                pd['area'] = area_list
                pd = pd[pd['area'] > area_limit]  
                
                
                # Parent Associations Single Process
                
                enclosed_i_list = []
            
                t0 = time.time()
                for i in tqdm(range(0,len(pd)),total=len(pd),desc='Calculating enclosed_i',disable=not output):
                    row = pd.iloc[i]
                    enclosed_i = make_point_enclosure_assoc(row,pd,img)
                    enclosed_i_list.append(enclosed_i)
                
                print('enclosed_i calculated! t='+str(time.time()-t0)+' s')
                
                pd['enclosed_i'] = enclosed_i_list
                
                               
        else:
            
            print('Calculating area with ',nproc,' processes')
            
            # Calculate area and enforce area limit Multi Process.
            
            t0 = time.time()
            print("Chunksize: ",len(pd)//nproc)
            with Pool(nproc) as p:
                area_list = list(p.imap(process_area, range(len(pd)),chunksize=len(pd)//nproc))
            
            print('Area calculated! t='+str(time.time()-t0)+' s')
        
            pd['area'] = area_list
            pd = pd[pd['area'] > area_limit]     # remove 1 pixel points
                    
            print(len(pd))
            
            ## Parent Associations Multi Process.
            print('Calculating enclosed_i with ',nproc,' processes')
            t0 = time.time()
            
            print("Chunksize: ",len(pd)//nproc)
            
            with Pool(nproc) as pool:
                enclosed_i_list = list(pool.imap(process_assoc, range(len(pd)),chunksize=len(pd)//nproc))
                
            print('enclosed_i calculated! t='+str(time.time()-t0)+' s')
            
            pd['enclosed_i'] = enclosed_i_list
            

       # print(pd)
        
        pd = correct_first_destruction(pd,output=not output) 
        
        print('Calculating parent_tags... ')
        
        t0_parent_tag = time.time()
        parent_tag_list = parent_tag_func_vectorized(pd)
        pd['parent_tag'] = parent_tag_list
        t1_parent_tag = time.time()
        print('parent_tag calculated! t='+str(t1_parent_tag-t0_parent_tag)+' s')
        
        #print(pd)
        #pd['parent_tag'] = pd.apply(lambda row: parent_tag_func(row,pd), axis=1)
        
        print('Assigning Class ...')
        t0_CLass = time.time()
        pd['Class'] = pd.apply(classify_single,axis=1)
        t1_Class = time.time()
        print('Class assigned! t='+str(t1_Class-t0_CLass)+' s')
        
        return pd
    
    else:
    
        return pd
    
    