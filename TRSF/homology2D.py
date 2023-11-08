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





def get_enclosing_mask(x, y, mask):
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



def calculate_area(row, img):

    mask = get_mask(row,img)
    area = np.sum(mask)
    return area




def get_mask(row, img):
    mask = np.zeros(img.shape)
    mask = np.logical_or(mask,np.logical_and(img <= row.Birth,img > row.Death))
    mask_enclosed = get_enclosing_mask(int(row.y1),int(row.x1),mask)
    return mask_enclosed




def compute_ph_components(img,local_bg,analysis_threshold_val,lifetime_limit=0,output=True,bg_map=False,area_limit=3):
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

        area_list = []
        for i in tqdm(range(0,len(pd)),total=len(pd),desc='Calculating area',disable=not output):
            area = calculate_area(pd.iloc[i],img)
            #print(area)
            area_list.append(area)
        print(len(pd))
        pd['area'] = area_list
        pd = pd[pd['area'] > area_limit]     # remove 1 pixel points
        
        print(len(pd))
        enclosed_i_list = []
        for index, row in tqdm(pd.iterrows(),total=len(pd),desc='Making point assoc',disable=not output):
            enclosed_i = make_point_enclosure_assoc(row,pd,img)
            enclosed_i_list.append(enclosed_i)
        pd['enclosed_i'] = enclosed_i_list

        # DROP ROW WITH AREA OF 1 if enclosed_i is not empty

        # DROP ROW WITH AREA OF 1
        
        pd = correct_first_destruction(pd,output=not output) 
        
        pd['parent_tag'] = pd.apply(lambda row: parent_tag_func(row,pd), axis=1)
        pd['Class'] = pd.apply(classify_single,axis=1)
        
    
        return pd
    else:
        return pd