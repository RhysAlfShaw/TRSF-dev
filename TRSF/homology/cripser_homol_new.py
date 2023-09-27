'''
author: Rhys Shaw
email: rhys.shaw@bristol.ac.uk
date: 16-09-2023
description: Alternative version of cripser homol processing.
'''

import cripser
import numpy as np
import pandas


def compute_ph_cripser(img,local_bg,sigma,maxdim=0,lifetime_limit=0,remove_1pixel=False,lower_cut_threshold=3):

    pd = cripser.computePH(-img,maxdim=0)
    pd = pandas.DataFrame(pd,columns=['dim','Birth','Death','x1','y1','z1','x2','y2','z2'],index=range(1,len(pd)+1))
    pd.drop(columns=['dim','z1','z2'],inplace=True)
    pd['lifetime'] = pd['Death'] - pd['Birth']
    pd['Birth'] = -pd['Birth'] # negative for max filtration
    pd['Death'] = -pd['Death'] # negative for max filtration
    pd['Death'] = np.where(pd['Death'] < local_bg*lower_cut_threshold, local_bg*lower_cut_threshold, pd['Death'])
    pd['lifetime'] = abs(pd['Death'] - pd['Birth'])
    
    pd = pd[pd['Birth']>local_bg*sigma]

    if lifetime_limit > 0:
        pd = pd[pd['lifetime'] > lifetime_limit]

    #pd = apply_confusion_limit(pd,confusion_limit=3)

    pd['dist'] = pd.apply(lambda row: remove_1pixel_points(row),axis=1)
    
    if remove_1pixel:
        pd = pd[pd['dist'] > 3]

    pd.sort_values(by='Birth',ascending=False,inplace=True,ignore_index=True)
    if len(pd) > 0:
        pd['enclosed_i'] = pd.apply(lambda row: make_point_enclosure_assoc(row,img,pd),axis=1)
        
        pd = correct_first_destruction(pd)

        #pd.drop(columns=['dist','enclosed_i'],inplace=True)
        
        pd['parent_tag'] = pd.apply(lambda row: parent_tag_func(row,pd), axis=1)
        pd['single'] = pd.apply(classify_single,axis=1)
        return pd
    else:
        return pd


def parent_tag_func(row,pd):
    # search for yourself in anouther points enclosed_i
    for i in range(len(pd)):
        if row.name in pd.iloc[i].enclosed_i:
            return pd.iloc[i].name
    return np.nan


def classify_single(row):
    if row.new_row == 0:
        if len(row.enclosed_i) == 0:
            if np.isnan(row.parent_tag):
                return 0 # no children, no parent.
            else:
                return 1 # no child has parent.
        else:
            return 2 # has children.
    else:
        return 3 # new row has children.
        

def apply_confusion_limit(pd,confusion_limit=3):
    '''
    Will remove all points that are within the confusion limit of each other. This is to ensure that later
    multiple gaussians are not fitted to the same source.

    Parameters
    ----------
    pd : pandas dataframe
        The persistence diagram.
    confusion_limit : float, optional
        The confusion limit. The default is 3.

    Returns
    -------
    pd : pandas dataframe
        The persistence diagram with points removed.

    Examples
    --------
    >>> from cripser_homol import apply_confusion_limit
    >>> pd = apply_confusion_limit(pd,confusion_limit=3)

    '''

    indexes_to_drop = []
    seperation = confusion_limit
    for k in range(0,len(pd)):
        point = pd.iloc[k]
        for i in range(0,len(pd)):
            if i == k:
                continue

            dist = np.sqrt((point.x1 - pd.iloc[i].x1)**2 + (point.y1 - pd.iloc[i].y1)**2)
            if dist <= seperation:
                if pd.iloc[i].lifetime < pd.iloc[k].lifetime:
                    indexes_to_drop.append(pd.index[k])
                else:
                    indexes_to_drop.append(pd.index[i])

    pd.drop(index=indexes_to_drop,inplace=True)

    return pd

def correct_first_destruction(pd):
    pd['new_row'] = 0
    for i in range(0,len(pd)):
        row = pd.iloc[i]
        enlosed_i = row['enclosed_i']
        # if enlosed_i is empty then skip
        if len(enlosed_i) > 1:
            # create a new row copied from the current row
            new_row = row.copy()
            new_row['Death'] = pd.iloc[enlosed_i[0]]['Death']
            
            # add new row to the dataframe
            pd = pd.append(new_row, ignore_index=True)

    return pd

def make_point_enclosure_assoc(row,img,pd):
    '''
    Returns a list of the indices of the points that are enclosed by the mask pd point.
    '''
    point = row
    mask = np.zeros(img.shape)
    mask = np.logical_or(mask,np.logical_and(img <= point.Birth,img > point.Death))
    #plt.imshow(mask)
    #plt.scatter(point.x1,point.y1)
    #plt.show()
    mask_enclosed = get_enclosing_mask(int(point.y1),int(point.x1),mask)
    
    # check if any brith points from other points are inside the mask
    
    encloses = []
    for i in range(len(pd)):
        point = pd.loc[i]
        if point.name == row.name:
            continue
        if mask_enclosed[int(point.x1),int(point.y1)]:
            # add column to pd    
            encloses.append(point.name)

    return encloses


def remove_1pixel_points(row):
    # create mask for each point
    dist = np.sqrt((row.x1-row.x2)**2 + (row.y1-row.y2)**2)
    return dist


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
