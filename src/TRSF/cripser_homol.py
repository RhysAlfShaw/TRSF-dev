'''
author: Rhys Shaw
email: rhys.shaw@bristol.ac.uk
date: 28-06-2023
description: This file contains the functions for the processing of cripser persistence diagrams.

'''


import cripser
import numpy as np
import matplotlib.pyplot as plt
import pandas




def compute_ph_cripser(img,local_bg,sigma,maxdim=0):
    '''
    Computed the persistence diagram of the image using cripser.
    Computes the cubical complex filtration of the image.
    Filters it by liftime and local background.
    Returns a pandas dataframe with the persistence diagram.

    
    Parameters
    ----------
    img : numpy array
        The image to compute the persistence diagram of.
    local_bg : float
        The local background.
    sigma : float   
        The standard deviation of the noise.
    maxdim : int, optional  
        The maximum dimension of the persistence diagram. The default is 0.

        
    Returns
    -------
    pd : pandas dataframe
        The persistence diagram.

        
    '''
    pd = cripser.computePH(img,maxdim=maxdim)
    pd = pandas.DataFrame(pd,columns=['dim','Birth','Death','x1','y1','z1','x2','y2','z2'],index=range(1,len(pd)+1))
    pd.drop(columns=['dim','z1','z2'])
    pd['lifetime'] = pd['Death'] - pd['Birth']
    pd.sort_values(by=['lifetime'],ascending=True,inplace=True)
    pd = pd[pd['lifetime'] > local_bg*sigma]

    return pd


def ph_precocessing(pd,img,local_bg,sigma):
    '''
    Basic Preprocessing of the persistence diagram. This includes:
        - Removing points with infinite death value.
        - Correcting the death value of points with children. (that cause more than one death).
        - Assigning a parent tag to each point.

    Parameters
    ----------
    pd : pandas dataframe
        the persistence diagram. (output of compute_ph_cripser)
    img : numpy array
        the image the ph was computed from. (used for mask creation during the process)
    local_bg : float
        the local background for the image above, used for death correction.
    sigma : float
        level above the noise specifed, used for death correction.    

    Returns
    -------
    pd : pandas dataframe
        the persistence diagram after the preprocessing.

    
    '''


    death_values = pd.Death.unique()
    value = death_values.max()
    del death_values

    pd = pd.apply(lambda row :alter_infinit_death(row,value,local_bg,sigma),axis=1)
    pd['encloses_i'] = pd.apply(lambda row: make_point_enclosure_assoc(row,img,pd),axis=1)
    pd['parent_tag'] = pd.apply(lambda row: assign_tag(row,pd),axis=1)
    pd['alt_Death'] = pd.apply(lambda row: death_correc(row,pd),axis=1)
    pd['alt_Death_x1'] = pd.apply(lambda row: alt_death_coord(row,'x',pd),axis=1)
    pd['alt_Death_y1'] = pd.apply(lambda row: alt_death_coord(row,'y',pd),axis=1)
    pd['len_enclosed'] = pd.apply(lambda row: len(row.encloses_i),axis=1)
    pd = create_new_row(pd)
    pd = pd.apply(alter_brith_and_death,axis=1)
    pd = pd.apply(lambda row: cap_death(row,local_bg,sigma),axis=1)
    pd.drop(columns=['alt_Death','alt_Death_x1','alt_Death_y1'],inplace=True)
    pd['encloses_i_len'] = pd.apply(lambda row: len(row.encloses_i),axis=1)
    pd = pd.dropna()
    pd['lifetime'] = pd.apply(lambda row: row.Birth - row.Death,axis=1)
    pd.drop(columns=['encloses_i_len','encloses_i'],inplace=True)
    return pd


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



def make_point_enclosure_assoc(row,img,pd):
    '''
    Returns a list of the indices of the points that are enclosed by the mask pd point.
    '''
    point = row
    mask = np.zeros(img.shape)
    mask = np.logical_or(mask,np.logical_and(img <= point.Birth,img > point.Death))
    mask_enclosed = get_enclosing_mask(int(point.y1),int(point.x1),mask)
    # check if any brith points from other points are inside the mask
    encloses = []
    for i in range(len(pd)):
        point = pd.iloc[i]
        if mask_enclosed[int(point.x1),int(point.y1)]:
            # add column to pd    
            encloses.append(point.name)
    return encloses



def death_correc(row,pd):
    '''
    Returns the corrected death value of the point.
    '''
    if len(row.encloses_i) == 1:
        return None
    else:
        # find highest birth point
        brith = []
        for i in row.encloses_i:
            found_row = pd[pd['parent_tag']==row.parent_tag].loc[i]
            if found_row.Birth == row.Birth:
                continue
            brith.append(-found_row.Birth)
        max_birth = np.max(brith)
        # identify the point with the highest birth point
        # remove self from list

        return max_birth



def alt_death_coord(row,coord,pd):
    '''
    alternative death coord to the new one.
    '''
    if row.alt_Death == None:
        return None
    if coord == 'x':
        target = row.alt_Death
        for i in row.encloses_i:
            found_row = pd[pd['parent_tag']==row.parent_tag].loc[i]
            if -found_row.Birth == target:
                return found_row.x1
    if coord =='y':
        target = row.alt_Death
        for i in row.encloses_i:
            found_row = pd[pd['parent_tag']==row.parent_tag].loc[i]
            if -found_row.Birth == target:
                return found_row.y1



def assign_tag(row,pd):
    '''
    assign a tag to the point based on which larger point it is enclosed by.
    '''
    row_index = row.name
    name_list = pd.index
    for i, list_i in enumerate(pd.encloses_i):
        list_name = name_list[i]
        # loop through list_i and check if row_index is in the list
        if row_index in list_i:
            return list_name  # return the name of the row that contains the list
    return row_index # if no match is found then return the row_index



def alter_infinit_death(row,value,local_bg,sigma):
    '''
    Correct the infinite death value.
    Recovers the brightness point.
    '''
    if row.Death == value:
        row.Death = local_bg*sigma
        row.lifetime = row.Death - row.Birth
    row.Birth = abs(row.Birth)
    row.Death = abs(row.Death)
    return row



def create_new_row(dataframe):
    '''
    Create new row altered data, to remove unesseary columns
    '''

    dataframe_copy = dataframe.copy()
    dataframe_copy = dataframe_copy.dropna(subset=['alt_Death'])
    dataframe_copy['Death'] = dataframe['alt_Death']
    dataframe_copy['x2'] = dataframe['alt_Death_x1']
    dataframe_copy['y2'] = dataframe['alt_Death_y1']
    dataframe = pandas.concat((dataframe,dataframe_copy))
    return dataframe



def alter_brith_and_death(row):
    '''
    makes the birth and death values positive. 
    Standard output from the algorithm is negative.
    '''
    row.Birth = np.abs(row.Birth)
    row.Death = np.abs(row.Death)
    return row



def cap_death(row,local_bg,sigma):
    '''
    Capt the lowest death value to the local background value * a significance threshold.
    '''
    if row.Death < local_bg*sigma:
        row.Death = local_bg*sigma
    return row    



def drop_row(row,max_len):
    '''
    Drop duplicate rows.
    '''
    if row.encloses_i_len > 2:
        if row.encloses_i_len == max_len:
            return row
        else:
            return None
    else:
        return row