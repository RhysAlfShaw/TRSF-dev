"""
author : Rhys Shaw
date : 27/06/2023
description : This file contains the functions for dealing with polygons.

"""

from skimage import measure
import numpy as np
from shapely.geometry import Polygon
from shapely.affinity import scale
from scipy.ndimage import convolve, generate_binary_structure
from shapely.ops import cascaded_union
from scipy import ndimage
from skimage import draw

def transform_mask_to_polygon(Collected_Region,img):
    point_list = []
    polygon = []
    for i in range(len(Collected_Region)):
        try:
            mask = Collected_Region[i].mask
            mask = mask.astype(int)
            # add the mask to the image
            # transform the mask to the correct size
            minr, minc, maxr, maxc = Collected_Region[i].bbox
            image_dims = np.shape(img)
            # add padding to the mask to make it the same size as the image
            mask = np.pad(mask,((minr,image_dims[0]-maxr),(minc,image_dims[1]-maxc)),'constant',constant_values=0)
            contours = measure.find_contours(mask, 0.5)
            largest_contour = np.flip(max(contours, key=len),axis=1)
            # change to ints
            largest_contour = largest_contour.astype(int)
            polygon.append(largest_contour)
        except:
            # these are all the point sources
            point_list.append(Collected_Region[i].centroid)
            pass

    return polygon, point_list


def polygons_to_mask(polygons, img_shape):
    masks = []
    # Create a mask for each polygon
    
    mask = np.zeros(img_shape, dtype=np.uint8)

    # Create a Shapely polygon from the input polygon
    shapely_polygon = polygons

    # Convert the Shapely polygon to a binary mask
    coords = np.array(shapely_polygon.exterior.coords)
    rr, cc = draw.polygon(coords[:, 1], coords[:, 0], img_shape)
    mask[rr, cc] = 1

    # Fill in the mask
    mask = ndimage.binary_fill_holes(mask).astype(int)
    masks.append(mask)

    return masks


def is_round_or_elliptical(polygon, round_threshold):
    # Scale the polygon to a circle/ellipse
    bounds = polygon.bounds
    diameter = max(bounds[2] - bounds[0], bounds[3] - bounds[1])
    scaled_polygon = scale(polygon, xfact=1.0/diameter, yfact=1.0/diameter)

    # Calculate the perimeter and area ratios
    perimeter_ratio = polygon.length / scaled_polygon.length
    area_ratio = polygon.area / scaled_polygon.area

    # Check if the ratios indicate a round or elliptical shape
    return perimeter_ratio <= round_threshold and area_ratio >= 1 - round_threshold



def split_polygons_r(polygons, round_threshold):
    round_polygons = []
    non_round_polygons = []
    for polygon in polygons:
        # Check if the polygon is round or elliptical
        if is_round_or_elliptical(polygon, round_threshold):
            round_polygons.append(polygon)
        else:
            non_round_polygons.append(polygon)
    return round_polygons, non_round_polygons


def reduce_polygons(polygons):
    areas = [polygon.area for polygon in polygons]
    largest_index = areas.index(max(areas))
    largest_polygon = polygons[largest_index]

    reduced_polygons = []
    for polygon in polygons:
        if polygon == largest_polygon or largest_polygon.contains(polygon):
            reduced_polygons.append(polygon)

    return reduced_polygons

def split_polygons(polygons):
    areas = [polygon.area for polygon in polygons]
    largest_index = areas.index(max(areas))
    largest_polygon = polygons[largest_index]

    inside_polygons = []
    outside_polygons = []
    for polygon in polygons:
        if polygon == largest_polygon or largest_polygon.contains(polygon):
            inside_polygons.append(polygon)
        else:
            outside_polygons.append(polygon)

    return inside_polygons, outside_polygons



def expand_polygon_region(component,img,bg,sigma):

    counter = 0
    convergence = False
    while convergence == False:

        counter = counter + 1
        neighborhood = generate_binary_structure(2, 2)
        dilated_mask = convolve(component, neighborhood, mode='constant', cval=0)
        inverted_mask2 = np.logical_not(component)

        # Obtain the non-overlapping parts of mask1
        non_overlapping = np.logical_and(inverted_mask2, dilated_mask)
        # differennce between two boolean arrays

        before = np.copy(non_overlapping)

        for i in range(non_overlapping.shape[0]):
            for j in range(non_overlapping.shape[1]):
                if non_overlapping[i,j] == False:
                    continue
                else:
                    temp_img = img*component
                    prop_pixel_val = img[i,j]
                    current_neightbour_vals = []
                    # left pixel 
                    current_neightbour_vals.append(temp_img[i,j-1])
                    # right pixel
                    current_neightbour_vals.append(temp_img[i,j+1])
                    # top pixel
                    current_neightbour_vals.append(temp_img[i-1,j])
                    # bottom pixel
                    current_neightbour_vals.append(temp_img[i+1,j])
                    # top left pixel
                    current_neightbour_vals.append(temp_img[i-1,j-1])
                    # top right pixel
                    current_neightbour_vals.append(temp_img[i-1,j+1])
                    # bottom left pixel
                    current_neightbour_vals.append(temp_img[i+1,j-1])
                    # bottom right pixel
                    current_neightbour_vals.append(temp_img[i+1,j+1])
            
                    del temp_img
                    # calculate the mean of the current_neightbour_vals
                    # remove zeros from the list
                    current_neightbour_vals = [x for x in current_neightbour_vals if x != 0]
                    mean_without = np.mean(current_neightbour_vals)
                    current_neightbour_vals.append(prop_pixel_val)
                    
                    if (mean_without>=prop_pixel_val) & (prop_pixel_val>=bg*sigma):
                        #print("pixel at index {},{} is kept".format(i,j))
                        #print("pixel_val = {} and mean without = {}".format(prop_pixel_val,mean_without))
                        # remove the pixel from the non_overlapping mask
                        non_overlapping[i,j] = True
                    else:
                        #print("pixel at index {},{} is removed".format(i,j))
                        #print("pixel_val = {} and mean without = {}".format(prop_pixel_val,mean_without))
                        non_overlapping[i,j] = False

        # change in the non_overlapping mask

        # Obtain the non-overlapping parts of mask1

        # checl if resid is all false

        #convergence = True
        if counter == 3:
            convergence = True
        component = np.logical_or(non_overlapping,component)
    return component
