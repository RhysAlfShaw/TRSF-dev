import numpy as np
import pandas
from scipy.ndimage import convolve, generate_binary_structure
from src.source_props.region_expansion import compute as cython_compute

def region_expansion_downhill(componet, Image, min_val, method=None, max_iter=3):
    counter = 0
    convergence = False
    while convergence == False:
        
        counter = counter + 1

        # Define the neighborhood structure for dilation
        neighborhood = generate_binary_structure(2, 2)

        dilated_mask = convolve(componet, neighborhood, mode='constant', cval=0)
        inverted_mask2 = np.logical_not(componet)

        # Obtain the non-overlapping parts of mask1
        non_overlapping = np.logical_and(inverted_mask2, dilated_mask)
        # differennce between two boolean arrays
        
        before = np.copy(non_overlapping)
        
        if method == "cython":
            non_overlapping = non_overlapping.astype(np.int32)
            componet = componet.astype(np.int32)
            Image = Image.astype(np.float32)
            non_overlapping = cython_compute(non_overlapping, Image, componet, min_val)
        else:
            non_overlapping = compute(non_overlapping, Image, componet, min_val)
        if counter == max_iter:
            convergence = True
        componet = np.logical_or(non_overlapping,componet)
    
    return componet

def compute(non_overlapping, Image, componet, min_val):
    for i in range(non_overlapping.shape[0]):
        for j in range(non_overlapping.shape[1]):
            if non_overlapping[i,j] == False:
                continue
            else:
                temp_img = Image*componet
                prop_pixel_val = Image[i,j]
                current_neightbour_vals = []

                # need to add a catch for image edges
                # left pixel 
                try:
                    current_neightbour_vals.append(temp_img[i,j-1])
                except:
                    continue
                # right pixel
                try:
                    current_neightbour_vals.append(temp_img[i,j+1])
                except:
                    continue
                # top pixel
                try:
                    current_neightbour_vals.append(temp_img[i-1,j])
                except:
                    continue
                # bottom pixel
                try:
                    current_neightbour_vals.append(temp_img[i+1,j])
                except:
                    continue
                # top left pixel
                try:
                    current_neightbour_vals.append(temp_img[i-1,j-1])
                except:
                    continue
                    # top right pixel
                try:
                    current_neightbour_vals.append(temp_img[i-1,j+1])
                except:
                    continue
                # bottom left pixel
                try:
                    current_neightbour_vals.append(temp_img[i+1,j-1])
                except:
                    continue
                # bottom right pixel
                try:
                    current_neightbour_vals.append(temp_img[i+1,j+1])
                except:
                    continue
                del temp_img
                # calculate the mean of the current_neightbour_vals
                # remove zeros from the list
                current_neightbour_vals = [x for x in current_neightbour_vals if x != 0]
                mean_without = np.mean(current_neightbour_vals)
                current_neightbour_vals.append(prop_pixel_val)
                
                if (mean_without>=prop_pixel_val) & (prop_pixel_val>=min_val):

                    #print("pixel at index {},{} is kept".format(i,j))
                    #print("pixel_val = {} and mean without = {}".format(prop_pixel_val,mean_without))
                    # remove the pixel from the non_overlapping mask
                    non_overlapping[i,j] = True
                else:
                    #print("pixel at index {},{} is removed".format(i,j))
                    #print("pixel_val = {} and mean without = {}".format(prop_pixel_val,mean_without))
                    non_overlapping[i,j] = False
    return non_overlapping
                

def find_touching_neighbors_between_masks(mask1, mask2):
    """
    Finds the locations of touching neighboring pixels between two binary masks.

    Arguments:
    mask1 -- first binary mask as a NumPy array
    mask2 -- second binary mask as a NumPy array

    Returns:
    touching_neighbors -- list of coordinates [(x1, y1), (x2, y2), ...] of touching neighboring pixels
    """

    # Define the neighborhood structure for connectivity
    neighborhood = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])

    # Perform a binary erosion to find the touching neighbors
    eroded_mask1 = np.logical_and.reduce(
        [mask1[i:i + 3, j:j + 3] for i in range(3) for j in range(3)] * neighborhood
    )

    eroded_mask2 = np.logical_and.reduce(
        [mask2[i:i + 3, j:j + 3] for i in range(3) for j in range(3)] * neighborhood
    )

    # Find the coordinates of the touching neighboring pixels
    touching_neighbors = np.argwhere(np.logical_and(eroded_mask1, eroded_mask2))

    return touching_neighbors.tolist()



def expand_polygon_with_mask(image, mask, iterations=1):
    """
    Expands a polygon represented by a binary mask by considering the intensity values of local neighboring pixels.

    Arguments:
    image -- input intensity image as a NumPy array
    mask -- binary mask representing the polygon
    iterations -- number of expansion iterations (default: 1)

    Returns:
    smoothed_mask -- expanded mask as a binary NumPy array
    """

    # Define the neighborhood structure for dilation
    neighborhood = generate_binary_structure(2, 2)

    # Perform polygon expansion
    for _ in range(iterations):
        # Dilate the mask
        dilated_mask = convolve(mask, neighborhood, mode='constant', cval=0)

        # Find the pixels that will enter the expanded polygon
        entering_pixels = np.logical_and(dilated_mask > 0, mask == 0)

        # Find the neighboring pixels already inside the polygon
        neighboring_pixels = convolve(mask.astype(np.uint8), neighborhood, mode='constant', cval=0) > 0

        # Find the pixels that meet the condition for expansion
        expansion_pixels = np.logical_and(entering_pixels, image < neighboring_pixels)

        # Update the mask by adding the expansion pixels
        mask = np.logical_or(mask, expansion_pixels)

    return mask


