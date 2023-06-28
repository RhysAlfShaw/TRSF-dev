"""
author: Rhys Shaw
email: rhys.shaw@bristol.ac.uk
date: 26-06-2023
description: This file contains the functions for the creating subimages of a larger image.

"""

def create_cutouts(image,size):
    # Create a list of the cutouts
    cutouts = []
    # Create a list of the coordinates of the cutouts
    coords = []
    # Iterate through the image
    for i in range(0,image.shape[0],size):
        for j in range(0,image.shape[1],size):
            # Create a cutout
            cutout = image[i:i+size,j:j+size]
            # Append the cutout to the list
            cutouts.append(cutout)
            # Append the coordinates to the list
            coords.append((i,j))
    return cutouts, coords



def tracible_create_cutouts(image, segment_size):
    """
    Cuts a large image into smaller segments with traceability.

    Args:
    - image: The large image as a numpy array.
    - segment_size: The size of each subimage segment (e.g., 256 for a 256x256 segment).

    Returns:
    - segments: A list of subimage segments.
    - index_map: A dictionary mapping the coordinates to the corresponding segment index.
    """
    height, width = image.shape
    segments = []
    index_map = {}

    for y in range(0, height, segment_size):
        for x in range(0, width, segment_size):
            segment = image[y:y+segment_size, x:x+segment_size]
            segment_index = len(segments)
            segments.append(segment)

            index_map[segment_index] = ((y, x), (y+segment_size-1, x+segment_size-1))

    return segments, index_map



def find_subimage_index_for_coord(index_map, segments, y_coord, x_coord):
    """
    Finds the index of the subimage segment that contains the given coordinate.

    Args:
    - index_map: A dictionary mapping the coordinates to the corresponding segment index.
    - segments: A list of subimage segments.
    - y_coord: The y coordinate of the pixel.
    - x_coord: The x coordinate of the pixel.

    Returns:
    - segment: The subimage segment that contains the given coordinate.
    """

    for segment_index, ((y_min, x_min), (y_max, x_max)) in index_map.items():
        if y_min <= y_coord <= y_max and x_min <= x_coord <= x_max:
            segment = segments[segment_index]
            break
    return segment