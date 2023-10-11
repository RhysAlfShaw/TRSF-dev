'''
from TRSF.homology.create_regions import create_regions, remove_low_flux_regions
import numpy as np
from skimage.measure import label, regionprops

def test_create_regions():
    # Create a binary image
    img = np.array([[0, 1, 0, 0],
                    [0, 1, 1, 0],
                    [0, 0, 0, 1],
                    [1, 1, 0, 0]], dtype=np.uint8)

    # Create an intensity image
    intensity_img = np.array([[10, 20, 30, 40],
                              [50, 60, 70, 80],
                              [90, 100, 110, 120],
                              [130, 140, 150, 160]], dtype=np.uint8)

    # Call the function to get the labels and regions
    labels, regions = create_regions(img, intensity_img)

    # Check if the labels and regions are of the correct types
    assert isinstance(labels, np.ndarray)
    assert isinstance(regions, list)

    # Check the shape of the labels image
    assert labels.shape == img.shape

    # Check the number of regions detected
    assert len(regions) == 2

    # Check the properties of each region
   
    assert regions[0].area == 4
   
    assert regions[0].mean_intensity == 67.5
    
    assert regions[1].area == 2
    assert regions[1].mean_intensity == 135
    
def test_remove_low_flux_regions():
    # Create a binary image
    img = np.array([[0, 1, 0, 0],
                    [0, 1, 1, 0],
                    [0, 0, 0, 1],
                    [1, 1, 0, 0]], dtype=np.uint8)

    # Create a list of regions
    class Region:
        def __init__(self, bbox):
            self.bbox = bbox

    regions = [Region((0, 0, 2, 2)), Region((1, 1, 3, 3)), Region((2, 3, 3, 4))]

    # Call the function to get the updated regions and mask
    updated_regions, mask = remove_low_flux_regions(img, regions, threshold=50)

    # Check if the updated regions and mask are of the correct types
    assert isinstance(updated_regions, np.ndarray)
    assert isinstance(mask, np.ndarray)

    # Check the number of remaining regions
    assert len(updated_regions) == 0

'''