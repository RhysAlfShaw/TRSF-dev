"""
author : Rhys Shaw
email : rhys.shaw@bristol.ac.uk
date : 27-06-2023
desc :  This file contains an object for storing a detection parameter set.

"""

class CustomRegionProps:
    def __init__(self, region_data):
        self._create_region(region_data)
    
    def _create_region(self, attributes):

        self.label = attributes['label']
        self.area = attributes['area']
        self.bbox = attributes['bbox']
        self.centroid = attributes['centroid']
        self.convex_area = attributes['convex_area']
        self.eccentricity = attributes['eccentricity']
        self.equivalent_diameter = attributes['equivalent_diameter']
        self.euler_number = attributes['euler_number']
        self.extent = attributes['extent']
        self.filled_area = attributes['filled_area']
        self.major_axis_length = attributes['major_axis_length']
        self.minor_axis_length = attributes['minor_axis_length']
        self.moments = attributes['moments']
        self.perimeter = attributes['perimeter']
        self.solidity = attributes['solidity']
        self.orientation = attributes['orientation']
        self.mask = attributes['mask']
        self.max_intensity = attributes['max_intensity']
        if self.major_axis_length == 0:
            self.major_axis_length = 1
        if self.minor_axis_length == 0:
            self.minor_axis_length = 1
    
    def __iter__(self):
        return iter(self.regions)

def dict_to_regionprops(region_data):
    region_list = []
    for region_attributes in region_data:
        region_list.append(CustomRegionProps(region_attributes))
    return region_list

def list_dict_to_regionprops(region_data):
    region_list = []
    for i in range(0,len(region_data)):
        region_list.append(CustomRegionProps(region_data[i]))
    return region_list
