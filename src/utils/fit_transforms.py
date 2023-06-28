"""
author : Rhys Shaw  
date : 27/06/2023
description : This file contains the functions for preprocessing the fitted gaussian catalogue.

"""
import math 

def reduce_angle(angle):
    reduced_angle = angle % (2*math.pi)  # Wrap angle within 0 to 2*pi range
    if reduced_angle > math.pi:  # Convert angle to the range -pi to pi
        reduced_angle -= 2*math.pi
    return reduced_angle

def transform_parameters(sigma_x, sigma_y, angle):
    # Transform the parameters to the range 0 to pi.
    # sigma_x and sigma_y are always positive.
    transformed_sigma_x = sigma_x.apply(lambda x: -x if x < 0 else x) 
    transformed_sigma_y = sigma_y.apply(lambda y: -y if y < 0 else y)
    transformed_angle = angle + math.pi
    return transformed_sigma_x, transformed_sigma_y, transformed_angle