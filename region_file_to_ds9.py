# read all facet_num.reg files from dir and save them as a ds9 region file.

import glob
import os
import numpy as np

from astropy.io import ascii
from astropy.table import Table
from astropy import units as u
from astropy.regions import read_ds9

PATH = '/data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/'

def main():
    # get all region files
    region_files = glob.glob(PATH + 'facet_*.reg')
    # remove any file names with polygons in them
    region_files = [f for f in region_files if 'polygons' not in f]
    print(region_files)
    # read all region files

    
    
        

if __name__ == "__main__":
    main()
    