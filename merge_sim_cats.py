import numpy as np
from astropy.table import Table, vstack
import glob
import os
import pdb

def find_files_in_directory(directory_path, file_name_structure):
 
    file_paths = glob.glob(os.path.join(directory_path, file_name_structure))

    return file_paths

PATH = '/data/typhon2/Rhys/data/SIMS/'

file_name_stucture_truth = "Simulation_Catalogue_Truth*.fits"
file_name_stucture_reco = "DRUIDsimulation_catalogue*.fits"

# in matching to prevent pile up on plane we will add a xmax*ns offset to x and y for both truth and reco
xmax = 5000
# merge Truth catalogues
found_files_truth = find_files_in_directory(PATH, file_name_stucture_truth)
# look for the combined catalogue and remove it
try:
    found_files_truth.remove(PATH+'Simulation_Catalogue_Truth_Combined.fits')
except ValueError:
    pass
catalogue_list = []
for file in found_files_truth:
    truth_catalogue = Table.read(file)
    # get the number which is the character before the .fits
    #pdb.set_trace()
    ns = int(file.split('/')[-1].split('.')[0][-1])
    truth_catalogue['x'] = truth_catalogue['x'] + ns*xmax 
    truth_catalogue['y'] = truth_catalogue['y'] + ns*xmax
    catalogue_list.append(truth_catalogue)
    
truth_catalogue = vstack(catalogue_list)

# write out truth catalogue
truth_catalogue.write(PATH+'Simulation_Catalogue_Truth_Combined.fits', overwrite=True)
    
# merge Reco catalogues
found_files_reco = find_files_in_directory(PATH, file_name_stucture_reco)
try:
    found_files_truth.remove(PATH+'DRUIDsimulation_catalogue_Combined.fits')
except ValueError:
    pass
catalogue_list = []
for file in found_files_reco:
    reco_catalogue = Table.read(file)
    # get the number which is the character before the .fits
    ns = int(file.split('/')[-1].split('.')[0][-1])
    reco_catalogue['Xc'] = reco_catalogue['Xc'] + ns*xmax 
    reco_catalogue['Yc'] = reco_catalogue['Yc'] + ns*xmax
    catalogue_list.append(reco_catalogue)
    
reco_catalogue = vstack(catalogue_list)

# write out reco catalogue
reco_catalogue.write(PATH+'DRUIDsimulation_catalogue_Combined.fits', overwrite=True)