'''
Command line interface for TRSF
'''

from .main import trsf
import sys

script_name = sys.argv[0]
arguments = sys.argv[1:]

print("Script name: {}".format(script_name))
print("Arguments: {}".format(arguments))

# arg 1 is the fits image path

img_path = arguments[0]
cutout_size = int(arguments[1])
sigma = int(arguments[2])
save_path = arguments[3]
save_type = arguments[4]

data = trsf(img_path,cutout_size,sigma)
print('Saving catalogue to {}'.format(save_path))
data.save_catalogue(save_path,save_type)