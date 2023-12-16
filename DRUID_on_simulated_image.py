from DRUID import sf
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.table import Table
#import argparse

#parser = argparse.ArgumentParser()

#parser.add_argument('--sn', type=int, default=1, help='simulation number.')
#parser.add_argument('--save_path', type=str, default='.', help='path to save the simulation to.')
#parser.add_argument('--detect_thresh', type=float, default=5, help='detection threshold for source finding.')
#parser.add_argument('--analysis_thresh', type=float, default=2, help='analysis threshold for source finding.')

#args = parser.parse_args()
# load the simulation
findmysource = sf(image=None, image_PATH='SImImage_TEST.fits', mode='Radio')
findmysource.set_background(detection_threshold=5,
                            analysis_threshold=2)

findmysource.phsf()
findmysource.source_characterising()

catalogue = findmysource.catalogue
image = findmysource.image

catalogue = Table.from_pandas(catalogue)
catalogue.write('DRUIDsimulation_catalogue_TEST.fits',overwrite=True)