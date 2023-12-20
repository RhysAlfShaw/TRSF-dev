from DRUID import sf
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.table import Table
import argparse
parser = argparse.ArgumentParser()

#parser.add_argument('--sn', type=int, default=1, help='simulation number.')
#parser.add_argument('--save_path', type=str, default='.', help='path to save the simulation to.')
#parser.add_argument('--detect_thresh', type=float, default=5, help='detection threshold for source finding.')
#parser.add_argument('--analysis_thresh', type=float, default=2, help='analysis threshold for source finding.')

# create the args object
args = parser.parse_args()
# load the simulation
findmysource = sf(image=None, image_PATH=str(args.save_path)+'Simulation_Image_radio'+str(args.sn)+'.fits', mode='Radio')
findmysource.set_background(detection_threshold=args.detect_thresh,
                            analysis_threshold=args.analysis_thresh)

findmysource.phsf()
findmysource.source_characterising()

catalogue = findmysource.catalogue
image = findmysource.image

catalogue = Table.from_pandas(catalogue)
catalogue.write(str(args.save_path)+'DRUIDsimulation_catalogue'+str(args.sn)+'.fits',overwrite=True)
