from TRSF import sf
# take system arguments
import argparse

from astropy.table import Table


parser = argparse.ArgumentParser(description='Process Lockman facet.')

parser.add_argument('--facet', type=str, default='0',
                    help='facet number')

parser.add_argument('--data_path', type=str, default='/data/lockman/',
                    help='data path')

parser.add_argument('--out_path', type=str, default='/data/lockman/',
                    help='output path')

args = parser.parse_args()

def main():
    findmysources = sf(image=None, image_PATH=args.data_path,
                   pb_PATH=None, mode='Radio',
                   cutup=True, cutup_size=1000,
                   output=False)

    findmysources.set_background(detection_threshold=5)
    findmysources.phsf()
    findmysources.source_characterising()
    findmysources.create_polygons_fast()
    findmysources.save_polygons_to_ds9(args.out_path + 'facet_' + args.facet + '_polygons.reg')
    catalogue = findmysources.catalogue
    catalogue_astropy_table = Table.from_pandas(catalogue)
    catalogue_astropy_table.write(args.out_path + 'facet_' + args.facet + '_catalogue.fits', overwrite=True)

if __name__ == "__main__":
    main()
    #python Lockman_facet_process.py --facet 01 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_01_beamcorrected_kernel25_mgain0p5_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/