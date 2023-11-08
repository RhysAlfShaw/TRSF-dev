#!/bin/bash

# Number of CPU cores available
NUM_CORES=$(nproc)

# Define a list of argument sets to pass to the Python script
ARGUMENTS=(
    "--facet 01 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_01_beamcorrected_kernel25_mgain0p5_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 02 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_02_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 03 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_03_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 04 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_04_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 05 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_05_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 06 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_06_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 07 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_07_kernel25-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 08 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_08_kernel25-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 09 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_09_kernel25-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 10 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_10_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 11 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_11_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 12 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_12_kernel25-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 13 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_13_kernel25-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 14 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_14_kernel25-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 15 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_15_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 16 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_16_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 17 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_17_kernel25-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 18 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_18_kernel25-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 19 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_19_kernel25-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 20 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_20_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 21 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_21_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 22 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_22_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 23 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_23_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 24 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_24_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
    "--facet 25 --data_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/facet_25_beamcorrected_kernel25_mgain0p6_pd8192_nmiter10_automask5-MFS-image-pb.fits.trimmed.scaled.fits --out_path /data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"
)

# Function to run the Python script with specified arguments
run_python_script() {
    python Lockman_facet_process.py $1
}


# Export the function so that it's accessible to parallel

export -f run_python_script

# Use parallel to run the script with different argument sets in parallel
for args in "${ARGUMENTS[@]}"; do
    echo "$args" | xargs -I {} -P $NUM_CORES -n 1 bash -c 'run_python_script "$@"' _ {} &
done
wait

echo "Combining facets catalogues..."

python merge_lockman_data.py

echo "All processes have completed."