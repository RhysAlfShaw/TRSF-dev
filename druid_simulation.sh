#!/bin/bash

# Define the functions to run
function run_create_simulated_data() {
    local script_args=("$@")
    python simulation_truth_catalogue.py "${script_args[@]}"
}

function run_source_find_in_simulated_data() {
    local script_args=("$@")
    python DRUID_on_simulated_image.py "${script_args[@]}"
}

export -f run_create_simulated_data
export -f run_source_find_in_simulated_data

# Define the parameter sets
param_sets=(
    "--N 1000 --xmax 5000 --ymax 5000 --sn 1 --seed 0 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 2 --seed 1 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 3 --seed 2 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 4 --seed 3 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 5 --seed 4 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 6 --seed 5 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 7 --seed 6 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 8 --seed 7 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 9 --seed 8 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 10 --seed 10 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 11 --seed 11 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 12 --seed 12 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 13 --seed 13 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 14 --seed 14 --save_path /data/typhon2/Rhys/data/SIMS/"
    "--N 1000 --xmax 5000 --ymax 5000 --sn 15 --seed 15 --save_path /data/typhon2/Rhys/data/SIMS/"
)

param_sets_2=(
    "--sn 1 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 2 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 3 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 4 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 5 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 6 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 7 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 8 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 9 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 10 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 11 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 12 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 13 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 14 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
    "--sn 15 --save_path /data/typhon2/Rhys/data/SIMS/ --detect_thresh 5 --analysis_thresh 2"
)

# Run the functions in parallel with the parameter sets
for param in "${param_sets[@]}"; do
    run_create_simulated_data $param &
done

wait

for param in "${param_sets_2[@]}"; do
    run_source_find_in_simulated_data $param &
done

# Wait for all background jobs to finish
wait

echo "Merging Catalogues"

python merge_sim_cats.py 

wait

echo "Job has finished. Mailing now.."
python /data/typhon2/email_script.py --rea "rhysalfshaw@gmail.com" --message "Simulated data has been created and source finding has been run."

echo "All functions have completed."