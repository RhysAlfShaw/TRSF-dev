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
    "--N 100 --xmax 500 --ymax 500 --sn 1 --seed 0 --save_path ."
    "--N 100 --xmax 500 --ymax 500 --sn 2 --seed 1 --save_path ."
)

param_sets_2=(
    "--sn 1 --save_path . --detect_thresh 5 --analysis_thresh 5"
    "--sn 2 --save_path . --detect_thresh 5 --analysis_thresh 5"
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

echo "All functions have completed."
