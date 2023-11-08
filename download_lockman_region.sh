#!/bin/bash
# list of file to download

files=(
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_01.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_02.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_03.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_04.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_05.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_06.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_07.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_08.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_09.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_10.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_11.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_12.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_13.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_14.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_15.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_16.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_17.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_18.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_19.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_20.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_21.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_22.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_23.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_24.reg"
    "https://lofar-surveys.org/public/Lockman-VLBI/facet_25.reg"
)

save_path="/data/typhon2/Rhys/data/LoFAR/Lockman_hole/facets/"

for file in "${files[@]}"
do
    wget -P $save_path $file
done
