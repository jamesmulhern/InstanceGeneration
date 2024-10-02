#!./bin/bash

# Requires that the get_utility.sh be run first

#echo "Generating Graphs"
#julia -q bin/generate.jl config_files/Winch_Get_Graphs.inst.yaml


for file in Winch_Sc3 Winch_Case_Study # Winch_Sc1 Winch_Sc2 
do
    echo "Running configuration file - ${file}"
    julia -q bin/generate.jl config_files/${file}.inst.yaml
done