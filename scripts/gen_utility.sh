#!./bin/bash
# Runs all of the utility computation

# Compute the utility datasets for the Scenarios
echo "Generating utility" 
for i in Scenario_1 Scenario_2 Scenario_3 All_Catagories
do
    file="${i}.utility.yaml"
    echo "Running file ${file}"
    julia -q bin/generate_utility.jl config_files/${file}
done