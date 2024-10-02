using Pkg
Pkg.activate(".")

using YAML
using DataFrames, CSV
using NearestNeighbors
using ArgParse

# Constants
const ParamDict = Dict{Symbol,Any}

#Includes
include("../src/gen_transit_graph.jl")
include("../src/DataFrameGraphs.jl")
include("../src/OptInstanceTools.jl")
include("../src/compute_vertex_weights.jl")
include("../src/modify_graphs.jl")
include("../src/opt_graph_funcs.jl")
include("../src/baseline_graph_funcs.jl")


#
# Arg Parse Setup
#
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "config_file"
            help = "Configuration YAML file.  This file defines the details of the instance"
            required = true
    end

    return parse_args(s)
end

#
# Save and Load Functions
#
function load_graph(data::Dict{Symbol, Any}, params::ParamDict)
    haskey(params, :folder) ? f = joinpath(params[:folder], params[:filename]) : f = params[:filename]
    rg = load_dataframegraph(f)
    return rg
end

function attempt_save(params::ParamDict, g::DataFrameGraph)
    if haskey(params, :filename)
        haskey(params, :folder) ? f = joinpath(params[:folder],params[:filename]) : f = params[:filename]
        save_dataframegraph(f, g)
    end
    return nothing
end

function save_instance(data::Dict{Symbol,Any}, params::ParamDict)
    opt_file = params[:filename]
    bg = data[Symbol(params[:baseline_key])]
    og = data[Symbol(params[:opt_key])]

    rename!(bg.nodes, :index => :idx)
    rename!(og.nodes, :index => :idx)
    select!(bg.nodes, [:idx, :weight, :x, :y])
    select!(bg.edges, [:src, :dst, :weight])
    select!(og.nodes, [:idx, :link_idx, :offset, :x, :y])
    select!(og.edges, [:src, :dst, :weight])
    save_opt_instance(opt_file, bg, og)
end

# 
# Graph generation functions 
#

function generate_osm_graph(data::Dict{Symbol,Any}, params::ParamDict)

    # Check if filename key exists
    if haskey(params, :filename) 
        fname = params[:filename]
        haskey(params, :folder) ? fname = joinpath(params[:folder], fname) : nothing
        retain = true
    else
        fname = "osm_graph"
        retain = false
    end

    # # Set Python Enviroment
    # run(Cmd(`/home/james/miniconda3/condabin/conda activate InstGen`))

    # Build command and run /usr/local/bin/python3
    dg_cmd = Cmd(`/home/james/miniconda3/bin/python3 src/get_osm_graph.py 
                 --file $(fname) 
                 --travel_speed $(params[:speed]) 
                 --place_str $(params[:place]) 
                 --CRS $(params[:crs])
                 --resolution $(params[:resolution]) 
                 --filters $(params[:filters])`)
    run(dg_cmd)

    # Reload graph
    dg = load_dataframegraph(fname)

    # Cleanup
    retain == false ? delete_dataframegraph(fname) : nothing

    return dg
end

function generate_transit_graph(data::Dict{Symbol, Any}, params::ParamDict)
    # Check if filename key exists
    if haskey(params, :filename) 
        haskey(params, :folder) ? fname = joinpath(params[:folder], params[:filename]) : fname = params[:filename]
        retain = true
    else
        fname = "transit_graph_temp"
        retain = false
    end

    generate_transit_graph(fname, params[:generation_params])
    tg = load_dataframegraph(fname)

    # Cleanup
    retain == false ? delete_dataframegraph(fname) : nothing

    return tg
end

#
# User Functions
#
function half_headway(headway::Float64)
    return headway / 2
end


function main()
    args = parse_commandline()
    file = args["config_file"]
    params = YAML.load_file(file, dicttype=ParamDict)

    data = Dict{Symbol,Any}()
    generators = params[:generators]
    for gen in generators
        print("Running: $(gen[:name]), Method: $(gen[:method]). ")
        val = getfield(Main, Symbol(gen[:method]))(data, gen[:params])

        if haskey(gen, :key_val)
            data[Symbol(gen[:key_val])] = val
            print("Data added under key :$(gen[:key_val]). ")
        end

        print("Complete\n")
    end
    return nothing
end


main()