using DataFrames
using Proj
import Distances: Haversine
using Graphs, MetaGraphsNext

using GTFSQuerys
using GTFSGraphs

#
# User Defined Types
#

#const coord_transform = Proj.Transformation("EPSG:4326","EPSG:26986")
const LabelTuple = NamedTuple{(:stop_id, :route_id, :direction_id), NTuple{3,String}}
Base.string(t::LabelTuple) = "$(t[:stop_id])_$(t[:route_id])_$(t[:direction_id])"


#
# Vertex and Edge Structs
#
mutable struct VertexData <: AbstractVertexData
    stop_id::String
    route_id::String
    direction_id::String
    stop_name::String
    stop_coord::CoordTuple
    position::XYTuple
    headway::Float64
    dwell_time::Float64
    exit_vertex::Bool
end

function VertexData() 
    return VertexData("", "", "", "", CoordTuple((NaN,NaN)), XYTuple((NaN,NaN)), NaN, NaN, false)
end

mutable struct EdgeData <: AbstractEdgeData
    weight::Float64
    shape_distance::Float64
    running_time::Float64
    s2s_time::Float64
    edge_type::String
end

function EdgeData()
    return EdgeData(NaN, NaN, NaN, NaN, "")
end

#
# Function to apply to graph
#

function set_vertex_params!(v::VertexData; kwargs...)
    df = kwargs[:df]

    #Set stop_name
    setfield!(v, :stop_name, df[1, :stop_name])

    # Set position
    setfield!(v, :position, XYTuple(coord_transform(v.stop_coord...)))

    # Set exit vertex flag
    loc_type = df[1, :location_type]
    parent = df[1, :parent_station]

    if loc_type == "2"
        setfield!(v, :exit_vertex, true)
    elseif loc_type == "0" && ismissing(parent)
        setfield!(v, :exit_vertex, true)
    else
        setfield!(v, :exit_vertex, false)
    end

    return nothing
end

function set_weight!(e::EdgeData; kwargs...)
    setfield!(e, :weight, e.s2s_time)
end

function link_weights!(e::EdgeData; kwargs...)

    src_df = kwargs[:src_df]
    dst_df = kwargs[:dst_df]

    fwd_transfer, _ = GTFSGraphs.transfer_headway(src_df, dst_df)
    setfield!(e, :weight, fwd_transfer)
    return nothing
end


function generate_transit_graph(filename::String, params::Dict{Symbol, Any})

    # Setup
    global coord_transform = Proj.Transformation("EPSG:4326",params[:crs])
    folder = "Data/boston_GTFS"
    folder = params[:GTFS_folder]
    q = GTFSQuery(folder, 
                cache_files="all", 
                #default_filters=get_date_filters("20240321", ["06:00:00","20:00:00"])
                default_filters=get_date_filters(params[:day], [params[:start_t],params[:end_t]])
                )

    gg = GTFSGraph{LabelTuple, VertexData, EdgeData}(
                    q, 
                    shape_params=[:shape_id, :route_pattern_id], 
                    distance_function=Haversine(6371000)
                )

    # Add Routes
    add!(gg, FilterDict(:municipality => x -> .==(x,[params[:municipality]])))
    add_headway_and_dwell!(gg)
    add_running_time!(gg)
    add_s2s_time!(gg)

    apply_vertex_function!(gg, [:stop_name, :location_type, :parent_station], set_vertex_params!)
    apply_edge_function!(gg, Vector{Symbol}(), set_weight!)


    #connect_vertexs!(gg, 65.0)
    #apply_edge_function!(gg, [:arrival_time, :departure_time], link_weights!, x -> ==(x.edge_type, "transfer"))

    save_graph(filename, gg)
end