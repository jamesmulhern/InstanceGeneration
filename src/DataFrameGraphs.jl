using DataFrames, CSV
using MetaGraphsNext
import Graphs.ne, Graphs.nv

struct DataFrameGraph
    nodes::DataFrame
    edges::DataFrame
end


# Build Graph Function from metagraph
"""
    DataFrameGraph(metagraph, node_func, edge_func)

Create a DataFrameGraph based on a provided MetaGraph.  The functions node_func and edge_func must map the
node and edge data from the MetaGraph to Dicts to be used for DataFrame construction. The function must be of 
the following types

    node_func(::MetaGraph.label_type, ::MetaGraph.vertex_data_type)::Dict
    edge_func(Tuple(::MetaGraph.label_type, ::MetaGraph.label_type), ::MetaGraph.edge_data_type)::Dict

Other configurations will cause errors or undefined behavior
"""
function DataFrameGraph(g::MetaGraph, node_func::Function, edge_func::Function)

    nodes = DataFrame()
    edges = DataFrame()

    for label in labels(g)
        push!(nodes, node_func(label, g[label]), cols=:union)
    end
    
    for e_label in edge_labels(g)
        push!(edges, edge_func(e_label, g[e_label[1], e_label[2]]), cols=:union)
    end

    return DataFrameGraph(nodes, edges)
end

#
# Load and Save
#

# Save Graph
function save_dataframegraph(file_name::String, g::DataFrameGraph)
    node_f = file_name * ".nodes.csv"
    edge_f = file_name * ".edges.csv"
    CSV.write(node_f, g.nodes)
    CSV.write(edge_f, g.edges)
    return [node_f, edge_f]
end

# Load Graph
function load_dataframegraph(file_name::String)
    # Check that files exist
    node_f = file_name * ".nodes.csv"
    edge_f = file_name * ".edges.csv"

    isfile(node_f) ? nothing : error(ArgumentError, "Input File $(node_f) does not exist")
    isfile(edge_f) ? nothing : error(ArgumentError, "Input File $(edge_f) does not exist")

    # Load
    nodes = DataFrame(CSV.File(node_f))
    edges = DataFrame(CSV.File(edge_f))

    # Create DataFramesGraph and return
    return DataFrameGraph(nodes, edges)
end

function delete_dataframegraph(f::String)
    rm(f*".nodes.csv")
    rm(f*".edges.csv")
end

#
# Methods
#

nv(g::DataFrameGraph) = nrow(g.nodes)
ne(g::DataFrameGraph) = nrow(g.edges)
