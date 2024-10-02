using DataFrames, CSV
using ZipFile



function save_opt_instance(filename::String, base_g::DataFrameGraph, route_g::DataFrameGraph)
    z = ZipFile.Writer(filename)
    
    # Write base graph
    f = ZipFile.addfile(z, "base.nodes.csv", method=ZipFile.Store)
    CSV.write(f, base_g.nodes)
    f = ZipFile.addfile(z, "base.edges.csv", method=ZipFile.Store)
    CSV.write(f, base_g.edges)

    # Write route graph
    f = ZipFile.addfile(z, "route.nodes.csv", method=ZipFile.Store)
    CSV.write(f, route_g.nodes)
    f = ZipFile.addfile(z, "route.edges.csv", method=ZipFile.Store)
    CSV.write(f, route_g.edges)

    close(z)
end

function read_opt_instance(filename::String)

    r = ZipFile.Reader(filename)

    names = [f.name for f in r.files]

    # Read base dataframes
    f_idx = findfirst(names .== "base.nodes.csv")
    base_nodes = CSV.read(r.files[f_idx], DataFrame)
    f_idx = findfirst(names .== "base.edges.csv")
    base_edges = CSV.read(r.files[f_idx], DataFrame)
    base  = DataFrameGraphs.DataFrameGraph(base_nodes,  base_edges)

    # Read route dataframes
    f_idx = findfirst(names .== "route.nodes.csv")
    route_nodes = CSV.read(r.files[f_idx], DataFrame)
    f_idx = findfirst(names .== "route.edges.csv")
    route_edges = CSV.read(r.files[f_idx], DataFrame)
    route = DataFrameGraphs.DataFrameGraph(route_nodes, route_edges)

    close(r)

    return base, route
end