# Allows for manual modification of a dataframe graph
function modify_graph!(data::Dict{Symbol,Any}, params::ParamDict)
    g = data[Symbol(params[:graph_key])]
    lookup_key = Symbol(params[:lookup_key])

    # Add Nodes
    haskey(params, :add_nodes) ? add_nodes!(g, lookup_key, params[:add_nodes]) : nothing
    haskey(params, :rmv_nodes) ? remove_nodes!(g, lookup_key, params[:rmv_nodes]) : nothing

    # Add and Remove Edges
    haskey(params, :rmv_edges) ? remove_edges!(g, lookup_key, params[:rmv_edges]) : nothing
    haskey(params, :add_edges) ? add_edges!(g, lookup_key, params[:add_edges]) : nothing


    attempt_save(params, g)
end

function add_nodes!(g::DataFrameGraph, lookup_key::Symbol, nodeList::Vector{Dict{Symbol, Any}})
    # Get node dict meta data
    fields = propertynames(g.nodes)
    types = eltype.(eachcol(g.nodes))
    
    # Loop over nodelist
    for node_dict in nodeList

        # Create row dict
        row = Dict{Symbol,Any}()

        # Add Data to row dict
        row[:index] = nrow(g.nodes) + 1
        for (field, type) in zip(fields, types)
            
            # skip if :index
            if field == :index
                continue
            end
            
            # Check if data provided and add
            if haskey(node_dict, field)
                row[field] = node_dict[field]
            else
                row[field] = type[]
            end
        end

        # Add to node dataframe
        append!(g.nodes, row)
    end

    return nothing
end

function remove_nodes!(g::DataFrameGraph, lookup_key::Symbol, nodeList::Vector{Dict{Symbol, Any}})
    warning("Function remove_nodes! is not implamented")
end


function add_edges!(g::DataFrameGraph, lookup_key::Symbol, edgeList::Vector{Dict{Symbol, Any}})
    # Get edge dict meta data
    fields = propertynames(g.edges)
    types = eltype.(eachcol(g.edges))

    for edge_dict in edgeList

        row = Dict{Symbol, Any}()
        row[:index] = nrow(g.edges) + 1
        row[:src] = g.nodes.index[findfirst(g.nodes[!, lookup_key] .== edge_dict[:src])]
        row[:dst] = g.nodes.index[findfirst(g.nodes[!, lookup_key] .== edge_dict[:dst])]

        for (field, type) in zip(fields, types)

            #Skip some values
            if in(field, Set([:index,:src,:dst]))
                continue
            end

            # Check if data provided and add
            if haskey(edge_dict, field)
                row[field] = edge_dict[field]
            else
                row[field] = type[]
            end
        end

        push!(g.edges, row)
    end

    return nothing
end

function remove_edges!(g::DataFrameGraph, lookup_key::Symbol, edgeList::Vector{Dict{Symbol, Any}})

    for edge_dict in edgeList
        src_idx = g.nodes.index[findfirst(g.nodes[!, lookup_key] .== edge_dict[:src])]
        dst_idx = g.nodes.index[findfirst(g.nodes[!, lookup_key] .== edge_dict[:dst])]

        row_idx = (g.edges.src .== src_idx) .&& (g.edges.dst .== dst_idx)
        if !any(row_idx)
            warning("No edges were identified for removal")
        end

        deleteat!(g.edges, row_idx)
    end

    # Reset Index
    g.edges.index = collect(1:nrow(g.edges))

    return nothing
end