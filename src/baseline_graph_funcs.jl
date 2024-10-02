
# Add the vertex weights to the walking graph based on the provided file
function add_vertex_weights_from_file!(data::Dict{Symbol, Any}, params::ParamDict)

    wg = data[Symbol(params[:data_key])]

    vertex_weight_file = params[:input_file]
    v_weights = DataFrame(CSV.File(vertex_weight_file))
    demandSum, scaledSum, localDemand = compute_vertex_weights(wg.nodes.x, wg.nodes.y, 
                                                               v_weights.x, v_weights.y, v_weights.utility)
    wg.nodes[!, :weight] = demandSum
    wg.nodes[!, :scaledWeight] = scaledSum
    wg.nodes[!, :localWeight] = localDemand

    attempt_save(params, wg)

    return nothing
end


# prunes the transit graph to only the nodes inside the region of interst, the nodes near the external 
# utility points and the links between the region and the exeranl points
function prune_transit!(data::Dict{Symbol,Any}, params::ParamDict)

    wg = data[Symbol(params[:walk_key])]
    tg = data[Symbol(params[:transit_key])]
    ext_u = DataFrame(CSV.File(params[:external_utility]))

    # Find transit nodes within 65 of a point in the walk graph
    tg.nodes.status = fill("", nv(tg))
    kdtree = KDTree(vcat(wg.nodes.x', wg.nodes.y'); leafsize = 10, reorder=true)
    _ , dist = nn(kdtree,vcat(tg.nodes.x',tg.nodes.y'))
    valid = dist .<= params[:dist_limit]
    tg.nodes.status[valid] .= "inside"

    # Find transit nodes nearest to exteranl points
    kdtree = KDTree(vcat(tg.nodes.x', tg.nodes.y'))
    idxs, dist = knn(kdtree, vcat(ext_u.x', ext_u.y'), 2)
    ext_u.tg_nodes = idxs
    ext_u.tg_dist = dist
    ext_u.near_v = zeros(Int64, nrow(ext_u))
    ext_u.near_d = fill(NaN, nrow(ext_u))

    # Mark edges that are inside region or enter region
    tg.edges.status = fill("", ne(tg))
    for row in eachrow(tg.edges)
        src = row.src
        dst = row.dst

        if tg.nodes.status[src] == "inside" && tg.nodes.status[dst] == "inside"
            row.status = "inside"
        end

        if tg.nodes.status[src] == "" && tg.nodes.status[dst] == "inside"
            row.status = "in"
        end

        if tg.nodes.status[src] == "inside" && tg.nodes.status[dst] == ""
            row.status = "out"
        end

    end

    # Find indexs of out edges
    idxs = findall(tg.edges.status .== "out")
    for idx in idxs
        if tg.edges.edge_type[idx] == "transfer"
            continue
        end

        edgeVec = Vector{Int64}()
        nodeVec = Vector{Int64}()
        checkVec = [idx]
        dir = tg.nodes.direction_id[tg.edges.src[idx]]

        while !isempty(checkVec)
            # Pop current edge and 
            edge = popfirst!(checkVec)
            src = tg.edges.src[edge]
            dst = tg.edges.dst[edge]

            # Add edge and dst to tracking vectors
            push!(edgeVec, edge)
            push!(nodeVec, dst)

            # Check if dst is in near list
            i = findfirst(vcat(ext_u.tg_nodes'...) .== dst)
            if !isnothing(i)
                ext_u.near_v[i[1]] = ext_u.tg_nodes[i[1]][i[2]]
                ext_u.near_d[i[1]] = ext_u.tg_dist[i[1]][i[2]]
                tg.nodes.status[nodeVec] .= "ext"
                tg.edges.status[edgeVec] .= "ext_link"
            end

            # Find outgoing edges and add to checkVec
            out_edges = findall(tg.edges.src .== dst)
            out_edges = out_edges[tg.nodes.direction_id[tg.edges.src[out_edges]] .== dir]
            out_edges = out_edges[.![e in edgeVec for e in out_edges]]
            push!(checkVec, out_edges...)
        end

    end
    
    # Add utility value
    tg.nodes.utility = zeros(Float64, nrow(tg.nodes))
    tg.nodes.utility[ext_u.near_v] = ext_u.util_sum

    # Remove edges that are not needed
    filter!(:status => n -> n in ["ext","inside"], tg.nodes)
    filter!(:status => n -> n in ["ext_link", "inside"], tg.edges)

    # Reset indexes
    tg.edges.index = collect(1:nrow(tg.edges))
    d = Dict( tg.nodes.index .=> collect(1:nrow(tg.nodes)))
    tg.edges.src = getindex.([d], tg.edges.src)
    tg.edges.dst = getindex.([d], tg.edges.dst)
    tg.nodes.index = collect(1:nrow(tg.nodes))

    attempt_save(params, tg)
    return nothing

end


# Merges the transit graph and walking graph into the baseline graph
function merge_graphs(data::Dict{Symbol, Any}, params::ParamDict)
    # Load graphs
    wg = data[Symbol(params[:walk_key])]
    tg = data[Symbol(params[:transit_key])]

    # Define columns
    nodes = DataFrame(index=Vector{Int}(), x=Vector{Float64}(), y=Vector{Float64}(), weight=Vector{Float64}(), type=Vector{String}())
    edges = DataFrame(index=Vector{Int}(), src=Vector{Int}(), dst=Vector{Int}(), weight=Vector{Float64}())

    w_nodes = collect(1:nv(wg))
    t_nodes = collect((nv(wg)+1):(nv(wg) + nv(tg)))
    w_map = Dict(wg.nodes.index .=> w_nodes)
    t_map = Dict(tg.nodes.index .=> t_nodes)
    t_map_rev = Dict(t_nodes .=> tg.nodes.index)

    # Add walk edges
    for n in eachrow(wg.nodes)
        push!(nodes, (w_map[n[:index]], n[:x], n[:y], n[:weight], "walk"))
    end
    for e in eachrow(wg.edges)
        l = nrow(edges)
        push!(edges, (l+1, w_map[e[:src]], w_map[e[:dst]], e[:weight]))
        push!(edges, (l+2, w_map[e[:dst]], w_map[e[:src]], e[:weight]))
    end

    # Add transit edges
    for n in eachrow(tg.nodes)
        push!(nodes, (t_map[n[:index]], n[:x], n[:y], n[:utility], "transit"))
    end
    for e in eachrow(tg.edges)
        l = nrow(edges)
        push!(edges, (l+1, t_map[e[:src]], t_map[e[:dst]], e[:weight]))
    end

    # Find nearest nodes (finds walk node nearest to transit nodes)
    kdtree = KDTree(vcat(nodes.x[w_nodes]', nodes.y[w_nodes]'); leafsize = 10, reorder=true)
    idxs, dist = nn(kdtree, vcat(nodes.x[t_nodes]', nodes.y[t_nodes]'))

    # Add linking edges
        # Only to vertexs inside the the walk network
        # add headway walk -> transit only

    walk_speed = params[:speed]

    for (i, idx) in pairs(idxs)
        t_idx = t_nodes[i]
        w_idx = w_nodes[idx]

        if tg.nodes.status[t_map_rev[t_idx]] != "inside" 
            continue
        end

        waiting_time = getfield(Main, Symbol(params[:waiting_time_func]))(tg.nodes.headway[t_map_rev[t_idx]])
        #waiting_time = tg.nodes.headway[t_map_rev[t_idx]] / 2

        l = nrow(edges)
        push!(edges, (l+1, t_idx, w_idx, dist[i] / walk_speed))                   # Transit exit edge
        push!(edges, (l+2, w_idx, t_idx, (dist[i] / walk_speed) + waiting_time))   # Enter transit edge

    end

    bg = DataFrameGraph(nodes, edges)
    attempt_save(params, bg)
    return bg
end

function convert_walk_to_baseline(data::Dict{Symbol, Any}, params::ParamDict)
    wg = data[Symbol(params[:walk_key])] 
    bg = wg

    bg.nodes.type = fill("walk", nrow(bg.nodes))
    attempt_save(params, bg)
    return bg
end
