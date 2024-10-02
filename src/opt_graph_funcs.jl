"""
    floyd_warshall_iterations!(dist_matrix)

Zeros the main diagal and exicutes the floyd_warshall APAP iterations in place. 
"""
function floyd_warshall_iterations!(dist::Matrix{Tf}, n::Ti) where {Ti<:Integer, Tf<:AbstractFloat}
    # Zero diag
    for i in 1:n
        dist[i,i] = zero(Tf)
    end

    # Main loops
    for k in 1:n
        for i in 1:n
            for j in 1:n
                if dist[j,i] > dist[k,i] + dist[j,k]
                    dist[j,i] = dist[k,i] + dist[j,k]
                end
            end
        end
    end
end

function dist(A::NTuple{2,Float64},B::NTuple{2,Float64})
    return sqrt((B[1]-A[1])^2 + (B[2]-A[2])^2)
end


# Labels all road vertexs as bus stop vertexs
function add_bus_stops!(data::Dict{Symbol,Any}, params::ParamDict)
    rg = data[Symbol(params[:road_key])]
    rg.nodes[!, :type] = fill("stop", nrow(rg.nodes))

    attempt_save(params, rg)
end

function add_bus_stops_spilt!(data::Dict{Symbol,Any}, params::ParamDict)
    rg = data[Symbol(params[:road_key])]
    rg.nodes[!, :type] = fill("intersection", nrow(rg.nodes))
    
    n = ne(rg)

    for i in 1:n
        src = rg.edges.src[i]
        dst = rg.edges.dst[i]

        mid = Dict{Symbol, Any}(:index=>nrow(rg.nodes) + 1)
        mid[:x] = (rg.nodes.x[src] + rg.nodes.x[dst])/2
        mid[:y] = (rg.nodes.y[src] + rg.nodes.y[dst])/2
        mid[:lat] = (rg.nodes.lat[src] + rg.nodes.lat[dst])/2
        mid[:lon] = (rg.nodes.lon[src] + rg.nodes.lon[dst])/2
        mid[:type] = "stop"

        push!(rg.nodes, mid, cols=:subset)

        e1 = Dict{Symbol, Any}(:index=>nrow(rg.edges) + 1)
        e1[:src] = src
        e1[:dst] = mid[:index]
        e1[:weight] = rg.edges.weight[i] / 2
        e1[:length] = rg.edges.length[i] / 2

        e2 = Dict{Symbol, Any}(:index=>nrow(rg.edges) + 2)
        e2[:src] = mid[:index]
        e2[:dst] = dst
        e2[:weight] = rg.edges.weight[i] / 2
        e2[:length] = rg.edges.length[i] / 2

        push!(rg.edges, e1)
        push!(rg.edges, e2)
    end
    
    attempt_save(params, rg)
    return nothing
end


# Uses trapizodal speed profile to generate fully connected transit graph
function generate_opt_graph_trap(data::Dict{Symbol, Any}, params::ParamDict)

    a_acc = params[:acc_rate] 
    a_dec = params[:dec_rate]
    v_max = params[:max_speed] #11.176 # 25mph, MA max speed limit for "thickly settled" (houses within 200ft)

    rg = data[Symbol(params[:road_key])]

    # Build blank dist matrix and time matrix
    n = nrow(rg.nodes)
    dist_mat = Array{Float64, 2}(undef, (n,n))
    fill!(dist_mat, Inf)
    t_mat = zeros(size(dist_mat))

    # Fill with edge data from road graph
    for row in eachrow(rg.edges)
        src = row.src
        dst = row.dst
        dist_mat[dst,src] = row.length
        dist_mat[src,dst] = row.length
    end
    
    # Compute shortest paths
    floyd_warshall_iterations!(dist_mat, n)
    

    # Iterate over indexes and compute the travel time based on trapizod accel profile
    dist_limit = (v_max^2/2) * (1/a_acc - 1/a_dec)
    for idx in eachindex(dist_mat)
        d = dist_mat[idx]
        if d > dist_limit
            t_mat[idx] = d/v_max + v_max/2/a_acc - v_max/2/a_dec
        else
            t_mat[idx] = sqrt(2*d*(1/a_acc - 1/a_dec))
        end
    end

    # Build dataframe with resulting edges
    # Limit to only nodes with type "stop" and edges between them
    nodes = rg.nodes[rg.nodes.type .== "stop", :]
    edges = DataFrame(index=Vector{Int}(), src=Vector{Int}(), dst=Vector{Int}(), length=Vector{Float64}(), weight=Vector{Float64}())
    for i in nodes.index
        for j in nodes.index
            i == j ? continue : nothing
            push!(edges, (nrow(edges) + 1, i, j, dist_mat[j,i], t_mat[j,i]))
        end
    end

    # Reset indexes
    d = Dict(rg.nodes.index .=> collect(1:nrow(rg.nodes)))
    rg.edges.src = getindex.([d], rg.edges.src)
    rg.edges.dst = getindex.([d], rg.edges.dst)
    rg.nodes.index = collect(1:nrow(rg.nodes))

    # Build opt Graph
    og = DataFrameGraph(nodes, edges)

    #haskey(params, :filename) ? save_dataframegraph(params[:filename], og) : nothing
    attempt_save(params, og)
    return og
end


# Creates a fully connected optimization graph using the all vertexs and the edge weights
function generate_opt_graph(data::Dict{Symbol, Any}, params::ParamDict)

    dg = data[Symbol(params[:road_key])]
    edges2 = copy(dg.edges)
    edges2.src = dg.edges.dst
    edges2.dst = dg.edges.src
    edges2.index .+= nrow(dg.edges)
    og = DataFrameGraph(dg.nodes, vcat(dg.edges,edges2))

    #haskey(params, :filename) ? save_dataframegraph(og, params[:filename]) : nothing
    attempt_save(params, og)
    return og
end


# Links the optimizaiton graph and the baseline graph together
function link_graphs!(data::Dict{Symbol, Any}, params::ParamDict)

    og = data[Symbol(params[:opt_key])]
    bg = data[Symbol(params[:baseline_key])]

    og.nodes[!, :link_idx] = zeros(Int64, nrow(og.nodes))
    og.nodes[!, :offset] = zeros(Float64, nrow(og.nodes))

    walk_nodes = bg.nodes.type .== "walk"
    kdtree = KDTree(vcat(bg.nodes.x[walk_nodes]', bg.nodes.y[walk_nodes]'); leafsize = 10, reorder=true)
    idxs, dist = nn(kdtree, vcat(og.nodes.x', og.nodes.y'))
    og.nodes.link_idx = idxs
    og.nodes.offset = dist ./ params[:speed]


    attempt_save(params, og)
end