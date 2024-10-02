using MetaGraphsNext
using Graphs
using DataFrames, CSV

include("../src/DataFrameGraph.jl")
include("../src/OptInstanceTools.jl")


function grid_layout(h_blocks, v_blocks; inside=1, buffer=1)
    b_w = 1 + inside
    start_idx = 1 + buffer

    width  = h_blocks*b_w + 1 + 2*buffer
    height = v_blocks*b_w + 1 + 2*buffer

    layout = ones(Int64, height, width)

    for i in 0:h_blocks
        idx = 1 + buffer + i*b_w
        layout[start_idx:end-buffer,idx] .= 2
    end

    for i in 0:v_blocks
        idx = 1 + buffer + i*b_w
        layout[idx,start_idx:end-buffer] .= 2
    end
    return layout
end

function generate_point_matrix(layout, spacing)
    dims = size(layout)
    x_p = collect(LinRange(0, spacing*(dims[2]-1), dims[2])) .- spacing*(dims[2]-1)/2
    X = x_p' .* ones(dims[1])
    y_p = collect(LinRange(0, spacing*(dims[1]-1), dims[1])) .- spacing*(dims[1]-1)/2
    Y = ones(dims[2])' .* y_p
    return X, Y
end

function Gaussian(A::Real, sigma::Real)::Function
    f(x,y,x_0,y_0) = A * exp(-((x-x_0)^2/(2*sigma^2) + (y-y_0)^2/(2*sigma^2)))
    return f
end

function make_adj_matrix(offsets, layout; type=1)
    dims = size(layout)
    nn = dims[1] * dims[2]
    w_adj = zeros(Int64, nn, nn)

    for i in 1:dims[1]
        for j in 1:dims[2]
            src = (j-1)*dims[1] + i
            layout[src] < type ? continue : nothing

            for p in offsets
                # Get coords of target node
                t_i = i + p[1]
                t_j = j + p[2]
                
                #Check limits
                t_i < 1 || t_i > dims[1] ? continue : nothing
                t_j < 1 || t_j > dims[2] ? continue : nothing

                dst = (t_j - 1) * dims[1] + t_i

                #Check if nodes is of right type
                layout[dst] < type ? continue : nothing

                w_adj[src,dst] = 1
            end
        end
    end

    return w_adj
end

function make_dist_matrix(adj, X, Y)
    dist = fill(NaN, size(adj))
    points = findall(adj .== 1)

    for p in points
        src = p[1]
        dst = p[2]

        dist[p] = sqrt((X[dst] - X[src])^2 + (Y[dst] - Y[src])^2)
    end

    return dist
end

function make_graph(adj, dist, X, Y, A)

    G = MetaGraph(
        DiGraph();  # underlying graph structure
        label_type = Int64,  # matrix index
        vertex_data_type = NTuple{3, Float64},  # Position and a value
        edge_data_type = Float64,  # edge weight
        weight_function=ed -> ed,
        default_weight=Inf,
        graph_data = "Graph"  # tag for the whole graph
    )

    nn = size(adj)[1]
    for i in 1:nn
        if any(adj[i,:] .== 1)
            G[i] = (X[i], Y[i], A[i])
        end
    end

    for i in 1:nn, j in 1:nn
        if adj[i,j] == 1
            G[i,j] = dist[i,j]
        end
    end

    return G
end

function generate_instance(folder::String, width::Int, height::Int, util_func::Function; spacing::Int=100, buffer::Int=1, inside::Int=1)

    filename = joinpath(folder, "grid_$(width)_$(height).zip")

    layout = grid_layout(width, height, buffer=buffer, inside=inside)
    X, Y = generate_point_matrix(layout, spacing)
    #display(layout)

    # Set utlitiy values (attractvieness)
    A = util_func(X,Y)
    

    four_p = [(1,0),(0,1),(-1,0),(0,-1)]
    eight_p = [(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1)]

    w_adj = make_adj_matrix(eight_p, layout; type=1)
    r_adj = make_adj_matrix(four_p, layout; type=2)


    w_d = make_dist_matrix(w_adj, X, Y)
    r_d = make_dist_matrix(r_adj, X, Y)

    Gw = make_graph(w_adj, w_d, X, Y, A)
    Gr = make_graph(r_adj, r_d, X, Y, zeros(Float64, size(A)))

    #
    # Generate Route Graph
    #

    # Make Edge Dataframe
    route_edges = DataFrame(src=[], dst=[], weight=[])
    for edge in edge_labels(Gr)
        push!(route_edges, (code_for(Gr, edge[1]), code_for(Gr, edge[2]), Gr[edge[1], edge[2]] / 12))
    end

    # Make node dataframe
    route_nodes = DataFrame(idx=[], link_idx=[], offset=[], x=[], y=[])
    for label in labels(Gr)
        push!(route_nodes, (code_for(Gr, label), label, 0, Gr[label][1], Gr[label][2]))
    end
    route = DataFrameGraph(route_nodes, route_edges)

    #
    # Genreate Base Graph 
    #
    function edge_func(edge::NTuple{2, Int64}, dist::Float64)
        return Dict(:src=>edge[1], :dst=>edge[2], :weight=>dist / 1.25)
    end

    function node_func(node::Int64, node_data::NTuple{3, Float64})
        return Dict(:idx=>node, :weight=>node_data[3], :x=>node_data[1], :y=>node_data[2])
    end
    base = DataFrameGraph(Gw, node_func, edge_func)

    # Save Output
    save_opt_instance(filename, base, route)
    return nothing
end

function main2()

    layout = grid_layout(6,2)
    spacing = 100
    X, Y = generate_point_matrix(layout, spacing)
    display(layout)

    g_func = Gaussian(200, 100)
    A = g_func.(X,Y,X,0)
    A += g_func.(X,Y,0,Y)

    four_p = [(1,0),(0,1),(-1,0),(0,-1)]
    eight_p = [(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1)]

    w_adj = make_adj_matrix(eight_p, layout; type=1)
    r_adj = make_adj_matrix(four_p, layout; type=2)


    w_d = make_dist_matrix(w_adj, X, Y)
    r_d = make_dist_matrix(r_adj, X, Y)

    Gw = make_graph(w_adj, w_d, X, Y, A)
    Gr = make_graph(r_adj, r_d, X, Y, zeros(Float64, size(A)))


    # Generate min paths for opt graph
    # fw_state = floyd_warshall_shortest_paths(Gr.graph, weights(Gr))
    # road_dist_mat = fw_state.dists

    # k = 20
    # nt = size(road_dist_mat)[1]

    # opt_dists = zeros(Float64, (nt,nt))
    # opt_adj = zeros(Int64, (nt,nt))

    # for i in 1:nt
    #     idxs = partialsortperm(road_dist_mat[i,:], 1:k+1)
    #     opt_adj[i, idxs] .= 1
    #     opt_adj[i,i] = 0
    #     opt_dists[i, idxs] = road_dist_mat[i,idxs]
    # end
    # display(opt_dists)
    # display(opt_adj)


    #
    # Generate Route Graph
    #

    # Make Edge Dataframe
    route_edges = DataFrame(src=[], dst=[], weight=[])
    # for i in 1:nt, j in 1:nt
    #     opt_adj[i,j] == 0 ? continue : nothing

    #     push!(route_edges, (i, j, opt_dists[i,j] / 12))
    # end 
    for edge in edge_labels(Gr)
        push!(route_edges, (code_for(Gr, edge[1]), code_for(Gr, edge[2]), Gr[edge[1], edge[2]] / 12))
    end

    # Make node dataframe
    route_nodes = DataFrame(idx=[], link_idx=[], offset=[])
    for label in labels(Gr)
        push!(route_nodes, (code_for(Gr, label), label, 0))
    end
    route = DataFrameGraph(route_nodes, route_edges)

    #
    # Genreate Base Graph 
    #
    function edge_func(edge::NTuple{2, Int64}, dist::Float64)
        return Dict(:src=>edge[1], :dst=>edge[2], :weight=>dist / 1.25)
    end

    function node_func(node::Int64, node_data::NTuple{3, Float64})
        return Dict(:idx=>node, :weight=>node_data[3])
    end
    base = DataFrameGraph(Gw, node_func, edge_func)

    # Save Output
    save_opt_instance("Opt_Testing/opt_data/zip_test.zip", base, route)

    # base2, route2 = read_opt_instance("Opt_Testing/opt_data/zip_test.zip")
    # display(base2)
    # display(route2)
end

function main()

    g_func = Gaussian(200, 100)
    util_func(X::Matrix{Float64}, Y::Matrix{Float64}) = g_func.(X,Y,X,0)
    #util_func(X::Matrix{Float64}, Y::Matrix{Float64}) = g_func.(X,Y,X,0) + g_func.(X,Y,0,Y)

    shapes = [(6,2), (12,4), (24,8)]
    #shapes = [(4,4), (8,8), (12,12)]

    for s in shapes
        generate_instance("Opt_Testing/opt_data", s[1], s[2], util_func)
    end

end

main()