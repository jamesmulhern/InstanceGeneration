using NearestNeighbors

function compute_vertex_weights(x::Vector{Float64}, y::Vector{Float64}, 
                                att_x::Vector{Float64}, att_y::Vector{Float64}, att_val::Vector{<:Real})
    
    @assert size(x) == size(y)
    @assert size(att_x) == size(att_y) == size(att_val)

    dist_limit = 500
    walk_speed = 1

    kdtree = KDTree(vcat(x', y'); leafsize = 10, reorder=true)
    idxs, dist = nn(kdtree,vcat(att_x',att_y'))

    # Init outputs
    np = length(x)
    demandSum = zeros(Float64, np)
    scaledSum = zeros(Float64, np)
    localDemand = zeros(Float64, np)

    for i in 1:np
        near_idxs = idxs .== i
        !any(near_idxs) ? continue : nothing

        # Check if max distance is violated
        max_dist = maximum(dist[near_idxs])
        if max_dist > dist_limit
            println("Point ($(x[i]),$(y[i])) has nearest neighbors over distance limit. Max distance: $(max_dist), Limit: $(dist_limit)")
        end

        # Compute sum and scaled sum of near demands
        demandSum[i] = sum(att_val[near_idxs])
        scaledSum[i] = sum(att_val[near_idxs] .* (walk_speed ./ dist[near_idxs]))

        # Compute local demand parameter
        local_idx = findall(near_idxs)
        for j in local_idx, k in local_idx
            j == k ? continue : nothing

            dist_val = sqrt((att_x[k]-att_x[j])^2 + (att_y[k]-att_y[j])^2)
            localDemand[i] += att_val[k] * att_val[j] * (walk_speed/dist_val)
        end
    end

    return demandSum, scaledSum, localDemand
end