using Pkg
Pkg.activate(".")

using YAML
using DataFrames, CSV
using NearestNeighbors
using Plots
using ArgParse

const ParamDict = Dict{Symbol, Any}

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "setup_file"
            help = "Setup YAML file.  This file defines the details of the utility computation"
            required = true
    end

    return parse_args(s)
end

function GaussianVolume(V::Real, sigma::Real)::Function
    A = V/2/pi/sigma^2
    f(x,y,x_0,y_0) = A * exp(-((x-x_0)^2/(2*sigma^2) + (y-y_0)^2/(2*sigma^2)))
    return f
end

function save_utility(points::DataFrame, params::ParamDict)
    # Build File name
    filename = ""
    segments = Vector{String}()

    haskey(params, :prefix) ? push!(segments, params[:prefix]) : nothing
    haskey(params, :label) ? push!(segments, params[:label]) : push!(segments, "Utility")
    haskey(params, :data) ? append!(segments, string.(params[:data])) : nothing

    for (i, s) in pairs(segments)
        i == 1 ? filename = s : filename = filename * "_" * s 
    end
    filename = filename * ".csv"

    # Trim points below utility_threshold
    if haskey(params, :utility_threshold) && params[:utility_threshold] != 0
        rows = points.utility .<= params[:utility_threshold] * maximum(points.utility)
        deleteat!(points,rows)
    end

    # Write output
    CSV.write(joinpath(params[:folder], filename), select(points, [:x, :y, :utility]))
    return nothing
end

function plot_utility(points::DataFrame, params::ParamDict)

    x_unq = unique(points.x)
    y_unq = unique(points.y)
    sort!(x_unq, rev=false)
    sort!(y_unq, rev=false)

    mat = fill(NaN, length(y_unq), length(x_unq))

    for (x, y, val) in zip(points.x, points.y, points.utility)
        x_idx = findfirst(x_unq .== x)
        y_idx = findfirst(y_unq .== y)
        mat[y_idx, x_idx] = val
    end

    # Create Color Bar
    if haskey(params, :color_bar)
        c_scale = cgrad(Symbol(params[:color_bar][:pallet]), rev=params[:color_bar][:rev], scale=Symbol(params[:color_bar][:scale]))
    else
        c_scale = cgrad(:roma, rev=true, scale=:exp)
    end
    fig = heatmap(x_unq, y_unq, mat, aspect_ratio=:equal, c=c_scale)

    # Add Title
    haskey(params,:title) ? plot!(fig, title=params[:title]) : nothing

    # Display
    display(fig)
end

function dist(A::NTuple{2,Float64},B::NTuple{2,Float64})
    return sqrt((B[1]-A[1])^2 + (B[2]-A[2])^2)
end

# Computes the external untility of the provided regions
function ext_utility(points::DataFrame, params::ParamDict)

    # Init Config dataframe
    config = DataFrame(type=Vector{String}(), num_points=Vector{Int}(), weight=Vector{Float64}(), util_per=Vector{Float64}())
    
    # Extract number of POI of each type and type data
    data_p = params[:data]
    for el in data_p[:datasets]
        name = split(el[:file], ".")[1]
        df = DataFrame(CSV.File(joinpath(data_p[:folder], el[:file])))
        push!(config, (name, nrow(df), el[:A], NaN))
    end

    # Compute utility per POI location
    config.util_per = config.weight .* params[:scale] ./ config.num_points

    display(config)

    # Find all exteranl POI files
    files = readdir(params[:input][:folder])
    deleteat!(files, findfirst(files .== params[:input][:points_file]))

    # Load exteranl points file
    ext_points = DataFrame(CSV.File(joinpath(params[:input][:folder], params[:input][:points_file])))
    select!(ext_points, [:id, :Location, :x, :y] .=> [:idx, :loc, :x, :y])
    ext_points.util_sum = zeros(Float64, nrow(ext_points))
    ext_points.num_POI = zeros(Int, nrow(ext_points))
    ext_points.dist_sum = zeros(Float64, nrow(ext_points))

    display(ext_points)


    for f in files

        # Determine the file location and the type
        type = "" 
        loc = ""
        e_idx = 0
        c_idx = 0
        for (i, t) in pairs(config.type)
            contains(f, t) ? (type = t; c_idx = i) : nothing
        end
        for (i, l) in pairs(ext_points.loc)
            contains(f, l) ? (loc = l; e_idx = i) : nothing
        end

        df = DataFrame(CSV.File(joinpath(params[:input][:folder], f)))

        ext_points.util_sum[e_idx] += config.util_per[c_idx] * nrow(df)

        ctr = (ext_points.x[e_idx], ext_points.y[e_idx])
        p = collect(zip(df.x,df.y))
        if !isempty(p)
            d = dist.([ctr], p)
            ext_points.num_POI[e_idx] += length(d)
            ext_points.dist_sum[e_idx] += params[:dist_factor] * sum(d)
        end

    end

    ext_points.dist_avg = ext_points.dist_sum ./ ext_points.num_POI
    display(ext_points)

    # Save Output
    CSV.write(joinpath(params[:output][:folder], params[:output][:file]), ext_points)


end

# Adds sparse data to the utility distrobution
function add_sparse!(x::Vector{Float64}, y::Vector{Float64}, params::ParamDict)
    val = zeros(Float64, size(x))

    # TODO: I can remove the volume form this section. it cancels out

    # Generate Smoothing Function
    smooth_func = getfield(Main, Symbol(params[:smooth_func]))(params[:volume], params[:spread])
    A_sum = 0

    # Compute POI contibutions
    for sub_d in params[:datasets]
        file = sub_d[:file]
        A = sub_d[:A]
        A_sum += A
        dest_df = DataFrame(CSV.File(joinpath(params[:folder], file)))
        n_dest = nrow(dest_df)

        println("A: $A, A_sum: $(A_sum), n_dest: $(n_dest)")

        for row in eachrow(dest_df)
            val += (A / n_dest) .* smooth_func.(x, y, row.x, row.y)
        end
    end

    val ./= (params[:volume]*A_sum)
    return val
end

# Adds dense data to the utility distrobution
function add_dense!(x::Vector{Float64}, y::Vector{Float64}, params::ParamDict)
    val = zeros(Float64, size(x))

    # Load Pop Data
    data = DataFrame(CSV.File(joinpath(params[:folder], params[:file])))
    data_sum = sum(data[!, Symbol(params[:data_field])])

    # Create kdtree for nearest neighbors
    kdtree = KDTree(vcat(data.x', data.y'); leafsize = 10, reorder=true)
    idxs, _ = nn(kdtree, vcat(x', y'))

    # Add data to the val vector, scalled so integral is 1
    val += (1 / data_sum) .* (data[idxs, Symbol(params[:data_field])] ./ params[:spacing]^2)
    return val
end

function main()
    args = parse_commandline()
    setup_file = args["setup_file"]
    #setup_file = "InstanceGeneration/config_files/All_Cat.utility.yaml"

    cfg = YAML.load_file(setup_file, dicttype=Dict{Symbol,Any})

    # Load Point Data
    points = DataFrame(CSV.File(joinpath(cfg[:points][:folder], cfg[:points][:file])))
    points[!,:utility] = zeros(Float64, nrow(points))

    # Loop over utility components
    for comp in cfg[:components]
        points[!, :utility] += getfield(Main, Symbol(comp[:method]))(points.x, points.y, comp[:params])
        points.utility .*= comp[:scale]
    end

    # Compute utility value from utility density. Increase by general scale value
    points.utility .*= cfg[:points][:gen_scale] * cfg[:points][:spacing]^2

    # Loop over actions
    for action in cfg[:actions]
        getfield(Main, Symbol(action[:method]))(points, action[:params])
    end

end

main()