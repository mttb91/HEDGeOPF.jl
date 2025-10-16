
###############################################################################
# Helper function
###############################################################################

"Create local path to destination folder and set it as current directory."
function _mkpath(settings::NamedTuple)

    path = joinpath(settings.PATH.output, first(split(settings.CASE.grid, ".")), settings.CASE.name)
    if !ispath(path)
        mkpath(path)
    else
        if !settings.CASE.append
            msg = "The destination folder must be unique. Rename the existing one to avoid overwriting."
            throw(DomainError(path, msg))
        end
    end
    cd(joinpath(pwd(), path))
    return nothing
end

"Create local path to folder if it does not exist already."
function _mkpath(path::String)

    if !ispath(path)
        mkpath(path)
    end
    return nothing
end

"Create folder if it does not exist already."
function _mkdir(path::String)

    if !isdir(path)
        mkdir(path)
    end
    return nothing
end

"Export a dictionary of DataFrames to Excel naming the sheets by the sorted dictionary keys."
function _to_xlsx(data::Dict{<:Any, _DF.DataFrame}, filepath::String)

    # Save static input information
    XLSX.openxlsx(filepath, mode="w") do xf
        for (i, key) in enumerate(sort!(collect(keys(data))))

            i == 1 ? XLSX.rename!(xf[i], String(key)) : XLSX.addsheet!(xf, String(key))
            XLSX.writetable!(xf[i], data[key])
        end
    end
end

###############################################################################
# I/O RNG
###############################################################################

"Freeze local RNG state and save it to file"
function save_rng(rng::_RND.AbstractRNG; filename::String = "rng_state.bin")

    open(joinpath(pwd(), filename), "w") do io
        serialize(io, rng)
    end
end

"Read RNG state if available else generate it"
function read_rng(settings::NamedTuple; filename::String = "rng_state.bin")

    if settings.CASE.append && isfile(joinpath(pwd(), filename))
        rng = open("rng_state.bin", "r") do io
            deserialize(io)
        end
    else
        rng = _RND.MersenneTwister(settings.CASE.baseseed)
    end
    return rng
end

###############################################################################
# I/O Polytope
###############################################################################

"Import polytope from .csv or .xlsx file"
function import_polytope(; path::Union{String, Nothing} = nothing, filename::String = "polytope")

    path = isnothing(path) ? dirname(pwd()) : path
    if isfile(joinpath(path, filename * ".xlsx"))
        polytope = _DF.DataFrame(XLSX.readtable(joinpath(path, filename * ".xlsx"), "polytope", infer_eltypes=true))
    elseif isfile(joinpath(path, filename * ".csv"))
        polytope = _DF.DataFrame(CSV.File(joinpath(path, filename * ".csv")))
    end
    b = polytope[!, end]
    A = _DF.select!(polytope, 1:(_DF.ncol(polytope) - 1))
    # Extract load indices
    cols = _DF.names(A)
    cols = [cols[contains.(cols, h)] for h in ["p", "q"]]
    ids = [parse.(Int, filter.(isdigit, h)) for h in cols]

    return (A = Matrix(A), b = b, ids = tuple(ids...))::PolyType
end

"Export polytope to .csv file"
function export_polytope(polytope::PolyType; filename::String = "polytope.csv")

    path = pwd()

    A, b, ids = polytope
    cols = reduce(vcat, [x .* string.(y) for (x, y) in zip(["p", "q"], ids)])
    push!(cols, "b")
    data = _DF.DataFrame(reduce(hcat, [A, b]), cols; copycols=false)
    ∉(filename, readdir(path)) && CSV.write(joinpath(path, filename), data)
    return nothing
end

###############################################################################
# I/O Power system
###############################################################################

"Extract data in COO from CSC matrix and save information to Dataframe"
function extract_coo_data(data::SparseMatrixCSC{<:Number, Int64})

    values = data.nzval
    colptr = data.colptr
    colidx = vcat([fill(j, d) for (j, d) in enumerate(diff(colptr))]...)
    return _DF.DataFrame(
        rows = data.rowval,
        cols = colidx,
        re = real.(values),
        im = imag.(values)
    )
end

"Export graph model based on PowerModels network in XLSX format"
function export_graph(model::_PM.AbstractPowerModel, topology::TopologyPerturbation; basename::String = "graph")

    path = joinpath(pwd(), "graph")
    _mkdir(path)

    # Update ref dictionary based on topology perturbation
    model = update_topology(model, topology)

    data = Dict{Symbol, _DF.DataFrame}()
    for element in Symbol.(["bus", "branch", "load", "shunt", "gen"])
        # Extract all features of each component table
        if !isempty(_PM.ref(model, element))
            vars = get_pm_key(model, element)
            filter!(x -> !in(x, ["index", "source_id"]), vars)
            data[element] = get_pm_value(model, element, sort!(vars), _DF.DataFrame)
        end
    end

    element = :gen
    vars = [:cost, :ncost]
    # Unfold cost vector into components
    _DF.select!(data[element], _DF.Not(vars),
        vars => _DF.ByRow((x, y) -> vcat(repeat([0.0], 3-y), x)) => [:c2, :c1, :c0]
    )
    _DF.select!(data[element], sort(names(data[element])))
    # Identify synchronous condensers (if any)
    ids = findall(x -> (x.pmax - x.pmin .== 0) && (x.qmax - x.qmin .!= 0), eachrow(data[element]))
    if !isempty(ids)
        data[:sync] = data[element][ids, [:gen_bus]]
    end

    # Add connection indices based on continuous bus mapping
    indices = calc_connection_indices(data)
    for (key, value) in indices
        data[key] = _DF.DataFrame(value, :auto)
    end
    # Add admittance matrices based on continuous bus mapping
    Y = calc_admittance_matrices(data, indices)
    for (key, value) in zip(keys(Y), Y)
        data[key] = extract_coo_data(value)
    end
    
    # Save graph data to xlsx
    filename = "$basename-$(string(topology.id)).xlsx"
    ∉(filename, readdir(path)) && _to_xlsx(data, joinpath(path, filename))

    return nothing
end

###############################################################################
# I/O Reproducibility
###############################################################################

"""
    generate_uid(path::String)

Generate unique identifier for AC-OPF instances by sorting them based on
total load active power, objective value and topology.
"""
function generate_uid(path::String)

    ls = _DF.DataFrame[]
    for file in filter(x -> contains(x, "info"), readdir(path))
        data = _DF.DataFrame(CSV.File(joinpath(path, file)))
        _DF.insertcols!(data, 1, :worker => parse(Int, file[6:end-4]))
        _DF.insertcols!(data, 2, :case => axes(data, 1))
        push!(ls, data)    
    end
    df = reduce(vcat, ls)
    # Sort AC-OPF instances based on total load active power and objective value
    _DF.sort!(df, [:pd_tot, :objective, :topology_id])
    _DF.insertcols!(df, 1, :uid => axes(df, 1))
    
    return df
end

###############################################################################
# Dataset split
###############################################################################

"""
    generate_split(setting::NamedTuple)

Split the AC-OPF dataset in folds for cross-validation, creating a reproducible mapping between
instance ID (in terms of worker, topology and case) and fold which it belongs to.
The mapping is saved as `map.csv` in the main dataset folder.

## Notes

The splitting strategy is automatically selected depending on whether a single or multiple topologies
are available:
- In case of single topology, splitting is performed as stratified in terms of total load active power.
Is is expected that AC-OPF instances are uniquely identified as sorted by this variable.
- In case of multiple topologies, these are distributed across the different folders, with each topology
beloning to a single fold.
"""
function generate_split(setting::NamedTuple)

    path = pwd()
    n_fold = setting.DATASET.num_folds
    n_quantile = setting.DATASET.num_quantiles
    rng = _RND.MersenneTwister(setting.CASE.baseseed)

    map = generate_uid(path)

    if length(unique(map.topology_id)) == 1
        folds = total_load_active_power_split(map, rng, n_fold, n_quantile)
    else
        folds = topology_split(map, rng, n_fold)
    end
    _DF.insertcols!(map, 2, :fold => folds)

    _DF.sort!(map, [:worker, :case, :topology_id])
    CSV.write(joinpath(path, "map.csv"), map)
    return nothing
end

"Split samples in CV folds for dataset with multiple topologies"
function topology_split(map::_DF.DataFrame, rng::_RND.AbstractRNG, n_fold::Int)

    _DF.sort!(map, :topology_id)

    ids = _RND.shuffle(rng, unique(map.topology_id))
    dim = length(copy(ids)) ÷ n_fold
    folds = zeros(Int, size(map, 1))
    for i in 1:n_fold

        num = 1:(i == n_fold ? length(ids) : dim)
        mask = in.(map.topology_id, Ref(ids[num]))
        folds[mask] .= i
        deleteat!(ids, num)
    end
    return folds
end

"Split samples in CV folds for dataset with single topology"
function total_load_active_power_split(map::_DF.DataFrame, rng::_RND.AbstractRNG, n_fold::Int, n_quantile::Int)

    _DF.sort!(map, :uid)
    @assert issorted(map.pd_tot)

    # Bin UID based on total load active power quantiles
    r = 1 / n_quantile
    quantiles = quantile(map.pd_tot, (0 + r):r:(1 - r); sorted=true)
    classes = searchsortedfirst.(Ref(quantiles), map.pd_tot)
    bins = Dict(i => findall(x -> x == i, classes) for i in 1:n_quantile)

    # Create folds with stratified sampling
    folds = zeros(Int, size(map, 1))
    for bin in values(bins)

        _RND.shuffle!(rng, bin)
        dim = length(bin) ÷ n_fold
        @assert dim > 0 "Not enough samples to create $n_fold folds with $n_quantile quantiles, reduce either."

        for (i, fold) in enumerate(_RND.shuffle(rng, 1:n_fold))
            num = 1:(i == n_fold ? length(bin) : dim)
            folds[bin[num]] .= fold
            deleteat!(bin, num)
        end
    end
    return folds
end
