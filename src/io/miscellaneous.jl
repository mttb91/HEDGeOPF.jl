
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

"Export a dictionary of DataFrame to zip folder of CSV files"
function _to_zip(data::Dict{<:Any, _DF.DataFrame}, path::String)

    _mkpath(dirname(path))

    _ZA.ZipWriter(path) do w
        for key in sort(collect(keys(data)))

            _ZA.zip_newfile(w, "$key.csv"; compress=true)
            io = IOBuffer()
            try
                CSV.write(io, data[key])
                write(w, take!(io))
            finally
                close(io)
            end
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
    if all(isreal.(values))
        ntp = (rows = data.rowval, cols = colidx, re = real.(values))
    else
        ntp = (rows = data.rowval, cols = colidx, re = real.(values), im = imag.(values))
    end
    return _DF.DataFrame(ntp)
end

"Export graph model based on PowerModels network in XLSX format"
function export_graph(model::_PM.AbstractPowerModel, topology::TopologyPerturbation)

    path = joinpath(pwd(), "graph")
    _mkdir(path)

    # Update ref dictionary based on topology perturbation
    model = update_topology(model, topology)

    data = Dict{Symbol, _DF.DataFrame}()
    elements = [:bus, :branch, :load, :shunt, :gen]
    for element in elements
        # Extract all features of each component table
        if !isempty(_PM.ref(model, element))
            vars = get_pm_key(model, element)
            filter!(x -> !in(x, ["source_id"]), vars)
            data[element] = get_pm_value(model, element, sort!(vars), _DF.DataFrame)
        end
    end

    element, vars = :gen, [:cost, :ncost]
    # Unfold cost vector into components
    _DF.select!(data[element], _DF.Not(vars),
        vars => _DF.ByRow((x, y) -> vcat(repeat([0.0], 3-y), x)) => [:c2, :c1, :c0]
    )
    _DF.select!(data[element], sort(names(data[element])))

    # Add admittance matrices based on continuous bus mapping
    indices = calc_connection_indices(data)
    Y = calc_admittance_matrices(data, indices)
    for (key, value) in zip(keys(Y), Y)
        data[Symbol("_$key")] = extract_coo_data(value)
    end
    # Add susceptance matrices for Fast Decoupled Power Flow
    B = calc_susceptance_matrices(data, indices)
    for (key, value) in zip(keys(B), B)
        data[Symbol("_$key")] = extract_coo_data(value)
    end

    # Save graph data to zip of CSVs
    filename = "$(topology.id).zip"
    ∉(filename, readdir(path)) && _to_zip(data, joinpath(path, filename))

    return nothing
end

###############################################################################
# I/O Reproducibility
###############################################################################

"""
    generate_uid(cleanup::Bool)

Generate unique identifier for AC-OPF instances by sorting them based on
total load active power, objective value and topology.
"""
function generate_uid(cleanup::Bool)

    path = pwd()
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
    # Resort to original order based on worker, case and topology
    _DF.sort!(df, [:worker, :case, :topology_id])
    CSV.write(joinpath(path, "map.csv"), df)

    if cleanup
        for file in filter(x -> contains(x, "info"), readdir(path))
            rm(joinpath(path, file); force=true)
        end
    end
    return nothing
end
