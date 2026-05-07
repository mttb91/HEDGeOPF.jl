
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

function _worker_from_info_filename(file::String)
    stem = splitext(file)[1]
    return parse(Int, split(stem, "-")[end])
end

function _worker_ids(path::String)
    files = filter(x -> startswith(x, "info-") && endswith(x, ".csv"), readdir(path))
    return sort(_worker_from_info_filename.(files))
end

function _quantile_bins_from_indices(map::_DF.DataFrame, ids::Vector{Int}, n_quantile::Int)

    bins = [Int[] for _ in 1:n_quantile]
    isempty(ids) && return bins

    values = map.pd_tot[ids]
    if n_quantile == 1
        append!(bins[1], ids)
        return bins
    end
    r = 1 / n_quantile
    q = quantile(values, (0 + r):r:(1 - r))
    classes = searchsortedfirst.(Ref(q), values)
    for (id, c) in zip(ids, classes)
        push!(bins[c], id)
    end

    return bins
end

function _sample_from_bins(bins::Vector{Vector{Int}}, target::Int, rng::_RND.AbstractRNG)

    available = sum(length, bins)
    target = min(target, available)
    target <= 0 && return Int[]

    local_bins = [copy(_RND.shuffle(rng, b)) for b in bins]
    selected = Int[]
    while length(selected) < target
        progressed = false
        for b in local_bins
            if !isempty(b)
                push!(selected, popfirst!(b))
                progressed = true
                length(selected) == target && break
            end
        end
        !progressed && break
    end
    return selected
end

function _assign_bins_to_folds!(
    folds::Vector{Int},
    bins::Vector{Vector{Int}},
    rng::_RND.AbstractRNG,
    n_fold::Int
)
    for bin in bins
        isempty(bin) && continue
        _RND.shuffle!(rng, bin)
        order = _RND.shuffle(rng, collect(1:n_fold))
        for (k, idx) in enumerate(bin)
            folds[idx] = order[mod1(k, n_fold)]
        end
    end
    return nothing
end

"Duplicate slack-active AC-OPF instances in place by appending rows to all worker csv files."
function duplicate_unfeasible_instances!(setting::NamedTuple)

    path = pwd()
    load_path = joinpath(path, "load")
    @assert isdir(load_path) "Missing folder `load` in dataset path."

    tol = setting.DATASET.separate_unfeasible.slack_tol
    for worker in _worker_ids(path)

        info_file = joinpath(path, "info-$worker.csv")
        info = _DF.DataFrame(CSV.File(info_file))
        if !hasproperty(info, :is_strict_feasible)
            _DF.insertcols!(info, :is_strict_feasible => trues(_DF.nrow(info)))
        end

        ids_candidate = findall(info.is_strict_feasible)
        isempty(ids_candidate) && continue

        pd_file = joinpath(load_path, "pd-$worker.csv")
        qd_file = joinpath(load_path, "qd-$worker.csv")
        psu_file = joinpath(load_path, "pd_slack_up-$worker.csv")
        psd_file = joinpath(load_path, "pd_slack_down-$worker.csv")
        qsu_file = joinpath(load_path, "qd_slack_up-$worker.csv")
        qsd_file = joinpath(load_path, "qd_slack_down-$worker.csv")

        required = [pd_file, qd_file, psu_file, psd_file, qsu_file, qsd_file]
        if !all(isfile.(required))
            @warn "Skipping worker $worker due to missing load/slack csv files."
            continue
        end

        pd = _DF.DataFrame(CSV.File(pd_file))
        qd = _DF.DataFrame(CSV.File(qd_file))
        pd_su = _DF.DataFrame(CSV.File(psu_file))
        pd_sd = _DF.DataFrame(CSV.File(psd_file))
        qd_su = _DF.DataFrame(CSV.File(qsu_file))
        qd_sd = _DF.DataFrame(CSV.File(qsd_file))

        n = _DF.nrow(info)
        @assert all(==(n), _DF.nrow.([pd, qd, pd_su, pd_sd, qd_su, qd_sd])) "Worker $worker has inconsistent row counts."

        active = falses(n)
        for i in ids_candidate
            active[i] = any(abs.(Float64.(collect(pd_su[i, :]))) .> tol) ||
                        any(abs.(Float64.(collect(pd_sd[i, :]))) .> tol) ||
                        any(abs.(Float64.(collect(qd_su[i, :]))) .> tol) ||
                        any(abs.(Float64.(collect(qd_sd[i, :]))) .> tol)
        end
        ids_dup = findall(active)
        isempty(ids_dup) && continue

        pd_unf = copy(pd[ids_dup, :])
        qd_unf = copy(qd[ids_dup, :])
        for (k, i) in enumerate(ids_dup)
            pd_unf[k, :] .= Float64.(collect(pd[i, :])) .+ Float64.(collect(pd_su[i, :])) .- Float64.(collect(pd_sd[i, :]))
            qd_unf[k, :] .= Float64.(collect(qd[i, :])) .+ Float64.(collect(qd_su[i, :])) .- Float64.(collect(qd_sd[i, :]))
        end
        pd_tot_unf = vec(sum(Matrix(pd_unf), dims=2))

        element_dirs = filter(x -> isdir(joinpath(path, x)) && x != "graph", readdir(path))
        for element in sort(element_dirs)
            folder = joinpath(path, element)
            files = filter(f -> endswith(f, "-$worker.csv"), readdir(folder))
            for file in files
                fpath = joinpath(folder, file)
                data = _DF.DataFrame(CSV.File(fpath))
                @assert _DF.nrow(data) == n "Worker $worker file $(fpath) is out-of-sync with info rows."

                data_dup = copy(data[ids_dup, :])
                if element == "load"
                    stem = splitext(file)[1]
                    var = replace(stem, "-$worker" => "")
                    if var == "pd"
                        data_dup = copy(pd_unf)
                    elseif var == "qd"
                        data_dup = copy(qd_unf)
                    end
                end

                append!(data, data_dup)
                CSV.write(fpath, data)
            end
        end

        info_dup = copy(info[ids_dup, :])
        info_dup.pd_tot .= pd_tot_unf
        info_dup.is_strict_feasible .= false
        append!(info, info_dup)
        CSV.write(info_file, info)
    end

    return nothing
end

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
        if !hasproperty(data, :is_strict_feasible)
            _DF.insertcols!(data, :is_strict_feasible => trues(_DF.nrow(data)))
        end
        _DF.insertcols!(data, 1, :worker => parse(Int, file[6:end-4]))
        _DF.insertcols!(data, 2, :case => axes(data, 1))
        push!(ls, data)    
    end
    df = reduce(vcat, ls)
    # Sort AC-OPF instances based on total load active power and objective value
    _DF.sort!(df, [:pd_tot, :objective, :topology_id, :is_strict_feasible])
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
