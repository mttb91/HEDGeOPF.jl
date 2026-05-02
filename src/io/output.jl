
"Convert string to SQL literal, escaping single quotes and backslashes."
function _sql_string(value::String)
    value = replace(value, "\\" => "/")
    return "'" * replace(value, "'" => "''") * "'"
end

"Create lookup dictionary for variable names in each component folder."
function _var_lookup()
    path = pwd()
    lookup = Dict{String, Vector{String}}()
    for comp in filter(isdir, readdir(path))
        # Keep only directories filled only with .csv files
        files = readdir(joinpath(path, comp))
        if !isempty(files) && all(endswith.(files, ".csv"))
            # Get all unique variable names in the component folder
            lookup[comp] = String.(unique(first.(split.(files, "-"))))
        end
    end
    return lookup
end

"Create nested lookup dictionary for cases in each worker and fold."
function _case_lookup(map::_DF.DataFrame)

    _DF.sort!(map, [:worker, :fold, :case, :topology_id])
    folds = sort(unique(map.fold))
    workers = sort(unique(map.worker))

    lookup = Dict(worker => Dict{Int, _DF.DataFrame}() for worker in workers)
    for (worker, value) in lookup
        map_worker = filter(:worker => ==(worker), map)
        for fold in folds
            map_fold = filter(:fold => ==(fold), map_worker)
            if !isempty(map_fold)
                value[fold] = _DF.select(map_fold, [:case, :uid, :topology_id])
            end
        end
    end
    return lookup
end

function _copy_topologies(out_dir::String)
    src = "graph"
    dst = joinpath(out_dir, src)
    _mkpath(dst)
    for file in filter(x -> endswith(x, ".zip"), sort(readdir(src)))
        cp(joinpath(src, file), joinpath(dst, file); force=true)
    end
    return nothing
end

"Remove all original dataset files"
function _cleanup(components::Vector{String}, cleanup:: Bool)
    if cleanup
        names = vcat(components, ["graph", "map.csv", "polytope.csv", "rng_state.bin"])
        for name in names
            rm(name; recursive=true)
        end
    end
    return nothing
end

function _write_parquet!(db::_DDB.DB, data::_DF.DataFrame, path::String)
    _mkpath(dirname(path))
    _DDB.register_data_frame(db, data, "temp")
    try
        sql = (
            "COPY (SELECT * FROM temp) " *
            "TO $(_sql_string(path)) (FORMAT PARQUET)"
        )
        _DDB.DBInterface.execute(db, sql)
    finally
        _DDB.unregister_data_frame(db, "temp")
    end
    return nothing
end

"Fuse cases across workers into each fold for a given variable."
function _combine_cases(var::String, component::String, cases::Dict{Int, Dict{Int, _DF.DataFrame}}, n_fold::Int)

    data = Dict(fold => _DF.DataFrame[] for fold in 1:n_fold)
    for worker in sort(collect(keys(cases)))
        filename = joinpath(component, "$var-$(worker).csv")
        df = _DF.DataFrame(CSV.File(
            filename, 
            header=1,
            types=Float32
            )
        )
        for fold in sort(collect(keys(cases[worker])))
            map = cases[worker][fold]
            chunk = df[map.case, :]
            _DF.insertcols!(chunk, 1, :uid => map.uid)
            _DF.insertcols!(chunk, 2, :topology_id => map.topology_id)
            push!(data[fold], chunk)
        end
    end
    return Dict(fold => reduce(_DF.vcat, data[fold]) for fold in 1:n_fold)
end


"""
    generate_split(setting::NamedTuple; dst::Union{String, Nothing} = nothing)

Split the AC-OPF dataset in folds for cross-validation and generate an `dst`
folder containing the split map and one parquet file per variable and fold.
Each parquet file has size S x (C + 2), where S is the number of AC-OPF instances
in the fold and C is the number of variables for the given component type. The
additional columns are:
- `uid`: unique identifier of the AC-OPF instance
- `topology_id`: identifier for the topology of the AC-OPF instance
"""
function generate_split(setting::NamedTuple; dst::Union{String, Nothing} = nothing)

    if isnothing(dst)
        dst = setting.DATASET.name
    end
    if !ispath(dst)
        _mkpath(dst)
    else
        msg = "The destination folder must be unique. Rename the existing one to avoid overwriting."
        throw(DomainError(dst, msg))
    end
    @assert isfile("map.csv") "The file map.csv does not exist in $(pwd()). Run `generate_uid` first."
    # Generate CV folds assignment for each sample in the map
    map = _DF.DataFrame(CSV.File("map.csv"))
    generate_cv_folds!(map, setting)

    # Generate vars and cases lookup dictionaries
    map_vars = _var_lookup()
    map_cases = _case_lookup(map)

    # Group cases by fold and save them in parquet format per variable
    db = _DDB.DB()
    _write_parquet!(db, map, joinpath(dst, "map.parquet"))
    for (component, vars) in map_vars
        for var in vars
            for (fold, data) in _combine_cases(var, component, map_cases, setting.DATASET.num_folds)
                _write_parquet!(db, data, joinpath(dst, component, "$(var)-$(fold).parquet"))
            end
        end
    end
     _DDB.DBInterface.close(db)
    
    _copy_topologies(dst)
    _cleanup(collect(keys(map_vars)), setting.DATASET.cleanup)
    return nothing
end

function generate_split(
    path::String;
    dst::Union{String, Nothing} = nothing,
    filename::String = "settings.yaml"
)

    setting = to_namedtuple(read_settings(joinpath(path, filename)));
    cd(joinpath(
        path,
        setting.PATH.output,
        first(split(setting.CASE.grid, ".")),
        setting.CASE.name
        )
    )

    generate_split(setting; dst=dst)
    return nothing
end
