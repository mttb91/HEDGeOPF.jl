
const DEFAULT_MAPPING = Dict{String, Dict{Symbol, Symbol}}(
    "var" => Dict(
        :pd => :pd,
        :qd => :qd,
        :p_fr => :pf,
        :p_to => :pt,
        :q_fr => :qf,
        :q_to => :qt,
        :p_gen => :pg,
        :q_gen => :qg,
        :va_bus => :va,
        :vm_bus => :vm,
    ),
    "comp" => Dict(
        :bus => :bus,
        :gen => :gen,
        :line => :branch,
        :load => :load,
    )
)

Base.@kwdef struct GlobalMap
    uid::Vector{Int} = Vector{Int}()
    fold::Vector{Int} = Vector{Int}()
    pd_tot::Vector{Float32} = Vector{Float32}()
    topology_id::Vector{Int} = Vector{Int}()
end

function _separate_numeric_from_string(s::AbstractString)
    left, right = split(s, ":"; limit=2)
    m = match(r"^([A-Za-z_]+)(\d+)$", left)
    m === nothing && error("Expected letters+digits before ':', got: $left (from $s)")
    text_part   = first(m.captures)
    number_part = parse(Int, last(m.captures))
    return (Symbol(text_part), number_part, Symbol(right))
end

function _get_column_map(code::Vector{Tuple{Symbol, Int, Symbol}})

    map_cols = Dict{Symbol, Dict{Symbol, Vector{Int}}}()
    for comp in unique(first.(code))
        map_cols[comp] = Dict{Symbol, Vector{Int}}()
        for var in unique(last.(filter(x -> first(x) == comp, code)))
            map_cols[comp][var] = findall(x -> first(x) == comp && last(x) == var, code)
        end
    end
    return map_cols
end

function _extract_info(rows::CSV.Rows)

    n_samples = sum(1 for _ in rows)

    cols = String.(rows.names)
    code = _separate_numeric_from_string.(cols)
    return (
        n_samples, 
        (names = cols, ids = getindex.(code, 2)),
        _get_column_map(code)
    )
end

function convert_dataset(path::String;
    n_samples::Union{Int, Nothing} = nothing,
    filename::String = "settings.yaml",
    dst::Union{String, Nothing} = nothing
)

    cd(path)

    setting = to_namedtuple(read_settings(filename));
    cd(joinpath(
        pwd(),
        setting.PATH.output,
        first(split(setting.CASE.grid, ".")),
        setting.CASE.name
    ))

    if isnothing(dst)
        dst = setting.DATASET.name
    end
    if !ispath(dst)
        _mkpath(dst)
    else
        msg = "The destination folder must be unique. Rename the existing one to avoid overwriting."
        throw(DomainError(dst, msg))
    end

    files = filter(x -> startswith(x, "pglib_opf_case"), readdir())
    if isempty(files)
        grid = first(split(setting.CASE.grid, "."))
        msg = replace("""
        No OPFLearn results .csv file found in the dataset path $(pwd()). You first need to
        separately generate the AC-OPF dataset with `OPFLearn` for $(grid) by using release
        0.1.2 and copy the result file in the folder at path $(pwd()) before running the
        conversion.
        """, "\n" => " ")
        error(msg)
    end
    if !isfile(joinpath(pwd(), "graph", "1.zip"))
        msg = replace("""
        No graph folder found in the dataset path $(pwd()). You first need to copy the graphs
        folder, generated with `HEDGeOPF` and containing the base topology file `1.zip` at
        at path $(pwd()) before running the conversion.
        """, "\n" => " ")
        error(msg)
    end
    path_src = joinpath(pwd(), first(files))

    if isnothing(n_samples)
        rows = CSV.Rows(path_src; header=1, types=Float32, reusebuffer=false)
    else
        rows = CSV.Rows(path_src; header=1, limit=n_samples, types=Float32, reusebuffer=false)
        msg = "Downsampling OPFLearn dataset by retaining only the first $(n_samples) samples."
        @warn msg
    end
    n_sample, cols, map_column = _extract_info(rows)
    chuck_size = floor(Int, n_sample / setting.DATASET.num_folds)
    @info "Splitting OPFLearn dataset into $(setting.DATASET.num_folds) folds."

    c = 0
    db = _DDB.DB()
    map = GlobalMap(topology_id = fill(1, n_sample))
    for (i, chunk) in enumerate(Base.Iterators.partition(rows, chuck_size))
        data = _DF.DataFrame(chunk)
        uids = collect(1:size(data, 1)) .+ c
        for (comp, vars) in map_column
            comp_dst = String(DEFAULT_MAPPING["comp"][comp])
            for (var, ids) in vars
                var_dst = get(DEFAULT_MAPPING["var"], var, nothing)
                if isnothing(var_dst)
                    continue
                end
                df = _DF.select(data, cols.names[ids]; copycols=true)

                if var_dst == :pd
                    append!(map.pd_tot, vec(sum(Matrix(df), dims=2)))
                end
                _DF.rename!(df, Symbol.(cols.ids[ids]))
                _DF.insertcols!(df, 1, :uid => uids)

                filepath = joinpath(dst, String(comp_dst), "$(var_dst)-$i.parquet")
                _write_parquet!(db, df, filepath)
            end
        end
        append!(map.uid, uids)
        append!(map.fold, fill(i, size(data, 1)))

        c += size(data, 1)
    end
    map = _DF.DataFrame(Dict(k => getfield(map, k) for k in propertynames(map)))
    _write_parquet!(db, map, joinpath(dst, "map.parquet"))
    _DDB.DBInterface.close(db)
    _copy_topologies(dst)
    return nothing
end
