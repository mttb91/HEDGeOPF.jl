
const DEFAULT_COMP = Dict{String, Vector{String}}(
    "branch" => ["pf", "qf", "pt", "qt"],
    "bus" => ["va", "vm"],
    "gen" => ["pg", "qg"],
    "load" => ["pd", "qd"],
)


function _decompress_h5_gz(
    path::AbstractString;
    overwrite::Bool = false,
    chunk_bytes::Int = 1 << 23,
)
    endswith(path, ".gz") || error("Expected a .gz file, got: $path")
    isfile(path) || error("Input file not found: $path")
    h5_path = replace(path, r"\.gz$" => "")
    if isfile(h5_path) && !overwrite
        return h5_path
    end
    mkpath(dirname(h5_path))

    _GZ.open(path, "r") do fin
        open(h5_path, "w") do fout
            buffer = Vector{UInt8}(undef, chunk_bytes)
            while !eof(fin)
                n_read = readbytes!(fin, buffer, chunk_bytes)
                n_read == 0 && break
                write(fout, view(buffer, 1:n_read))
            end
        end
    end
    return h5_path
end

function _read_h5(
    path::AbstractString,
    cols::Vector{String};
    rows::Union{Nothing, Vector{Int}} = nothing,
)
    HDF5.h5open(path, "r") do file

        out = Dict{String, Any}()
        for col in cols
            haskey(file, col) || error("Variable $(col) not found in H5 file $(path)")
            data = file[col]
            ndim = ndims(data)

            if isnothing(rows)
                value = read(data)
            else
                if ndim == 1
                    value = read(data)[rows]
                elseif ndim == 2
                    value = read(data)[:, rows]
                else
                    error("Expected at most 2D datasets in H5 file, got $(ndim)D for variable $(col) in file $(path)")
                end
            end
            if ndim == 2
                value = permutedims(value)
            end
            out[col] = value
        end
        return out
    end
end

function _read_graph()

    filepath = joinpath(pwd(), "graph", "1.zip")
    if !isfile(filepath)
        msg = replace("""
        No graph folder found in the dataset path $(pwd()). You first need to copy the graphs
        folder, generated with `HEDGeOPF` and containing the base topology file `1.zip` at
        at path $(pwd()) before running the conversion.
        """, "\n" => " ")
        error(msg)
    end

    data = Dict{String, _DF.DataFrame}()
    for comp in keys(DEFAULT_COMP)
        bytes = _ZA.zip_readentry(_ZA.ZipReader(read(filepath)), "$comp.csv")
        data[comp] = CSV.read(IOBuffer(bytes), _DF.DataFrame)
    end
    return data
end

function build_map(path::AbstractString)
    path_meta = _decompress_h5_gz(joinpath(path, "ACOPF", "meta.h5.gz"))
    path_input = _decompress_h5_gz(joinpath(path, "input.h5.gz"))
    
    data = Dict{Symbol, Any}(
        :pd_tot => vec(sum(first(values(_read_h5(path_input, ["pd"]))), dims=2)),
        :objective => first(values(_read_h5(path_meta, ["primal_objective_value"]))),
    )
    data[:row] = 1:length(data[:pd_tot])

    data = _DF.DataFrame(data)
    _DF.sort!(data, [:pd_tot, :objective])
    data[!, :uid] = 1:size(data, 1)
    data[!, :topology_id] .= 1
    return data
end


function convert_dataset(path::String;
    folder::String = "train",
    n_samples::Union{Int, Nothing} = nothing,
    filename::String = "settings.yaml",
    dst::Union{String, Nothing} = nothing
)

    cd(path)

    setting = read_settings(filename);
    if folder == "train"
        value = isnothing(n_samples) ? setting["DATASET"]["num_samples"] : n_samples
    elseif folder == "test"
        value = nothing
    else
        error("Invalid folder name: $folder. Expected 'train' or 'test'.")
    end
    setting["DATASET"]["num_samples"] = value
    setting = to_namedtuple(setting)
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

    graph = _read_graph()
    files = filter(x -> startswith(x, "input"), readdir(folder))
    if isempty(files)
        grid = first(split(setting.CASE.grid, "."))
        msg = replace("""
        No PGLearn results .h5.gz file found in the dataset path $(pwd()). You first need to
        separately download the AC-OPF dataset with `PGLearn` for $(grid) and copy the result
        file in the folder at path $(pwd()) before running the conversion.
        """, "\n" => " ")
        error(msg)
    end
    path_src = joinpath(pwd(), folder)

    map = build_map(path_src)
    generate_cv_folds!(map, setting)
    @info "Splitting PGLearn dataset into $(setting.DATASET.num_folds) folds."

    db = _DDB.DB()
    for fold in unique(map.fold)
        mask = map.fold .== fold
        uids = map.uid[mask]
        rows = map.row[mask]

        for (comp, vars) in DEFAULT_COMP
            filename = comp == "load" ? "input" : "ACOPF//primal"
            filepath = _decompress_h5_gz(joinpath(path_src, "$filename.h5.gz"))
            
            for (var, data) in _read_h5(filepath, vars; rows=rows)
                @assert size(data, 2) == size(graph[comp], 1) "Mismatch in number of $(comp)."
                df = _DF.DataFrame(data, Symbol.(graph[comp].index))
                _DF.insertcols!(df, 1, :uid => uids)

                filepath = joinpath(dst, comp, "$(var)-$fold.parquet")
                _write_parquet(db, df, filepath)
            end
        end
    end
    _write_parquet(db, map, joinpath(dst, "map.parquet"))
    _DDB.DBInterface.close(db)
    _copy_topologies(dst)
    return nothing
end