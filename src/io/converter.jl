
const DEFAULT_MAPPING = Dict(
    :load   => ([:grid, :nodes, :load], [:pd, :qd]),
    :bus    => ([:solution, :nodes, :bus], [:va, :vm]),
    :gen    => ([:solution, :nodes, :generator], [:pg, :qg]),
    :branch => ([:solution, :edges, [:ac_line, :transformer]], [:pt, :qt, :pf, :qf])
)

struct GridNodes
    bus::Vector{Vector{Float32}}
    generator::Vector{Vector{Float32}}
    load::Vector{Vector{Float32}}
    shunt::Vector{Vector{Float32}}
end

struct EdgeWithFeatures
    senders::Vector{Int}
    receivers::Vector{Int}
    features::Vector{Vector{Float32}}
end

struct LinkEdge
    senders::Vector{Int}
    receivers::Vector{Int}
end

struct GridEdges
    ac_line::EdgeWithFeatures
    transformer::EdgeWithFeatures
    generator_link::LinkEdge
    load_link::LinkEdge
    shunt_link::LinkEdge
end

struct Grid
    nodes::GridNodes
    edges::GridEdges
    context::Vector{Vector{Vector{Float32}}}
end

struct SolutionNodes
    bus::Vector{Vector{Float32}}
    generator::Vector{Vector{Float32}}
end

struct SolutionEdges
    ac_line::EdgeWithFeatures
    transformer::EdgeWithFeatures
end

struct Solution
    nodes::SolutionNodes
    edges::SolutionEdges
end

struct Metadata
    objective::Float64
end

struct Root
    grid::Grid
    solution::Solution
    metadata::Metadata
end


function _read_json(path::AbstractString; type::Type=Root)
    return JSON.parsefile(path, type)
end

function _reorder_index(
    example::Root;
    labels::Vector{Symbol} = [:ac_line, :transformer],
    feature_names::Vector{Symbol} = [:f_bus, :t_bus, :transformer],
)

    # Read branch ordering from HEDGeOPF graph data
    path_zip = joinpath(pwd(), "graph", "1.zip")
    if !isfile(path_zip)
        msg = "Zip file for graph topology not found at $(path_zip). You need to manually place the
        original `1.zip` file generated with `HEDGeOPF` for the same case in the `graph` folder at
        path $(pwd())."
        error(msg)
    end
    bytes = _ZA.zip_readentry(_ZA.ZipReader(read(path_zip)), "branch.csv")
    data = CSV.read(IOBuffer(bytes), _DF.DataFrame)
    _DF.select!(data, feature_names)
    ids_branch_ref = Int.(Matrix(data))

    # Read branch ordering from one sample of OPFData (obtained by stacking ac lines and transformers)
    blocks = Matrix{Int}[]
    for label in labels
        senders = foldl(getfield, vcat(label, :senders); init=example.grid.edges)
        receivers = foldl(getfield, vcat(label, :receivers); init=example.grid.edges)
        @assert length(senders) == length(receivers)

        value = Matrix{Int}(undef, length(senders), 3)
        # Original indexing is python based (0-based), revert to Julia format
        value[:, 1] = senders .+ 1
        value[:, 2] = receivers .+ 1
        value[:, 3] .= label == first(labels) ? 0 : 1
        push!(blocks, value)
    end
    ids_branch = reduce(vcat, blocks)
    @assert size(ids_branch, 1) == size(ids_branch_ref, 1)

    # Build lookup from order in OPFData to original graph ordering
    order = Vector{Int}()
    for row in eachrow(unique(ids_branch_ref, dims=1))
        idx = findall(r -> all(r .== row), eachrow(ids_branch))
        if isempty(idx)
            error("Branch row $(row) from graph not found in OPFData sample edges.")
        end
        append!(order, idx)
    end
    @assert isequal(ids_branch[order, :], ids_branch_ref)
    return order
end

function _get_dims(
    example::Root;
    mapping::Dict{Symbol, Tuple{Vector, Vector{Symbol}}} = DEFAULT_MAPPING
)
    dims = Dict{Symbol, Int}()
    for (comp, keys) in mapping
        src_keys, _ = keys
        if comp != :branch
            value = foldl(getfield, src_keys, init=example)
        else
            value = Vector{Float32}[]
            for branch_type in last(src_keys)
                path = vcat(first(src_keys, 2), [branch_type], [:features])
                append!(value, foldl(getfield, path, init=example))
            end
        end
        value = reduce(hcat, value)
        dims[comp] = size(value, 2)
    end
    return dims
end

function _init_store(
    dims::Dict{Symbol, Int};
    mapping::Dict{Symbol, Tuple{Vector, Vector{Symbol}}} = DEFAULT_MAPPING
)
    store = Dict{Symbol, Dict{Symbol, Matrix{Float32}}}()
    for (comp, spec) in mapping
        vars = last(spec)
        inner = Dict{Symbol, Matrix{Float32}}()
        for v in vars
            inner[v] = fill(Float32(NaN), dims[comp], dims[:samples])
        end
        store[comp] = inner
    end
    store[:metadata] = Dict(
        :uid => fill(Float32(NaN), dims[:samples], 1),
        :fold => fill(Float32(NaN), dims[:samples], 1),
        :pd_tot => fill(Float32(NaN), dims[:samples], 1),
        :objective => fill(Float32(NaN), dims[:samples], 1),
        :topology_id => fill(Float32(1.0), dims[:samples], 1)
    )
    return store
end

function _add_sample!(
    store::Dict{Symbol, Dict{Symbol, Matrix{Float32}}},
    sample::Root,
    col::Int,
    index::Vector{Int};
    mapping::Dict{Symbol, Tuple{Vector, Vector{Symbol}}} = DEFAULT_MAPPING
)

    for (comp, (src_keys, dst_keys)) in mapping

        if comp != :branch
            value = foldl(getfield, src_keys, init=sample)
        else
            value = Vector{Float32}[]
            for branch_type in last(src_keys)
                path = vcat(first(src_keys, 2), [branch_type], [:features])
                append!(value, foldl(getfield, path, init=sample))
            end
        end
        # Reorder branch rows to match HEDGeOPF ordering
        if comp == :branch
            value = value[index]
        end
        for (i, key) in enumerate(dst_keys)
            store[comp][key][:, col] .= getindex.(value, i)
        end
    end

    key = :metadata
    store[key][:pd_tot][col, :] .= sum(store[:load][:pd][:, col])
    store[key][:objective][col, :] .= Float32(sample.metadata.objective)

    return nothing
end

function process_group(
    paths::NamedTuple{(:src, :dst), Tuple{String, String}},
    index::Vector{Int},
    dims::Dict{Symbol, Int},
    db::_DDB.DB
)

    files = first.(split.(filter(f -> startswith(f, "example_"), readdir(paths.src)), "."))
    uids = parse.(Int, last.(split.(files, "_")))
    uid_min = minimum(uids)
    @assert isequal(sort(uids), collect(uid_min:uid_min + length(files) - 1))

    dims[:samples] = length(files)
    store = _init_store(dims)
    store[:metadata][:fold] .= Float32(parse(Int, split(basename(paths.src), "_")[end]) + 1)

    Base.Threads.@threads for i in eachindex(files)

        col = uids[i] - uid_min + 1
        store[:metadata][:uid][col, :] .= Float32(uids[i] + 1)
        sample = _read_json(joinpath(paths.src, files[i] * ".json"));
        _add_sample!(store, sample, col, index)
    end

    map = _DF.DataFrame(
        Dict(k => in(k, [:pd_tot, :objective]) ? vec(v) : Int.(vec(v))
        for (k, v) in pop!(store, :metadata))
    )
    
    for (comp, value) in store
        for var in collect(keys(value))
            name = "$(var)-$(first(map.fold)).parquet"
            array = _DF.DataFrame(permutedims(pop!(value, var)), Symbol.(1:dims[comp]))
            _DF.insertcols!(array, 1, :uid => map.uid)

            _write_parquet!(db, array, joinpath(paths.dst, String(comp), name))
        end
    end
    return map
end

function convert_dataset(path::String;
    filename::String = "settings.yaml",
    dst::Union{String, Nothing} = nothing
)

    cd(path)

    setting = read_settings(filename);
    network = instantiate_network(setting);
    setting = to_namedtuple(setting);
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

    groups = filter(x -> startswith(x, "group_"), readdir())
    if isempty(groups)
        grid = first(split(setting.CASE.grid, "."))
        msg = replace("""
        No OPFDataset group folders found in the dataset path $(pwd()). You first need to
        separately download and unzip every group of OPFDataset for $(grid) (e.g., by using
        the `OPFDataset` dataset in Pytorch Geometric) and copy them in the folder at path
        $(pwd()) before running the conversion.
        """, "\n" => " ")
        error(msg)
    end

    db = _DDB.DB()
    map = _DF.DataFrame[]
    dims = Dict{Symbol, Int}()
    index = Vector{Int}()
    for (i, group) in enumerate(groups)
        paths = (src = joinpath(pwd(), group), dst = dst)
        if i == 1
            example = _read_json(joinpath(paths.src, first(readdir(paths.src))))
            index = _reorder_index(example)
            dims = _get_dims(example)
        end
        push!(map, process_group(paths, index, dims, db))
    end
    map = reduce(vcat, map)
    _write_parquet!(db, map, joinpath(dst, "map.parquet"))
    _DDB.DBInterface.close(db)
    return nothing
end