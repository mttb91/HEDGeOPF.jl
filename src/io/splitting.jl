
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
belonging to a single fold.
"""
function generate_split(setting::NamedTuple)

    path = pwd()
    n_fold = setting.DATASET.num_folds
    n_quantile = setting.DATASET.num_quantiles
    rng = _RND.MersenneTwister(setting.CASE.baseseed)

    filename = joinpath(path, "map.csv")
    @assert isfile(filename) "The file map.csv does not exist in $path. Run generate_uid(path) first."
    map = _DF.DataFrame(CSV.File(filename))

    if length(unique(map.topology_id)) == 1
        folds = total_load_active_power_split(map, rng, n_fold, n_quantile)
    else
        folds = topology_split(map, rng, n_fold)
    end
    _DF.insertcols!(map, 2, :fold => folds)

    _DF.sort!(map, [:worker, :case, :topology_id])
    CSV.write(filename, map)
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
