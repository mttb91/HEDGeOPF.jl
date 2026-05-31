
"""
    generate_cv_folds!(map:: _DF.DataFrame, setting::NamedTuple)

Generate cross-validation folds for AC-OPF instances based on selected splitting strategy and insert
the fold assignment as new column "fold" in the `map` dataframe.

The splitting strategy is automatically selected depending on whether a single or multiple topologies
are available:
- **In case of single topology**, splitting is performed as stratified in terms of total load active power.
Is is expected that AC-OPF instances are uniquely identified as sorted by this variable.
- **In case of multiple topologies**, these are distributed across the different folders, with each topology
belonging to a single fold.

"""
function generate_cv_folds!(map:: _DF.DataFrame, setting::NamedTuple)

    n_fold = setting.DATASET.num_folds
    n_sample = setting.DATASET.num_samples
    rng = _RND.MersenneTwister(setting.CASE.baseseed)
    @assert n_sample <= size(map, 1) "Number of samples to split cannot be greater than the total number of available samples."

    if length(unique(map.topology_id)) == 1
        folds = total_load_active_power_split(
            map,
            rng,
            n_fold,
            setting.DATASET.num_quantiles;
            n_sample = n_sample
        )
    else
        folds = topology_split(map, rng, n_fold; n_sample = n_sample)
    end
    _DF.insertcols!(map, 2, :fold => folds)
    # Remove unassigned instances from split-specific map in any
    # HACK: it is assumed that unassigned instances have fold value equal to -1
    filter!(:fold => !=(-1), map)

    @assert all(map.fold .> 0) "Valid fold values are integer-based and greater than 0."
    return nothing
end

"Split samples in CV folds for dataset with multiple topologies"
function topology_split(
    map::_DF.DataFrame,
    rng::_RND.AbstractRNG,
    n_fold::Int;
    n_sample::Union{Nothing, Int} = nothing
)
    _DF.sort!(map, :topology_id)

    if isnothing(n_sample)
        n_sample = size(map, 1)
    end

    ids = _RND.shuffle(rng, unique(map.topology_id))
    topo_idx = Dict(id => findall(==(id), map.topology_id) for id in ids)
    dim = n_sample ÷ n_fold
    rem = n_sample % n_fold
    folds = fill(-1, size(map, 1))
    for i in 1:n_fold
        target = dim + (i <= rem ? 1 : 0)
        assigned = 0

        while assigned < target
            id = popfirst!(ids)
            idx = topo_idx[id]
            num = target - assigned

            if length(idx) <= num
                folds[idx] .= i
                assigned += length(idx)
            else
                _RND.shuffle!(rng, idx)
                folds[idx[1:num]] .= i
                assigned += num
            end
        end
    end
    return folds
end

"Split samples in CV folds for dataset with single topology"
function total_load_active_power_split(
    map::_DF.DataFrame,
    rng::_RND.AbstractRNG,
    n_fold::Int,
    n_quantile::Int;
    n_sample::Union{Nothing, Int} = nothing
)
    _DF.sort!(map, :uid)
    @assert issorted(map.pd_tot)

    # Limit bin size if n_samples is provided
    bin_size = isnothing(n_sample) ? typemax(Int) : n_sample ÷ n_quantile

    # Bin UID based on total load active power quantiles
    r = 1 / n_quantile
    quantiles = quantile(map.pd_tot, (0 + r):r:(1 - r); sorted=true)
    classes = searchsortedfirst.(Ref(quantiles), map.pd_tot)
    bins = Dict(i => findall(x -> x == i, classes) for i in 1:n_quantile)

    # Create folds with stratified sampling
    folds = fill(-1, size(map, 1))
    for bin in values(bins)

        _RND.shuffle!(rng, bin)
        budget = min(length(bin), bin_size)
        dim = budget ÷ n_fold
        @assert dim > 0 "Not enough samples to create $n_fold folds with $n_quantile quantiles, reduce either."

        for (i, fold) in enumerate(_RND.shuffle(rng, 1:n_fold))
            num = 1:(i == n_fold ? budget - dim * (n_fold - 1) : dim)
            folds[bin[num]] .= fold
            deleteat!(bin, num)
        end
    end
    return folds
end
