
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

Fold assignment supports downsampling (i.e., limiting the number of AC-OPF instances used in
the dataset split). Instance to be removed are assigned fold value equal to -1 and filtered
out from the split-specific map.
"""
function generate_cv_folds!(map::_DF.DataFrame, setting::NamedTuple)

    n_fold = setting.DATASET.num_folds
    n_sample = setting.DATASET.num_samples
    rng = _RND.MersenneTwister(setting.CASE.baseseed)
    if isnothing(n_sample)
        n_sample = size(map, 1)
    else
        msg = "Available sample number is less than requested one for dataset split."
        n_sample > size(map, 1) && throw(ArgumentError(msg))
    end
    has_single_topology = length(unique(map.topology_id)) == 1

    if has_single_topology
        folds = total_load_active_power_split(
            map,
            rng,
            n_fold,
            n_sample,
            setting.DATASET.num_quantiles,
        )
    else
        folds = topology_split(map, rng, n_fold, n_sample)
    end
    _DF.insertcols!(map, 2, :fold => folds)
    # Remove unassigned instances from split-specific map in any
    # HACK: it is assumed that unassigned instances have fold value equal to -1
    filter!(:fold => !=(-1), map)

    _check_split(map, has_single_topology, n_sample)
    return nothing
end


"Check if fold assignment is valid and consistent"
function _check_split(map::_DF.DataFrame, has_single_topology::Bool, n_sample::Int)

    if !allunique(map.uid)
        error("Split map contains duplicated uid values after fold assignment.")
    end
    if any(map.fold .<= 0)
        error("Valid fold values are integer-based and greater than 0.")
    end

    # Verify that each topology appears at most in a single fold
    if !has_single_topology
        topo_fold = Dict{Int, Int}()
        for (id, fold) in zip(map.topology_id, map.fold)
            if haskey(topo_fold, id)
                if topo_fold[id] != fold
                    error("Topology $(id) appears in multiple folds: $(topo_fold[id]) and $(fold).")
                end
            else
                topo_fold[id] = fold
            end
        end
    end
    if size(map, 1) < n_sample
        msg = "$(size(map, 1)) AC-OPF instances are allocated to dataset split out of $(n_sample) requested."
        @warn msg
    end
    return nothing
end


"Split samples in CV folds for dataset with multiple topologies"
function topology_split(
    map::_DF.DataFrame,
    rng::_RND.AbstractRNG,
    n_fold::Int,
    n_sample::Int
)
    _DF.sort!(map, :topology_id)

    topo_ids = _RND.shuffle(rng, unique(map.topology_id))
    topo_row = Dict(id => findall(==(id), map.topology_id) for id in topo_ids)

    dim = n_sample ÷ n_fold
    rem = n_sample % n_fold
    folds = fill(-1, size(map, 1))
    residuals = Dict{Int, Int}(i => dim + (i <= rem ? 1 : 0) for i in 1:n_fold)
    for i in keys(residuals)

        while !isempty(topo_ids)
            id = popfirst!(topo_ids)
            rows = topo_row[id]
            # Assign topology to fold if all samples fit in remaining budget
            if length(rows) <= residuals[i]
                folds[rows] .= i
                residuals[i] -= length(rows)
            else
                pushfirst!(topo_ids, id)
                break
            end
        end
    end
    while !isempty(topo_ids)
        # Sort folds by largest residuals
        for i in first.(sort(collect(residuals), by = last, rev = true))
            if residuals[i] <= 0
                continue
            end
            if isempty(topo_ids)
                break
            end
            rows = topo_row[popfirst!(topo_ids)]
            num = min(length(rows), residuals[i])
            _RND.shuffle!(rng, rows)
            folds[rows[1:num]] .= i
            residuals[i] -= num
        end
        if maximum(values(residuals)) <= 0
            break
        end
    end

    return folds
end

"Split samples in CV folds for dataset with single topology"
function total_load_active_power_split(
    map::_DF.DataFrame,
    rng::_RND.AbstractRNG,
    n_fold::Int,
    n_sample::Int,
    n_quantile::Int,
)
    _DF.sort!(map, :uid)
    if !issorted(map.pd_tot)
        error("UIDs should be sorted by total load active power for stratified splitting")
    end

    # Limit bin size if n_samples is provided
    bin_size = n_sample ÷ n_quantile
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
        if dim <= 0
            msg = (
                "Not enough samples in quantile bin to create $n_fold folds with " *
                "$n_quantile quantiles. Reduce either or increase the dataset size."
            )
            error(msg)
        end

        for (i, fold) in enumerate(_RND.shuffle(rng, 1:n_fold))
            num = 1:(i == n_fold ? budget - dim * (n_fold - 1) : dim)
            folds[bin[num]] .= fold
            deleteat!(bin, num)
        end
    end
    return folds
end
