
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
    n_quantile = setting.DATASET.num_quantiles
    split_kwargs = setting.DATASET.separate_unfeasible
    rng = _RND.MersenneTwister(setting.CASE.baseseed)

    if length(unique(map.topology_id)) == 1
        if split_kwargs.active
            ratio = effective_unfeasible_ratio(
                map,
                split_kwargs.unfeasible_ratio,
                split_kwargs.total_size
            )
            folds = mixed_total_load_active_power_split(
                map,
                rng,
                n_fold,
                n_quantile,
                ratio;
                total_size = split_kwargs.total_size
            )
        else
            folds = total_load_active_power_split(map, rng, n_fold, n_quantile)
        end
    else
        folds = topology_split(map, rng, n_fold)
    end
    _DF.insertcols!(map, 2, :fold => folds)
    # Remove unassigned instances from split-specific map in any
    # HACK: it is assumed that unassigned instances have fold value equal to -1
    filter!(:fold => !=(-1), map)

    @assert all(map.fold .> 0) "Valid fold values are integer-based and greater than 0."
    return nothing
end

"Return effective unfeasible ratio constrained by available data."
function effective_unfeasible_ratio(map::_DF.DataFrame, desired_ratio::Float64, total_size::Int)

    n_target = total_size < 0 ? _DF.nrow(map) : min(total_size, _DF.nrow(map))
    n_unf = count(.!map.is_strict_feasible)
    achievable = iszero(n_target) ? 0.0 : min(1.0, n_unf / n_target)
    ratio = min(desired_ratio, achievable)

    if desired_ratio > achievable
        @warn "Desired unfeasible ratio $(desired_ratio) is not achievable with available data. Using $(ratio)."
    end
    return ratio
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

"Split fixed-topology dataset with feasible/unfeasible strata while preserving target ratio globally."
function mixed_total_load_active_power_split(
    map::_DF.DataFrame,
    rng::_RND.AbstractRNG,
    n_fold::Int,
    n_quantile::Int,
    unfeasible_ratio::Float64;
    total_size::Int = -1
)

    @assert hasproperty(map, :is_strict_feasible) "Column `is_strict_feasible` is required for mixed split."

    ids_f = findall(map.is_strict_feasible)
    ids_u = findall(.!map.is_strict_feasible)

    if total_size < 0
        n_f_target = length(ids_f)
        if iszero(n_f_target)
            n_u_target = length(ids_u)
        elseif unfeasible_ratio <= 0.0
            n_u_target = 0
        else
            n_u_target = min(length(ids_u), round(Int, unfeasible_ratio * n_f_target / (1 - unfeasible_ratio)))
        end
    else
        n_target = min(total_size, _DF.nrow(map))
        n_u_target = min(length(ids_u), round(Int, unfeasible_ratio * n_target))
        n_f_target = min(length(ids_f), n_target - n_u_target)

        n_selected = n_u_target + n_f_target
        if n_selected < n_target
            rem = n_target - n_selected
            add_f = min(rem, max(length(ids_f) - n_f_target, 0))
            n_f_target += add_f
            rem -= add_f
            add_u = min(rem, max(length(ids_u) - n_u_target, 0))
            n_u_target += add_u
        end
    end

    bins_f = _quantile_bins_from_indices(map, ids_f, n_quantile)
    bins_u = _quantile_bins_from_indices(map, ids_u, n_quantile)

    n_selected = n_f_target + n_u_target
    if n_selected > 0
        realized_ratio = n_u_target / n_selected
        if !isapprox(realized_ratio, unfeasible_ratio; atol=1e-12, rtol=0.0)
            @warn "Requested unfeasible ratio of $(unfeasible_ratio) could not be matched exactly after availability/rounding fallback. Using $(realized_ratio) with $(n_u_target) unfeasible and $(n_f_target) feasible samples."
        end
    end

    selected_f = _sample_from_bins(bins_f, n_f_target, rng)
    selected_u = _sample_from_bins(bins_u, n_u_target, rng)

    folds = fill(-1, _DF.nrow(map))
    bins_f_sel = _quantile_bins_from_indices(map, selected_f, n_quantile)
    bins_u_sel = _quantile_bins_from_indices(map, selected_u, n_quantile)
    _assign_bins_to_folds!(folds, bins_f_sel, rng, n_fold)
    _assign_bins_to_folds!(folds, bins_u_sel, rng, n_fold)

    return folds
end
