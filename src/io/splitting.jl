
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

The setting `DATASET.unfeasible_ratio` controls the dataset ratio (in %) between unfeasible and feasible
instances. If it is zero, unfeasible instances are removed and the split is performed on the feasible
ones only. Otherwise, it is supported only for multiple topologies, see [`topology_split`](@ref).
"""
function generate_cv_folds!(map::_DF.DataFrame, setting::NamedTuple)

    n_fold = setting.DATASET.num_folds
    n_sample = setting.DATASET.num_samples
    ratio = setting.DATASET.unfeasible_ratio
    rng = _RND.MersenneTwister(setting.CASE.baseseed)
    has_single_topology = length(unique(map.topology_id)) == 1

    if ratio == 0
        # Keep only feasible instances, which are no longer paired
        filter!(:is_strictly_feasible => identity, map)
        map.pair .= -1
    end
    n_target = something(n_sample, size(map, 1))
    if n_target > size(map, 1)
        msg = "Available sample number is less than requested one for dataset split."
        throw(ArgumentError(msg))
    end

    if has_single_topology
        if ratio == 0
            folds = total_load_active_power_split(map, rng, n_fold, n_target, setting.DATASET.num_quantiles)
        else
            throw(ArgumentError("Option `unfeasible_ratio` must be zero for a single topology."))
        end
    else
        if ratio == 0
            folds = topology_split(map, rng, n_fold, n_target)
        else
            folds = topology_split(map, rng, n_fold, n_sample, ratio)
            n_target = 0 # the feasible sample budget is checked in `topology_split`
        end
    end
    _DF.insertcols!(map, 2, :fold => folds)
    # Remove unassigned instances from split-specific map in any
    # HACK: it is assumed that unassigned instances have fold value equal to -1
    filter!(:fold => !=(-1), map)

    _check_split(map, has_single_topology, n_target)
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


"Select `m` out of `n` indices with a random start and a fixed stride."
function _systematic(rng::_RND.AbstractRNG, n::Int, m::Int)
    m >= n && return collect(1:n)
    stride = n / m
    return floor.(Int, rand(rng) * stride .+ (0:(m - 1)) .* stride) .+ 1
end

"""
    topology_split(map, rng, n_fold, n_sample, ratio)

Split samples in CV folds for dataset with multiple topologies, controlling the dataset ratio `ratio` (in %)
between unfeasible and feasible instances. `nothing` keeps all unfeasible instances of the selected
topologies. `n_sample` is the total number of feasible instances, with `nothing` using as many as fit.
It is reduced, with a warning, if `ratio` is above the natural one and `n_sample` cannot be met.
"""
function topology_split(
    map::_DF.DataFrame,
    rng::_RND.AbstractRNG,
    n_fold::Int,
    n_sample::Union{Int, Nothing},
    ratio::Union{Float64, Nothing}
)
    _DF.sort!(map, :topology_id)
    topo = _DF.combine(_DF.groupby(map, :topology_id), :is_strictly_feasible => sum => :nf, _DF.nrow => :n)
    n_topo = _DF.nrow(topo)
    if n_topo < n_fold
        throw(ArgumentError("At least $n_fold topologies are required to create $n_fold folds."))
    end
    nf, nu = topo.nf, topo.n .- topo.nf
    rho = isnothing(ratio) ? 0.0 : ratio / 100

    # Select topologies spread along their unfeasible ratio and distribute them across folds in blocks of
    # `n_fold`. A budget is met with the smallest number of complete blocks, while without a budget all
    # topologies are used and the last incomplete block goes to the folds with fewest feasible instances
    order = sortperm(tuple.(nu ./ nf, rand(rng, n_topo)))
    fold_topo = zeros(Int, n_topo)
    F, U, cap = zeros(Int, n_fold), zeros(Int, n_fold), zeros(Int, n_fold)
    for m in (isnothing(n_sample) ? (n_topo:n_fold:n_topo) : (n_fold:n_fold:(n_fold * (n_topo ÷ n_fold))))
        sel = order[_systematic(rng, n_topo, m)]
        pos = sortperm(((0:(m - 1)) .÷ n_fold) .+ rand(rng, m))
        full = m - m % n_fold
        fold_topo .= 0
        fold_topo[sel[pos[1:full]]] .= mod1.(1:full, n_fold)
        fold_topo[sel[pos[(full + 1):m]]] .= sortperm((fold_topo .== (1:n_fold)')' * nf)[1:(m - full)]
        onehot = fold_topo .== (1:n_fold)'
        F, U = onehot' * nf, onehot' * nu
        cap = isnothing(ratio) ? F : min.(F, floor.(Int, U ./ rho .+ 1e-9))
        sum(cap) >= something(n_sample, typemax(Int)) && break
    end
    count(>(0), fold_topo) ÷ n_fold < 20 && @warn "Fewer than 20 topologies per fold: fold balance may be poor."

    # Number of feasible instances in the split
    n_split = min(something(n_sample, sum(cap)), sum(cap))
    short = something(n_sample, n_split) - n_split
    if short > 0
        if !(!isnothing(ratio) && rho > sum(nu) / sum(nf)) && short > max(1e-3 * n_sample, 1)
            throw(ArgumentError("Requested $n_sample feasible instances, but at most $n_split can be split in $n_fold folds."))
        end
        @warn "Requested $n_sample feasible instances cannot be met: $n_split are used."
    end
    f = copy(cap)
    for _ in 1:(sum(f) - n_split)
        f[argmax(f)] -= 1 # cut the largest folds first
    end
    u = isnothing(ratio) ? U : min.(floor.(Int, rho .* f .+ 1e-9), U)

    # Select instances within each fold, preserving their distribution in total load active power
    fold_row = fold_topo[searchsortedfirst.(Ref(topo.topology_id), map.topology_id)]
    folds = fill(-1, size(map, 1))
    for k in 1:n_fold, (feasible, n) in ((true, f[k]), (false, u[k]))
        rows = findall((fold_row .== k) .& (map.is_strictly_feasible .== feasible))
        folds[rows[sortperm(map.pd_tot[rows])][_systematic(rng, length(rows), n)]] .= k
    end

    # Unpair instances whose pair partner was removed
    kept = _DS.counter(filter(>(0), map.pair[folds .> 0]))
    map.pair .= ifelse.(getindex.(Ref(kept), map.pair) .== 2, map.pair, -1)
    return folds
end
