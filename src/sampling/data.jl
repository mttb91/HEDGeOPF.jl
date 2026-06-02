
###############################################################################
# Topology perturbation
###############################################################################

"""
    perturb_topology(pm::_PM.AbstractPowerModel, rng::_RND.AbstractRNG, k::Int;
        ids_bus_ref_valid::Vector{Int} = Vector{Int}()
    )

Generate a single topology perturbation by removing up to `k` branches at random and by
keeping only perturbations that result in a single connected component. Returns:
- `ids_branch`: the ids of removed branches
- `ids_bus`: the ids of disconnected buses (empty in case islanding is not modelled)
- `ids_ref`: the ids of the reference buses, one for each connected component
- `ids_gen`: the ids of disconnected generators (empty in case islanding is not modelled)
"""
function perturb_topology(pm::_PM.AbstractPowerModel, rng::_RND.AbstractRNG, k::Int;
    ids_bus_ref_valid::Vector{Int} = Vector{Int}()
)

    ks = 1:k
    nbus = length(_PM.ref(pm, :bus))
    edges = get_pm_value(pm, :branch, ["index", "f_bus", "t_bus"], Array{Any, 2})
    dims = 1:size(edges, 1)

    ids_branch = Vector{Int}()
    ids_ref, ids_bus, ids_gen = Vector{Int}(), Vector{Int}(), Vector{Int}()
    if !iszero(k)
        while true
            # Generate ids of branches to be removed
            n = _RND.rand(rng, ks)
            ids = sort(_RND.shuffle(rng, dims)[1:n])
            ids_branch = edges[ids, 1]
            # Get largest and all minor connected components
            lcc, ccs = connected_components(edges[:, end-1:end], ids, nbus)

            if isempty(ccs)
                push!(ids_ref, define_ref_bus(pm, lcc, ids_bus_ref_valid))
                break
            end
        end
    end

    return (
        ids_branch = ids_branch,
        ids_bus = ids_bus,
        ids_ref = ids_ref,
        ids_gen = ids_gen
    )
end

"""
    perturb_generation(pm::_PM.AbstractPowerModel, rng::_RND.AbstractRNG, k::Int)

Generate a single generation perturbation by dropping at most `k` generators at random as
long as at least an eligible reference generator candidate is retained in the system `pm`.
Return:
- `ids_gen_faulted`: the ids of faulted generators
- `ids_bus_ref_valid`: the ids of remaining eligible reference bus candidates
"""
function perturb_generation(pm::_PM.AbstractPowerModel, rng::_RND.AbstractRNG, k::Int)

    # Define candidate reference buses for the given model
    data = define_candidate_ref_buses(pm)
    ids_gen = data.index
    mask_valid = data.mask_ref .& data.mask_unshared
    ids_gen_valid = data.index[mask_valid]
    ids_bus_ref_valid = data.gen_bus[mask_valid]

    if iszero(k)
        ids = Vector{Int}()
    else
        ks = 0:k
        while true
            # Generate ids of generator to be removed
            n = _RND.rand(rng, ks)
            if !iszero(n)
                ids = _RND.shuffle(rng, ids_gen)[1:n]
                sort!(ids)
            else
                ids = Vector{Int}()
            end
            # Keep pertrubation only if eligible reference generators are retained
            if !isempty(setdiff(ids_gen_valid, ids))
                break
            end
        end
    end
    ids_pos = findall(in(ids), data.index)
    return (
        ids_gen_faulted = ids,
        ids_bus_ref_valid = setdiff(ids_bus_ref_valid, data.gen_bus[ids_pos])
    )
end

"Wrapper generator for combined topology perturbations of different types"
function Base.iterate(gen::TopologyPerturbationGenerator, state::Nothing = nothing)
    setting = gen.setting
    gen.state > ceil(Int, setting.TOPOLOGY.num_topo * 1.2) + 1 && return nothing

    pm = gen.model
    if gen.state === 1
        # Generate intact topology at first iteration
        limits = vec(sum(get_pm_value(pm, :gen, ["pmin", "pmax"], Array{Any, 2}), dims=1))
        perturbation = TopologyPerturbation(
            id = gen.state,
            pg_tot_bounds = limits
        )
    else
        p_gen = perturb_generation(pm, gen.rng, setting.TOPOLOGY.k_gen)
        p_branch = perturb_topology(pm, gen.rng, setting.TOPOLOGY.k_branch;
            ids_bus_ref_valid = p_gen.ids_bus_ref_valid
        )

        # Add here additional perturbations

        # Compute limits in total generator active power for given topology
        ids_gen = sort(collect(_PM.ids(pm, :gen)))
        ids_gen_remove = Set(sort(vcat(p_branch.ids_gen, p_gen.ids_gen_faulted)))
        limits = get_pm_value(pm, :gen, ["pmin", "pmax"], Array{Any, 2};
            mask = findall(x -> !in(x, ids_gen_remove), ids_gen)
        )

        perturbation = TopologyPerturbation(
            gen.state,
            p_branch.ids_branch,
            p_branch.ids_bus,
            p_branch.ids_ref,
            p_branch.ids_gen,
            p_gen.ids_gen_faulted,
            vec(sum(limits, dims=1))
        )
    end
    gen.state += 1
    return (perturbation, nothing)
end

"""
    generate_topologies(generator::TopologyPerturbationGenerator, num::Int)

Generate `num` topology perturbations with generator `generator`, returning, alongside the
list of topologies, the `mapping` grouping topologies by bounds in total load active power.
"""
function generate_topologies(generator::TopologyPerturbationGenerator, num::Int)

    mapping = Dict{Vector{Float64}, Vector{Int}}()
    topologies = Vector{TopologyPerturbation}()
    for (i, t) in enumerate(generator)
        push!(topologies, t)
        key = t.pg_tot_bounds
        if !haskey(mapping, key)
            mapping[key] = [t.id]
        else
            push!(mapping[key], t.id)
        end
        if i >= num
            break
        end
    end
    return topologies, mapping
end

###############################################################################
# Load perturbation
###############################################################################

"""
    record_convergence!(
        convergence::Dict{Vector{Float64}, _DF.DataFrame},
        info::_DF.DataFrame,
        mapping::Dict{Vector{Float64}, Vector{Int}}
    )

Record OPF convergence for different topologies, grouping results by bounds in
total load active power.
"""
function record_convergence!(
    convergence::Dict{Vector{Float64}, _DF.DataFrame},
    info::_DF.DataFrame,
    mapping::Dict{Vector{Float64}, Vector{Int}}
)
    for (b, ids) in mapping
        mask = in.(Int.(info.id), Ref(ids))
        if !haskey(convergence, b)
            convergence[b] = info[mask, [:pd_tot, :status]]
        else
            convergence[b] = vcat(convergence[b], info[mask, [:pd_tot, :status]])
        end
    end
    return nothing
end

"""
    update_pd_distributions!(
        distributions::Dict{Vector{Float64}, _DIST.Distribution},
        convergence::Dict{Vector{Float64}, _DF.DataFrame}
    )

Update the uniform distributions in total load active power based on convergence
results, choosing for each distribution the minimum support between the distribution
own extrema and the global bounds refined based on convergence.
"""
function update_pdtot_distributions!(
    distributions::Dict{Vector{Float64}, _DIST.Distribution},
    convergence::Dict{Vector{Float64}, _DF.DataFrame}
)
    data = reduce(vcat, values(convergence))
    delta = first(diff([extrema(data.pd_tot)...])) * 0.01
    bounds_new = [extrema(data[Bool.(data.status), :pd_tot])...]
    bounds_new += [-delta, delta]

    for bounds in keys(distributions)
        bmin = max(minimum(bounds), minimum(bounds_new))
        bmax = min(maximum(bounds), maximum(bounds_new))
        distributions[bounds] = _DIST.Uniform(bmin, bmax)
    end
    return nothing
end

"Estimate the remaining number of input OPF samples to be generated for each topology"
function estimate_sample_number(
    counter::ConvergenceCounter,
    convergence::Dict{Vector{Float64}, _DF.DataFrame},
    mapping::Dict{Vector{Float64}, Vector{Int}}
)
    idx = minimum(reduce(vcat, values(mapping))) - 1
    rate = zeros(length(counter.n_sample))
    for (b, ids) in mapping
        rate[ids .- idx] .= sum(convergence[b].status) / size(convergence[b], 1)
    end
    rate = vec(minimum(hcat(rate, counter.n_converged ./ counter.n_iter), dims=2))
    target = (counter.n_sample - counter.n_converged) ./ rate

    # Do not generate additional samples for topologies with zero convergence
    ids_zero = findall(isinf, target)
    target[ids_zero] .= -1.0
    n_sample_new = repeat([maximum(counter.n_sample)], length(ids_zero))

    return ceil.(Int, target), n_sample_new
end

"""
    generate_load_samples!(
        samples::Vector{InputSample},
        polytope::PolyType,
        map_dist::Dict{_DIST.Distribution, Vector{Int}},
        rng::_RND.AbstractRNG,
        settings::NamedTuple
    )

Generate active and reactive load power perturbations for each topology specified
in `sample`. Sampling is performed in the convex load `polytope` uniformly
in total load active power, where different bounds are applied on total active power
depending on the sample topology as specified in `map_dist`.
"""
function generate_load_samples!(
    samples::Vector{InputSample},
    polytope::PolyType,
    map_dist::Dict{_DIST.Distribution, Vector{Int}},
    rng::_RND.AbstractRNG,
    settings::NamedTuple
)

    # Get load polytope
    A0, b0, (ids_pd, ids_qd) = polytope
    ncon = length(ids_pd) * 2 + length(ids_qd) * 2 + length(findall(in(ids_qd), ids_pd)) * 2 + 2

    # Generate total load active power samples
    pd_tot = Vector{Float64}[]
    for (dist, ids_topo) in map_dist
        for id in ids_topo
            n_sample = count(x -> x.topology.id == id, samples)
            push!(pd_tot, _RND.rand(rng, dist, n_sample))
        end
    end
    ids = reduce(vcat, collect(values(map_dist)))
    pd_tot = pd_tot[sortperm(ids)]
    pd_tot = reduce(vcat, pd_tot)
    ranges = getfield.(getfield.(samples, :topology), :pg_tot_bounds)
    bounds = hcat(
        min.(pd_tot .+ 0.01, last.(ranges)),
        max.(pd_tot .- 0.01, first.(ranges))
    )
    bounds = bounds .* [1 -1]
    # Sample polytope uniformly in total load active power
    model = chebyshev_model(A0, b0, settings);
    loads = sample_polytope_uniformly(
        permutedims(bounds),
        convert(_SA.SharedMatrix{Float64}, A0),
        b0, model, 1, ncon, rng
    )
    # Relax total active power bounds for OPF simulation
    delta = max(0.01, (b0[ncon-1] + b0[ncon]) * 0.001)
    bounds = permutedims(pd_tot .+ [-delta delta])

    @assert length(samples) === length(loads)
    for (i, (sample, load)) in enumerate(zip(samples, loads))
        sample.pd_tot = bounds[:, i]
        sample.data["load"] = Dict(
            "pd" => InputData(load[1:length(ids_pd)], ids_pd), 
            "qd" => InputData(load[(length(ids_pd)+1):end], ids_qd)
        )
    end
    bounds = nothing
    loads = nothing
end

###############################################################################
# Wrapper functions
###############################################################################

"Instantiate `n_sample` samples for each topology perturbation in `topologies`"
function instantiate_input_samples(topologies::Vector{TopologyPerturbation}, n_sample::Vector{Int})

    ls = Vector{InputSample}[]
    for (i, t) in enumerate(topologies)
        push!(ls, [InputSample(topology=t) for _ in 1:n_sample[i]])
    end
    return reduce(vcat, ls)
end

"Wrapper function to generate batch of input sample (e.g. load or gen cost perturbations) for a set of topology perturbations"
function generate_input_samples!(
    samples::Vector{InputSample},
    distributions::Dict{String, Dict{Vector{Float64}, _DIST.Distribution}},
    polytope::PolyType,
    rng::_RND.AbstractRNG,
    settings::NamedTuple
)
    ids = getfield.(getfield.(samples, :topology), :id)
    topologies = [samples[findfirst(==(v), ids)].topology for v in unique(ids)]

    var = "load"
    # Extract extrema from existing total load active power distributions
    if !isempty(distributions[var])
        bounds = reduce(hcat, collect.(extrema.(values(distributions[var]))))
        bounds_max = [minimum(bounds[1, :]), maximum(bounds[2, :])]
    else
        bounds_max = Vector{Float64}()
    end 

    map_dist = Dict{_DIST.Distribution, Vector{Int}}()
    for t in topologies
        key = t.pg_tot_bounds
        a, b = copy(key)
        # Use already-refined distribution bounds if available
        if !isempty(bounds_max)
            a = max(key[1], bounds_max[1])
            b = min(key[2], bounds_max[2])
        end
        if !haskey(distributions[var], key)
            # Add new distribution if missing
            distributions[var][key] = _DIST.Uniform(a, b)
        end
        dist = distributions[var][key]
        if !haskey(map_dist, dist)
            map_dist[dist] = [t.id]
        else
            push!(map_dist[dist], t.id)
        end
    end
    generate_load_samples!(samples, polytope, map_dist, rng, settings);

    # Add here other functions to generate input samples for additional variables (e.g. generator cost)
    # Order of input sample generation for different variables matters for reproducibility

    return nothing
end