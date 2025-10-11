
###############################################################################
# Topology perturbation
###############################################################################

"""
    perturbe_topology(pm::_PM.AbstractPowerModel, k::Int, rng::_RND.AbstractRNG;
        ids_gen_faulted::Vector{Int} = Vector{Int}()
    )

Generate a single topology perturbation by removing up to `k` branches, defining and recording
a reference bus for every resulting island. If the island has no self-balancing capabilities 
(i.e. either generation or consumption is missing), the indices of isolated buses and generators
are recorded as well. 
"""
function perturbe_topology(pm::_PM.AbstractPowerModel, k::Int, rng::_RND.AbstractRNG;
    ids_gen_faulted::Vector{Int} = Vector{Int}()
)

    nbus = length(_PM.ref(pm, :bus))
    ngen = length(_PM.ref(pm, :gen))
    edges = get_pm_value(pm, :branch, ["f_bus", "t_bus"], Array{Any, 2})
    dim = size(edges, 1)

    bus_load = vec(get_pm_value(pm, :load, ["load_bus"], Array{Any, 2}))
    ids_gen_active = deleteat!(collect(1:ngen), ids_gen_faulted)
    data = get_pm_value(pm, :gen, ["gen_bus", "pmin", "pmax", "qmin", "qmax"], _DF.DataFrame; mask = ids_gen_active)
    bus_gen = data.gen_bus
    # Keep only generators with both active and reactive power support as reference bus candidates
    mask = .!iszero.(data.pmax .- data.pmin) .& .!iszero(data.qmax .- data.qmin)
    bus_gen_ref = bus_gen[mask]

    ids_ref, ids_bus, ids_gen = Vector{Int}(), Vector{Int}(), Vector{Int}()

    # Generate ids of branches to be removed
    n = _RND.rand(rng, 1:k)
    ids_branch = _RND.shuffle(rng, 1:dim)[1:n]
    sort!(ids_branch)
    # Get largest and all minor connected components
    lcc, ccs = connected_components(edges, ids_branch, nbus)
    # Define reference bus for largest connected component if missing
    push!(ids_ref, define_ref_bus(pm, lcc, bus_gen_ref))
    
    for cc in ccs
        sort!(cc)
        push!(ids_ref, define_ref_bus(pm, cc, bus_gen_ref))
        # Record the nodes and generators belonging to islands without self-balancing capabilities
        if !(any(in.(cc, Ref(bus_gen_ref))) && any(in.(cc, Ref(bus_load))))

            append!(ids_bus, cc)
            append!(ids_gen, ids_gen_active[findall(in.(bus_gen, Ref(cc)))])
        end
    end
    return (
        ids_branch = ids_branch,
        ids_bus = sort(ids_bus),
        ids_ref = sort(ids_ref),
        ids_gen = sort(ids_gen)
    )
end

"Generate a single generation perturbation by dropping at most one generator"
function perturbe_generation(pm::_PM.AbstractPowerModel, rng::_RND.AbstractRNG)

    dim = length(_PM.ref(pm, :gen))
    n = _RND.rand(rng, 0:1)
    if !iszero(n)
        ids = _RND.shuffle(rng, 1:dim)[1:n]
        sort!(ids)
    else
        ids = Vector{Int}()
    end
    return ids
end

"Wrapper generator for combined topology perturbations of different types"
function Base.iterate(gen::TopologyPerturbationGenerator, state::Nothing = nothing)
    setting = gen.setting.TOPOLOGY
    gen.state > setting.num_topo + 1 && return nothing

    pm = gen.model
    if gen.state === 1
        # Generate intact topology at first iteration
        limits = vec(sum(get_pm_value(pm, :gen, ["pmin", "pmax"], Array{Any, 2}), dims=1))
        perturbation = TopologyPerturbation(
            id = gen.state,
            pg_tot_bounds = limits
        )
    else
        if setting.generation
            ids_gen_faulted = perturbe_generation(pm, gen.rng)
        else
            ids_gen_faulted = Vector{Int}()
        end
        ids = perturbe_topology(pm, setting.k, gen.rng; ids_gen_faulted = ids_gen_faulted)

        # Add here additional perturbations

        ngen = length(_PM.ref(pm, :gen))
        mask = deleteat!(collect(1:ngen), sort(vcat(ids.ids_gen, ids_gen_faulted)))
        limits = vec(sum(get_pm_value(pm, :gen, ["pmin", "pmax"], Array{Any, 2}; mask = mask), dims=1))

        perturbation = TopologyPerturbation(
            gen.state,
            ids.ids_branch,
            ids.ids_bus,
            ids.ids_ref,
            ids.ids_gen,
            ids_gen_faulted,
            limits
        )
    end
    gen.state += 1
    return (perturbation, nothing)
end

"""
    generate_topologies(generator::TopologyPerturbationGenerator, dim::Int)

Generate topology perturbations, returning for each a number of input samples equal to `dim`
and the `mapping` grouping topologies by bounds in total load active power.
"""
function generate_topologies(generator::TopologyPerturbationGenerator, dim::Int)

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
        if i >= dim
            break
        end
    end
    return topologies, mapping
end

###############################################################################
# Load perturbation
###############################################################################

"Estimate required number of input samples per batch, accounting for convergence rate"
function estimate_count(data::Vector{Int64}, n_batches::Int, n_items::Int)

    if iszero(data[end])
        target = data[begin]
    else
        rate = data[end-1] / data[end]
        target = (data[begin] - data[end-1]) / rate
    end
    return repeat([ceil(Int, (target ÷ n_items) / n_batches)], n_batches)
end

"Instantiate uniform distribution in total load active power"
function dist_uniform(polytope::PolyType, settings::NamedTuple)

    # Get total load active power bounds for the polytope
    # HACK: THIS DOES NOT ACCOUNT FOR LOADS WITH FIXED POWER
    pd_max = load_power_bound(polytope, settings)
    pd_min = load_power_bound(polytope, settings; upper=false)

    return _DIST.Uniform(pd_min, pd_max)
end

"Instantiate nonparametric distribution on convergence obtained through Kernel Density Estimation"
function dist_nonparametric(data::Vector{Matrix{Float64}}; epsilon::Float64 = 1E-6)

    # Get total active power values for converged cases
    # HACK: CONVERGENCE MUST BE THE SECOND COLUMN
    data = reduce(vcat, data)
    data = data[.!iszero.(data[:, 2]), 1]

    grid = kde_lscv(data; boundary=extrema(data), npoints=2^15)
    # Inverse probability density function
    pdf = 1 ./ (grid.density .+ epsilon)
    # Derive probability mass
    x = (grid.x[1:end-1] .+ grid.x[2:end]) ./ 2
    mass = ((pdf[1:end-1] .+ pdf[2:end]) ./ 2) .* diff(grid.x)
    mass ./= sum(mass)

    return _DIST.DiscreteNonParametric(x, mass)
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
    bounds = hcat(
        min.(pd_tot .+ 0.01, b0[ncon-1]),
        max.(pd_tot .- 0.01, -b0[ncon])
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