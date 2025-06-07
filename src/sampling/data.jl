
Base.@kwdef mutable struct InputSample
    data::Vector{Float64}
    ids::Union{Vector{Int}, Nothing} = nothing
end

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

    return _DS.Uniform(pd_min, pd_max)
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

    return _DS.DiscreteNonParametric(x, mass)
end

"Generate batch of input load samples for OPF simulation by sampling the polytope uniformly in total load active power"
function generate_load_samples!(samples::Vector{Dict{String, <:Any}}, dist::_DS.Distribution, polytope::PolyType, n_samples::Tuple{Int, Int}, rng::_RND.AbstractRNG, settings::NamedTuple)

    # Get polytope
    A0, b0, ids = polytope
    ids_pd, ids_qd = ids
    nvar = length(ids_pd)
    ncon = length(ids_pd) * 2 + length(ids_qd) * 2 + length(findall(in(ids_qd), ids_pd)) * 2 + 2

    # Polytope total active power bounds
    # HACK: THIS DOES NOT ACCOUNT FOR LOADS WITH FIXED POWER
    pd_tot_max = load_power_bound(polytope, settings)
    pd_tot_min = load_power_bound(polytope, settings; upper=false)
    # Generate total load active power samples
    pd_tot = rand(rng, dist, first(n_samples))
    bounds = hcat(min.(pd_tot .+ 0.01, pd_tot_max), max.(pd_tot .- 0.01, pd_tot_min))
    bounds = bounds .* [1 -1]
    # Instantiate model for chebyshev center
    model = chebyshev_model(A0, b0, settings);
    # Sample polytope uniformly in total load active power
    loads = sample_polytope_uniformly(permutedims(bounds), convert(_SA.SharedMatrix{Float64}, A0), b0, model, last(n_samples), ncon, rng)
    # Relax total active power bounds
    delta = max(0.01, (b0[ncon-1] + b0[ncon]) * 0.001)
    bounds = pd_tot .+ [-delta delta]
    # Re-organize data structure
    for (i, load) in enumerate(loads)
        j = ((i - 1) ÷ last(n_samples)) + 1
        value = Dict(
            "info" => bounds[j, :],
            "load" => Dict("pd" => InputSample(load[1:nvar], ids_pd), "qd" => InputSample(load[(nvar+1):end], ids_qd))        
        )
        if isempty(get(samples, i, Dict()))
            push!(samples, value)
        else
            merge!(samples[i], value)
        end
    end
    loads=nothing
end

"Wrapper function to generate batch of input samples for OPF simulation"
function generate_input_samples(dist::_DS.Distribution, polytope::PolyType, n_samples::Tuple{Int, Int}, rng::_RND.AbstractRNG, settings::NamedTuple)

    samples = Vector{Dict{String, <:Any}}()
    generate_load_samples!(samples, dist, polytope, n_samples, rng, settings);

    # Add here other functions to generate input samples for additional variables (e.g. generator cost)
    # Order of input sample generation for different variables matters for reproducibility

    return samples
end