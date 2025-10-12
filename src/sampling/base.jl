
"Generate batch of OPF instances uniformly in total load active power for multiple topology perturbations of the same power system"
function generate_batch!(
    distributions::Dict{String, Dict{Vector{Float64}, _DIST.Distribution}},
    convergence::Dict{Vector{Float64}, _DF.DataFrame},
    generator::TopologyPerturbationGenerator,
    model::_PM.AbstractPowerModel,
    polytope::PolyType,
    n_sample::Vector{Int},
    rng::_RND.AbstractRNG,
    setting::NamedTuple
)

    dim = length(n_sample)
    database = instantiate_database(model, setting.MODEL.duals)

    topologies, mapping = generate_topologies(generator, dim)
    samples = instantiate_input_samples(topologies, n_sample)
    generate_input_samples!(samples, distributions, polytope, rng, setting)

    counter = ConvergenceCounter(
        getfield.(topologies, :id),
        n_sample,
        zeros(Int, dim),
        zeros(Int, dim)
    )
    info = generate_opf_samples!(counter, model, database, samples)
    record_convergence!(convergence, info, mapping)

    n_sample = estimate_sample_number(counter, convergence, mapping)
    if any(n_sample .> 0)
        samples = instantiate_input_samples(topologies, n_sample)
        generate_input_samples!(samples, distributions, polytope, rng, setting)
        info = generate_opf_samples!(counter, model, database, samples)
        record_convergence!(convergence, info, mapping)
    end
    export_graph.(Ref(model), topologies)
end

function generate_opf_instances(model::_PM.AbstractPowerModel, polytope::PolyType, rng::_RND.AbstractRNG, setting::NamedTuple)

    # Batch samples
    n_topo = setting.TOPOLOGY.num_topo + 1
    n_batch = setting.CASE.num_batches
    n_topo_batch = div(n_topo, n_batch)
    if iszero(n_topo_batch)
        n_sample = ceil(Int, setting.CASE.num_samples / n_batch)
        n_samples = repeat([[n_sample]], n_batch)
    else
        n_sample = ceil(Int, setting.CASE.num_samples / n_topo)
        repeats = repeat([n_topo_batch], n_batch)
        repeats[end] += n_topo - sum(repeats)
        n_samples = repeat.([[n_sample]], repeats)
    end

    # Instantiate topology generator
    generator = TopologyPerturbationGenerator(model=model, rng=rng, setting=setting)
    # Instantiate relevant dictionaries
    convergence = Dict{Vector{Float64}, _DF.DataFrame}()
    distributions = Dict("load" => Dict{Vector{Float64}, _DIST.Distribution}())

    for n_sample in n_samples

        generate_batch!(
            distributions,
            convergence,
            generator,
            model,
            polytope,
            n_sample, rng, setting
        )
        if iszero(n_topo_batch)
            generator = TopologyPerturbationGenerator(
                model=model,
                rng=rng,
                setting=setting
            )
        end
        update_pdtot_distributions!(distributions["load"], convergence)
    end
end

"Generate OPF dataset by sampling load instances uniformly in total load active power"
function generate_dataset(path::String; filename::String = "settings.yaml")

    t1 = time()
    cd(path)

    setting = read_settings(filename);
    network = instantiate_network(setting);
    setting = update_settings(setting);
    rng = read_rng(setting)
    init_workers!(setting)

    model = instantiate_model(network, "ACP", setting)
    if setting.CASE.append
        polytope = import_polytope()
    else
        polytope = instantiate_polytope(model);
        export_polytope(polytope)
    end
    generate_opf_instances(model, polytope, rng, setting)
    
    save_rng(rng)
    if setting.CASE.uid
        generate_uid()
    end
    _DC.nprocs() > 1 && _DC.rmprocs(_DC.workers());
    println("Total elapsed time: $(round((time() - t1)/3600, digits=2)) [h]")
    return nothing
end
