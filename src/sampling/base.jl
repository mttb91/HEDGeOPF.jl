
"Generate OPF instances uniformly in total load active power for one power system configuration"
function generate_opf_instances(model::_PM.AbstractPowerModel, polytope::PolyType, n_samples::Int,
    config::String, rng::_RND.AbstractRNG, settings::NamedTuple
)   
    n_items = settings.CASE.num_items
    counter = [n_samples, 0, 0];
    convergence = Matrix{Float64}[]
    database = instantiate_database(model, settings.MODEL.duals);
    # Generate OPF instances in batches, monitoring convergence
    for num in estimate_count(counter, settings.CASE.num_batches, n_items)
        
        if isempty(convergence)
            dist = dist_uniform(polytope, settings);
        else
            dist = dist_nonparametric(convergence);
        end
        samples = generate_input_samples(dist, polytope, (num, n_items), rng, settings);
        push!(convergence, generate_opf_samples!(counter, model, database, samples, config));
    end
    # Generate one additional batch if target number of feasible samples is not reached
    num = first(estimate_count(counter, 1, n_items))
    if num > 0
        dist = dist_nonparametric(convergence);
        samples = generate_input_samples(dist, polytope, (num, n_items), rng, settings);
        push!(convergence, generate_opf_samples!(counter, model, database, samples, config));
    end
    return nothing
end

"Generate OPF dataset by sampling load instances uniformly in total load active power"
function generate_dataset(path::String; filename::String = "settings.yaml")

    t1 = time()
    cd(path)

    setting = read_settings(filename);
    network = instantiate_network(setting);
    setting = update_settings(setting);
    init_workers!(setting)

    # Generate model for base power system configuration
    models = Dict{String, Any}("C0" => instantiate_model(network, "ACP", setting));
    n_samples = ceil(Int, setting.CASE.num_samples / length(models))

    # Generate OPF instance for each power system configuration (only one for now)
    rng = read_rng(setting)
    for (config, model) in models
        if config == "C0"
            if setting.CASE.append
                polytope = import_polytope()
            else
                polytope = instantiate_polytope(model);
                export_polytope(polytope)
            end
        end
        generate_opf_instances(model, polytope, n_samples, config, rng, setting)
        export_graph(model, config)
    end
    
    save_rng(rng)
    if setting.CASE.uid
        generate_uid()
    end
    _DC.nprocs() > 1 && _DC.rmprocs(_DC.workers());
    println("Total elapsed time: $(round((time() - t1)/3600, digits=2)) [h]")
    return nothing
end
