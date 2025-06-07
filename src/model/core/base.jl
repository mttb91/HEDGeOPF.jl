
"Instantiate and pre-process power network dictionary in `PowerModels` format."
function instantiate_network(settings::Dict{String, <:Any})

    # Generate the powermodels network
    network = _PM.parse_file(joinpath(pwd(), settings["PATH"]["input"], settings["CASE"]["grid"]))
    # Reset bus index
    update_bus_ids!(network)
    # Add active and reactive load power bounds
    set_pm_value!(network, "load", compute_load_power_bounds(network, settings)...)
    # Add curtailment cost for loads
    set_pm_value!(network, "load", ["curt_cost"], settings["MODEL"]["voll"] * network["baseMVA"])

    return network
end

"Instantiate PowerModels model with optimizer defined by `settings`"
function instantiate_model(network::Dict{String, <:Any}, type::String, settings::NamedTuple)

    build_method = settings.MODEL.duals ? build_opf_slack_dual : build_opf_slack
    type = getfield(_PM, Symbol("$(type)PowerModel"))
    model = _PM.instantiate_model(network, type, build_method; setting = Dict("output" => Dict()))
    # Specify solver options
    solver = getfield(getfield(Main, Symbol(settings.SOLVER.nlp)), :Optimizer)
    solver = _PM.optimizer_with_attributes(solver,
        [string(k) => v for (k, v) in pairs(settings.SOLVER.nlp_options)]...
    )
    JuMP.set_optimizer(model.model, solver)

    return model
end

"Solve the OPF problem for a given input sample, recording input/output variables"
function solve_model!(pm::_PM.AbstractPowerModel, db::NamedTuple, sample::Dict{String, <:Any})

    # Update and solve the optimisation model
    update_model!(pm, sample)
    results = _PM.optimize_model!(pm)
    results["input"] = sample
    # Extract relevant data to the database
    extract_data!(db, results, pm)

    results = nothing
    return extract_info(pm)
end
