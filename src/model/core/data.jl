
###############################################################################
# Network modification function
###############################################################################

function update_bus_ids!(network::Dict{String, <:Any})

    # Mapping from original sorted bus index to new one
    ids_bus = get_pm_value(network, "bus", ["index"], Array{Any, 2})
    mapping = Dict(vec(ids_bus) .=> axes(ids_bus, 1))
    
    _PM.update_bus_ids!(network, mapping)
    return nothing
end

function compute_load_power_bounds(network::Dict{String, <:Any}, settings::Dict{String, <:Any})
    setting = settings["SAMPLING"]
    pd_vars, pd_bounds = compute_load_active_bounds(network, setting)
    qd_vars, qd_bounds = compute_load_reactive_bounds(network, setting)

    return vcat(pd_vars, qd_vars), hcat(pd_bounds, qd_bounds)
end

function compute_load_active_bounds(network::Dict{String, <:Any}, settings::Dict{String, <:Any})

    ref, vars, vars_new = "load", ["pd"], ["pmin", "pmax"]

    # Get relevant settings
    delta_pd = settings["delta_pd"] / 100
    # Compute minimum and maximum load active power bounds
    data = repeat(get_pm_value(network, ref, vars, Array{Any, 2}), 1, 2) .* [1 - delta_pd 1 + delta_pd]
    return vars_new, sort(data, dims=2)
end

function compute_load_reactive_bounds(network::Dict{String, <:Any}, settings::Dict{String, <:Any})

    ref, vars, vars_new = "load", ["pd", "qd"], ["qmin", "qmax", "qp_ratio_min", "qp_ratio_max"]

    # Get settings
    min_pf = settings["min_pf"]
    max_pf = settings["max_pf"]
    delta_pf = settings["delta_pf"]
    delta_qd = settings["delta_qd"] / 100
    # Get the base loading setpoints
    data = get_pm_value(network, ref, vars, _DF.DataFrame)
    # Compute minimum and maximum load reactive power bounds
    bounds = sort(repeat(data.qd, 1, 2) .* [1 - delta_qd 1 + delta_qd], dims=2)
    # Compute minimum power factor reduced by pre-defined threshold accounting for reactive power sign
    _DF.transform!(data, vars => ((x, y) -> cos.(atan.(./(y, x)))) => "pf")
    _DF.transform!(data,
        last(_DF.names(data), 2)
        => ((x, y) -> copysign.(maximum.(vcat.(y .- delta_pf, min_pf)), x))
        => "pf_min"
    )
    # Compute reactive/active power minimum and maximum ratios accounting for reactive power sign
    ratio_min = copysign.(tan.(acos.(maximum.(vcat.(data.pf, max_pf)))), data.pf_min)
    ratios = sort(reduce(hcat, [ratio_min, tan.(acos.(data.pf_min))]), dims=2)

    return vars_new, hcat(bounds, ratios)
end

###############################################################################
# Optimisation model updates
###############################################################################

"Extract from `model` specific convergence information required for dataset generation"
function extract_info(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default)

    feasible = JuMP.is_solved_and_feasible(pm.model)
    if feasible
        pd_tot = sum(JuMP.value.(values(_PM.var(pm, nw, :pd))))
    else
        pd_tot = sum(JuMP.fix_value.(_PM.var(pm, nw, :pd_fix)).data)
    end
    return [pd_tot, Float64(feasible)]
end


"""
    update_load!(pm::_PM.AbstractPowerModel, sample::Dict{String, <:Any}, bounds::AbstractVector{Float64};
        nw::Int=_PM.nw_id_default
    )

Update in-place the OPF `model` by:
- fixing the input active and reactive load variables `pd_fix` and `qd_fix` to the new setpoints in `samples` 
- replace the lower and upper bounds in the constraint on total load active power

## Notes

The 
The load `sample` may contain values for a subset of the loads in the power system. Once the load variables
are fixed to the new setpoints, `sample` is updated to record the active and reactive power values of every
load in the system.

"""
function update_load_power!(pm::_PM.AbstractPowerModel, sample::Dict{String, <:Any};
    nw::Int=_PM.nw_id_default
)
    ref = "load"
    # Update active and reactive load power setpoints
    for var in ["pd", "qd"]
        label = Symbol(var * "_fix")
        value = sample[ref][var]
        JuMP.fix.(_PM.var(pm, nw, label)[value.ids], value.data)
        delete!(sample[ref], var)
    end
    # Update total load active power bounds
    JuMP.set_normalized_rhs.(pm.model[:pd_tot], pop!(sample, "info") .* [-1, 1])
    return nothing
end


"""
    update_gen_cost!(model::_PM.AbstractPowerModel, sample::Dict{String, <:Any})

Update in-place the linear and quadratic generator cost coefficients in the OPF `model`
objective function based on the input `sample`.

"""
function update_gen_cost!(pm::_PM.AbstractPowerModel, sample::Dict{String, <:Any};
    nw::Int=_PM.nw_id_default
)
    ref, var = "gen", :pg
    ids = sort(first(_PM.var(pm, nw, var).axes))
    for key in filter(x -> in(x, ["c1", "c2"]), keys(sample[ref]))
        value = sample[ref][key]
        mask = isnothing(value.ids) ? ids : value.ids
        # Linear or quadratic coefficients
        vars = Vector{JuMP.VariableRef}[]
        if key == "c1"
            push!(vars, (_PM.var(pm, nw, var)[mask]).data)
        elseif key == "c2"
            append!(vars, [(_PM.var(pm, nw, var)[mask]).data, (_PM.var(pm, nw, var)[mask]).data])
        end
        !isempty(vars) && JuMP.set_objective_coefficient.(pm.model, vars..., value.data)
    end
    return nothing
end


"Wrapper function to update PowerModels optimisation `model` based on the given input `sample`"
function update_model!(pm::_PM.AbstractPowerModel, sample::Dict{String, <:Any})

    if haskey(sample, "load")
        update_load_power!(pm, sample)
    end
    if haskey(sample, "gen")
        update_gen_cost!(pm, sample)
    end

    for (key, value) in sample
        isempty(value) && delete!(sample, key)
    end
    return nothing
end
