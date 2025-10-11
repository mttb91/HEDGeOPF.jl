
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

"Return the reference bus for a subgraph of model `pm` defined by `buses`, defining it if missing"
function define_ref_bus(pm::_PM.AbstractPowerModel, buses::Vector{Int}, bus_gen::Vector{Int})

    ref, vars = :bus, ["bus_i", "bus_type"]
    data = get_pm_value(pm, ref, vars, Array{Any, 2}; mask=buses)
    if all(data[:, 2] .!= 3)

        # Select a generator with both active and reactive support as reference bus if possible
        bus_gen = bus_gen[in.(bus_gen, Ref(buses))]
        if isempty(bus_gen)
            bus_ref = first(buses)
        else
            bus_ref = first(bus_gen)
        end
    else
        bus_ref = first(data[data[:, 2] .=== 3, 1])
    end
    return bus_ref
end

"Update in-place the `ref` dictionary of a PowerModel model `pm` based on `topology` perturbation"
function update_topology!(pm::_PM.AbstractPowerModel, topology::TopologyPerturbation)

    t = topology
    pm = deepcopy(pm)
    !isempty(t.ids_ref) && set_pm_value!(pm, :bus, ["bus_type"], 3; mask = t.ids_ref)
    !isempty(t.ids_bus) && set_pm_value!(pm, :bus, ["is_connected"], 0; mask = topology.ids_bus)
    !isempty(t.ids_gen) && set_pm_value!(pm, :gen, ["gen_status"], 0; mask = topology.ids_gen)
    !isempty(t.ids_branch) && set_pm_value!(pm, :branch, ["br_status"], 0; mask = topology.ids_branch)
    !isempty(t.ids_gen_faulted) && set_pm_value!(pm, :gen, ["pmin", "pmax", "qmin", "qmax"], 0.0; mask = topology.ids_gen_faulted)
    if !isempty(t.ids_bus)
        for key in ["load", "shunt"]
            if isempty(_PM.ref(pm, Symbol(key)))
                continue
            end
            bus = get_pm_value(pm, Symbol(key), [key * "_bus"], Array{Any, 2})
            ids = findall(x -> in(x, t.ids_bus), vec(bus))
            !isempty(ids) && set_pm_value!(pm, Symbol(key), ["status"], 0; mask = ids)
        end
    end
    return nothing
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
    fix_theta_ref!(pm::_PM.AbstractPowerModel, ids_fixed::Vector{Int};
        nw::Int=_PM.nw_id_default
    )

Update in-place the reference bus(es) of model `pm` by:
- unfixing the voltage angle for the references buses that are not in `ids_fixed`
- fixing to zero the voltage angle for the new reference buses among `ids_fixed`

"""
function fix_theta_ref!(pm::_PM.AbstractPowerModel, ids_fixed::Vector{Int};
    nw::Int=_PM.nw_id_default
)
    var = :va
    ref = _PM.var(pm, nw, var)
    ids = sort(first(ref.axes))
    ids_fixed_old = ids[JuMP.is_fixed.(ref[ids].data)]
    
    if !isequal(ids_fixed, ids_fixed_old)

        JuMP.unfix.(_PM.var(pm, nw, var, setdiff(ids_fixed_old, ids_fixed)))
        JuMP.fix.(_PM.var(pm, nw, var, setdiff(ids_fixed, ids_fixed_old)), 0.0; force = true)
    end
    return nothing
end


"""
    fix_to_zero!(pm::_PM.AbstractPowerModel, ids_fixed::Vector{Int}, element::Symbol, var::Symbol, props::Vector{String};
        nw::Int=_PM.nw_id_default
    )

Fix to zero the variable `var` of model `pm` at indices `ids_fixed`, re-defining bounds for initially fixed variables.

## Notes

Bounds values are derived from the `ref` dictionary of the PowerModels model, at `props` of `element`.
It is assumed that variable `var` is bounded at both sides. If a single bound is provided under `props`, the
variable is considered as symmetrically bounded.
"""
function fix_to_zero!(pm::_PM.AbstractPowerModel, ids_fixed::Vector{Int}, element::Symbol, var::Symbol, props::Vector{String};
    nw::Int=_PM.nw_id_default
)
    ref = _PM.var(pm, nw, var)
    ids = sort(first(ref.axes))
    if isa(first(ids), Tuple)
        ids_fixed = ids_fixed * 2
        ids_fixed = sort(vcat(ids_fixed .- 1, ids_fixed))
    end
    ids_fixed_new = ids[ids_fixed]
    ids_fixed_old = ids[JuMP.is_fixed.(ref[ids].data)]
    
    if !isequal(ids_fixed_new, ids_fixed_old)

        mask = setdiff(ids_fixed_old, ids_fixed_new)
        if !isempty(mask)
            if isa(first(ids), Tuple)
                bounds = get_pm_value(pm, element, props, Array{Any, 2}; mask=first.(mask)[1:2:end])
                bounds = repeat(bounds; inner=(2, 1))
            else
                bounds = get_pm_value(pm, element, props, Array{Any, 2}; mask=mask)
            end

            JuMP.unfix.(_PM.var(pm, nw, var, mask))
            if length(props) === 1
                JuMP.set_lower_bound.(_PM.var(pm, nw, var, mask), -bounds[:, 1])
                JuMP.set_upper_bound.(_PM.var(pm, nw, var, mask), bounds[:, 1])
            elseif length(props) === 2
                @assert any(diff(bounds, dims=2) .>= 0)
                JuMP.set_lower_bound.(_PM.var(pm, nw, var, mask), bounds[:, 1])
                JuMP.set_upper_bound.(_PM.var(pm, nw, var, mask), bounds[:, 2])
            end
        end
        mask = setdiff(ids_fixed_new, ids_fixed_old)
        JuMP.fix.(_PM.var(pm, nw, var, mask), 0.0; force = true)
    end    
    return nothing
end


"""
    update_branch_status!(pm::_PM.AbstractPowerModel, ids_faulted::Vector{Int};
        nw::Int=_PM.nw_id_default
    )

Set to zero the multiplicative parameter in the active/reactive power flow equation at branches `ids_faulted` .
"""
function update_branch_status!(pm::_PM.AbstractPowerModel, ids_faulted::Vector{Int};
    nw::Int=_PM.nw_id_default
)

    var = :br_status
    ref = _PM.var(pm, nw, var)
    ids = sort(first(ref.axes))
    ids_faulted_old = ids[iszero.(JuMP.parameter_value.(ref[ids].data))]

    is_changed = false
    if !isequal(ids_faulted, ids_faulted_old)
        mask = setdiff(ids_faulted_old, ids_faulted)
        JuMP.set_parameter_value.(_PM.var(pm, nw, var, mask), 1.0)
        mask = setdiff(ids_faulted, ids_faulted_old)
        JuMP.set_parameter_value.(_PM.var(pm, nw, var, ids_faulted), 0.0)
        is_changed = true
    end
    return is_changed
end


"""
    update_topology!(pm::_PM.AbstractPowerModel, ids_branch::Vector{Int}, ids_bus::Vector{Int}, ids_ref::Vector{Int}, ids_gen::Vector{Int})

Modify in-place the topology of model `pm` by fixing to zero:
- the voltage angle of new reference buses specified by `ids_ref`
- the voltage magnitude of buses `ids_bus`
- the value of power flow equation at indices `ids_branch`
- the active and reactive power of generators at indices `ids_gen`

## Notes

Variable found initially fixed are unfixed and bounded with original limits.
No variable is removed from model `pm`.

"""
function update_topology!(pm::_PM.AbstractPowerModel, ids_branch::Vector{Int}, ids_bus::Vector{Int}, ids_ref::Vector{Int}, ids_gen::Vector{Int})

    is_changed = update_branch_status!(pm, ids_branch)
    if is_changed
        fix_theta_ref!(pm, ids_ref)
        fix_to_zero!(pm, ids_bus, :bus, :vm, ["vmin", "vmax"])
    end
    fix_to_zero!(pm, ids_gen, :gen, :pg, ["pmin", "pmax"])
    fix_to_zero!(pm, ids_gen, :gen, :qg, ["qmin", "qmax"])

    return nothing
end


"""
    update_load_power!(pm::_PM.AbstractPowerModel, sample::Dict{String, Dict{String, InputData}}, pd_tot::Vector{Float64};
        nw::Int=_PM.nw_id_default
    )

Update in-place model `pm` by:
- fixing the input active and reactive load variables `pd_fix` and `qd_fix` to the new setpoints in `samples` 
- replace the lower and upper bounds in the constraint on total load active power

"""
function update_load_power!(pm::_PM.AbstractPowerModel, sample::Dict{String, Dict{String, InputData}}, pd_tot::Vector{Float64};
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
    JuMP.set_normalized_rhs.(pm.model[:pd_tot], pd_tot .* [-1, 1])
    return nothing
end


"""
    update_gen_cost!(model::_PM.AbstractPowerModel, sample::Dict{String, <:Any})

Update in-place the linear and quadratic generator cost coefficients in the OPF `model`
objective function based on the input `sample`.

"""
function update_gen_cost!(pm::_PM.AbstractPowerModel, sample::Dict{String, Dict{String, InputData}};
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
function update_model!(pm::_PM.AbstractPowerModel, sample::InputSample)

    t = sample.topology
    update_topology!(pm, t.ids_branch, t.ids_bus, t.ids_ref, sort(vcat(t.ids_gen, t.ids_gen_faulted)))

    if haskey(sample.data, "load")
        update_load_power!(pm, sample.data, sample.pd_tot)
    end
    if haskey(sample.data, "gen")
        update_gen_cost!(pm, sample.data)
    end

    for (key, value) in sample.data
        isempty(value) && delete!(sample.data, key)
    end
    return nothing
end
