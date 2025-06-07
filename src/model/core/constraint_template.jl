
"Formulation-agnostic constraint on active and reactive nodal power balance"
function constraint_power_balance_slack(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)

    bus = _PM.ref(pm, nw, :bus, i)
    bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = _PM.ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_sw = _PM.ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)
    bus_storage = _PM.ref(pm, nw, :bus_storage, i)

    bus_pd = Dict(k => _PM.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PM.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    return constraint_power_balance_slack(pm, nw, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
end

"Formulation-agnostic constraint on reactive power setpoint of load with slack"
function constraint_load_fixed_and_slack_active_power(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)

    pmin = _PM.ref(pm, nw, :load, i, "pmin")
    pmax = _PM.ref(pm, nw, :load, i, "pmax")

    return constraint_load_fixed_and_slack_active_power(pm, nw, i, pmin, pmax)
end

"Formulation-agnostic constraint on reactive power setpoint of load with slack"
function constraint_load_fixed_and_slack_reactive_power(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)

    qmin = _PM.ref(pm, nw, :load, i, "qmin")
    qmax = _PM.ref(pm, nw, :load, i, "qmax")
    qp_ratio_min = _PM.ref(pm, nw, :load, i, "qp_ratio_min")
    qp_ratio_max = _PM.ref(pm, nw, :load, i, "qp_ratio_max")

    return constraint_load_fixed_and_slack_reactive_power(pm, nw, i, qmin, qmax, qp_ratio_min, qp_ratio_max)
end

"Formulation-agnostic constraint on load total active power"
function constraint_load_total_fixed_and_slack_active_power(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default)

    pg_min = sum(k["pmin"] for k in values(_PM.ref(pm, nw, :gen)))
    pg_max = sum(k["pmax"] for k in values(_PM.ref(pm, nw, :gen)))

    return constraint_load_total_fixed_and_slack_active_power(pm, nw, pg_min, pg_max)
end

function constraint_thermal_limit_from(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    constraint_thermal_limit_from(pm, nw, f_idx)
end

function constraint_thermal_limit_to(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    constraint_thermal_limit_to(pm, nw, t_idx)
end