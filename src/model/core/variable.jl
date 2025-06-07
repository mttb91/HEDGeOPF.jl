
function variable_load_power(pm::_PM.AbstractPowerModel; kwargs...)
    variable_load_fixed_power(pm; kwargs...)
    variable_load_slack_power(pm; kwargs...)
    expression_load_power(pm; kwargs...)
end

function variable_load_fixed_power(pm::_PM.AbstractPowerModel; kwargs...)
    variable_load_fixed_power_real(pm; kwargs...)
    variable_load_fixed_power_imaginary(pm; kwargs...)
end

function variable_load_fixed_power_real(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, report::Bool=false)

    var = _PM.var(pm, nw)[:pd_fix] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)],
        base_name="$(nw)_pd_fix",
    )

    for (i, load) in _PM.ref(pm, nw, :load)
        JuMP.fix(var[i], load["pd"])
    end

    report && _PM.sol_component_value(pm, nw, :load, :pd_fix, _PM.ids(pm, nw, :load), var)
end

function variable_load_fixed_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, report::Bool=false)

    var = _PM.var(pm, nw)[:qd_fix] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)],
        base_name="$(nw)_qd_fix"
    )

    for (i, load) in _PM.ref(pm, nw, :load)
        JuMP.fix(var[i], load["qd"])
    end

    report && _PM.sol_component_value(pm, nw, :load, :qd_fix, _PM.ids(pm, nw, :load), var)
end

"variables for active and reactive load power slack"
function variable_load_slack_power(pm::_PM.AbstractPowerModel; kwargs...)
    variable_load_slack_power_real(pm; kwargs...)
    variable_load_slack_power_imaginary(pm; kwargs...)
end

function variable_load_slack_power_real(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, report::Bool=true)

    var1 = _PM.var(pm, nw)[:pd_slack_up] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_pd_slack_up",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "pd_slack_up_start", 0.0)
    )
    var2 = _PM.var(pm, nw)[:pd_slack_down] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_pd_slack_down",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "pd_slack_down_start", 0.0)
    )

    baseMVA = _PM.ref(pm, :baseMVA)
    for (i, load) in _PM.ref(pm, nw, :load)
        if (load["pmax"] - load["pmin"]) <= 1E-3 / baseMVA
            JuMP.fix(var1[i], 0.0)
            JuMP.fix(var2[i], 0.0)
        else
            JuMP.set_lower_bound(var1[i], 0.0)
            JuMP.set_lower_bound(var2[i], 0.0)
        end
    end

    report && _PM.sol_component_value(pm, nw, :load, :pd_slack_up, _PM.ids(pm, nw, :load), var1)
    report && _PM.sol_component_value(pm, nw, :load, :pd_slack_down, _PM.ids(pm, nw, :load), var2)
end

function variable_load_slack_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, report::Bool=true)

    var1 = _PM.var(pm, nw)[:qd_slack_up] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_qd_slack_up",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "q_slack_up_start", 0.0)
    )
    var2 = _PM.var(pm, nw)[:qd_slack_down] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_qd_slack_down",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "qd_slack_down_start", 0.0)
    )

    baseMVA = _PM.ref(pm, :baseMVA)
    for (i, load) in _PM.ref(pm, nw, :load)
        if (load["qmax"] - load["qmin"]) <= 1E-3 / baseMVA
            JuMP.fix(var1[i], 0.0)
            JuMP.fix(var2[i], 0.0)
        else
            JuMP.set_lower_bound(var1[i], 0.0)
            JuMP.set_lower_bound(var2[i], 0.0)
        end
    end

    report && _PM.sol_component_value(pm, nw, :load, :qd_slack_up, _PM.ids(pm, nw, :load), var1)
    report && _PM.sol_component_value(pm, nw, :load, :qd_slack_down, _PM.ids(pm, nw, :load), var2)
end

function variable_branch_power_apparent_squared(pm::_PM.AbstractACPModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)

    var = _PM.var(pm, nw)[:s] = JuMP.@variable(pm.model,
        [(l,i,j) in _PM.ref(pm, nw, :arcs)],
        base_name="$(nw)_s",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "s_start")
    )

    if bounded
        for (l, branch) in _PM.ref(pm, nw, :branch)
            idx = (l, branch["f_bus"], branch["t_bus"])
            JuMP.set_upper_bound(var[idx], branch["rate_a"]^2)
            idx = (l, branch["t_bus"], branch["f_bus"])
            JuMP.set_upper_bound(var[idx], branch["rate_a"]^2)
        end
    end
    report && _PM.sol_component_value_edge(pm, nw, :branch, :sf, :st, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), var)
end
