
function constraint_power_balance_slack(pm::_PM.AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = _PM.var(pm, n, :vm, i)
    p    = _PM.get(_PM.var(pm, n),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = _PM.get(_PM.var(pm, n),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = _PM.get(_PM.var(pm, n),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = _PM.get(_PM.var(pm, n),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = _PM.get(_PM.var(pm, n),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = _PM.get(_PM.var(pm, n),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = _PM.get(_PM.var(pm, n),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = _PM.get(_PM.var(pm, n),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = _PM.get(_PM.var(pm, n), :p_dc, Dict()); _PM._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = _PM.get(_PM.var(pm, n), :q_dc, Dict()); _PM._check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    pd_fix        = _PM.get(_PM.var(pm, n), :pd_fix, Dict()); _PM._check_var_keys(pd_fix, bus_pd, "active power", "load")
    qd_fix        = _PM.get(_PM.var(pm, n), :qd_fix, Dict()); _PM._check_var_keys(qd_fix, bus_pd, "reactive power", "load")
    pd_slack_up   = _PM.get(_PM.var(pm, n), :pd_slack_up, Dict()); _PM._check_var_keys(pd_slack_up, bus_pd, "active power slack", "load")
    qd_slack_up   = _PM.get(_PM.var(pm, n), :qd_slack_up, Dict()); _PM._check_var_keys(qd_slack_up, bus_qd, "reactive power slack", "load")
    pd_slack_down = _PM.get(_PM.var(pm, n), :pd_slack_down, Dict()); _PM._check_var_keys(pd_slack_down, bus_pd, "active power slack", "load")
    qd_slack_down = _PM.get(_PM.var(pm, n), :qd_slack_down, Dict()); _PM._check_var_keys(qd_slack_down, bus_qd, "reactive power slack", "load")

    JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd_fix[d] - pd_slack_up[d] + pd_slack_down[d] for (d,pd) in bus_pd)
        - sum(gs for (i,gs) in bus_gs)*vm^2
    )

    JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd_fix[d] - qd_slack_up[d] + qd_slack_down[d] for (d,qd) in bus_qd)
        + sum(bs for (i,bs) in bus_bs)*vm^2
    )

    return nothing
end

function constraint_load_fixed_and_slack_reactive_power(pm::_PM.AbstractACPModel, n::Int, i::Int, qmin::Float64, qmax::Float64, qp_ratio_min::Float64, qp_ratio_max::Float64)

    pd_fix = _PM.var(pm, n, :pd_fix, i)
    qd_fix = _PM.var(pm, n, :qd_fix, i)
    pd_slack_up = _PM.var(pm, n, :pd_slack_up, i)
    qd_slack_up = _PM.var(pm, n, :qd_slack_up, i)
    pd_slack_down = _PM.var(pm, n, :pd_slack_down, i)
    qd_slack_down = _PM.var(pm, n, :qd_slack_down, i)

    if !all(JuMP.is_fixed.([qd_fix, qd_slack_up, qd_slack_down]))

        qd = qd_fix + qd_slack_down - qd_slack_up
        JuMP.@constraint(pm.model, qmin <= qd <= qmax)

        if !all(JuMP.is_fixed.([pd_fix, pd_slack_up, pd_slack_down]))
            pd = pd_fix + pd_slack_down - pd_slack_up
            if JuMP.fix_value(pd_fix) < 0
                pd *= -1.0
            end
            JuMP.@constraint(pm.model, qd <= pd * qp_ratio_max)
            JuMP.@constraint(pm.model, qd >= pd * qp_ratio_min)
        end
    end
    return nothing
end

function constraint_thermal_limit_from(pm::_PM.AbstractACPModel, n::Int, f_idx)
    p_fr = _PM.var(pm, n, :p, f_idx)
    q_fr = _PM.var(pm, n, :q, f_idx)
    s_fr = _PM.var(pm, n, :s, f_idx)

    JuMP.@constraint(pm.model, s_fr == p_fr^2 + q_fr^2)
end

function constraint_thermal_limit_to(pm::_PM.AbstractACPModel, n::Int, t_idx)
    p_to = _PM.var(pm, n, :p, t_idx)
    q_to = _PM.var(pm, n, :q, t_idx)
    s_to = _PM.var(pm, n, :s, t_idx)

    JuMP.@constraint(pm.model, s_to == p_to^2 + q_to^2)
end