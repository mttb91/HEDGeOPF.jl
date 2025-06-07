
 "Create a named reference empty container in `PowerModels` constraint dictionary"
function _reference_constraint(pm::_PM.AbstractPowerModel, n::Int, component::Symbol, name::Symbol; cols::Union{Int, Nothing} = nothing)
 
    if !haskey(_PM.con(pm, n), name)
        dims = []
        push!(dims, _PM.ids(pm, n, component))
        !isnothing(cols) && push!(dims, 1:cols)
        _PM.con(pm, n)[name] = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(undef, dims...)
    end
    return nothing
end

function constraint_load_fixed_and_slack_active_power(pm::_PM.AbstractPowerModel, n::Int, i::Int, pmin::Float64, pmax::Float64)

    pd_fix = _PM.var(pm, n, :pd_fix, i)
    pd_slack_up = _PM.var(pm, n, :pd_slack_up, i)
    pd_slack_down = _PM.var(pm, n, :pd_slack_down, i)

    if !all(JuMP.is_fixed.([pd_fix, pd_slack_up, pd_slack_down]))
        
        pd = pd_fix + pd_slack_down - pd_slack_up
        JuMP.@constraint(pm.model, pmin <= pd <= pmax)
    end
    return nothing
end

function constraint_load_total_fixed_and_slack_active_power(pm::_PM.AbstractPowerModel, n::Int, pg_min::Float64, pg_max::Float64)

    pd = _PM.var(pm, n, :pd_fix) .+_PM.var(pm, n, :pd_slack_down) .-_PM.var(pm, n, :pd_slack_up)

    JuMP.@constraint(pm.model, pd_tot[i = 1:2], [-sum(pd), sum(pd)][i] <= [-pg_min, pg_max][i])
end

function expression_load_power(pm::_PM.AbstractPowerModel; kwargs...)

    expression_load_power_real(pm; kwargs...)
    expression_load_power_imaginary(pm; kwargs...)
end

function expression_load_power_real(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, report::Bool=true)

    pd = _PM.var(pm, nw)[:pd] = Dict{Int, Any}()

    for (i, load) in _PM.ref(pm, nw, :load)
        pd[i] = JuMP.@expression(pm.model, _PM.var(pm, nw, :pd_fix, i) + _PM.var(pm, nw, :pd_slack_down, i) - _PM.var(pm, nw, :pd_slack_up, i))
    end
    report && _PM.sol_component_value(pm, nw, :load, :pd, _PM.ids(pm, nw, :load), pd)
end

function expression_load_power_imaginary(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, report::Bool=true)

    qd = _PM.var(pm, nw)[:qd] = Dict{Int, Any}()

    for (i, load) in _PM.ref(pm, nw, :load)
        qd[i] = JuMP.@expression(pm.model, _PM.var(pm, nw, :qd_fix, i) + _PM.var(pm, nw, :qd_slack_down, i) - _PM.var(pm, nw, :qd_slack_up, i))
    end
    report && _PM.sol_component_value(pm, nw, :load, :qd, _PM.ids(pm, nw, :load), qd)
end