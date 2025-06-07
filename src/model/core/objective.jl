
"Minimize sum of active power generation cost and active and reactive power load slack (upward and downward)"
function objective_min_fuel_and_slack_cost(pm::_PM.AbstractPowerModel; kwargs...)

    expression_pg_cost(pm; kwargs...)
    expression_slack_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum(_PM.var(pm, n, :pg_cost, i) for (i, gen) in nw_ref[:gen])
            + sum(_PM.var(pm, n, :slack_cost, i) for (i, bus) in nw_ref[:load])
        for (n, nw_ref) in _PM.nws(pm)
        )
    )
end

"Formulation-agnostic expression for generator active power cost"
function expression_pg_cost(pm::_PM.AbstractPowerModel)

    for (n, nw_ref) in _PM.nws(pm)
        pg_cost = _PM.var(pm, n)[:pg_cost] = Dict{Int, Any}()

        for (i, gen) in nw_ref[:gen]
            x = _PM.var(pm, n, :pg, i)

            if gen["model"] == 2
                cost_terms = reverse(gen["cost"])

                if length(cost_terms) == 0
                    pg_cost[i] = 0.0
                elseif length(cost_terms) == 1
                    pg_cost[i] = cost_terms[1]
                elseif length(cost_terms) == 2
                    pg_cost[i] = JuMP.@expression(pm.model, cost_terms[1] + cost_terms[2]*x)
                elseif length(cost_terms) == 3
                    pg_cost[i] = JuMP.@expression(pm.model, cost_terms[1] + cost_terms[2]*x + cost_terms[3]*x^2)
                end
            else
                Memento.error(_LOGGER, "Only cost model of type 2 is supported at this time for generator $i")
            end
        end
    end
end

"Formulation-agnostic expression for nodal slack power cost"
function expression_slack_cost(pm::_PM.AbstractPowerModel; ratio::Float64 = 0.2)

    for (n, nw_ref) in _PM.nws(pm)
        cost = _PM.var(pm, n)[:slack_cost] = Dict{Int, Any}()

        for (i, load) in nw_ref[:load]
            pd_slack_up   = _PM.var(pm, n, :pd_slack_up, i)
            qd_slack_up   = _PM.var(pm, n, :qd_slack_up, i)
            pd_slack_down = _PM.var(pm, n, :pd_slack_down, i)
            qd_slack_down = _PM.var(pm, n, :qd_slack_down, i)
            
            cost[i] = JuMP.@expression(pm.model, (pd_slack_up + (pd_slack_down + qd_slack_up + qd_slack_down) * ratio) .* load["curt_cost"])
        end
    end
end
