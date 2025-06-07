
function _build_opf_slack(pm::_PM.AbstractPowerModel)

    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)    

    variable_load_power(pm)

    objective_min_fuel_and_slack_cost(pm)

    _PM.constraint_model_voltage(pm)

    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in _PM.ids(pm, :bus)
        constraint_power_balance_slack(pm, i)
    end

    for i in _PM.ids(pm, :branch)
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        _PM.constraint_voltage_angle_difference(pm, i)
    end

    for i in _PM.ids(pm, :dcline)
        _PM.constraint_dcline_power_losses(pm, i)
    end

    for i in _PM.ids(pm, :load)
        constraint_load_fixed_and_slack_active_power(pm, i)
        constraint_load_fixed_and_slack_reactive_power(pm, i)
    end
    constraint_load_total_fixed_and_slack_active_power(pm)
end

"Build optimal power flow model with active and reactive load slack variables"
function build_opf_slack(pm::_PM.AbstractPowerModel)

    _build_opf_slack(pm)

    for i in _PM.ids(pm, :branch)
        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end
end

"Build optimal power flow model with active and reactive load slack variable and dual values reporting"
function build_opf_slack_dual(pm::_PM.AbstractPowerModel)

    _build_opf_slack(pm)
    variable_branch_power_apparent_squared(pm)

    for i in _PM.ids(pm, :branch)
        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end
end
