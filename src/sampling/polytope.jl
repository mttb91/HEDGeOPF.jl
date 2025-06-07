
###############################################################################
# LP optimisation problems on convex polytopes
###############################################################################

"Check is a point is within convex polytope"
function _isin(point::AbstractVector{Float64}, polytope::Tuple)
    A, b = polytope
    return all(A * point .<= b)
end

"Generic JuMP LP model"
function _lp_model(settings::NamedTuple)

    solver = getfield(getfield(Main, Symbol(settings.SOLVER.lp)), :Optimizer)
    solver = JuMP.optimizer_with_attributes(solver,
        [string(k) => v for (k, v) in pairs(settings.SOLVER.lp_options)]...
    )
    model = JuMP.Model(solver)
    JuMP.set_silent(model)
    return model
end

"LP mathematical model of polytope"
function _polytope_model(polytope::PolyType, settings::NamedTuple)

    A, b, ids = polytope

    model = _lp_model(settings)
    JuMP.@variable(model, x[i = axes(A, 2)])
    JuMP.@constraint(model, A * x .<= b)

    return model, x, length.(ids)
end

"Check if polytope intersection is empty"
function _isempty(polytope::PolyType, settings::NamedTuple)

    model, _, _ = _polytope_model(polytope, settings)
    JuMP.@objective(model, Min, 1)
    JuMP.optimize!(model);

    return !JuMP.is_solved_and_feasible(model; allow_local=false)
end

"Find upper/lower limits for active/reactive power of a given load in the polytope"
function load_power_bound(polytope::PolyType, idx::Int, settings::NamedTuple; var::String = "pd", upper::Bool = true)

    model, x, dims = _polytope_model(polytope, settings)
    # Objective
    idx += ifelse(contains(var, "pd"), 0, first(dims))
    if upper
        JuMP.@objective(model, Max, x[idx])
    else
        JuMP.@objective(model, Min, x[idx])
    end
    JuMP.optimize!(model);

    return JuMP.objective_value(model)
end

"Find upper/lower limit for total active/reactive load power in the polytope"
function load_power_bound(polytope::PolyType, settings::NamedTuple; var::String = "pd", upper::Bool = true)

    model, x, dims = _polytope_model(polytope, settings)
    # Objective
    idx = ifelse(contains(var, "pd"), 1:first(dims), first(dims) + 1:sum(dims))
    if upper
        JuMP.@objective(model, Max, sum(x[idx]))
    else
        JuMP.@objective(model, Min, sum(x[idx]))
    end
    JuMP.optimize!(model);

    return JuMP.objective_value(model)
end

"Model for chebyshev centre computation"
function chebyshev_model(A::AbstractMatrix{Float64}, b::AbstractVector{Float64}, settings::NamedTuple;
    epsilon::Float64=1E-6
)

    model = _lp_model(settings)
    JuMP.@variable(model, x[i = axes(A, 2)])
    JuMP.@variable(model, r >= epsilon)
    JuMP.@constraint(model, con, A * x + norm.(eachcol(permutedims(A))) * r .<= b)
    JuMP.@objective(model, Max, r)

    return model, x, con
end

"Find chebyshev centre of instantiated polytope model"
function chebyshev_centre(model, x)

    JuMP.optimize!(model)

    if JuMP.is_solved_and_feasible(model)
        centre = JuMP.value.(x)
    else
        centre = nothing
    end
    return centre 

end

###############################################################################
# Polytope generation
###############################################################################

"""
    instantiate_polytope(model::_PM.AbstractPowerModel)

Initialise a closed convex H-polytope in the form `Ax <= b`, defining matrices `A` and `b`.

## Notes

The following load constraints are accounted for each load in this exact order:

- maximum active power as `p <= pmax`
- minimum active power as `-p <= -pmin`
- maximum reactive power as `q <= qmax`
- maximum reactive power as `-q <= -qmin`
- maximum reactive power as `-p * rmax + q <= 0`
- minimum reactive power as `p * rmin - q <= 0`
- maximum total active power consumption as `sum(p) <= sum(pg_max)`
- minimum total active power consumption as `-sum(p) <= -sum(pg_min)`

Fixed loads (equal lower and upper bounds) are not accounted for in the
polytope definition. Similarly, the p-q coupling constraint is not
formulated for loads with decoupled active and reactive power.
"""
function instantiate_polytope(model::_PM.AbstractPowerModel)

    # Extract relevant network data
    vars = ["pmin", "pmax", "qmin", "qmax", "qp_ratio_min", "qp_ratio_max"]
    data = get_pm_value(model, :load, vars, _DF.DataFrame)
    # Keep only variables with nonzero range
    pd = filter(row -> !iszero(row.pmax - row.pmin), data)[!, first(vars, 2)]
    qd = filter(row -> !iszero(row.qmax - row.qmin), data)[!, last(first(vars, 4), 2)]
    pd_res = sum(ifelse.(iszero.(data.pmax - data.pmin), data.pmax, 0.0))
    # Masking for removal of coupling constraints
    pd_diag = ifelse.(iszero.(data.pmax - data.pmin), data.pmax, copysign.(1.0, (data.pmin .+ data.pmax) ./ 2))
    qd_diag = ifelse.(iszero.(data.qmax - data.qmin), data.qmax, 1.0)
    pd_mask = .!iszero.(pd_diag)
    qd_mask = .!iszero.(qd_diag)
    nvar_pd = _DF.nrow(pd)
    nvar_qd = _DF.nrow(qd)
    ids = (findall(pd_mask), findall(qd_mask))

    # Maximum and minimum load active power constraint
    pd_max = hcat(I(nvar_pd), zeros(nvar_pd, nvar_qd))
    pd_min = hcat(-I(nvar_pd), zeros(nvar_pd, nvar_qd))
    # Maximum and minimum load reactive power constraint
    qd_max = hcat(zeros(nvar_qd, nvar_pd), I(nvar_qd)) 
    qd_min = hcat(zeros(nvar_qd, nvar_pd), -I(nvar_qd))
    # Maximum and minimum load reactive/active power ratio constraint
    mask = pd_mask .& qd_mask
    pd_diag = diagm(pd_diag)
    qd_diag = diagm(qd_diag)[:, qd_mask]
    qd_ratio_max = hcat((-pd_diag .* data[!, "qp_ratio_max"])[:, pd_mask], qd_diag)
    qd_ratio_max = qd_ratio_max[mask, :]
    qd_ratio_min = hcat((pd_diag .* data[!, "qp_ratio_min"])[:, pd_mask], -qd_diag)
    qd_ratio_min = qd_ratio_min[mask, :]
    # Max total active demand constraint
    pd_max_tot = hcat(ones(1, nvar_pd), zeros(1, nvar_qd))
    # Min total active demand constraint
    pd_min_tot = hcat(-ones(1, nvar_pd), zeros(1, nvar_qd))

    A = vcat(pd_max, pd_min, qd_max, qd_min, qd_ratio_max, qd_ratio_min, pd_max_tot, pd_min_tot)
    b = vcat(pd.pmax, -pd.pmin, qd.qmax, -qd.qmin,
        zeros(sum(mask)),
        zeros(sum(mask)),
        sum(get_pm_value(model, :gen, ["pmax"], Array{Any, 2})) - pd_res,
        -sum(get_pm_value(model, :gen, ["pmin"], Array{Any, 2})) + pd_res
    )

	return (A = A, b = b, ids = ids)::PolyType
end
