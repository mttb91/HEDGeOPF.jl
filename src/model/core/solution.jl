
Base.@kwdef mutable struct InfoEntry
    termination_status::Vector{Int} = Int[]
    pd_tot::Vector{Float64} = Float64[]
    objective::Vector{Float64} = Float64[]
    solve_time::Vector{Float64} = Float64[]
end

Base.@kwdef mutable struct NodeEntry{T}
    var::Symbol
    ids::Vector{Int}
    data::Vector{Vector{T}} = Vector{T}[]
    mask_lb::Union{Vector{Bool}, Nothing} = nothing
    dual_lb::Union{Vector{Vector{T}}, Nothing} = nothing
    mask_ub::Union{Vector{Bool}, Nothing} = nothing
    dual_ub::Union{Vector{Vector{T}}, Nothing} = nothing
end

Base.@kwdef mutable struct EdgeEntry{T}
    var::Symbol
    ids::Vector{Tuple{Int, Int, Int}}
    data::Vector{Vector{T}} = Vector{T}[]
    mask_lb::Union{Vector{Bool}, Nothing} = nothing
    dual_lb::Union{Vector{Vector{T}}, Nothing} = nothing
    mask_ub::Union{Vector{Bool}, Nothing} = nothing
    dual_ub::Union{Vector{Vector{T}}, Nothing} = nothing
end

function _init_duals!(fields::Dict{Symbol, <:Any}, value::JuMP.Containers.DenseAxisArray)

    for (label, sense) in zip(["lb", "ub"], ["lower", "upper"])
        # Check if bound exist for given variable
        f = getfield(JuMP, Symbol("has_$(sense)_bound"))
        mask = (f.(value[fields[:ids]])).data
        if any(mask)
            merge!(fields, Dict(
                Symbol("mask_$label") => mask, 
                Symbol("dual_$label") => Vector{Float32}[])
            )
        end
    end
end

"Instantiate struct to record optimisation results of a given node variable `var` in a PowerModels `model`"
function _init_fields(value::JuMP.Containers.DenseAxisArray, var::Symbol, record_duals::Bool)
    # Initialize fields for primal variable results
    fields = Dict(:var => var, :ids => sort(first(value.axes), by=x -> x[1]))
    if record_duals
        _init_duals!(fields, value)
    end
    return NamedTuple(fields)
end

function _init_fields(value::JuMP.Containers.DenseAxisArray, var::Symbol, ids::Vector{Tuple{Int, Int, Int}}, record_duals::Bool)
    # Initialize fields for primal variable results
    fields = Dict(:var => var, :ids => sort(ids, by=x -> x[1]))
    if record_duals
        _init_duals!(fields, value)
    end
    return NamedTuple(fields)
end

function _init_fields(value::Dict{Int, <:Any}, var::Symbol, record_duals::Bool)
    return (var = var, ids = sort(collect(keys(value)), by=x -> x[1])) 
end

function _init_fields(value::Dict{Int, <:Any}, var::Symbol, ids::Vector{Tuple{Int, Int, Int}}, record_duals::Bool)
    return (var = var, ids = sort(ids, by=x -> x[1])) 
end

function init_node_fields(model::_PM.AbstractPowerModel, var::Symbol, record_duals::Bool)
    return _init_fields(_PM.var(model, var), var, record_duals)
end

function init_edge_fields(model::_PM.AbstractPowerModel, var::String, name::String, record_duals::Bool)

    suffix_var = contains(var, "_dc") ? "_dc" : ""
    # HACK: all branch/dcline variable names are assumed to end with "f" or "t"
    if name[end] == 'f'
        suffix_ref = "from"
    elseif name[end] == 't'
        suffix_ref = "to"
    end
    var = Symbol(var)
    value = _PM.var(model, var)
    if isa(first(value), JuMP.VariableRef)
        ids = _PM.ref(model, Symbol(:arcs_, suffix_ref, suffix_var))
    else
        ids = collect(keys(value))
    end
    return _init_fields(value, var, ids, record_duals)
end

"Empty the data collection database"
function empty_database!(db::NamedTuple)

    T = first(typeof(first(values(db.bus))).parameters)

    for (key, data) in pairs(db) 
        if key == :check
            continue
        elseif key == :info
            for field in fieldnames(typeof(data))
                empty!(getfield(data, field))
            end
        else
            for entry in values(data)
                stype = typeof(entry)
                for field in filter(x -> Vector{Vector{T}} <: fieldtype(stype, x), fieldnames(stype))
                    collection = getfield(entry, field)
                    if !isnothing(collection)
                        empty!(collection)
                    end
                end
            end
        end
    end
    return nothing
end

"Instantiate an empty database to collect all non-constant variables OPF `model` variables"
function instantiate_database(model::_PM.AbstractPowerModel, record_duals::Bool)

    # Add names for control variables of the optimisation problem
    results = _PM.optimize_model!(model)["solution"]
    filter!(x -> isa(x.second, Dict) & !isempty(x.second), results)
    # Map PowerModels variable labels to JuMP variable names
    labels = [get_pm_key(results, k) for k in keys(results)]
    vars = [k for (k, v) in _PM.var(model) if !isempty(v)]
    
    db = Dict{Symbol, Any}()
    # Initialize check to verify if input variables are recorded
    db[:check] = [false]
    # Initialize data collector with general information
    db[:info] = InfoEntry()

    for (element, names) in zip(Symbol.(keys(results)), labels)
        list = []
        for name in names
            # Initialize an empty struct for the variable to be recorded
            if in(element, [:branch, :dcline])
                var = name[1:end-1] * (element == :branch ? "" : "_dc")
                in(var, String.(vars)) || continue
                value = EdgeEntry{Float32}(; init_edge_fields(model, var, name, record_duals)...)
            else
                var = Symbol(name)
                in(var, vars) || continue
                value = NodeEntry{Float32}(; init_node_fields(model, var, record_duals)...)
            end
            push!(list, name => value)
        end
        # Merge dictionaries, appending nested dicts in case of common keys
        if !haskey(db, element)
            merge!(db, Dict(element => Dict(list...)))
        else
            mergewith!(merge, db, Dict(element => Dict(list...)))
        end
    end
    return NamedTuple(db)
end

"""
    extract_data!(db::NamedTuple, results::Dict{String, <:Any}, model::_PM.AbstractPowerModel)

Extract data from output dictionary `results`, saving it to the OPF instance database `db`.

## Notes

The `results` dictionary is organized with two main keys:
- `input` : contains all input data of the OPF simulation that varies across OPF instances
- `output`: contains the `PowerModels` OPF results dictionary
"""
function extract_data!(db::NamedTuple, results::Dict{String, <:Any}, model::_PM.AbstractPowerModel)

    if !first(db.check)
        # Initialize database for non-constant input data
        for (element, value) in results["input"]
            element = Symbol(element)
            if in(element, [:branch, :dcline])
                list = [var => EdgeEntry{Float32}(var = Symbol(var), ids = value[var].ids) for var in keys(value)]
            else
                list = [var => NodeEntry{Float32}(var = Symbol(var), ids = value[var].ids) for var in keys(value)]
            end
            if haskey(db, element)
                mergewith!(merge, getfield(db, element), Dict(list...))
            else
                @error "It is not possible to add field $element to the database being a NamedTuple."
            end
        end
        setindex!(db.check, true, 1)
    end
    
    solved = JuMP.is_solved_and_feasible(model.model)
    if solved
        # Extract global OPF solution information
        data = Dict(
            :termination_status => Int(solved),
            :pd_tot => sum(vec(get_pm_value(results["solution"], "load", ["pd"], Array{Any, 2}))),
            [Symbol(h) => results[h] for h in ["objective", "solve_time"]]...
        )
        for field in fieldnames(typeof(db.info))
            push!(getfield(db.info, field), data[field])
        end
        # Extract not-constant input OPF data that exist for the given model
        for key in (h for h in keys(db) if haskey(results["input"], string(h)))
            for (var, value) in db[key]
                if haskey(results["input"][string(key)], var)
                    push!(getfield(value, :data), Float32.(results["input"][string(key)][var].data))
                end
            end
        end
        # Extract OPF primal and (if requested) dual solution results
        for key in filter(k -> in(string(k), keys(results["solution"])), keys(db))
            for (var, value) in filter(p -> !haskey(get(results["input"], string(key), Dict()), p.first), db[key])
                push!(getfield(value, :data),
                    Float32.(vec(get_pm_value(results["solution"], string(key), [var], Array{Any, 2})))
                )
                for sense in ["lower", "upper"]
                    suffix = "_$(first(sense))b"
                    mask = getfield(value, Symbol("mask$suffix"))
                    if !isnothing(mask)
                        push!(getfield(value, Symbol("dual$suffix")),
                            Float32.(get_variable_dual(model, value.var, sense; ids=value.ids[mask]))
                        )
                    end
                end
            end
        end
    end
    return nothing
end

"Update input load sample if load slack variables are active at the optimum"
function relax_load!(sample::Dict{String, <:Any}, results::Dict{String, <:Any})

    ref = "load"
    for var in ["pd", "qd"]
        data = get_pm_value(results["solution"], ref, var * "_slack_" .* ["up", "down"], Array{Any, 2})
        sample[ref][var].data .+= vec(diff(data, dims=2))
    end
    return nothing
end

"""
    get_variable_dual(model::_PM.AbstractPowerModel, var::Symbol, sense::Symbol;
        ids::Vector{<:Any} = Vector()
    )

Retrieve the dual value(s) of a given variable `var` bound, specified by `sense`, at the
indices where the bound exists. Optionally, it is possible to select a subset of the
variable dual values through the `ids` keywork. By default, all dual values are returned
as ordered by sorted variable index.                      
"""
function get_variable_dual(model::_PM.AbstractPowerModel, var::Symbol, sense::String;
    ids::Vector{<:Any} = Vector()
)
    ref = _PM.var(model, var)
    # Set the ids equal to the variable indices
    ids = isempty(ids) ? sort!(first(ref.axes), by=x -> x[1]) : ids
    # Get masking of variable indices for which the the given bound exist
    f = getfield(JuMP, Symbol("has_$(sense)_bound"))
    mask = f.(ref[ids]).data
    # Get the dual value for the given variable bound
    f = getfield(JuMP, Symbol("$(uppercasefirst(sense))BoundRef"))
    # Extract the dual value only for variables whose bound exists
    value = zeros(Float64, length(ids))
    setindex!(value, JuMP.dual.(f.(ref[ids[mask]])).data, axes(ids, 1)[mask])

    return value
end
