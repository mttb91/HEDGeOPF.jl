
###############################################################################
# Common helper functions
###############################################################################

function _check_vector_lengths(data::Any, mask::Vector)
    len = size(data, 1)
    msg = "The length of the masking array ($(size(mask, 1))) does not match the number of values ($len)."
    ((len == 1) || (len == size(mask, 1))) || throw(DimensionMismatch(msg))    
end

function _check_keys(data::Dict{T, <:Any}, element::T, vars::Vector{String}) where T
    # Get all properties that do not still exist for the given element
    vars = filter(x -> !haskey(first(values(data[element])), x), vars)
    if !isempty(vars)
        msg = "Values have been assigned to the following newly added $element property(s): $(join(vars, ", "))"
        @warn msg
    end
end

function _sort_keys!(keys::Vector{String})
    # Sort string keys by their numeric value
    @. keys = keys[$sortperm(parse(Int64, keys))]
    return nothing
end

function _sort_keys!(keys::Vector{Int64})
    # Sort string keys by their numeric value
    sort!(keys)
    return nothing
end

"Convert a dictionary of nested dictionaries, each with the same keys, to an array of the given type."
function _nested_dict2arr(dict::Dict{T, Any}, keys::Vector{T}, vars::Vector{String}, type::Type{Matrix{Any}}) where T

    # Initialize array to be returned
    dest = type(undef, length(keys), length(vars))
    # Extract a vector of values for each specified nested dictionary key
    @. dest = getindex.([dict[key] for key = keys], $permutedims(vars))
    # Create deepcopy of specific column if required
    cols = eachindex(vars)[isa.(view(dest, begin, :), Vector)]
    @views dest[:, cols] = deepcopy(dest[:, cols])
    # Convert to proper datatype if concrete
    type = typejoin(typeof.(first(eachrow(dest)))...)
    if isconcretetype(type)
        dest = convert(Array{type, 2}, dest)
    end

    return dest
end

function _nested_dict2arr(dict::Dict{T, Any}, keys::Vector{T}, vars::Vector{String}, type::Type{_DF.DataFrame}) where T

    # Extract a vector of values for each specified nested dictionary key
    list = []
    for var = vars
        push!(list, [dict[key][var] for key = keys])
    end
    # Generate the dataframe 
    dest = type(list, vars; copycols=false)
    # Create deepcopy of specific column if required
    cols = vars[isa.(first.(eachcol(dest)), Vector)]
    dest[!, cols] = deepcopy(dest[!, cols])

    return dest
end

function _get_pm_value(dict::Dict{T, Any}, vars::Vector{String}, mask::Vector{Int}, type::DataType) where T

    # Sort the nested dictionary keys
    nids = collect(keys(dict))
    _sort_keys!(nids)
    # Keep only masked keys
    if !isempty(mask)
        keepat!(nids, mask)
    end
    # Extract value from nested dictionaries
    return _nested_dict2arr(dict, nids, vars, type)
end


###############################################################################
# Power system network dictionary helper functions
###############################################################################

"Get all property names of an existing `element` in a Powermodels network dictionary."
function get_pm_key(model::Dict{String, <:Any}, element::String)

    # Get the first dictionary of values for the given element
    data = first(values(model[element]))
    return collect(keys(data))
end

function get_pm_key(model::_PM.AbstractPowerModel, element::Symbol)

    # Get the first dictionary of values for the given element
    data = first(values(_PM.ref(model, element)))
    return collect(keys(data))
end

"""
    get_pm_value(model::Dict{String, <:Any}, element::String, vars::Vector{String}, type::DataType;
        mask::Vector{Int} = Vector{Int}()
    )
    
Get the values of existing property(s) of a given `element` from a Powermodels `network` or `ref` dictionary.

## Notes

-   The argument `type` defines the type of the returned object and accepts only two possible values,
namely `Array{Any, 2}` and `DataFrame`.
-   By default, values are returned as ordered by sorted `element` index
unless a `mask` of element indices is specified. 

"""
function get_pm_value(model::Dict{String, <:Any}, element::String, vars::Vector{String}, type::DataType;
    mask::Vector{Int} = Vector{Int}()
)
    data = model[element]
    return _get_pm_value(data, vars, mask, type)
end

function get_pm_value(model::_PM.AbstractPowerModel, element::Symbol, vars::Vector{String}, type::DataType;
    mask::Vector{Int} = Vector{Int}()
)
    data = _PM.ref(model, element)
    return _get_pm_value(data, vars, mask, type)
end


"""
    set_pm_value!(model::Dict{String, <:Any}, element::String, vars::Vector{String}, value::Any;
        mask::Vector{Int} = Vector{Int}()
    )

Set value(s) for existing or new property name(s) `vars` of a given `element` in a Powermodels `network` or `ref` dictionary.

## Notes

Three value assignment methods are supported: element-by-element, element-wise and masked assignment. For the case of
element-by-element or masked assignment, it is assumed that values are provided as ordered by sorted element index.
It is not possible to use different masking if multiple properties are provided.
 
"""
function set_pm_value!(model::Dict{String, <:Any}, element::String, vars::Vector{String}, value::Any;
    mask::Vector{Int} = Vector{Int}()
)
    # Convert value to proper format if necessary
    if length(vars) != 1
        value = isa(value, Vector{<:Real}) ? permutedims(value) : value
    end
    # Sort the nested dictionary keys numerically
    nids = collect(keys(model[element]))
    _sort_keys!(nids)
    # Keep only masked keys
    if !isempty(mask)
        _check_vector_lengths(value, mask)
        keepat!(nids, mask)
    end
    # Give a warning is values are assigned to non-existing properties
    _check_keys(model, element, vars)
    # Assign the values to the given variables
    setindex!.([model[element][key] for key = nids], value, permutedims(vars))

    return nothing
end

function set_pm_value!(model::_PM.AbstractPowerModel, element::Symbol, vars::Vector{String}, value::Any;
    mask::Vector{Int} = Vector{Int}()
)
    # Convert value to proper format if necessary
    if length(vars) != 1
        value = isa(value, Vector{<:Real}) ? permutedims(value) : value
    end
    # Sort the nested dictionary keys numerically
    nids = collect(keys(_PM.ref(model, element)))
    _sort_keys!(nids)
    # Keep only masked keys
    if !isempty(mask)
        _check_vector_lengths(value, mask)
        keepat!(nids, mask)
    end
    # Assign the values to the given variables
    setindex!.([_PM.ref(model, element, key) for key = nids], value, permutedims(vars))

    return nothing
end
