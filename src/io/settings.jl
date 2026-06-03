
"Convert a dictionary of nested dictionaries to a NamedTuple."
function to_namedtuple(data::Dict)
    return (; zip(Symbol.(keys(data)), to_namedtuple.(values(data)))...)
end
to_namedtuple(data::Any) = data

"Read a YAML configuration file into a dictionary of nested dictionaries and pre-preprocess it"
function read_settings(filename::String)

    # Read YAML file to nested dictionary
    data =  YAML.load_file(joinpath(pwd(), filename); dicttype=Dict{String, Any})
    # Check on settings values
    check_settings(data)

    return data
end

"Convert settings `dict` to `NamedTuple` and copy file to destination folder"
function update_settings(settings::Dict{String, <:Any}; filename::String = "settings.yaml")
    # Convert dictionary to namedtuple
    settings = to_namedtuple(settings)
    save_settings(settings; filename=filename)
    return settings
end

"Create destination folder and copy settings .yaml file to it"
function save_settings(settings::NamedTuple; filename::String = "settings.yaml")

    if isfile(abspath(filename))
        base_path = abspath(filename)
    else
        msg = "The settings file does not exist in $(pwd()). Please provide a valid path"
        throw(DomainError(filename, msg))
    end
    _mkpath(settings)
    cp(base_path, joinpath(pwd(), "settings.yaml"); force=true)
    return nothing
end

function check_settings(settings::Dict{String, <:Any})

    value = settings["CASE"]
    for key in ["num_samples", "num_batches", "num_items", "baseseed"]
        @assert isa(value[key], Int) "Option `$key` must be of Int type."
    end

    value, key = settings["PARALLEL"], "cpu_ratio"
    @assert value[key] > 0 "Option `$key` cannot be zero because serial dataset generation does not work. "

    value = settings["SAMPLING"]
    for key in ["delta_pd", "delta_qd", "delta_pf", "max_pf", "min_pf"]
        @assert isa(value[key], Float64) "Option `$key` must be of Float64 type."
    end
    
    key = "delta_pd"
    @assert value[key] > 0 "Invalid value for option `$key`: $(value[key]). It must be positive."
    @assert value[key] <= 100 (
        """
        Invalid value for option `$key`: $(value[key])%. This value represents the percentage range of variation 
        in active power around the nominal loading scenario. It must not exceed 100%, as the active power of a 
        single load cannot be unrestricted in sign. To create a negative (positive) load at a bus that has already 
        a positive (negative) one, add a new load to the power system dictionary or .m file with negative (positive) 
        nominal active power. 
        """
    )
    for key in ["delta_pf", "min_pf", "max_pf"]
        @assert (value[key] > 0) & (value[key] <= 1) "Invalid value for option `$key`: $(value[key]). It must be in range (0, 1]."
    end
    @assert value["max_pf"] > value["min_pf"] "The value of option `max_pf` must greater that `min_pf`."
    return nothing
end