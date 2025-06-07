
###############################################################################
# Helper function
###############################################################################

"Create local path to destination folder and set it as current directory."
function _mkpath(settings::NamedTuple)

    path = joinpath(settings.PATH.output, first(split(settings.CASE.grid, ".")), settings.CASE.name)
    if !ispath(path)
        mkpath(path)
    else
        if !settings.CASE.append
            msg = "The destination folder must be unique. Rename the existing one to avoid overwriting."
            throw(DomainError(path, msg))
        end
    end
    cd(joinpath(pwd(), path))
    return nothing
end

"Create local path to folder if it does not exist already."
function _mkpath(path::String)

    if !ispath(path)
        mkpath(path)
    end
    return nothing
end

"Create folder if it does not exist already."
function _mkdir(path::String)

    if !isdir(path)
        mkdir(path)
    end
    return nothing
end

"Export a dictionary of DataFrames to Excel naming the sheets by the sorted dictionary keys."
function _to_xlsx(data::Dict{<:Any, _DF.DataFrame}, filepath::String)

    # Save static input information
    XLSX.openxlsx(filepath, mode="w") do xf
        for (i, key) in enumerate(sort!(collect(keys(data))))

            i == 1 ? XLSX.rename!(xf[i], String(key)) : XLSX.addsheet!(xf, String(key))
            XLSX.writetable!(xf[i], data[key])
        end
    end
end

###############################################################################
# I/O RNG
###############################################################################

"Freeze local RNG state and save it to file"
function save_rng(rng::_RND.AbstractRNG; filename::String = "rng_state.bin")

    open(joinpath(pwd(), filename), "w") do io
        serialize(io, rng)
    end
end

"Read RNG state if available else generate it"
function read_rng(settings::NamedTuple; filename::String = "rng_state.bin")

    if settings.CASE.append && isfile(joinpath(pwd(), filename))
        rng = open("rng_state.bin", "r") do io
            deserialize(io)
        end
    else
        rng = _RND.MersenneTwister(settings.CASE.baseseed)
    end
    return rng
end

###############################################################################
# I/O Polytope
###############################################################################

"Import polytope from .csv or .xlsx file"
function import_polytope(; path::Union{String, Nothing} = nothing, filename::String = "polytope")

    path = isnothing(path) ? dirname(pwd()) : path
    if isfile(joinpath(path, filename * ".xlsx"))
        polytope = _DF.DataFrame(XLSX.readtable(joinpath(path, filename * ".xlsx"), "polytope", infer_eltypes=true))
    elseif isfile(joinpath(path, filename * ".csv"))
        polytope = _DF.DataFrame(CSV.File(joinpath(path, filename * ".csv")))
    end
    b = polytope[!, end]
    A = _DF.select!(polytope, 1:(_DF.ncol(polytope) - 1))
    # Extract load indices
    cols = _DF.names(A)
    cols = [cols[contains.(cols, h)] for h in ["p", "q"]]
    ids = [parse.(Int, filter.(isdigit, h)) for h in cols]

    return (A = Matrix(A), b = b, ids = tuple(ids...))::PolyType
end

"Export polytope to .csv file"
function export_polytope(polytope::PolyType; filename::String = "polytope.csv")

    path = pwd()

    A, b, ids = polytope
    cols = reduce(vcat, [x .* string.(y) for (x, y) in zip(["p", "q"], ids)])
    push!(cols, "b")
    data = _DF.DataFrame(reduce(hcat, [A, b]), cols; copycols=false)
    ∉(filename, readdir(path)) && CSV.write(joinpath(path, filename), data)
    return nothing
end

###############################################################################
# I/O Power system
###############################################################################

"Export graph model based on PowerModels network in XLSX format"
function export_graph(model::_PM.AbstractPowerModel, config::String; filename::String = "graph.xlsx")

    path = joinpath(pwd(), config)
    
    data = Dict{Symbol, _DF.DataFrame}()
    for element in Symbol.(["bus", "branch", "load", "shunt", "gen"])
        # Extract all features of each component table
        if !isempty(_PM.ref(model, element))
            vars = get_pm_key(model, element)
            filter!(x -> !in(x, ["index", "source_id"]), vars)
            data[element] = get_pm_value(model, element, sort!(vars), _DF.DataFrame)
        end
    end

    element = :gen
    vars = [:cost, :ncost]
    # Unfold cost vector into components
    _DF.select!(data[element], _DF.Not(vars),
        vars => _DF.ByRow((x, y) -> vcat(repeat([0.0], 3-y), x)) => [:c2, :c1, :c0]
    )
    _DF.select!(data[element], sort(names(data[element])))
    # Identify synchronous condensers (if any)
    ids = findall(x -> x.pmax - x.pmin .== 0, eachrow(data[element]))
    if !isempty(ids)
        data[:sync] = data[element][ids, [:gen_bus]]
    end

    # Save graph data to xlsx
    ∉(filename, readdir(path)) && _to_xlsx(data, joinpath(path, filename))

    return nothing
end

###############################################################################
# I/O Reproducibility
###############################################################################

"Generate unique identifier for AC-OPF instances"
function generate_uid()

    ls_main = _DF.DataFrame[]
    configs = filter(entry -> isdir(entry), readdir())
    for config in configs
        
        path = joinpath(pwd(), config)
        ls_sub = _DF.DataFrame[]
        for file in filter(x -> contains(x, "info"), readdir(path))
            data = _DF.DataFrame(CSV.File(joinpath(path, file)))
            _DF.insertcols!(data, 1, :worker => parse(Int, file[6:end-4]))
            _DF.insertcols!(data, 2, :case => axes(data, 1))
            push!(ls_sub, data)    
        end
        df = reduce(vcat, ls_sub)
        _DF.insertcols!(df, 1, :config => parse(Int, config[2:end]))
        push!(ls_main, df)   
    end
    df = reduce(vcat, ls_main)
    # Sort AC-OPF instances based on total load active power and objective value
    _DF.sort!(df, [:pd_tot, :objective, :config])
    _DF.insertcols!(df, 1, :uid => axes(df, 1))
    _DF.sort!(df, [:config, :worker, :case])

    CSV.write("map.csv", df)
    return nothing
end