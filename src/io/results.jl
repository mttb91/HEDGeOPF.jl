
"Export data to .csv file at the pre-defined path folder, emptying the database."
function export_data!(db::NamedTuple, config::String, worker::Int)

    fields = [:data, :dual_lb, :dual_ub]
    for (element, entry) in filter(x -> !in(x.first, [:check, :info]), pairs(db))

        element = string(element)
        _mkpath(joinpath(config, element))
        for (var, value) in entry
            # Check which variables are recorded (among primals and duals)
            mask = .!isnothing.(getfield.(Ref(value), fields))
            for (n, f) in zip(["", "_min", "_max"][mask], fields[mask])

                isempty(getfield(value, f)) && continue
                # Convert to dataframe
                cols = isequal(f, :data) ? value.ids : value.ids[getfield(value, Symbol(replace(string(f), "dual" => "mask")))]
                cols = isa(first(cols), Tuple) ? first.(cols) : cols
                data = _DF.DataFrame(reduce(hcat, getfield(value, f))', string.(cols); copycols=false)
                # Save data
                filename = joinpath(pwd(), config, element, "$(var)$(n)-$(worker-1).csv")
                CSV.write(filename, data; append = isfile(filename) ? true : false)
            end
        end
    end
    # Write general information to file
    data = _DF.DataFrame([f => getfield(db.info, f) for f in fieldnames(typeof(db.info))]; copycols=false)
    filename = joinpath(pwd(), config, "info-$(worker-1).csv")
    CSV.write(filename, data; append = isfile(filename) ? true : false)
    # Empty data collector
    empty_database!(db)

    return nothing
end

"Export to file a batch of OPF instances (i.e. input-solution pairs)"
function export_batch!(db::NamedTuple, worker::Int, config::String;
    threshold::Int = 10,
    final_batch::Bool = false)

    dim = length(db.info.termination_status)
    test_a = dim >= threshold && !final_batch
    test_b = final_batch && !iszero(dim)

    if test_a | test_b
        export_data!(db, config, worker)
    end
    return nothing
end
