
@testset "test OPF results database initialization" begin
    
    network = _PM.parse_file("data/grids/pglib_opf_case5_pjm.m")
    optimizer = _PM.optimizer_with_attributes(Ipopt.Optimizer,
        "print_level" => 0)
    pm = _PM.instantiate_model(network, _PM.ACPPowerModel, _PM.build_opf)
    JuMP.set_optimizer(pm.model, optimizer)
   
    @testset "initialization base case" begin
        db = instantiate_database(pm, false)

        @test isa(db, NamedTuple)
        @test all(haskey.(Ref(db), [:branch, :bus, :gen, :info]))
        @test isa(db.info, InfoEntry)
        @test isa(db.branch, Dict{String, EdgeEntry{Float32}})
        @test isa(db.gen, Dict{String, NodeEntry{Float32}})
        @test isequal(db.gen["pg"].ids, sort(parse.(Int, keys(network["gen"]))))
        @test isequal(db.branch["pf"].ids, Tuple.(eachrow(_PM.component_table(network, "branch", ["f_bus", "t_bus"]))))
        @test all(reduce(vcat, [isnothing.(getfield.(values(db.gen), h)) for h in [:mask_lb, :mask_ub]]))
        @test all(reduce(vcat, [isnothing.(getfield.(values(db.bus), h)) for h in [:mask_lb, :mask_ub]]))
        @test all(reduce(vcat, [isnothing.(getfield.(values(db.branch), h)) for h in [:mask_lb, :mask_ub]]))
    end

    @testset "initialization with dual variable recording" begin
        db = instantiate_database(pm, true)

        @test all(isnothing.(getfield.(Ref(db.bus["va"]), [:mask_lb, :mask_ub])))
        @test all(isnothing.(getfield.(Ref(db.gen["pg_cost"]), [:mask_lb, :mask_ub])))
        @test all(isa.(getfield.(Ref(db.gen["pg"]), [:mask_lb, :mask_ub]), Vector{Bool}))
        @test all(isa.(getfield.(Ref(db.gen["qg"]), [:mask_lb, :mask_ub]), Vector{Bool}))
        @test all(isa.(getfield.(Ref(db.bus["vm"]), [:mask_lb, :mask_ub]), Vector{Bool}))
        @test all(reduce(vcat, [isa.(getfield.(values(db.branch), h), Vector{Bool}) for h in [:mask_lb, :mask_ub]]))
    end

    @testset "initialization with altered topology" begin
        net = deepcopy(network)
        id = "2"
        net["gen"][id]["gen_status"] = 0
        delete!(net["branch"], id)
        m = _PM.instantiate_model(net, _PM.ACPPowerModel, _PM.build_opf)
        JuMP.set_optimizer(m.model, optimizer)
        db = instantiate_database(m, true)

        @test !isequal(db.gen["pg"].ids, sort(parse.(Int, keys(net["gen"]))))
        @test !in(parse(Int, id), db.gen["pg"].ids)
        @test !in((2, 1, 4), db.branch["pf"].ids)
        @test !in((2, 4, 1), db.branch["pt"].ids)
        @test isequal(sum(db.gen["pg"].mask_lb), length(net["gen"]) - 1)
        @test isequal(sum(db.branch["pf"].mask_lb), length(net["branch"]))
    end
end

@testset "test recording and clearing of OPF results in database" begin

    setup = deepcopy(SETUP)
    pm = deepcopy(DATA["5_pjm"]["pm"]);
    d = get_pm_value(pm, :load, ["pd", "qd"], _DF.DataFrame)
    db = instantiate_database(pm, false)
    res = _PM.optimize_model!(pm)
    res["input"] = Dict{String, Any}()
    res["id"] = 1
    sol = res["solution"]

    @testset "recording of primal OPF inputs and results" begin

        extract_data!(db, res, pm)

        @test !any(e.locked for entries in (db.bus, db.gen, db.load, db.branch) for e in values(entries))
        @test isequal(first(db.info.termination_status), 1)
        @test isapprox(first(db.info.pd_tot), sum(d.pd); atol = 1e-6)
        @test isapprox(first(db.info.objective), 17552.0; atol = 1e0)
        @test all(.!isempty.(getfield.(values(db.bus), :data)))
        @test all(.!isempty.(getfield.(values(db.gen), :data)))
        @test all(.!isempty.(getfield.(values(db.branch), :data)))
        @test all(isapprox.(first(db.load["pd"].data), d.pd; atol = 1e-06))
        @test all(isapprox.(first(db.gen["pg"].data), _PM.component_table(sol, "gen", "pg")[:, end]; atol = 1e-06))
        @test all(isapprox.(first(db.branch["pf"].data), _PM.component_table(sol, "branch", "pf")[:, end]; atol = 1e-06))
    end

    @testset "emptying the database" begin
        
        extract_data!(db, res, pm)

        @test all(.!isempty.(getfield.(values(db.bus), :data)))
        @test all(.!isempty.(getfield.(values(db.gen), :data)))
        @test all(.!isempty.(getfield.(values(db.branch), :data)))

        empty_database!(db)

        @test all(isempty.(getfield.(values(db.bus), :data)))
        @test all(isempty.(getfield.(values(db.gen), :data)))
        @test all(isempty.(getfield.(values(db.branch), :data)))
        @test all(.!isempty.(getfield.(values(db.bus), :ids)))
        @test all(.!isempty.(getfield.(values(db.gen), :ids)))
        @test all(.!isempty.(getfield.(values(db.branch), :ids)))
    end

    @testset "record of primal inputs/outputs for non-strictly feasible OPF instance" begin

        pm = DATA["5_pjm"]["pm"]
        ids = sort(collect(_PM.ids(pm, :load)))
        # Loads at their upper bound exceed the total generation capacity, hence require slack activation
        pmax = vec(get_pm_value(pm, :load, ["pmax"], Array{Any, 2}))
        slack_vars = ["pd_slack_up", "pd_slack_down", "qd_slack_up", "qd_slack_down"]

        load_values(results, var) = vec(get_pm_value(results["solution"], "load", [var], Array{Any, 2}))
        # Setpoints are sampled only for the first `n` loads, the others remain at their nominal value
        function solve_setpoints(pd, n = length(ids))
            model = deepcopy(pm)
            sample = Dict("load" => Dict(
                "pd" => InputData(pd[1:n], ids[1:n]), "qd" => InputData(d.qd[1:n], ids[1:n])
            ))
            update_load_power!(model, sample, total_load_active_power_bounds(model))
            res = _PM.optimize_model!(model)
            res["input"], res["id"] = sample, 1
            return model, res
        end
        model, res = solve_setpoints(pmax)

        @testset "preparation of unfeasible instance" begin

            results = deepcopy(res)
            process_unfeasible_instance!(results, model)

            @test results["is_strictly_feasible"] === false
            @test results["objective"] == res["objective"]
            @test results["solution"]["load"] == res["solution"]["load"]
            @test results["input"]["load"]["pd"].data == pmax
            @test results["input"]["load"]["qd"].data == d.qd

            results = deepcopy(res)
            delete!(results["input"]["load"], "qd")
            @test_throws ArgumentError process_unfeasible_instance!(results, model)
        end

        @testset "preparation of feasible instance" begin

            results = deepcopy(res)
            process_feasible_instance!(results, model)
            cost = sum(JuMP.value, values(_PM.var(model, 0, :slack_cost)))

            @test results["is_strictly_feasible"] === true
            @test isempty(results["input"]["load"])
            @test all(all(iszero(load[var]) for var in slack_vars) for load in values(results["solution"]["load"]))
            @test load_values(results, "pd") == load_values(res, "pd")
            @test cost > 0
            @test isapprox(results["objective"], res["objective"] - cost)
            @test isapprox(results["objective"], sum(JuMP.value, values(_PM.var(model, 0, :pg_cost))))
        end

        @testset "recording of unfeasible and feasible instances" begin

            db = instantiate_database(deepcopy(pm), false)
            extract_instance!(db, deepcopy(res), model)
            cost = sum(JuMP.value, values(_PM.var(model, 0, :slack_cost)))

            @test db.info.is_strictly_feasible == [false, true]
            @test db.info.termination_status == [1, 1]
            @test isapprox(db.info.pd_tot, [sum(pmax), sum(load_values(res, "pd"))]; atol = 1e-5)
            @test isapprox(db.info.objective, res["objective"] .- [0, cost])
            # Unfeasible instance: sampled setpoints and optimal slack; feasible instance: corrected setpoints, no slack
            @test all(isapprox.(db.load["pd"].data, [Float32.(pmax), Float32.(load_values(res, "pd"))]; atol = 1e-6))
            @test all(isapprox.(db.load["qd"].data, [Float32.(d.qd), Float32.(load_values(res, "qd"))]; atol = 1e-6))
            @test any(>(1e-4), db.load["pd_slack_up"].data[1])
            for var in slack_vars
                @test db.load[var].data[1] == Float32.(load_values(res, var))
                @test all(iszero, db.load[var].data[2])
            end
            # Everything else is recorded identically in both instances, with the same row order
            @test all(length(entry.data) == 2 for entry in values(db.load))
            for element in [:bus, :gen, :branch], entry in values(db[element])
                @test length(entry.data) == 2
                @test entry.data[1] == entry.data[2]
            end
        end

        @testset "recording of unfeasible instance with loads not sampled" begin

            # The last load is not sampled, its nominal setpoint must still be recorded in the unfeasible instance
            model, res = solve_setpoints(pmax, length(ids) - 1)
            db = instantiate_database(deepcopy(pm), false)
            extract_instance!(db, res, model)

            @test db.info.is_strictly_feasible == [false, true]
            @test db.load["pd"].data[1] == Float32.([pmax[1:end-1]; d.pd[end]])
            @test db.load["qd"].data[1] == Float32.(d.qd)
            @test isapprox(first(db.info.pd_tot), sum(pmax[1:end-1]) + d.pd[end]; atol = 1e-5)
        end

        @testset "recording of strictly feasible instance" begin

            model, res = solve_setpoints(d.pd)
            db = instantiate_database(deepcopy(pm), false)
            extract_instance!(db, res, model)

            @test db.info.is_strictly_feasible == [true]
            @test isapprox(db.info.pd_tot, [sum(load_values(res, "pd"))]; atol = 1e-5)
            @test isapprox(first(db.info.objective), 17552.0; atol = 1e0)
            for var in slack_vars
                @test all(iszero, first(db.load[var].data))
            end
        end

        @testset "no instance recorded if the model is not solved" begin

            model = deepcopy(pm)
            db = instantiate_database(deepcopy(pm), false)
            extract_instance!(db, Dict{String, Any}("input" => Dict{String, Any}(), "id" => 1), model)

            @test !JuMP.is_solved_and_feasible(model.model)
            @test isempty(db.info.pd_tot)
        end

        @testset "recording of solutions from input samples" begin

            for (pd, is_feasible) in [(pmax, [false, true]), (d.pd, [true])]
                model = deepcopy(pm)
                db = instantiate_database(model, false)
                topology = TopologyPerturbation(id = 1, pd_tot_bounds = total_load_active_power_bounds(model))
                sample = InputSample(
                    topology = topology,
                    pd_tot = topology.pd_tot_bounds,
                    data = Dict("load" => Dict("pd" => InputData(pd, ids), "qd" => InputData(d.qd, ids)))
                )
                out = solve_model!(model, db, sample)

                @test db.info.is_strictly_feasible == is_feasible
                @test isapprox(first(db.info.pd_tot), sum(pd); atol = 1e-3)
                # A single solve is counted, with the load total of the feasible instance
                @test out[[1, 3]] == [1.0, 1.0]
                @test isapprox(out[2], last(db.info.pd_tot); atol = 1e-3)
            end
        end
    end
end

@testset "test initialisation of database entries for input data" begin

    pm = deepcopy(DATA["5_pjm"]["pm"])
    ids_load = sort(collect(_PM.ids(pm, :load)))
    ids_gen = sort(collect(_PM.ids(pm, :gen)))
    pd0 = vec(get_pm_value(pm, :load, ["pd"], Array{Any, 2}))
    pos = [1, length(ids_load)]

    @testset "new entry of an input variable" begin
        entry = HEDGeOPF._new_input_entry(pm, :gen, "c1")

        @test isa(entry, NodeEntry{Float32})
        @test isequal(entry.var, :c1)
        @test isequal(entry.ids, ids_gen)
        @test !entry.locked
        @test isempty(entry.data)

        # Inactive components are not part of the entry ids
        net = deepcopy(DATA["5_pjm"]["net"])
        net["gen"]["2"]["gen_status"] = 0
        m = instantiate_model(net, "ACP", to_namedtuple(SETUP))
        @test isequal(HEDGeOPF._new_input_entry(m, :gen, "c1").ids, setdiff(ids_gen, 2))

        @test_throws ArgumentError HEDGeOPF._new_input_entry(pm, :branch, "rate_a")
        @test_throws ArgumentError HEDGeOPF._new_input_entry(pm, :dcline, "rate_a")
    end

    @testset "input at a subset of ids, with default values from the model" begin
        entry = NodeEntry{Float32}(var = :pd, ids = ids_load)
        init_input_entry!(entry, InputData(pd0[pos] .+ 0.1, ids_load[pos]), pm, :load)

        @test entry.locked
        @test isequal(entry.pos, pos)
        @test isa(entry.base, Vector{Float32})
        @test isequal(entry.base, Float32.(pd0))
        @test isempty(entry.data)
    end

    @testset "input at a subset of non-contiguous ids, with zero as default value" begin
        ids = [1, 3, 4, 5]
        entry = NodeEntry{Float32}(var = :c1, ids = ids)

        @test_logs (:warn, r"not a model property") init_input_entry!(entry, InputData([1.0, 2.0], [3, 5]), pm, :gen)
        @test entry.locked
        @test isequal(entry.pos, [2, 4])
        @test isequal(entry.base, zeros(Float32, length(ids)))
    end

    @testset "input at all ids needs no default value" begin
        entry = NodeEntry{Float32}(var = :c1, ids = ids_gen)

        @test_logs init_input_entry!(entry, InputData(ones(length(ids_gen)), nothing), pm, :gen)
        @test entry.locked
        @test isnothing(entry.base)
        @test isequal(entry.pos, eachindex(ids_gen))
    end

    @testset "input sharing its name with a solution variable recorded with duals" begin
        db = instantiate_database(pm, true)
        entry = db.gen["pg"]
        masks = (entry.mask_lb, entry.mask_ub)

        @test !entry.locked
        @test_logs (:warn, r"dual values are recorded") init_input_entry!(entry, InputData(ones(2), ids_gen[1:2]), pm, :gen)
        @test entry.locked
        @test isequal((entry.mask_lb, entry.mask_ub), masks)
        @test isequal(entry.ids, ids_gen)

        # No warning if duals are not recorded
        entry = instantiate_database(pm, false).gen["pg"]
        @test_logs init_input_entry!(entry, InputData(ones(2), ids_gen[1:2]), pm, :gen)
        @test entry.locked
    end

    @testset "invalid input" begin
        entry = NodeEntry{Float32}(var = :pd, ids = ids_load)

        @test_throws ArgumentError init_input_entry!(
            entry, InputData(ones(2), [first(ids_load), maximum(ids_load) + 1]), pm, :load
        )
        @test_throws DimensionMismatch init_input_entry!(entry, InputData(ones(1), ids_load[pos]), pm, :load)
        @test_throws DimensionMismatch init_input_entry!(entry, InputData(ones(2), nothing), pm, :load)
        @test !entry.locked
    end
end

@testset "test recording of input data in the database" begin

    pm = deepcopy(DATA["5_pjm"]["pm"])
    ids_load = sort(collect(_PM.ids(pm, :load)))
    ids_gen = sort(collect(_PM.ids(pm, :gen)))
    pd0 = vec(get_pm_value(pm, :load, ["pd"], Array{Any, 2}))
    res = _PM.optimize_model!(pm)
    res["id"] = 1
    sol_pd = vec(get_pm_value(res["solution"], "load", ["pd"], Array{Any, 2}))
    # Input defined on a subset of the loads, not coinciding with the solution
    pos = [1, length(ids_load)]
    vals = pd0[pos] .+ 0.123
    expected = Float32.(pd0)
    expected[pos] .= vals

    @testset "subset input is completed over all element ids and takes precedence over solution" begin
        db = instantiate_database(pm, false)
        entry = db.load["pd"]
        res["input"] = Dict("load" => Dict("pd" => InputData(vals, ids_load[pos])))
        extract_data!(db, res, pm)

        @test db.load["pd"] === entry
        @test entry.locked
        @test isequal(length(entry.data), 1)
        @test isapprox(first(entry.data), expected)
        @test !isapprox(first(entry.data)[pos], sol_pd[pos]; atol = 1e-3)
    end

    @testset "input at all ids is recorded as it is" begin
        db = instantiate_database(pm, false)
        vals_all = pd0 .+ 0.456
        res["input"] = Dict("load" => Dict("pd" => InputData(vals_all, nothing)))
        extract_data!(db, res, pm)

        @test isapprox(first(db.load["pd"].data), vals_all)
    end

    @testset "same variable from input in some instances and from solution in others" begin
        db = instantiate_database(pm, false)
        entry = db.load["pd"]
        input = Dict("load" => Dict("pd" => InputData(vals, ids_load[pos])))
        for in_dict in (Dict{String, Any}(), input, Dict{String, Any}())
            res["input"] = in_dict
            extract_data!(db, res, pm)
        end

        @test isequal(length(entry.data), 3)
        @test all(isequal.(length.(entry.data), length(ids_load)))
        @test isapprox(entry.data[1], sol_pd; atol = 1e-6)
        @test isapprox(entry.data[2], expected)
        @test isapprox(entry.data[3], sol_pd; atol = 1e-6)
    end

    @testset "input-only variable is added to the database" begin
        db = instantiate_database(pm, false)
        res["input"] = Dict("gen" => Dict("c1" => InputData([1.0, 2.0], ids_gen[[2, 4]])))
        Logging.with_logger(Logging.NullLogger()) do
            extract_data!(db, res, pm)
            extract_data!(db, res, pm)
        end

        @test haskey(db.gen, "c1")
        @test isequal(db.gen["c1"].ids, ids_gen)
        @test isequal(length(db.gen["c1"].data), 2)
        @test all(isequal.(db.gen["c1"].data, Ref(Float32[0, 1, 0, 2, 0])))
    end

    @testset "input of an element that is not in the database" begin
        db = instantiate_database(pm, false)
        res["input"] = Dict("shunt" => Dict("gs" => InputData(vals, nothing)))

        @test_throws ArgumentError extract_data!(db, res, pm)
    end
end
