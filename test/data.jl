
@testset "test functions to read/write nested PowerModels dictionaries" begin
    
    net = _PM.parse_file("data/grids/pglib_opf_case5_pjm.m")
    pm = _PM.instantiate_model(net, _PM.ACPPowerModel, _PM.build_opf)
    element, vars = "gen", ["pmin", "pmax"]

    @testset "read component properties from a PowerModels network or `ref` dictionary" begin
        data1 = get_pm_value(net, element, vars, Array{Any, 2})
        data2 = get_pm_value(net, element, vars, _DF.DataFrame)
        data3 = get_pm_value(pm, Symbol(element), vars, Array{Any, 2})
        data4 = get_pm_value(pm, Symbol(element), vars, _DF.DataFrame)
        
        @test all(isa.([data1, data3], Matrix))
        @test all(isa.([data2, data4], _DF.DataFrame))
        @test data1 == Matrix(data2) == data3 == Matrix(data4)
        @test data1 == _PM.component_table(net, element, vars)[:, 2:end]
    end

    @testset "modify component properties in a PowerModels network or `ref` dictionary" begin
        net1, net2, net3 = deepcopy(net), deepcopy(net), deepcopy(net);
        value = rand(length(net[element]), 2)
        ids = [1, 3, 4]
        set_pm_value!(net1, element, vars, first(value))
        set_pm_value!(net2, element, vars, value)
        set_pm_value!(net3, element, vars, value[ids, :]; mask=ids)

        @test all(isequal.(_PM.component_table(net1, element, vars)[:, 2:end], first(value)))
        @test isequal(_PM.component_table(net2, element, vars)[:, 2:end], value)
        @test isequal(_PM.component_table(net3, element, vars)[ids, 2:end], value[ids, :])
    end
end

@testset "test functions to modify in-place PowerModels optimisation model" begin
    
    setup = deepcopy(SETUP)
    setup["CASE"]["grid"] = "pglib_opf_case5_pjm.m"
    pm = instantiate_model(instantiate_network(setup), "ACP", to_namedtuple(setup))
    res = _PM.optimize_model!(pm)
    # Relevant input data
    baseMVA = _PM.ref(pm, 0, :baseMVA)
    load = get_pm_value(pm, :load, ["pd", "qd"], Array{Any, 2})
    cl = vec(first.(get_pm_value(pm, :gen, ["cost"], Array{Any, 2})))
    # Relevant output data
    pg = vec(get_pm_value(res["solution"], "gen", ["pg"], Array{Any, 2}))    

    @testset "load active/reactive power" begin
        m = deepcopy(pm)
        load1 = copy(load)
        id_pd, id_qd = 1, 3
        load1[id_pd, 1] += 100 / baseMVA
        load1[id_qd, 2] -= 50 / baseMVA
        sample = Dict(
            "info" => sum(load1[:, 1]) * [0.98, 1.02], 
            "load" => Dict(
                "pd" => InputSample(load1[[id_pd], 1], [id_pd]),
                "qd" => InputSample(load1[[id_qd], 2], [id_qd])
            )
        )
        update_model!(m, sample)
        res1 = _PM.optimize_model!(m)
        ids = sort(first(_PM.var(m, 0, :pd_fix).axes))

        @test string(res1["termination_status"]) == "LOCALLY_SOLVED"
        @test !isequal(load, load1)
        @test isequal(load1, reduce(hcat, [JuMP.value.(_PM.var(m, 0, h, ids).data) for h in [:pd_fix, :qd_fix]]))
        @test sum(getindex.(values(res1["solution"]["gen"]), "pg")) > sum(load1[:, 1])
        @test sum(getindex.(values(res1["solution"]["gen"]), "pg")) > sum(getindex.(values(res["solution"]["gen"]), "pg"))
        @test sum(getindex.(values(res1["solution"]["gen"]), "qg")) < sum(getindex.(values(res["solution"]["gen"]), "qg"))
    end

    @testset "generator linear cost" begin
        m = deepcopy(pm)
        cl1 = copy(cl)
        ids = [2, 4]
        cl1[ids] = cl[reverse(ids)]
        sample = Dict(
            "gen" => Dict("c1" => InputSample(data = cl1))
        )
        update_model!(m, sample)
        res1 = _PM.optimize_model!(m)
        pg1 = vec(get_pm_value(res1["solution"], "gen", ["pg"], Array{Any, 2}))

        @test string(res1["termination_status"]) == "LOCALLY_SOLVED"
        @test !isequal(cl, cl1)
        @test isapprox(pg[last(ids)], 0.0; atol = 1e-6)
        @test isapprox(pg1[first(ids)], 0.0; atol = 1e-6)
    end
end