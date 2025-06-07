
@testset "test OPF formulation with load slack variables" begin

    setup = deepcopy(SETUP)
    
    @testset "5-bus case" begin
        setup["CASE"]["grid"] = "pglib_opf_case5_pjm.m"
        res = _PM.optimize_model!(instantiate_model(instantiate_network(setup), "ACP", to_namedtuple(setup)));
        sol = res["solution"]

        @test string(res["termination_status"]) == "LOCALLY_SOLVED"
        @test isapprox(res["objective"], 17552.0; atol = 1e0)
        @test all(get.(values(sol["load"]), "pd_slack_up", NaN) .< 1e-6)
        @test all(get.(values(sol["load"]), "qd_slack_up", NaN) .< 1e-6)
        @test all(get.(values(sol["load"]), "pd_slack_down", NaN) .< 1e-6)
        @test all(get.(values(sol["load"]), "qd_slack_down", NaN) .< 1e-6)
    end

    @testset "5-bus case with quadratic branch apparent power variables" begin
        setup["MODEL"]["duals"] = true
        setup["CASE"]["grid"] = "pglib_opf_case5_pjm.m"
        res = _PM.optimize_model!(instantiate_model(instantiate_network(setup), "ACP", to_namedtuple(setup)));
        sol = res["solution"]
    
        @test string(res["termination_status"]) == "LOCALLY_SOLVED"
        @test isapprox(res["objective"], 17552.0; atol = 1e0)
        @test all(get.(values(sol["load"]), "pd_slack_up", NaN) .< 1e-6)
        @test all(get.(values(sol["load"]), "qd_slack_up", NaN) .< 1e-6)
        @test all(get.(values(sol["load"]), "pd_slack_down", NaN) .< 1e-6)
        @test all(get.(values(sol["load"]), "qd_slack_down", NaN) .< 1e-6)

        setup["MODEL"]["duals"] = false
    end

    @testset "5-bus case with active load slack variables" begin
        setup["CASE"]["grid"] = "pglib_opf_case5_pjm.m"
        net = instantiate_network(setup);
        net["gen"]["3"]["gen_status"] = 0
        net["load"]["1"]["pd"] = 400.0 / net["baseMVA"]
        net["load"]["2"]["pd"] = 500.0 / net["baseMVA"]
        net["load"]["3"]["pd"] = 0.0
        net["load"]["3"]["qd"] = 0.0
        res = _PM.optimize_model!(instantiate_model(net, "ACP", to_namedtuple(setup)));
        sol = res["solution"]

        @test string(res["termination_status"]) == "LOCALLY_SOLVED"
        @test isapprox(res["objective"], 478477.8; atol = 1e0)
        @test any(get.(values(sol["load"]), "pd_slack_up", NaN) .> 1e-6)
        @test all(get.(values(sol["load"]), "pd_slack_down", NaN) .< 1e-6)
        @test isapprox(sqrt(sum(sol["branch"]["1"][h]^2 for h in ["pf", "qf"])), net["branch"]["1"]["rate_a"], atol = 1e-6)
    end

    @testset "5-bus case with active load slack and quadratic branch apparent power variables" begin
        setup["MODEL"]["duals"] = true
        setup["CASE"]["grid"] = "pglib_opf_case5_pjm.m"
        net = instantiate_network(setup);
        net["gen"]["3"]["gen_status"] = 0
        net["load"]["1"]["pd"] = 400.0 / net["baseMVA"]
        net["load"]["2"]["pd"] = 500.0 / net["baseMVA"]
        net["load"]["3"]["pd"] = 0.0
        net["load"]["3"]["qd"] = 0.0
        res = _PM.optimize_model!(instantiate_model(net, "ACP", to_namedtuple(setup)));
        sol = res["solution"]

        @test string(res["termination_status"]) == "LOCALLY_SOLVED"
        @test isapprox(res["objective"], 478477.8; atol = 1e0)
        @test any(get.(values(sol["load"]), "pd_slack_up", NaN) .> 1e-6)
        @test all(get.(values(sol["load"]), "pd_slack_down", NaN) .< 1e-6)
        @test isapprox(sqrt(sum(sol["branch"]["1"][h]^2 for h in ["pf", "qf"])), net["branch"]["1"]["rate_a"], atol = 1e-6)
        @test isapprox(sqrt(sol["branch"]["1"]["sf"]), net["branch"]["1"]["rate_a"], atol = 1e-6)

        setup["MODEL"]["duals"] = false
    end    

    @testset "30-bus case" begin
        setup["CASE"]["grid"] = "pglib_opf_case30_ieee.m"
        res = _PM.optimize_model!(instantiate_model(instantiate_network(setup), "ACP", to_namedtuple(setup)));
        sol = res["solution"]

        @test string(res["termination_status"]) == "LOCALLY_SOLVED"
        @test isapprox(res["objective"], 8208.5; atol = 1e0)
        @test all(get.(values(sol["load"]), "pd_slack_up", NaN) .< 1e-6)
        @test all(get.(values(sol["load"]), "qd_slack_up", NaN) .< 1e-6)
        @test all(get.(values(sol["load"]), "pd_slack_down", NaN) .< 1e-6)
        @test all(get.(values(sol["load"]), "qd_slack_down", NaN) .< 1e-6)
    end
end