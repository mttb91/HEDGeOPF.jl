
@testset "test YAML file reading and content" begin
    
    setup = read_settings("data/settings.yaml")

    @testset "Invalid configuration file" begin
        @test_throws AssertionError read_settings("data/settings - error(1).yaml")
        @test_throws KeyError read_settings("data/settings - error(2).yaml")
    end

    @testset "Valid configuration file" begin

        @test isa(setup, Dict{String, Any})
        @test all(haskey.(Ref(setup), ["CASE", "MODEL", "PARALLEL", "PATH", "SAMPLING", "SOLVER"]))
        @test all(haskey.(Ref(setup["CASE"]), ["append", "baseseed", "grid", "name", "num_batches", "num_items", "num_samples", "uid"]))
        @test all(haskey.(Ref(setup["MODEL"]), ["duals", "voll"]))
        @test all(haskey.(Ref(setup["PARALLEL"]), ["cpu_ratio"]))
        @test all(haskey.(Ref(setup["PATH"]), ["input", "output"]))
        @test all(haskey.(Ref(setup["SAMPLING"]), ["delta_pd", "delta_qd", "delta_pf", "max_pf", "min_pf"]))
        @test all(haskey.(Ref(setup["SOLVER"]), ["lp", "nlp", "lp_options", "nlp_options"]))
    end

    for (key, value) in zip(
        ["delta_pd", "delta_qd", "delta_pf", "max_pf", "min_pf"],
        [100.0, 100.0, 0.05, 0.99, 0.01]
    )
        setup["SAMPLING"][key] = value
    end
    setup["PATH"]["input"] = "data/grids"
    setup["MODEL"]["voll"] = 3000.0
    setup["SOLVER"]["lp"] = "HiGHS"
    setup["SOLVER"]["nlp"] = "Ipopt"
    setup["SOLVER"]["lp_options"] = Dict("solver" => "ipm")
    setup["SOLVER"]["nlp_options"] = Dict(
        "tol" => 1e-8,
        "max_cpu_time" => 1000.0,
        "print_level" => 0
    )

    return setup
end