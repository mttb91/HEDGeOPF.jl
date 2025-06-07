
@testset "test OPF results database initialization" begin
    
    network = _PM.parse_file("data/grids/pglib_opf_case5_pjm.m")
    optimizer = _PM.optimizer_with_attributes(Ipopt.Optimizer,
        "print_level" => 0)
    pm = _PM.instantiate_model(network, _PM.ACPPowerModel, _PM.build_opf)
    JuMP.set_optimizer(pm.model, optimizer)
   
    @testset "initialization base case" begin
        db = instantiate_database(pm, false)

        @test isa(db, NamedTuple)
        @test all(haskey.(Ref(db), [:branch, :bus, :check, :gen, :info]))
        @test isa(db.check, Vector{Bool})
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
    pm = instantiate_model(instantiate_network(setup), "ACP", to_namedtuple(setup));
    d = get_pm_value(pm, :load, ["pd", "qd"], _DF.DataFrame)
    db = instantiate_database(pm, false)
    res = _PM.optimize_model!(pm)
    res["input"] = Dict{String, Any}()
    sol = res["solution"]

    @testset "recording of primal OPF inputs and results" begin

        extract_data!(db, res, pm)

        @test first(db.check)
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
end