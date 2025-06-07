
@testset "test generic polytope functions" begin
    
    setup = to_namedtuple(deepcopy(SETUP))
    n = 3
    p = (
        A = vcat(diagm(ones(n)), diagm(-ones(n))),
        b = vcat(ones(n), ones(n)),
        ids = (collect(1:n), collect(1:0))
    )::PolyType
    p1 = (A = p.A, b = copy(p.b), ids = p.ids)::PolyType
    p1.b[n+1] = -2

    @testset "empty polytope" begin

        @test !_isempty(p, setup)
        @test _isempty(p1, setup)
    end

    @testset "chebyshev centre" begin

        m, x, _ = chebyshev_model(p.A, p.b, setup)
        m1, x1, _ = chebyshev_model(p1.A, p1.b, setup)
        c = chebyshev_centre(m, x)
        c1 = chebyshev_centre(m1, x1)

        @test all(iszero.(c))
        @test all(isnothing.(c1))
    end

    @testset "variable bound" begin
        
        id = 1
        @test isequal(load_power_bound(p, id, setup; upper=true), 1.0)
        @test isequal(load_power_bound(p, id, setup; upper=false), -1.0)
        @test isequal(load_power_bound(p, setup; upper=true), 1.0 * n)
        @test isequal(load_power_bound(p, setup; upper=false), -1.0 * n)
    end
end

@testset "test custom polytope function" begin

    setup = deepcopy(SETUP)
    net = instantiate_network(setup)
    nt = to_namedtuple(setup)
    pm = instantiate_model(net, "ACP", nt)
    gen = get_pm_value(pm, :gen, ["pmin", "pmax"], _DF.DataFrame)
    load = get_pm_value(pm, :load, ["pd", "pmin", "pmax",  "qd", "qmin", "qmax", "qp_ratio_min", "qp_ratio_max"], _DF.DataFrame)
    nvar = size(load, 1) * 2

    @testset "polytope generation" begin

        p = instantiate_polytope(pm)
        id = 2

        @test isa(p, PolyType)
        @test all(!iszero, Matrix(load[!, ["pd", "qd"]]))
        @test isequal(size(p.A, 2), nvar)
        @test isequal(size(p.A, 1), nvar * 3 + 2)
        @test isequal(size(p.A, 2), sum(length.(p.ids)))
        @test !_isempty(p, nt)
        @test isequal(load_power_bound(p, id, nt; var="pd", upper=true), load[id, "pmax"])
        @test isequal(load_power_bound(p, id, nt; var="pd", upper=false), load[id, "pmin"])
        @test isequal(load_power_bound(p, id, nt; var="qd", upper=true), load[id, "qmax"])
        @test isequal(load_power_bound(p, id, nt; var="qd", upper=false), load[id, "qmin"])
        @test isequal(load_power_bound(p, nt; var="pd", upper=true), min(sum(load.pmax), sum(gen.pmax)))
        @test isequal(load_power_bound(p, nt; var="pd", upper=false), max(sum(load.pmin), sum(gen.pmin)))
    end

    @testset "polytope generation with fixed loads" begin

        net1 = deepcopy(net)
        ids = (1, 2)
        set_pm_value!(net1, "load", ["pd", "pmin", "pmax"], 0.0; mask=[first(ids)])
        set_pm_value!(net1, "load", ["qd", "qmin", "qmax"], 0.0; mask=[last(ids)])
        p = instantiate_polytope(instantiate_model(net1, "ACP", nt))

        @test isa(p, PolyType)
        @test any(iszero, get_pm_value(net1, "load", ["pd", "qd"], Array{Any, 2}))
        @test isequal(size(p.A, 2), sum(length.(p.ids)))
        @test isequal(size(p.A, 1), (sum(length.(p.ids)) + length(intersect(p.ids...))) * 2 + 2)
        @test !in(first(ids), first(p.ids))
        @test !in(last(ids), last(p.ids))
        @test !_isempty(p, nt)
    end
end
