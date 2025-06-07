
@testset "test R functions" begin
    
    @testset "use RCall" begin
        x = 10
        y = convert(typeof(x), 
            RCall.R"""
            y <- $x + 5
            """
        )
        @test isequal(y, 15)
    end

    @testset "import Volesti R package" begin
        output = false
        try
            RCall.R"""
            library(volesti)
            """
            output = true
        catch
        end
        @test output
    end

    @testset "use Volesti R package" begin
        setup = to_namedtuple(deepcopy(SETUP))
        n, s = 3, 10
        p = (
            A = vcat(diagm(ones(n)), diagm(-ones(n))),
            b = vcat(ones(n), ones(n)),
            ids = (collect(1:n), collect(1:0))
        )::PolyType
        x0 = zeros(n)
        samples = sample_polytope(p.A, p.b, x0, s)

        @test isequal(size(samples), (n, s))
        @test all(all(p.A * sample .<= p.b) for sample in eachcol(samples))
    end
end