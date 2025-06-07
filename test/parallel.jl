
@testset "test remote worker initialization" begin
    
    setup = deepcopy(SETUP)
    n_cpu = ceil(Int, Sys.CPU_THREADS * 0.25)
    n_cpu = n_cpu < 2 ? 2 : min(n_cpu, 4)
    setup["PARALLEL"]["cpu_ratio"] = (n_cpu / Sys.CPU_THREADS) * 100
    nt = to_namedtuple(setup)
    init_workers!(nt)

    @test isequal(_DC.nprocs(), n_cpu + 1)
    @test all(_DC.remotecall_fetch.(pwd, _DC.workers()) .== pwd())
    @test all(_DC.remotecall_fetch.(() -> isdefined(Main, Symbol(nt.SOLVER.lp)), _DC.workers()))

    _DC.nprocs() > 1 && _DC.rmprocs(_DC.workers())
end