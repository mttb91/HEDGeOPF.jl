
# Helper: build a synthetic map DataFrame
function _make_map(; n::Int = 60, n_topo::Int = 3, seed::Int = 0)
    rng = _RND.MersenneTwister(seed)
    topology_id = repeat(1:n_topo, outer = ceil(Int, n / n_topo))[1:n]
    pd_tot      = sort(0.5 .+ _RND.rand(rng, Float64, n))
    objective   = _RND.rand(rng, Float64, n)
    worker      = ones(Int, n)
    case        = 1:n
    return _DF.DataFrame(;
        uid = 1:n, worker, case, topology_id, pd_tot, objective,
        is_strictly_feasible = trues(n), pair = fill(-1, n)
    )
end

# Helper: build a minimal setting NamedTuple for generate_cv_folds!
function _make_setting(;
    n_fold::Int    = 3,
    n_sample       = nothing,
    n_quantile::Int = 4,
    seed::Int      = 0,
    unfeasible_ratio = 0.0
)
    return (
        CASE    = (baseseed = seed,),
        DATASET = (
            num_folds        = n_fold,
            num_samples      = n_sample,
            num_quantiles    = n_quantile,
            unfeasible_ratio = unfeasible_ratio,
        ),
    )
end

@testset "test dataset splitting functions" begin

    @testset "topology_split" begin

        @testset "each topology assigned to exactly one fold" begin
            map     = _make_map(n = 60, n_topo = 6)
            rng     = _RND.MersenneTwister(0)
            folds   = topology_split(map, rng, 3, _DF.nrow(map))

            assigned = map[folds .!= -1, :]
            assigned[!, :fold] = folds[folds .!= -1]

            # every fold index that appears is in 1..n_fold
            @test all(in.(unique(assigned.fold), Ref(1:3)))
            # no topology_id appears in more than one fold
            check = Bool[]
            for t in unique(assigned.topology_id)
                push!(check, length(unique(assigned.fold[assigned.topology_id .== t])) == 1)
            end
            @test all(check)
        end

        @testset "n_sample budget is respected" begin
            n_sample = 30
            map    = _make_map(n = 60, n_topo = 6)
            rng    = _RND.MersenneTwister(0)
            folds  = topology_split(map, rng, 3, n_sample)

            @test count(!=(-1), folds) <= n_sample
        end

        @testset "reproducibility with fixed seed" begin
            map1 = _make_map(n = 60, n_topo = 6)
            map2 = _make_map(n = 60, n_topo = 6)
            rng1 = _RND.MersenneTwister(42)
            rng2 = _RND.MersenneTwister(42)

            f1 = topology_split(map1, rng1, 3, _DF.nrow(map1))
            f2 = topology_split(map2, rng2, 3, _DF.nrow(map2))
            @test f1 == f2
        end

        @testset "insufficient topologies handled without hard failure" begin
            # 2 topology IDs, 3 folds required
            map = _make_map(n = 20, n_topo = 2)
            rng = _RND.MersenneTwister(0)
            folds = topology_split(map, rng, 3, _DF.nrow(map))
            @test any(folds .!= -1)
        end
    end

    @testset "total_load_active_power_split" begin

        @testset "samples distributed across all folds" begin
            map   = _make_map(n = 80, n_topo = 1)
            rng   = _RND.MersenneTwister(0)
            folds = total_load_active_power_split(map, rng, 4, _DF.nrow(map), 4)

            @test sort(unique(folds[folds .!= -1])) == 1:4
        end

        @testset "stratification: each fold receives samples from multiple quantile classes" begin
            map    = _make_map(n = 160, n_topo = 1)
            rng    = _RND.MersenneTwister(0)
            n_fold = 4
            n_q    = 4
            folds  = total_load_active_power_split(map, rng, n_fold, _DF.nrow(map), n_q)

            r        = 1 / n_q
            quantiles = Statistics.quantile(map.pd_tot, (0 + r):r:(1 - r); sorted = true)
            classes   = searchsortedfirst.(Ref(quantiles), map.pd_tot)

            n_classes_per_fold = [length(unique(classes[folds .== f])) for f in 1:n_fold]
            @test all(>(1), n_classes_per_fold)
        end

        @testset "n_sample budget respected (total assigned ≤ n_sample)" begin
            n_sample = 40
            map   = _make_map(n = 160, n_topo = 1)
            rng   = _RND.MersenneTwister(0)
            folds = total_load_active_power_split(map, rng, 4, n_sample, 4)

            @test count(!=(-1), folds) <= n_sample
        end

        @testset "reproducibility with fixed seed" begin
            map1 = _make_map(n = 80, n_topo = 1)
            map2 = _make_map(n = 80, n_topo = 1)
            rng1 = _RND.MersenneTwister(7)
            rng2 = _RND.MersenneTwister(7)

            f1 = total_load_active_power_split(map1, rng1, 3, _DF.nrow(map1), 4)
            f2 = total_load_active_power_split(map2, rng2, 3, _DF.nrow(map2), 4)
            @test f1 == f2
        end

        @testset "issorted(pd_tot) precondition guard fires" begin
            map = _make_map(n = 40, n_topo = 1)
            # Keep uid order unchanged but break pd_tot monotonicity.
            tmp = map.pd_tot[1]
            map.pd_tot[1] = map.pd_tot[2]
            map.pd_tot[2] = tmp
            rng = _RND.MersenneTwister(0)
            @test_throws ErrorException total_load_active_power_split(map, rng, 3, _DF.nrow(map), 4)
        end

        @testset "empty-bin guard fires for constant pd_tot" begin
            map = _make_map(n = 40, n_topo = 1)
            map[!, :pd_tot] .= 1.0               # all identical → some quantile bins empty
            HEDGeOPF._DF.sort!(map, :uid)
            rng = _RND.MersenneTwister(0)
            @test_throws ErrorException total_load_active_power_split(map, rng, 3, _DF.nrow(map), 4)
        end

        @testset "dim > 0 guard fires when n_sample too small for n_fold and n_quantile" begin
            # n_sample=2 → bin_size = 2÷4 = 0 → dim = 0÷3 = 0
            map = _make_map(n = 40, n_topo = 1)
            rng = _RND.MersenneTwister(0)
            @test_throws ErrorException total_load_active_power_split(map, rng, 3, 2, 4)
        end
    end

    @testset "generate_cv_folds!" begin

        @testset "single topology: fold column inserted, all folds valid, no -1 rows" begin
            map     = _make_map(n = 60, n_topo = 1)
            setting = _make_setting(n_fold = 3, n_quantile = 4)
            generate_cv_folds!(map, setting)

            @test _DF.hasproperty(map, :fold)
            @test all(map.fold .> 0)
            @test all(in.(unique(map.fold), Ref(1:3)))
        end

        @testset "multiple topologies: each topology_id in exactly one fold" begin
            map     = _make_map(n = 60, n_topo = 6)
            setting = _make_setting(n_fold = 3)
            generate_cv_folds!(map, setting)

            @test _DF.hasproperty(map, :fold)
            @test all(map.fold .> 0)

            n_folds_per_topo = [length(unique(map.fold[map.topology_id .== t]))
                                for t in unique(map.topology_id)]
            @test all(==(1), n_folds_per_topo)
        end

        @testset "downsampling via n_sample (single topology)" begin
            n = 60
            n_sample = 30
            map     = _make_map(n = n, n_topo = 1)
            setting = _make_setting(n_fold = 3, n_quantile = 4, n_sample = n_sample)
            generate_cv_folds!(map, setting)

            @test _DF.nrow(map) <= n_sample
            @test _DF.nrow(map) < n
        end

        @testset "downsampling via n_sample (multiple topologies)" begin
            n = 60
            n_sample = 30
            map     = _make_map(n = n, n_topo = 6)
            setting = _make_setting(n_fold = 3, n_sample = n_sample)
            generate_cv_folds!(map, setting)

            @test _DF.nrow(map) <= n_sample
            @test _DF.nrow(map) < n
        end

        @testset "n_sample > nrow assertion fires" begin
            map     = _make_map(n = 20)
            setting = _make_setting(n_sample = 9999)
            @test_throws ArgumentError generate_cv_folds!(map, setting)
        end
    end

    @testset "topology_split with unfeasible ratio" begin

        quiet(f) = Logging.with_logger(f, Logging.NullLogger())

        # Map with `nf` feasible instances per topology, a varying number of which has an unfeasible twin
        function _make_pair_map(; n_topo::Int = 40, nf::Int = 20, seed::Int = 0)
            map = _make_map(n = n_topo * nf, n_topo = n_topo, seed = seed)
            n_pair = _RND.rand(_RND.MersenneTwister(seed), 2:(nf ÷ 2 + 2), n_topo)
            twins = map[(map.uid .- 1) .÷ n_topo .< n_pair[map.topology_id], :]
            unfeasible = copy(twins)
            unfeasible.uid .= n_topo * nf .+ (1:_DF.nrow(twins))
            unfeasible.is_strictly_feasible .= false
            unfeasible.pair .= unfeasible.uid
            map.pair[twins.uid] .= unfeasible.uid
            return vcat(map, unfeasible)
        end

        @testset "ratio below natural: totals, ratio cap, consistent pairs" begin
            map = _make_pair_map()
            quiet() do
                generate_cv_folds!(map, _make_setting(n_fold = 4, n_sample = 600, unfeasible_ratio = 20.0))
            end
            @test count(map.is_strictly_feasible) == 600
            @test count(.!map.is_strictly_feasible) <= 120
            for g in _DF.groupby(map, :fold)
                @test count(.!g.is_strictly_feasible) <= 0.2 * count(g.is_strictly_feasible)
            end
            for g in _DF.groupby(filter(:pair => >(0), map), :pair)
                @test _DF.nrow(g) == 2 && length(unique(g.fold)) == 1 && count(g.is_strictly_feasible) == 1
            end
            @test all(length(unique(map.fold[map.topology_id .== t])) == 1 for t in unique(map.topology_id))
        end

        @testset "zero ratio: feasible instances only, original strategy" begin
            map = _make_pair_map()
            quiet() do
                generate_cv_folds!(map, _make_setting(n_fold = 4, n_sample = 600, unfeasible_ratio = 0.0))
            end
            @test all(map.is_strictly_feasible) && all(map.pair .== -1) && _DF.nrow(map) <= 600
        end

        @testset "natural ratio: unfeasible instances of selected topologies are kept" begin
            map = _make_pair_map()
            quiet() do
                generate_cv_folds!(map, _make_setting(n_fold = 4, unfeasible_ratio = nothing))
            end
            @test count(map.is_strictly_feasible) == 800
            @test count(.!map.is_strictly_feasible) == count(.!_make_pair_map().is_strictly_feasible)
        end

        @testset "no budget keeps all topologies, also when they are not a multiple of the folds" begin
            map = _make_pair_map(n_topo = 43)
            quiet() do
                generate_cv_folds!(map, _make_setting(n_fold = 4, unfeasible_ratio = nothing))
            end
            @test length(unique(map.topology_id)) == 43
            @test count(map.is_strictly_feasible) == 43 * 20
            n_topo = _DF.combine(_DF.groupby(map, :fold), :topology_id => (x -> length(unique(x))) => :n).n
            @test maximum(n_topo) - minimum(n_topo) <= 1
        end

        @testset "ratio above natural reduces the feasible sample number" begin
            map = _make_pair_map()
            setting = _make_setting(n_fold = 4, n_sample = 700, unfeasible_ratio = 100.0)
            @test_logs (:warn, r"cannot be met") match_mode = :any generate_cv_folds!(map, setting)
            @test count(map.is_strictly_feasible) < 700
            @test count(.!map.is_strictly_feasible) <= count(map.is_strictly_feasible)
        end

        @testset "binding sample budget that cannot be met" begin
            map = _make_pair_map(n_topo = 41) # 40 topologies fit in 4 folds: 800 of 820 feasible instances
            setting = _make_setting(n_fold = 4, n_sample = 820, unfeasible_ratio = 10.0)
            err = quiet() do
                @test_throws ArgumentError generate_cv_folds!(map, setting)
            end
            @test occursin("800", err.value.msg)
        end

        @testset "reproducibility with fixed seed" begin
            folds(seed) = quiet() do
                map = _make_pair_map()
                generate_cv_folds!(map, _make_setting(n_fold = 4, n_sample = 600, unfeasible_ratio = 20.0, seed = seed))
                return map.fold
            end
            @test folds(3) == folds(3)
            @test folds(3) != folds(4)
        end

        @testset "single topology: ratio above zero is not supported" begin
            map = _make_map(n = 60, n_topo = 1)
            @test_throws ArgumentError generate_cv_folds!(map, _make_setting(unfeasible_ratio = 10.0))
        end
    end
end