
@testset "test perturbation functions" begin

    setup = deepcopy(SETUP)

    @testset "test generation perturbation" begin

        k_min = 0
        k_max = 2
        setup["CASE"]["grid"] = "pglib_opf_case5_pjm.m"
        pm = deepcopy(DATA["5_pjm"]["pm"]);
        data = define_candidate_ref_buses(pm)
        ids_gen = data.index
        ids_bus_ref_valid = data.gen_bus[data.mask_ref .& data.mask_unshared]

        @testset "output contract" begin

            rng = HEDGeOPF._RND.MersenneTwister(1234)
            seq = [perturb_generation(pm, rng, k_min, k_max) for _ in 1:100]

            @test all(in.(length.(getfield.(seq, :ids_gen_faulted)), Ref(collect(k_min:k_max))))
            @test all(in.(unique(reduce(vcat, getfield.(seq, :ids_gen_faulted))), Ref(ids_gen)))
            @test all(in.(unique(reduce(vcat, getfield.(seq, :ids_bus_ref_valid))), Ref(ids_bus_ref_valid)))
        end

        @testset "reproducibility with fixed seed" begin
            n = 30
            rng1 = HEDGeOPF._RND.MersenneTwister(1234)
            rng2 = HEDGeOPF._RND.MersenneTwister(1234)

            seq1 = [perturb_generation(pm, rng1, k_min, k_max) for _ in 1:n]
            seq2 = [perturb_generation(pm, rng2, k_min, k_max) for _ in 1:n]

            @test seq1 == seq2
        end

        @testset "perturb generation with too many faulted generators" begin

            # All generators are faulted, so no eligible reference bus can be retained
            rng = HEDGeOPF._RND.MersenneTwister(1234)
            k_min = 5
            k_max = 5

            err = @test_throws ArgumentError perturb_generation(pm, rng, k_min, k_max)
            @test occursin("Unable to generate a valid generation perturbation", err.value.msg)
        end
    end

    @testset "test topology perturbation" begin

        setup["CASE"]["grid"] = "pglib_opf_case14_ieee.m"
        settings = to_namedtuple(setup)
        pm = deepcopy(DATA["14_ieee"]["pm"]);

        # Identify candidate reference buses
        data = define_candidate_ref_buses(pm)
        ids_gen = data.index
        ids_bus_ref_valid = data.gen_bus[data.mask_ref .& data.mask_unshared]

        @testset "perturb topology without islanding" begin

            k_min = 1
            k_max = 2
            rng1 = HEDGeOPF._RND.MersenneTwister(1234)
            rng2 = HEDGeOPF._RND.MersenneTwister(1234)
            ids_all_bus = sort(collect(_PM.ids(pm, :bus)))
            ids_all_branch = sort(collect(_PM.ids(pm, :branch)))

            topos1 = [perturb_topology(pm, rng1, k_min, k_max) for _ in 1:100]
            topos2 = [perturb_topology(pm, rng2, k_min, k_max) for _ in 1:100]
            # Test reproducibility with fixed seed
            @test topos1 == topos2
            # Test output contract
            @test all(isempty.(getfield.(topos1, :ids_bus)))
            @test all(isempty.(getfield.(topos1, :ids_gen)))
            @test all(length.(getfield.(topos1, :ids_ref)) .== 1)
            @test all(.>=(length.(getfield.(topos1, :ids_branch)), k_min))
            @test all(.<=(length.(getfield.(topos1, :ids_branch)), k_max))
            # Reference bus cannot change if islanding is not allowed
            @test all(first.(getfield.(topos1, :ids_ref)) .== 1)
            # Test consistency of returned ids with the model
            ids_branch_removed = unique(reduce(vcat, getfield.(topos1, :ids_branch)))
            @test all(in.(ids_branch_removed, Ref(ids_all_branch)))
        end

        @testset "perturb topology with faulted generator" begin

            # Fault at reference bus generator
            ids_bus_ref_valid_red = setdiff(ids_bus_ref_valid, [1])

            k_min = 1
            k_max = 3
            rng1 = HEDGeOPF._RND.MersenneTwister(1234)
            topos1 = [
                perturb_topology(pm, rng1, k_min, k_max; ids_bus_ref_valid=ids_bus_ref_valid_red)
                for _ in 1:100
            ]

            # Reference bus must be moved to other eligible generator (gen 2 in 14_ieee case)
            @test all(first.(getfield.(topos1, :ids_ref)) .== 2)
        end

        @testset "perturb topology with too many faulted branches" begin

            # All branches are faulted, so islanding is unavoidable
            k_min = 20
            k_max = 20
            rng1 = HEDGeOPF._RND.MersenneTwister(1234)
            err = @test_throws ArgumentError perturb_topology(pm, rng1, k_min, k_max)
            @test occursin("Unable to generate a valid topology perturbation", err.value.msg)
        end

        @testset "topology generator iteration contract" begin
            merge!(
                setup["TOPOLOGY"],
                Dict("k_min_gen" => 0, "k_max_gen" => 0, "k_min_branch" => 1, "k_max_branch" => 6, "num_topo" => 100)
            )
            settings = to_namedtuple(setup)
            rng = HEDGeOPF._RND.MersenneTwister(1234);
            gen = TopologyPerturbationGenerator(model = pm, rng = rng, setting = settings);

            # First iterate return always the original unperturbed topology
            t = first(iterate(gen))
            @test t.id == 1
            @test isempty(t.ids_branch)
            @test isempty(t.ids_bus)
            @test isempty(t.ids_ref)
            @test isempty(t.ids_gen)
            @test isempty(t.ids_gen_faulted)

            bounds_tot = vec(sum(get_pm_value(pm, :gen, ["pmin", "pmax"], Array{Any, 2}), dims = 1))
            @test t.pg_tot_bounds == bounds_tot

            t = first(iterate(gen))
            @test t.id == 2
            @test !isempty(t.ids_branch)
            @test !isempty(t.ids_ref)
        end

        @testset "topology generator iteration with fixed topology" begin
            merge!(
                setup["TOPOLOGY"],
                Dict("k_min_gen" => 0, "k_max_gen" => 2, "k_min_branch" => 1, "k_max_branch" => 6, "num_topo" => 0)
            )
            settings = to_namedtuple(setup)
            rng = HEDGeOPF._RND.MersenneTwister(1234);
            gen = TopologyPerturbationGenerator(model = pm, rng = rng, setting = settings);

            ts = []
            for t in gen
                push!(ts, t)
            end
            @test length(ts) == 1
            @test ts[1].id == 1
        end

        @testset "topology generator iteration with generation perturbation" begin
            merge!(
                setup["TOPOLOGY"],
                Dict("k_min_gen" => 0, "k_max_gen" => 2, "k_min_branch" => 1, "k_max_branch" => 6, "num_topo" => 1000)
            )
            settings = to_namedtuple(setup)
            rng = HEDGeOPF._RND.MersenneTwister(1234);
            gen = TopologyPerturbationGenerator(model = pm, rng = rng, setting = settings);
            bounds_pg = get_pm_value(pm, :gen, ["pmin", "pmax"], Array{Any, 2})

            check1, check2 = Bool[], Bool[]
            for t in gen
                bounds_pg_perturbed = vec(sum(bounds_pg[setdiff(ids_gen, t.ids_gen_faulted), :], dims=1))
                push!(check1, isequal(t.pg_tot_bounds, bounds_pg_perturbed))
                # Generator perturbations that result in no eligible ref bus should be discarded
                push!(check2, bounds_pg_perturbed[end] > 0.0)
            end
            @test all(check1)
            @test all(check2)
        end

        @testset "wrapper to generate topologies" begin
            merge!(
                setup["TOPOLOGY"], 
                Dict("k_min_gen" => 0, "k_max_gen" => 2, "k_min_branch" => 1, "k_max_branch" => 6, "num_topo" => 500)
            )
            settings = to_namedtuple(setup)
            rng = HEDGeOPF._RND.MersenneTwister(1234)
            gen = TopologyPerturbationGenerator(model = pm, rng = rng, setting = settings);

            num = 100
            topologies, mapping = generate_topologies(gen, num)

            # Test 1: correct number of topologies returned
            @test length(topologies) == num
            # Test 2: all topologies collected (including intact topology id=1)
            ids_topo = sort(getfield.(topologies, :id))
            @test ids_topo == sort(collect(1:num))
            # Test 3: all topologies appear in mapping
            ids_map = sort(reduce(vcat, collect(values(mapping))))
            @test ids_map == ids_topo
            # Test 4: each topology's bounds matches mapping key
            check = Bool[]
            lookup = Dict(getfield.(topologies, :id) .=> getfield.(topologies, :pg_tot_bounds))
            for (bounds_key, ids_in_group) in mapping
                for id in ids_in_group
                    push!(check, lookup[id] == bounds_key)
                end
            end
            @test all(check)
        end
    end
end
