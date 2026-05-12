
@testset "test functions to modify/process PowerModels model in and out-of-place" begin

    setup = deepcopy(SETUP)

    @testset "out-of-place modifications" begin

        setup["CASE"]["grid"] = "pglib_opf_case14_ieee.m"
        settings = to_namedtuple(setup)
        pm = deepcopy(DATA["14_ieee"]["pm"]);
        buses = sort(collect(_PM.ids(pm, :bus)))
        bus_gen = vec(get_pm_value(pm, :gen, ["gen_bus"], Array{Any, 2}))
        bus_type = vec(get_pm_value(pm, :bus, ["bus_type"], Array{Any, 2}))
        bus_ref_old = findfirst(==(3), bus_type)
        baseMVA = _PM.ref(pm, 0, :baseMVA)

        @testset "identify generators sharing buses" begin

            # No generator shares bus in 14_ieee case
            @test all(is_node_shared(pm, bus_gen) .== false)
            # Generators 1 and 2 share bus in modified 14_ieee case
            net = deepcopy(DATA["14_ieee"]["net"]);
            net["gen"]["2"]["gen_bus"] = 1
            pm_new = instantiate_model(net, "ACP", settings);
            bus_gen_new = vec(get_pm_value(pm_new, :gen, ["gen_bus"], Array{Any, 2}))
            @test isequal(is_node_shared(pm_new, bus_gen_new), [true, true, false, false, false])
        end

        @testset "identify candidate reference buses" begin
            
            data = define_candidate_ref_buses(pm)
            ids_bus_ref_valid = data.gen_bus[data.mask_ref .& data.mask_unshared]
            # In 14_ieee only first two generators are eligible
            @test isequal(ids_bus_ref_valid, [1, 2])
            # Generators [1, 2] are faulted -> no valid reference bus
            data = define_candidate_ref_buses(pm; ids_gen_faulted = [1, 2])
            ids_bus_ref_valid = data.gen_bus[data.mask_ref .& data.mask_unshared]
            @test isempty(ids_bus_ref_valid)
        end

        @testset "define reference bus for a subgraph" begin

            @testset "existing reference bus in subgraph is retained" begin
                data = define_candidate_ref_buses(pm)
                bus_gen_valid = data.gen_bus[data.mask_ref .& data.mask_unshared]
                @test define_ref_bus(pm, buses, bus_gen_valid) == bus_ref_old
            end
            @testset "existing reference bus in subgraph is faulted -> picks first candidate gen bus" begin
                data = define_candidate_ref_buses(pm; ids_gen_faulted = [1])
                bus_gen_valid = data.gen_bus[data.mask_ref .& data.mask_unshared]
                @test define_ref_bus(pm, buses, bus_gen_valid) == first(bus_gen_valid)
            end
            @testset "no reference bus in subgraph -> picks first candidate gen bus" begin
                buses_sel = setdiff(buses, [bus_ref_old])
                @test define_ref_bus(pm, buses_sel, bus_gen) == first(intersect(bus_gen, buses_sel))
            end
            @testset "no reference bus and no gen bus in subgraph -> falls back to first bus" begin
                buses_sel = setdiff(buses, bus_gen)
                @test define_ref_bus(pm, buses_sel, Int[]) == first(buses_sel)
            end
        end

        @testset "generate new reference model aligned with topology perturbation" begin

            get_sol = sol -> Dict(
                var => get_pm_value(sol, comp, [var], Array{Any, 2})
                for (var, comp) in zip(["va", "vm", "pg", "qg"], ["bus", "bus", "gen", "gen"])
            )

            pg_bounds = vec(sum(get_pm_value(pm, :gen, ["pmin", "pmax"], Array{Any, 2}), dims=1))
            t = TopologyPerturbation(
                id = 1,
                ids_branch = [8, 9, 10],
                ids_bus = collect(6:14),
                ids_ref = [2, 6],
                ids_gen = [4, 5],
                ids_gen_faulted = [1],
                pg_tot_bounds = pg_bounds
            )
            pm_new = update_topology(pm, t);

            @testset "bus type consistency" begin
                bus_type = vec(get_pm_value(pm_new, :bus, ["bus_type"], Array{Any, 2}))

                @test isequal(findall(x -> x == 3, bus_type), t.ids_ref)
                # Only bus 3 has a generator that is still active and non-slack
                @test isequal(findall(x -> x == 2, bus_type), [3])
                # All remaining bus are PQ-type (load only or de-energized)
                @test isequal(findall(x -> x == 1, bus_type), setdiff(buses, [t.ids_ref; 3]))
            end

            @testset "policy equivalence" begin
                # Modify voltage magnitude bounds of disconnected buses
                mask = vec(get_pm_value(pm_new, :bus, ["is_connected"], Array{Any, 2}))
                set_pm_value!(pm_new, :bus, ["vmin", "vmax"], 0.0; mask = findall(x -> x == 0, mask))
                pm_new = instantiate_model(pm_new.data, "ACP", settings);

                sol = _PM.optimize_model!(pm_new)["solution"]
                sol_dict = get_sol(sol)

                pm_t = deepcopy(pm);
                update_topology!(pm_t, t.ids_branch, t.ids_bus, t.ids_ref, sort(vcat(t.ids_gen, t.ids_gen_faulted)))
                sol_t = _PM.optimize_model!(pm_t)["solution"]
                sol_t_dict = get_sol(sol_t)

                # Deactivated components are not present in `sol_dict``
                @test all(isapprox.(sol_dict["vm"], sol_t_dict["vm"]; atol = 1e-8))
                @test all(isapprox.(
                    sol_dict["va"][setdiff(buses, t.ids_bus)],
                    sol_t_dict["va"][setdiff(buses, t.ids_bus)]; atol = 1e-8)
                )
                for var in ["pg", "qg"]
                    @test all(isapprox.(
                        sol_dict[var],
                        sol_t_dict[var][setdiff(collect(1:length(sol_t_dict["pg"])), t.ids_gen)];
                        atol = 1e-6 / baseMVA)
                    )
                end
            end
        end
    end

    @testset "in-place modifications" begin

        setup["CASE"]["grid"] = "pglib_opf_case5_pjm.m"
        pm = deepcopy(DATA["5_pjm"]["pm"]);
        ids_bus = sort(collect(_PM.ids(pm, :bus)))
        baseMVA = _PM.ref(pm, 0, :baseMVA)

        @testset "update reference bus(es)" begin
            pm_new = deepcopy(pm);
            var = _PM.var(pm_new, 0, :va)
            # Test 1: reference bus does not change
            ids_fixed = [4]
            fix_theta_ref!(pm_new, ids_fixed)
            @test JuMP.is_fixed.(_PM.var(pm, 0, :va)) == JuMP.is_fixed.(_PM.var(pm_new, 0, :va))
            # Test 2: multiple assignment with partial change
            ids_fixed = [1, 4]
            fix_theta_ref!(pm_new, ids_fixed)
            @test all(JuMP.is_fixed.(var)[ids_fixed])
            @test !any(JuMP.is_fixed.(var)[setdiff(ids_bus, ids_fixed)])
            @test JuMP.is_fixed.(_PM.var(pm, 0, :va)) != JuMP.is_fixed.(_PM.var(pm_new, 0, :va))
            # Test 3: complete change of reference bus
            ids_fixed = [3]
            fix_theta_ref!(pm_new, ids_fixed)
            res = _PM.optimize_model!(pm_new)
            @test JuMP.is_fixed.(var)[first(ids_fixed)]
            @test isapprox(JuMP.value.(var)[first(ids_fixed)], 0.0; atol = 1e-8)
            @test !isequal(pm, pm_new)
            # Test 4: if no reference index is provided, current reference bus is retained
            fix_theta_ref!(pm_new, Int[])
            @test JuMP.is_fixed.(var)[first(ids_fixed)]
            @test !any(JuMP.is_fixed.(var)[setdiff(ids_bus, ids_fixed)])
        end

        @testset "fix variable to zero" begin
            pm_new = deepcopy(pm);
            pg_bounds = get_pm_value(pm_new, :gen, ["pmin", "pmax"], Array{Any, 2})
            nm = :pg

            # Case 1: fix to zero in-place with no variable previously fixed
            ids_fixed = [3]
            ids_bounded = setdiff(collect(1:size(pg_bounds, 1)), ids_fixed)
            var = _PM.var(pm_new, 0, nm)
            @test all(JuMP.has_lower_bound.(var))
            @test all(JuMP.has_upper_bound.(var))
            fix_to_zero!(pm_new, ids_fixed, :gen, nm, ["pmin", "pmax"])
            @test JuMP.is_fixed.(var)[first(ids_fixed)]
            @test all(JuMP.has_lower_bound.(var)[ids_bounded])
            @test all(JuMP.has_upper_bound.(var)[ids_bounded])
            # Case 2: fix to zero in-place while reconstructing bounds for other variable
            ids_fixed_new = [1, 2]
            ids_bounded_new = setdiff(collect(1:size(pg_bounds, 1)), ids_fixed_new)
            fix_to_zero!(pm_new, ids_fixed_new, :gen, nm, ["pmin", "pmax"])
            @test all(JuMP.is_fixed.(var)[ids_fixed_new])
            @test isequal(JuMP.lower_bound.(var[ids_bounded_new]).data, pg_bounds[ids_bounded_new, 1])
            @test isequal(JuMP.upper_bound.(var[ids_bounded_new]).data, pg_bounds[ids_bounded_new, 2])
            # Ref dictionary is never modified
            @test pg_bounds == get_pm_value(pm_new, :gen, ["pmin", "pmax"], Array{Any, 2})
            # A single bound is considered as symmetrical around zero
            fix_to_zero!(pm_new, ids_fixed, :gen, nm, ["pmax"])
            ids_bounded_new = setdiff(ids_bounded_new, ids_fixed)
            # The fixed variable is now bounded with pmax as upper bound and -pmax as lower bound
            @test isequal(JuMP.lower_bound.(var[ids_fixed_new]).data, -pg_bounds[ids_fixed_new, 2])
            @test isequal(JuMP.upper_bound.(var[ids_fixed_new]).data, pg_bounds[ids_fixed_new, 2])
            # Variables previously bounded are not modified
            @test isequal(JuMP.lower_bound.(var[ids_bounded_new]).data, pg_bounds[ids_bounded_new, 1])
            @test isequal(JuMP.upper_bound.(var[ids_bounded_new]).data, pg_bounds[ids_bounded_new, 2])
            # Case 3: error for unsupported variable
            @test_throws AssertionError fix_to_zero!(pm_new, ids_fixed, :branch, :pf, ["rate_a"])
        end

        @testset "update branch status parameter" begin
            pm_new = deepcopy(pm);
            ids_branch = sort(collect(_PM.ids(pm_new, :branch)))

            nm, nw = :br_status, 0
            var = _PM.var(pm_new, nw, nm)
            # Case 1: initial state is fully intact
            @test all(isone.(JuMP.parameter_value.(var.data)))
            # Case 2: single branch is out-of-service
            ids_faulted = [1]
            ids_active = setdiff(ids_branch, ids_faulted)
            updated = update_branch_status!(pm_new, ids_faulted)
            values = get_pm_value(
                _PM.optimize_model!(pm_new)["solution"],
                "branch", ["pf", "pt", "qf", "qt"], Array{Any, 2}
            )
            @test updated
            @test iszero(JuMP.parameter_value(var[first(ids_faulted)]))
            @test all(isone.(JuMP.parameter_value.(var[ids_active])))
            @test all(isapprox.(values[first(ids_faulted), :], 0.0; atol = 1e-6 / baseMVA))
            @test all(any(.!isapprox.(values[ids_active, :], 0.0; atol = 1e-6 / baseMVA) , dims=2))
            # Case 2: other topology change with same faulted branch -> no update
            updated = update_branch_status!(pm_new, ids_faulted)
            @test !updated
            # Case 3: new faulted branches while previous one is fixed
            ids_faulted = [2, 5]
            ids_active = setdiff(ids_branch, ids_faulted)
            updated = update_branch_status!(pm_new, ids_faulted)
            values = get_pm_value(
                _PM.optimize_model!(pm_new)["solution"],
                "branch", ["pf", "pt", "qf", "qt"], Array{Any, 2}
            )
            @test updated
            @test all(iszero.(JuMP.parameter_value.(var[ids_faulted])))
            @test all(isone.(JuMP.parameter_value.(var[ids_active])))
            @test all(isapprox.(values[first(ids_faulted), :], 0.0; atol = 1e-6 / baseMVA))
            @test all(any(.!isapprox.(values[ids_active, :], 0.0; atol = 1e-6 / baseMVA) , dims=2))
        end

        @testset "update topology logic" begin

            get_sol = sol -> Dict(
                var => get_pm_value(sol, comp, [var], Array{Any, 2})
                for (var, comp) in zip(["va", "vm", "pg", "pf"], ["bus", "bus", "gen", "branch"])
            )

            setup["CASE"]["grid"] = "pglib_opf_case14_ieee.m"
            settings = to_namedtuple(setup)
            pm = deepcopy(DATA["14_ieee"]["pm"]);
            baseMVA = _PM.ref(pm, 0, :baseMVA)
            sol_ref = _PM.optimize_model!(pm)["solution"]
            sol_ref_dict = get_sol(sol_ref)
            pm_new = deepcopy(pm);

            @testset "intact topology (same reference bus) -> no update" begin
                ids_ref = [1]
                ids_branch, ids_bus, ids_gen = Int[], Int[], Int[]
                update_topology!(pm_new, ids_branch, ids_bus, ids_ref, ids_gen)
                sol = _PM.optimize_model!(pm_new)["solution"]
                sol_dict = get_sol(sol)
                for var in keys(sol_ref_dict)
                    @test all(isapprox.(sol_dict[var], sol_ref_dict[var]; atol = 1e-6 / baseMVA))
                end
            end
            @testset "intact topology with faulted generator and new reference bus" begin
                ids_gen = [1]
                ids_ref = [2]
                ids_branch, ids_bus = Int[], Int[]
                update_topology!(pm_new, ids_branch, ids_bus, ids_ref, ids_gen)
                sol = _PM.optimize_model!(pm_new)["solution"]
                sol_dict = get_sol(sol)
                @test isapprox(sol_dict["va"][first(ids_ref)], 0.0; atol = 1e-8)
                @test isapprox(sol_dict["pg"][first(ids_gen)], 0.0; atol = 1e-6 / baseMVA)
            end
            @testset "faulty topology resulting in single connected component with original slack" begin
                ids_ref_old = 2
                ids_ref = [1]
                ids_branch = [1, 4]
                ids_gen, ids_bus = Int[], Int[]
                update_topology!(pm_new, ids_branch, ids_bus, ids_ref, ids_gen)
                sol = _PM.optimize_model!(pm_new)["solution"]
                sol_dict = get_sol(sol)
                @test isapprox(sol_dict["va"][first(ids_ref)], 0.0; atol = 1e-8)
                @test !iszero(sol_dict["va"][ids_ref_old])
                @test all(isapprox.(sol_dict["pf"][ids_branch], 0.0; atol = 1e-6 / baseMVA))
            end
            @testset "faulty topology resulting in two connected components, one of which unfeasible" begin
                ids_branch = [8, 9, 10]
                ids_bus = collect(6:14)
                ids_gen = [4, 5]
                ids_ref = [1, 6]
                update_topology!(pm_new, ids_branch, ids_bus, ids_ref, ids_gen)
                sol = _PM.optimize_model!(pm_new)["solution"]
                sol_dict = get_sol(sol)
                @test JuMP.is_solved_and_feasible(pm_new.model)
                @test all(isapprox.(sol_dict["vm"][ids_bus], 0.0; atol = 1e-8))
                @test all(isapprox.(sol_dict["pf"][ids_branch], 0.0; atol = 1e-6 / baseMVA))
                @test all(isapprox.(sol_dict["pf"][11:end], 0.0; atol = 1e-6 / baseMVA))
            end
        end

        @testset "update power system setpoints" begin

            get_sol = sol -> Dict(
                var => get_pm_value(sol, comp, [var], Array{Any, 2})
                for (var, comp) in zip(["pg", "qg", "pd", "qd"], ["gen", "gen", "load", "load"])
            )

            setup["CASE"]["grid"] = "pglib_opf_case5_pjm.m"
            settings = to_namedtuple(setup)
            pm = deepcopy(DATA["5_pjm"]["pm"]);
            baseMVA = _PM.ref(pm, 0, :baseMVA)
            sol_ref = _PM.optimize_model!(pm)["solution"]
            sol_ref_dict = get_sol(sol_ref)

            # Relevant input data
            sd = get_pm_value(pm, :load, ["pd", "qd"], Array{Any, 2})
            c1 = vec(first.(get_pm_value(pm, :gen, ["cost"], Array{Any, 2})))
            bounds_pg = vec(sum(get_pm_value(pm, :gen, ["pmin", "pmax"], Array{Any, 2}), dims = 1))

            @testset "update load active/reactive power" begin

                pm_new = deepcopy(pm);
                sd_new = copy(sd)
                id_pd, id_qd = 1, 3
                sd_new[id_pd, 1] += 100 / baseMVA
                sd_new[id_qd, 2] -= 50 / baseMVA
                sample = Dict(
                    "load" => Dict(
                        "pd" => InputData(sd_new[[id_pd], 1], [id_pd]),
                        "qd" => InputData(sd_new[[id_qd], 2], [id_qd])
                    )
                )
                update_load_power!(pm_new, sample, bounds_pg)
                sol = _PM.optimize_model!(pm_new)["solution"]
                sol_dict = get_sol(sol)
                ids = sort(first(_PM.var(pm_new, 0, :pd_fix).axes))

                @test JuMP.is_solved_and_feasible(pm_new.model)
                @test !isequal(sd, sd_new)
                @test isequal(sd_new, reduce(hcat, [JuMP.value.(_PM.var(pm_new, 0, h, ids).data) for h in [:pd_fix, :qd_fix]]))
                @test sum(sol_dict["pg"]) > sum(sd_new[:, 1])
                @test sum(sol_dict["pg"]) > sum(sol_ref_dict["pg"])
                @test sum(sol_dict["qg"]) < sum(sol_ref_dict["qg"])
            end

            @testset "update generator linear cost" begin
                pm_new = deepcopy(pm);
                c1_new = copy(c1)
                # Case 1: permute costs to make gen 2 out-of-market 
                ids = [2, 4]
                c1_new[ids] = c1[reverse(ids)]
                sample = Dict(
                    "gen" => Dict("c1" => InputData(c1_new[ids], ids))
                )
                update_gen_cost!(pm_new, sample)
                sol = _PM.optimize_model!(pm_new)["solution"]
                sol_dict = get_sol(sol)

                @test JuMP.is_solved_and_feasible(pm_new.model)
                @test !isequal(c1, c1_new)
                @test isapprox(sol_ref_dict["pg"][last(ids)], 0.0; atol = 1e-6)
                @test isapprox(sol_dict["pg"][first(ids)], 0.0; atol = 1e-6)
                @test sol_ref_dict["pg"][first(ids)] > 0.0
                @test sol_dict["pg"][last(ids)] > 0.0
            end
        end
    end
end
