
@testset "test functions to read, write and process PowerModels nested dictionaries" begin

    setup = deepcopy(SETUP)
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

    @testset "compute load active bounds" begin
        net1 = deepcopy(net)
        net1["load"]["1"]["pd"] *= -1

        vars, bounds = compute_load_active_bounds(net1, setup["SAMPLING"])
        pd = vec(get_pm_value(net1, "load", ["pd"], Array{Any, 2}))
        delta = setup["SAMPLING"]["delta_pd"] / 100
        expected = sort(repeat(pd, 1, 2) .* [1 - delta 1 + delta], dims = 2)

        @test vars == ["pmin", "pmax"]
        @test size(bounds) == (length(pd), 2)
        @test all(bounds[:, 1] .<= bounds[:, 2])
        @test isapprox(bounds, expected; atol = 1e-12)
    end

    @testset "compute load reactive bounds" begin
        id_pd, id_qd = 1, 2
        net1 = deepcopy(net)
        net1["load"][string(id_pd)]["pd"] *= -1
        net1["load"][string(id_qd)]["qd"] *= -1

        vars, bounds = compute_load_reactive_bounds(net1, setup["SAMPLING"])
        pd = vec(get_pm_value(net1, "load", ["pd"], Array{Any, 2}))
        qd = vec(get_pm_value(net1, "load", ["qd"], Array{Any, 2}))
        delta = setup["SAMPLING"]["delta_qd"] / 100
        expected_q = sort(repeat(qd, 1, 2) .* [1 - delta 1 + delta], dims = 2)

        @test vars == ["qmin", "qmax", "qp_ratio_min", "qp_ratio_max"]
        @test size(bounds) == (length(qd), 4)
        @test !isequal(sign(maximum(bounds[id_pd, 3:4])), sign(pd[id_pd]))
        @test all(isequal.(sign.(maximum(bounds[:, 3:4], dims=2)), sign.(qd)))
        @test all(bounds[:, 1] .<= bounds[:, 2])
        @test all(bounds[:, 3] .<= bounds[:, 4])
        @test isapprox(bounds[:, 1:2], expected_q; atol = 1e-12)
        @test all(isfinite.(bounds))
    end

    @testset "detect synchronous condensers" begin
        id_sync = 3
        net1 = deepcopy(net)
        net1["gen"][string(id_sync)]["pmax"] = 0.0
        set_pm_value!(net1, "gen", ["pmax", "qmax", "qmin"], 0.0; mask = [4])
        gen_ids = vec(get_pm_value(net1, "gen", ["index"], Array{Any, 2}))
        vars, mask = is_synchronous(net1)

        @test vars == ["is_synchronous"]
        @test isa(mask, BitVector) || isa(mask, Vector{Bool})
        @test all(isequal.(gen_ids[mask], [id_sync]))
        @test all(isequal.(gen_ids[.!mask], setdiff(gen_ids, [id_sync])))
    end

    @testset "derivation of power system matrices" begin
        
        net_new = deepcopy(net)
        pm_new = deepcopy(pm)
        Ybus_pm = _PM.calc_admittance_matrix(net_new).matrix

        buses = vec(get_pm_value(pm_new, :bus, ["bus_i"], Array{Any, 2}))
        n_bus = length(buses)

        @testset "admittance matrices with intact system" begin

            data = Dict{Symbol, _DF.DataFrame}()
            for key in keys(filter(p->isa(p.second, Dict{Int, Any}) && !isempty(p.second), _PM.ref(pm_new)))
                vars = get_pm_key(pm_new, key)
                data[key] = get_pm_value(pm_new, key, vars, _DF.DataFrame)
            end
            bus_to_idx = buses .=> 1:length(buses)
        
            if haskey(data, :shunt)
                data_shunt = data[:shunt][!, ["gs", "bs"]]
            else
                data_shunt = _DF.DataFrame()
            end
            data_branch = data[:branch][!, ["b_fr", "b_to", "br_r", "br_x", "g_fr", "g_to", "tap", "shift", "br_status"]]

            indices = HEDGeOPF._calc_connection_indices(data, bus_to_idx, :_index)
            Ybus, Yf, Yt = HEDGeOPF._calc_admittance_matrices(n_bus, data_branch, data_shunt, indices)
            @test isapprox(Ybus, Ybus_pm; rtol=1e-10, atol=1e-10)
            @test isequal(size(Yf, 1), sum(data_branch[!, "br_status"]))
            @test isequal(Ybus, transpose(Ybus))
        end
        @testset "admittance matrices with topology outage" begin
            
            id = 2
            net_new["branch"][string(id)]["br_status"] = 0
            set_pm_value!(pm_new, :branch, ["br_status"], 0; mask = [id])
            br_status = vec(get_pm_value(pm_new, :branch, ["br_status"], Array{Any, 2}))
            Ybus_pm = _PM.calc_admittance_matrix(net_new).matrix

            indices = calc_connection_indices(pm_new)
            Ybus, Yf, Yt = calc_admittance_matrices(pm_new, indices)
            @test isapprox(Ybus, Ybus_pm; rtol=1e-10, atol=1e-10)
            @test isequal(length(Ybus.nzval), sum(br_status) * 2 + n_bus)
            @test isequal(size(Yf, 1), sum(br_status))
            @test isequal(Ybus, transpose(Ybus))
        end
    end

    data = Dict{String, Dict{String, Any}}()
    setup = deepcopy(SETUP)
    for case in ["5_pjm", "14_ieee", "30_ieee"]
        setup["CASE"]["grid"] = "pglib_opf_case$(case).m"
        data[case] = Dict("net" => instantiate_network(setup))
        data[case]["pm"] = instantiate_model(data[case]["net"], "ACP", to_namedtuple(setup))
    end
    return data
end

