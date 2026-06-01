
@testset "test output writing functions" begin

    @testset "_sql_string" begin

        @testset "backslashes converted to forward slashes" begin
            @test HEDGeOPF._sql_string("C:\\Users\\data") == "'C:/Users/data'"
        end

        @testset "single quotes escaped" begin
            @test HEDGeOPF._sql_string("it's a test") == "'it''s a test'"
        end

        @testset "backslash and single quote combined" begin
            @test HEDGeOPF._sql_string("C:\\path\\it's") == "'C:/path/it''s'"
        end
    end

    @testset "_var_lookup" begin

        @testset "detects components and variable names from csv-only dirs" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath(joinpath(tmp, "bus"))
                    mkpath(joinpath(tmp, "gen"))
                    touch(joinpath(tmp, "bus", "pg-1.csv"))
                    touch(joinpath(tmp, "bus", "pg-2.csv"))
                    touch(joinpath(tmp, "bus", "vm-1.csv"))
                    touch(joinpath(tmp, "gen", "qg-1.csv"))

                    result = HEDGeOPF._var_lookup()

                    @test haskey(result, "bus")
                    @test haskey(result, "gen")
                    @test sort(result["bus"]) == ["pg", "vm"]
                    @test result["gen"] == ["qg"]
                end
            end
        end

        @testset "ignores dirs that contain non-csv files" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath(joinpath(tmp, "bus"))
                    touch(joinpath(tmp, "bus", "pg-1.csv"))
                    touch(joinpath(tmp, "bus", "extra.zip"))    # impure dir

                    result = HEDGeOPF._var_lookup()
                    @test !haskey(result, "bus")
                end
            end
        end

        @testset "ignores empty directories" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath(joinpath(tmp, "empty_comp"))

                    result = HEDGeOPF._var_lookup()
                    @test !haskey(result, "empty_comp")
                end
            end
        end
    end

    @testset "_case_lookup" begin

        # Build a map with 2 workers, folds 1 and 2
        function _make_case_map()
            return HEDGeOPF._DF.DataFrame(;
                uid         = 1:8,
                worker      = [1, 1, 1, 1, 2, 2, 2, 2],
                case        = [1, 2, 3, 4, 1, 2, 3, 4],
                fold        = [1, 1, 2, 2, 1, 2, 2, 1],
                topology_id = ones(Int, 8)
            )
        end

        @testset "correct nesting: worker → fold → DataFrame with required columns" begin
            map    = _make_case_map()
            result = HEDGeOPF._case_lookup(map)

            @test haskey(result, 1) && haskey(result, 2)
            @test haskey(result[1], 1) && haskey(result[1], 2)
            @test haskey(result[2], 1) && haskey(result[2], 2)

            required_cols = [:case, :uid, :topology_id]
            col_present = [HEDGeOPF._DF.hasproperty(result[w][f], c)
                        for w in [1, 2], f in [1, 2], c in required_cols]
            @test all(col_present)
        end

        @testset "fold absent for worker when no rows match" begin
            map = HEDGeOPF._DF.DataFrame(;
                uid         = 1:4,
                worker      = [1, 1, 2, 2],
                case        = [1, 2, 1, 2],
                fold        = [1, 1, 1, 1],  # worker 2 has no fold 2
                topology_id = ones(Int, 4)
            )
            result = HEDGeOPF._case_lookup(map)

            @test haskey(result[1], 1)
            @test !haskey(result[2], 2)
        end

        @testset "uid values match the original map" begin
            map    = _make_case_map()
            result = HEDGeOPF._case_lookup(map)

            # Worker 1, fold 1 should contain uid 1 and 2
            uids_w1_f1 = sort(result[1][1].uid)
            @test uids_w1_f1 == [1, 2]
        end
    end

    @testset "_combine_cases" begin

        # Write a minimal Float32 CSV and return the cases Dict
        function _setup_combine(tmp; n_rows = 5, n_fold = 2)
            comp = "bus"
            var  = "vm"
            mkpath(joinpath(tmp, comp))

            # Worker 1 CSV: 5 rows, 2 float columns
            df = HEDGeOPF._DF.DataFrame(;
                v1 = Float32.(1:n_rows),
                v2 = Float32.(n_rows+1:2*n_rows)
            )
            CSV.write(joinpath(tmp, comp, "$(var)-1.csv"), df)

            # cases: worker 1, fold 1 → rows 1-3, fold 2 → rows 4-5
            cases = Dict(
                1 => Dict(
                    1 => HEDGeOPF._DF.DataFrame(; case = [1, 2, 3], uid = [1, 2, 3], topology_id = ones(Int, 3)),
                    2 => HEDGeOPF._DF.DataFrame(; case = [4, 5],    uid = [4, 5],    topology_id = ones(Int, 2)),
                )
            )
            return var, comp, cases
        end

        @testset "happy path: uid prepended, Float32 columns, correct row counts" begin
            mktempdir() do tmp
                cd(tmp) do
                    var, comp, cases = _setup_combine(tmp)
                    result = HEDGeOPF._combine_cases(var, comp, cases, 2)

                    @test haskey(result, 1) && haskey(result, 2)
                    @test HEDGeOPF._DF.nrow(result[1]) == 3
                    @test HEDGeOPF._DF.nrow(result[2]) == 2
                    @test HEDGeOPF._DF.hasproperty(result[1], :uid)
                    # uid is the first column
                    @test HEDGeOPF._DF.names(result[1])[1] == "uid"
                    # data columns are Float32
                    @test eltype(result[1].v1) == Float32
                    # uid values are unique within each fold
                    @test allunique(result[1].uid)
                    @test allunique(result[2].uid)
                end
            end
        end

        @testset "multiple workers merged into same fold, uid uniqueness preserved" begin
            mktempdir() do tmp
                cd(tmp) do
                    comp = "bus"
                    var  = "vm"
                    mkpath(joinpath(tmp, comp))

                    # Two workers, 4 rows each
                    for w in 1:2
                        df = HEDGeOPF._DF.DataFrame(;
                            v1 = Float32.(((w-1)*4+1):((w-1)*4+4))
                        )
                        CSV.write(joinpath(tmp, comp, "$(var)-$(w).csv"), df)
                    end

                    # Both workers contribute to both folds, with non-overlapping uids
                    cases = Dict(
                        1 => Dict(
                            1 => HEDGeOPF._DF.DataFrame(; case = [1, 2], uid = [1, 2],   topology_id = ones(Int, 2)),
                            2 => HEDGeOPF._DF.DataFrame(; case = [3, 4], uid = [3, 4],   topology_id = ones(Int, 2)),
                        ),
                        2 => Dict(
                            1 => HEDGeOPF._DF.DataFrame(; case = [1, 2], uid = [5, 6],   topology_id = ones(Int, 2)),
                            2 => HEDGeOPF._DF.DataFrame(; case = [3, 4], uid = [7, 8],   topology_id = ones(Int, 2)),
                        ),
                    )

                    result = HEDGeOPF._combine_cases(var, comp, cases, 2)

                    @test HEDGeOPF._DF.nrow(result[1]) == 4   # 2 from worker 1 + 2 from worker 2
                    @test HEDGeOPF._DF.nrow(result[2]) == 4
                    @test allunique(result[1].uid)
                    @test allunique(result[2].uid)
                    # all expected uids present
                    @test sort(result[1].uid) == [1, 2, 5, 6]
                    @test sort(result[2].uid) == [3, 4, 7, 8]
                end
            end
        end

        @testset "bounds check fires for out-of-range case index" begin
            mktempdir() do tmp
                cd(tmp) do
                    var, comp, _ = _setup_combine(tmp)
                    bad_cases = Dict(
                        1 => Dict(
                            1 => HEDGeOPF._DF.DataFrame(; case = [6], uid = [1], topology_id = [1]),  # row 6 out of 5
                            2 => HEDGeOPF._DF.DataFrame(; case = [1], uid = [2], topology_id = [1]),
                        )
                    )
                    @test_throws ErrorException HEDGeOPF._combine_cases(var, comp, bad_cases, 2)
                end
            end
        end

        @testset "empty-fold guard fires when no worker contributes to a fold" begin
            mktempdir() do tmp
                cd(tmp) do
                    var, comp, _ = _setup_combine(tmp)
                    # Only fold 1 has data; fold 2 is completely absent from cases
                    cases_no_fold2 = Dict(
                        1 => Dict(
                            1 => HEDGeOPF._DF.DataFrame(; case = [1, 2], uid = [1, 2], topology_id = ones(Int, 2)),
                        )
                    )
                    err = @test_throws ErrorException HEDGeOPF._combine_cases(var, comp, cases_no_fold2, 2)
                    @test occursin("fold 2", err.value.msg)
                end
            end
        end
    end

    @testset "generate_uid" begin

        # Write two synthetic info CSVs that generate_uid expects
        function _write_info_csvs(tmp)
            for worker in 1:2
                df = HEDGeOPF._DF.DataFrame(;
                    pd_tot      = Float64.(worker:2:(worker + 4)),
                    objective   = ones(Float64, 3),
                    topology_id = ones(Int, 3)
                )
                CSV.write(joinpath(tmp, "info-$(worker).csv"), df)
            end
        end

        @testset "output contract: map.csv written with correct structure" begin
            mktempdir() do tmp
                cd(tmp) do
                    _write_info_csvs(tmp)
                    generate_uid(false)

                    @test isfile(joinpath(tmp, "map.csv"))
                    map = HEDGeOPF._DF.DataFrame(CSV.File(joinpath(tmp, "map.csv")))

                    # Required columns exist
                    for col in [:uid, :worker, :case, :pd_tot, :objective, :topology_id]
                        @test HEDGeOPF._DF.hasproperty(map, col)
                    end
                    # UIDs are unique and span 1:nrow
                    @test allunique(map.uid)
                    @test sort(map.uid) == collect(axes(map, 1))
                    # 2 workers × 3 rows each
                    @test HEDGeOPF._DF.nrow(map) == 6
                end
            end
        end

        @testset "uid order matches sort by pd_tot, objective, topology_id" begin
            mktempdir() do tmp
                cd(tmp) do
                    _write_info_csvs(tmp)
                    generate_uid(false)

                    map = HEDGeOPF._DF.DataFrame(CSV.File(joinpath(tmp, "map.csv")))
                    HEDGeOPF._DF.sort!(map, [:pd_tot, :objective, :topology_id])
                    @test map.uid == collect(axes(map, 1))
                end
            end
        end

        @testset "cleanup=true removes info files, keeps map.csv" begin
            mktempdir() do tmp
                cd(tmp) do
                    _write_info_csvs(tmp)
                    generate_uid(true)

                    @test isfile(joinpath(tmp, "map.csv"))
                    info_files = filter(x -> contains(x, "info"), readdir(tmp))
                    @test isempty(info_files)
                end
            end
        end

        @testset "cleanup=false keeps info files" begin
            mktempdir() do tmp
                cd(tmp) do
                    _write_info_csvs(tmp)
                    generate_uid(false)

                    info_files = filter(x -> contains(x, "info"), readdir(tmp))
                    @test length(info_files) == 2
                end
            end
        end
    end
end