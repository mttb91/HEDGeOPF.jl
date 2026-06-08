
@testset "test graph connected components function" begin

    ref_edges = Int[
        1 2
        2 3
        3 4
        4 5
        5 6
        3 7
    ]
    normalize_component(c::Vector{Int}) = sort(copy(c))
    normalize_partition(parts::Vector{Vector{Int}}) = sort(normalize_component.(parts); by = first)

    @testset "intact graph has one connected component" begin
        edges = ref_edges[begin:end-2, :]
        lcc, ccs = connected_components(edges, Int[], 5)

        @test sort(lcc) == [1, 2, 3, 4, 5]
        @test isempty(ccs)
    end

    @testset "single edge removal creates one minor island" begin
        edges = ref_edges[begin:end-1, :]
        edges[end, 1] = 3
        lcc, ccs = connected_components(edges, [2], 6)

        @test sort(lcc) == [3, 4, 5, 6]
        @test normalize_partition(ccs) == [[1, 2]]
    end

    @testset "multiple removals create multiple islands" begin
        lcc, ccs = connected_components(ref_edges, [2, 4], 7)

        @test sort(lcc) == [3, 4, 7]
        @test normalize_partition(ccs) == [[1, 2], [5, 6]]
    end

    @testset "removed edge ids are order- and duplicate-invariant" begin
        lcc1, ccs1 = connected_components(ref_edges, [2, 4], 7)
        lcc2, ccs2 = connected_components(ref_edges, [4, 2, 2], 7)

        @test sort(lcc1) == sort(lcc2)
        @test normalize_partition(ccs1) == normalize_partition(ccs2)
    end
end
