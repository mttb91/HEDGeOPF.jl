
"Identify the largest and minor connected components in graph defined by `edges`, with edges at `ids` being removed"
function connected_components(edges::Matrix{Int}, ids::Vector{Int}, n::Int)

    ds = _DS.IntDisjointSets(n)

    ids = Set(ids)
    for i in axes(edges, 1)
        if !in(i, ids)
            _DS.union!(ds, edges[i, 1], edges[i, 2])
        end
    end

    components = Dict{Int, Vector{Int}}()
    for i in 1:n
        r = _DS.find_root!(ds, i)
        if !haskey(components, r)
            components[r] = Int[]
        end
        push!(components[r], i)
    end
    ccs = collect(values(components))
    lcc = popat!(ccs, argmax(length.(ccs)))
    return lcc, ccs
end
