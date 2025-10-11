
Base.@kwdef mutable struct TopologyPerturbationGenerator
    model::_PM.AbstractPowerModel
    rng::_RND.AbstractRNG
    setting::NamedTuple
    state::Int = 1
end

Base.@kwdef struct TopologyPerturbation
    id::Int
    ids_branch::Vector{Int} = Int[]
    ids_bus::Vector{Int} = Int[]
    ids_ref::Vector{Int} = Int[]
    ids_gen::Vector{Int} = Int[]
    ids_gen_faulted::Vector{Int} = Int[]
    pg_tot_bounds::Vector{Float64}
end

Base.@kwdef mutable struct InputData
    data::Vector{Float64}
    ids::Union{Vector{Int}, Nothing} = nothing
end

Base.@kwdef mutable struct InputSample
    topology::TopologyPerturbation
    pd_tot::Vector{Float64} = Float64[]
    data::Dict{String, Dict{String, InputData}} = Dict{String, Dict{String, InputData}}()
end

const PolyType = NamedTuple{
    (:A, :b, :ids),
    Tuple{Matrix{Float64}, Vector{Float64}, Tuple{Vector{Int}, Vector{Int}}}
}
