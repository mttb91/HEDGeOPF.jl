
"Sample uniformly from convex polytope using `R` library `volesti`"
function _sample_polytope(A::AbstractMatrix{Float64}, b::AbstractVector{Float64}, n_samples::Int,
    walk::String,
    walk_length::Int,
    nburns::Int,
    x0::AbstractVector{Float64},
    seed::Int32)

    return convert(Matrix{Float64},
        RCall.R"""
        library(volesti)
        P = Hpolytope(A = $A, b = $b);
        sample_points(P, n = $n_samples, random_walk = list(
            "seed" = $seed,
            "walk" = $walk, 
            "nburns" = $nburns, 
            "walk_length" = $walk_length,
            "starting_point" = $x0
            )
        )
        """
    )
end

"Sample uniformly in convex polytope `Ax <= b` with an already-defined Chebyshev center by calling `R` package `volesti`"
function sample_polytope(A::AbstractMatrix{Float64}, b::AbstractVector{Float64}, x0::AbstractVector{Float64}, n_samples::Int;
    nburns::Int = 0,
    walk::String = "CDHR",
    walk_length::Union{Int, Nothing} = nothing,
    seed::Int32 = _RND.rand(Int32)
    )
    # Set walk length as `10 + 10*size(A, 2)` if no value is provided
    walk_length = isnothing(walk_length) ? 10 + 10 * size(A, 2) : walk_length

    return _sample_polytope(A, b, n_samples, walk, walk_length, nburns, x0, seed)
end