using HEDGeOPF
using LinearAlgebra
using Test

import HiGHS
import Ipopt
import HEDGeOPF: _DC, _PM, _DF, JuMP, RCall
import HEDGeOPF: _isempty

@testset "HEDGeOPF" begin

    _PM.silence()

    global SETUP = include("settings.jl")

    # OPF model

    global DATA = include("network.jl")
    include("model.jl")
    include("opf.jl")
    include("solution.jl")

    # Graph

    include("graph.jl")

    # Polytope

    include("polytope.jl")
    include("rcall.jl")

    # Sampling

    include("parallel.jl")
    include("sampling.jl")
    
end
