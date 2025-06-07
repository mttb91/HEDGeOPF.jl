using HEDGeOPF
using LinearAlgebra
using Test

import HiGHS
import Ipopt
import HEDGeOPF: _DC, _PM, _DF, JuMP, RCall
import HEDGeOPF: _isempty

@testset "HEDGeOPF" begin

    global SETUP = include("settings.jl")

    # OPF model

    include("data.jl")
    include("opf.jl")
    include("solution.jl")

    # Polytope

    include("polytope.jl")
    include("rcall.jl")

    # Sampling

    include("parallel.jl")
    
end
