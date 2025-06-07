module HEDGeOPF

using LinearAlgebra

import PowerModels as _PM
import JuMP

import Random as _RND
import DataFrames as _DF
import Distributions as _DS
import KernelDensity: kde_lscv
import Serialization: serialize, deserialize
import CSV
import XLSX
import YAML
import RCall

import Distributed as _DC
import SharedArrays as _SA
import ProgressMeter as _PMT

function __init__()
    BLAS.set_num_threads(1)
end

const PolyType = NamedTuple{
    (:A, :b, :ids),
    Tuple{Matrix{Float64}, Vector{Float64}, Tuple{Vector{Int}, Vector{Int}}}
}

## OPF problem

include("model/core/base.jl")
include("model/core/constraint_template.jl")
include("model/core/constraint.jl")
include("model/core/data_basic.jl")
include("model/core/data.jl")
include("model/core/objective.jl")
include("model/core/solution.jl")
include("model/core/variable.jl")

include("model/form/acp.jl")
include("model/prob/opf.jl")

# Sampling

include("sampling/base.jl")
include("sampling/data.jl")
include("sampling/parallel.jl")
include("sampling/polytope.jl")
include("sampling/volesti.jl")

## I/O

include("io/miscellaneous.jl")
include("io/results.jl")
include("io/settings.jl")

_PM.silence()

## Exports

# Exports everything except internal symbols, which are defined as those whose name
# starts with an underscore. If you don't want all of these symbols in your environment,
# then use `import` instead of `using`. Do not add custom-defined symbols to this exclude
# list. Instead, rename them with an underscore.

const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__, all=true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS || startswith(sym_string, "_") || startswith(sym_string, "@_")
        continue
    end
    if !(Base.isidentifier(sym) || (startswith(sym_string, "@") &&
         Base.isidentifier(sym_string[2:end])))
       continue
    end
    @eval export $sym
end


end
