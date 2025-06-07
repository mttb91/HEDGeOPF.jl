# HEDGeOPF.jl Documentation

```@meta
CurrentModule = HEDGeOPF
```

## Overview

[HEDGeOPF.jl](https://github.com/mttb91/HEDGeOPF.jl) is a Julia package to generate datasets of AC Optimal Power Flow (AC-OPF) instances for standardized training and testing of Neural Networks (NNs) that learn to approximate this problem. It implements a methodology that delivers high-quality datasets without compromising on efficiency and scalability. Please refer to our [publication](https://ieeexplore.ieee.org/abstract/document/10761586) for details on the methodology, altough this has further evolved in the code with respect to the paper. An extended version will be made available soon.

## Installation

The package HEDGeOPF relies on the R package [volesti](https://www.rdocumentation.org/packages/volesti/) for polytope sampling via [RCall.jl](https://github.com/JuliaInterop/RCall.jl). Therefore, these external dependencies must be installed for HEDGeOPF to fully work. Currently, this can be done manually through the following steps:

1. Perform a system-wide installation of R from [CRAN](https://cran.r-project.org/). Is is recommended to install a R ≥ 4.4 release. In particular, the latest version of volesti is build with R 4.4.3.
2. Open an R session and run the following command to install the latest release of volesti:

```r
install.packages("volesti", repos = "https://cloud.r-project.org")
```

It is now possible to install HEDGeOPF using the Julia package manager.

```julia
] add https://github.com/mttb91/HEDGeOPF.jl.git
```

As explained [here](https://juliainterop.github.io/RCall.jl/stable/installation/#Customizing-the-R-installation-using-R_HOME), in order for RCall to locate the R installation, run

```julia
] build RCall
```

Running HEDGeOPF requires two solvers: one for Linear Programming (LP) problems and the other for nonlinear (NLP) ones. The open-source solvers HiGHS and Ipopt are recommended for the first and second problem class, respectively, and can be installed via the package manager with

```julia
] add HiGHS, Ipopt
```

Test that the package works by running

```julia
] test HEDGeOPF
```

HEDGeOPF's tests verify, among other things, if RCall connection to R is established correctly and if the R package volesti is installed and fully working.
