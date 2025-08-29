# Quick Start Guide

Once the installation of HEDGeOPF is successfully completed, it is possible to generate an AC-OPF dataset of `N` samples for a given power system test case by following few steps. Start by placing at a path of your choosing:

* the network data file in .m format (e.g., `pglib_opf_case5_pjm.m`)
* the configuration YAML file [`settings.yaml`](https://github.com/mttb91/HEDGeOPF.jl/blob/main/examples/settings.yaml) that controls the simulation.

Open the configuration YAML file and set entries `grid` and `num_samples` under section `CASE` to the grid file name as string (e.g., `"pglib_opf_case5_pjm.m"`) and the desired number of AC-OPF samples `N` (e.g., 10000), respectively. If needed, change entry `cpu_ratio` under section `PARALLEL` to the desired percentage of CPU threads relative to `Sys.CPU_THREADS` that are used for distributed computing.

```yaml
CASE:
    grid : "pglib_opf_case5_pjm.m"
    ...
    num_samples : 10000

...

PARALLEL:
    cpu_ratio : 50.0
```

Save and close the YAML file. Then, in Julia run

```julia
using HEDGeOPF
import HiGHS, Ipopt

path = "/path/of/your/choosing"
generate_dataset(path; filename = "settings.yaml")
```
