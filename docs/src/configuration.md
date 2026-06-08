# Configuration Options

The generation of an AC-OPF dataset is completely controlled through an input YAML configuration file [`settings.yaml`](https://github.com/mttb91/HEDGeOPF.jl/blob/main/examples/settings.yaml). As exemplified below, the configuration options of this file are organized into sections and are explained in detail hereafter.

```yaml
CASE:
    name : "example"
    grid : "pglib_opf_case5_pjm.m"
    ...

TOPOLOGY:
    k_min_branch : 1
    k_max_branch : 5
    ...

SAMPLING:
    delta_pd : 100.0
    ...
```

## [`CASE` Options](@id case-options)

This section defines global parameters that control how the AC-OPF dataset is generated, with a few options deserving a more in-depth characterization.

| Key           |    Type    | Description                                                                  |
|---------------|:----------:|------------------------------------------------------------------------------|
| `name`        | `String`   | Name of the folder where the dataset is saved                                |
| `grid`        | `String`   | Full name of the power system grid (MATPOWER format `.m` file)               |
| `append`      | `Bool`     | Whether to append new results to an existing dataset                         |
| `baseseed`    | `Int64`    | Random number generator seed to control reproducibility                      |
| `num_items`   | `Int64`    | Number of load samples generated for a single total load active power level  |
| `num_samples` | `Int64`    | Total number of AC-OPF samples to be generated                               |
| `num_batches` | `Int64`    | Number of batches in which the total samples are processed                   |

### Behaviour of `append`

When `append` is set to `false`, generation fails if the target dataset folder `name` already exists. When, instead, this option is `true`, new AC-OPF instances are appended to existing files. In this case, sampling resumes from the previously saved Random Number Generator (RNG) state, which is stored as `rng_state.bin` in the results folder. See [Extending Dataset](@ref "Extending Dataset") for a practical example.

!!! warning
    If `append` is `true` and `rng_state.bin` is missing, **change the RNG seed** `baseseed` to avoid regenerating the same dataset.

!!! warning
    If extension is run on a different machine, keep the same number of distributed workers used in the previous run. See [`PARALLEL` Options](@ref "parallel-options") to understand how to control the number of parallel workers.

### Behaviour of `num_batches`

Option `num_batches` splits generation into multiple batches to reduce memory usage when `num_samples` is large.

!!! warning
    **Keep batch size above 2000-2500 samples**. After processing the first batch, the sampling region in total load active power is trimmed to exclude areas at the distribution extrema where the AC-OPF never converges. If the batch size is too small to approximate uniform sampling well, regions that are actually feasible for the AC-OPF may be incorrectly removed.

### Behaviour of `num_items`

As detailed in the reference [publication](https://arxiv.org/abs/2508.19083), load sampling is performed uniformly in terms of total load active power by slicing the convex polytope and sampling uniformly from it. The option `num_items` controls how many load samples are generated from a single polytope slice.

!!! note
    Sampling from a polytope slice requires computing first its Chebyshev center. Increasing `num_items` implies computing fewer Chebyshev centers, reducing runtime.

## [`DATASET` Options](@id dataset-options)

The following options control how the NN-ready dataset is processed and split into reproducible cross-validation folds. See [Regenerating Dataset Splits](@ref "Regenerating Dataset Splits") for a practical example.

| Key             | Type                     | Description                                                                                   |
|-----------------|:------------------------:|-----------------------------------------------------------------------------------------------|
| `name`          | `String`                 | Folder where the NN-ready dataset is saved                                                    |
| `cleanup`       | `Bool`                   | Whether to keep only the NN-ready dataset and `settings.yaml`                                 |
| `num_folds`     | `Int64`                  | Number of cross-validation folds for dataset partitioning                                     |
| `num_samples`   | `Int64` or `Nothing`     | Number of samples in the NN-ready dataset (`nothing` means use all)                           |
| `num_quantiles` | `Int64`                  | Number of total load active power quantiles for OPF samples binning                           |

!!! note
    Option `cleanup` should be set to `false` if the user may need to extend the dataset or generate a new NN-ready split from the same data source.

!!! note
    Currently, users have limited control over the splitting strategy. This is selected internally based on whether the dataset includes a single topology or multiple topology perturbations. In the first case, AC-OPF instances are partitioned into folds stratified by total load active power, with samples binned into `num_quantiles` quantiles. In the second case, topology perturbations are distributed randomly across folds.

## [`TOPOLOGY` Options](@id topology-options)

These options control how the topology is perturbed. Currently, only the status of branches (AC lines and transformers) and dispatchable generators can be perturbed independently. The minimum and maximum number of generators and/or branches randomly removed are controlled for each topology. The total number of generated perturbations (as composition, not product, of generator and branch removals) is controlled by `num_topo`. The final number of topologies is `num_topo` + 1, since the original, intact topology is also included.

| Key            |   Type    | Unit | Description                                                               |
|----------------|:---------:|:----:|---------------------------------------------------------------------------|
| `k_min_branch` | `Int64`   | [1]  | Minimum number of branches removed per topology perturbation              |
| `k_max_branch` | `Int64`   | [1]  | Maximum number of branches removed per topology perturbation              |
| `k_min_gen`    | `Int64`   | [1]  | Minimum number of generators removed per topology perturbation            |
| `k_max_gen`    | `Int64`   | [1]  | Maximum number of generators removed per topology perturbation            |
| `num_topo`     | `Int64`   | [1]  | Total number of perturbations to be generated                             |

## [`SAMPLING` Options](@id sampling-options)

These options control how the input sampling space is created. Currently, it consists solely of a convex polytope defined in terms of load active and reactive power variables. Approximately `CASE.num_samples/(num_topo + 1)` load setpoints are generated for each topology perturbation by sampling uniformly in total load active power.

| Key            |   Type    | Unit | Description                                                               |
|----------------|:---------:|:----:|---------------------------------------------------------------------------|
| `delta_pd`     | `Float64` | [%]  | Percentage variation in load active power around the nominal values       |
| `delta_qd`     | `Float64` | [%]  | Percentage variation in load reactive power around the nominal values     |
| `delta_pf`     | `Float64` | [1]  | Maximum reduction in power factor w.r.t. the nominal absolute value       |
| `max_pf`       | `Float64` | [1]  | Maximum allowable load power factor                                       |
| `min_pf`       | `Float64` | [1]  | Minimum allowable load power factor                                       |

!!! note
    The value `delta_pd` cannot exceed 100% since the active power of a load cannot be unrestricted in sign to preserve the sign relation between active and reactive power. To create a load with negative (positive) active power at a bus that already has a positive (negative) one, add a new load to the power system dictionary or .m file with negative (positive) nominal active power.

## [`PARALLEL` Options](@id parallel-options)

The `PARALLEL` section controls distributed computing settings to be applied throughout the simulation.

| Key          | Type      | Unit | Description                                                         |
|--------------|:---------:|:----:|---------------------------------------------------------------------|
| `cpu_ratio`  | `Float64` | [%]  | Percentage of CPU threads to use w.r.t. `Sys.CPU_THREADS` count     |

## [`MODEL` Options](@id model-options)

These options are related to the modified PowerModels OPF model with slack variables for active and reactive load power.

| Key          | Type      | Unit     | Description                                                    |
|--------------|:---------:|:--------:|----------------------------------------------------------------|
| `duals`      | `Bool`    | –        | Whether to record dual values for every primal AC-OPF variable |
| `voll`       | `Float64` | [€/MWh]  | Value Of Lost Load                                             |

When `duals` is set to `true`, by design HEDGeOPF automatically retrieves the dual values, if available, of every `JuMP.VariableRef` defined in the model and, yet, does not look for any `JuMP.ConstraintRef` object. This choice stems from the fact that PowerModels mainly employs anonymous, non-containerized JuMP constraints in model definition, making it difficult to retrieve their (dual) values when inspecting results.

!!! warning
    To record the dual value of branch apparent power, the OPF model is modified when `duals` is set to `true` by adding variables for the **square of the branch apparent power** at the from and to buses. When accessing and using the primal values for branch apparent power in the dataset results, the user should remember that they are squared.

## [`PATH` Options](@id path-options)

These define input and output file paths relative to the input configuration YAML file (`basepath`). Overall, the absolute path of the grid file is `basepath/PATH.input/`. Similarly, the absolute path of the dataset is composed as `basepath/PATH.output/CASE.grid/CASE.name/`, with the NN-ready dataset being located at `basepath/PATH.output/CASE.grid/CASE.name/DATASET.name`.

| Key           | Type      | Description                                                       |
|---------------|:---------:|-------------------------------------------------------------------|
| `input`       | `String`  | Relative path to the folder containing the grid file              |
| `output`      | `String`  | Relative path to the folder where the AC-OPF dataset is saved     |

## [`SOLVER` Options](@id solver-options)

The `SOLVER` section specifies the LP and NLP solvers employed in HEDGeOPF.

| Key           | Type      | Description                                                       |
|---------------|:---------:|-------------------------------------------------------------------|
| `lp`          | `String`  | Julia package name of linear programming solver                   |
| `nlp`         | `String`  | Julia package name of nonlinear programming solver                |
| `lp_options`  | `Pair`    | Key-value pairs of options for the linear programming solver      |
| `nlp_options` | `Pair`    | Key-value pairs of options for the nonlinear programming solver   |

The following example shows how `lp_options` and `nlp_options` can be specified in the YAML file.

```yaml
SOLVER:
    lp          : "HiGHS"
    lp_options  :
        solver            : "ipm"  
    nlp         : "Ipopt"
    nlp_options :
        max_cpu_time      : 1000.0
        mumps_mem_percent : 10
        print_level       : 0
```

HEDGeOPF is designed to be independent of the specific LP and NLP solver choice (as long as the selected solvers support the variables and constraints used in the optimization models). Users can control this by installing the relevant Julia packages and updating the `SOLVER` section accordingly. Solver options, such as those of [HiGHS](https://ergo-code.github.io/HiGHS/dev/options/definitions/) and [Ipopt](https://coin-or.github.io/Ipopt/OPTIONS.html), are typically available in the solver documentation.

!!! note
    Currently only HiGHS and Ipopt solvers have been tested with HEDGeOPF, with Ipopt being equipped with the default sequential MUMPS linear solver.
