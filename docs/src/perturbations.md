# Perturbation Modes

This page describes, at a semantic level, which power system input data can be perturbed when generating an AC-OPF dataset with HEDGeOPF and how different perturbation modes interact with each other.

Currently, the following perturbation types are supported:

* Faults at AC lines and transformers
* Faults at dispatchable generators
* Variation of load active and reactive power setpoints

Faulted components are not removed from the model or from the AC-OPF results. Instead, they remain represented explicitly, with optimal generator setpoints and/or branch flows being constrained to zero by model limits/status settings. Faulted components are also signalled at the level of static power system data (see the [Dataset Format](@ref "Graph Level")). This, however, is intentionally handled differently for branches and generators to reflect different representation needs in downstream NN models for the AC-OPF task.

## Topology Perturbations

Topology can be perturbed by deactivating randomly selected AC lines and/or transformers. The number of branches placed out of service varies across topology perturbations and is drawn at random between lower and upper bounds specified in the [`TOPOLOGY` Options](@ref "topology-options") of the YAML configuration file.

!!! note
    Only topology perturbations that result in a single connected component containing all buses are considered valid and are processed in the AC-OPF dataset generation pipeline.

!!! note
    Faulted branches are signalled at the level of static power system data by setting `br_status` to zero in the `branch.csv` file.

## Generator Perturbations

Dispatchable generators can also be deactivated at random. At each perturbation, the number of faulted generators is randomly selected between lower and upper bounds specified in the [`TOPOLOGY` Options](@ref "topology-options") of the YAML configuration file.

!!! note
    Only generator perturbations that leave at least one bus connected to a generator with nonzero active and reactive power capability are considered valid, so that a feasible slack/reference bus can be assigned.

!!! note
    Faulted generators are signalled at the level of static data by setting generator active and reactive power limits (i.e., entries `pmin`, `pmax`, `qmin` and `qmax`) to zero in the `gen.csv` file. In addition, the type of the connection bus (entry `bus_type` in the `bus.csv` file) is switched from PV (or slack) to PQ if the bus does not host any other generator with nonzero capability. If needed, a new slack bus is then selected among eligible candidates.

## Load Perturbations

As detailed in the reference [publication](https://arxiv.org/abs/2508.19083), load setpoints in active and reactive power are sampled from a convex polytope. This is bounded by:

* per-load lower and upper limits on active and reactive power;
* per-load two-sided constraints on the coupling between active and reactive power;
* minimum and maximum limits on total load active power, which depend on the total installed active generator capacity;

where the per-load limits can be controlled (in a basic, limited way) through the [`SAMPLING` Options](@ref "sampling-options") of the YAML configuration file.

Generating a single load perturbation requires the following three sequential steps:

1. Draw a sample from a uniform distribution in total load active power.
2. Slice the convex polytope by intersecting it with a hyperplane defined by the sampled total active load.
3. Sample the polytope slice uniformly at random to generate active and reactive power load setpoints.

## Workflow

Altogether, the generation of perturbations is organized into two nested loops based on the following heuristic logic:

* The outer level handles generation of *structural* perturbations. These are secondary (i.e., less frequent) sources of variability in the power system that may alter the convex hyperspace for load sampling.
* At the inner level, multiple *operational* perturbations are sampled for each structural one. Operational perturbations are primary sources of variation in power system operations (e.g., changes in load setpoints, generator costs, etc.).

Within this framework, in HEDGeOPF input samples to the AC-OPF problem are obtained by:

1. Generating a structural perturbation that combines branch and generator outages (assuming that both perturbation modes are enabled).
2. Updating the load-sampling feasible region, if needed, to reflect any change in the installed generation capability induced by the structural perturbation.
3. Sampling the resulting polytope uniformly in total load active power, generating multiple load perturbations for the same structural one.

This workflow is repeated for multiple structural perturbations.
