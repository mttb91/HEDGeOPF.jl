# HEDGeOPF

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mttb91.github.io/HEDGeOPF.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mttb91.github.io/HEDGeOPF.jl/dev/)
[![Build Status](https://github.com/mttb91/HEDGeOPF.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/mttb91/HEDGeOPF.jl/actions/workflows/ci.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mttb91/HEDGeOPF.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mttb91/HEDGeOPF.jl)

HEDGeOPF.jl is a Julia package for generating high-quality, synthetic datasets of AC Optimal Power Flow (AC-OPF) instances, designed to support standardization in training and testing of Neural Networks (NNs) that approximate this problem.

It works by sampling the unreduced input load convex space uniformly in terms of total load active power and by generating AC-OPF instances using a modified [PowerModels](https://github.com/lanl-ansi/PowerModels.jl) formulation with load slack variables. The package ensures efficiency and scalability by:

* modifying the AC-OPF model in-place with [JuMP](https://github.com/jump-dev/JuMP.jl),
* leveraging extensively Julia's distributed computing capabilities,
* using R's [volesti](https://github.com/GeomScale/volesti) package via [RCall](https://juliainterop.github.io/RCall.jl/stable/) for uniform sampling in high-dimensional convex polytopes.

## Documentation

The package [documentation](https://mttb91.github.io/HEDGeOPF.jl/dev/) provides useful information, including [installation](https://mttb91.github.io/HEDGeOPF.jl/dev/#Installation) and [quick-start](https://mttb91.github.io/HEDGeOPF.jl/dev/quickstartguide/) guides.

## Development

HEDGeOPF.jl is research-grade software and is constantly being improved and extended. If you have suggestions for improvement, please contact us via the Issues page on the repository.

## Acknowledgements

This work has been partly financed by the Research Fund for the Italian Electrical System under the Three-Year Research Plan 2022–2024 (DM MITE n. 337, 15.09.2022), in compliance with the Decree of April 16th, 2018, and by the EU funds Next-GenerationEU (Piano Nazionale di Ripresa e Resilienza (PNRR) – Missione 4 Componente 2, Linea d’Investimento 3.3 — D.M. 352 09/04/2022.

## Citing HEDGeOPF

If you find HEDGeOPF useful in your work, we kindly request that you cite the following [preprint](https://arxiv.org/abs/2508.19083), which extends and supersedes our earlier [publication](https://ieeexplore.ieee.org/abstract/document/10761586). The new preprint includes a thorough comparison with other open-source AC-OPF dataset generation methods based on newly proposed quality metrics, alongside a detailed explanation of the methodology.

```bibtex
@misc{hedgeopf2025,
  author={Matteo Baù and Luca Perbellini and Samuele Grillo},
  title={A Principled Framework to Evaluate Quality of AC-OPF Datasets for Machine Learning: Benchmarking a Novel, Scalable Generation Method},
  year={2025},
  eprint={2508.19083},
  archivePrefix={arXiv},
  primaryClass={eess.SY},
  url={https://arxiv.org/abs/2508.19083}, 
}
```

```bibtex
@inproceedings{10761586,
  author={Perbellini, Luca and Baù, Matteo and Grillo, Samuele},
  booktitle={2024 IEEE 8th Forum on Research and Technologies for Society and Industry Innovation (RTSI)}, 
  title={An Efficient and Scalable Algorithm for the Creation of Representative Synthetic AC-OPF Datasets}, 
  year={2024},
  pages={653-658},
  doi={10.1109/RTSI61910.2024.10761586}
}
```

## License

This code is provided under a [BSD 3-Clause License](/LICENSE.md).
