
# MSPD: Multiscale Stochastic Particle Method for Diatomic Gas Flows

This repository extends the SPARTACUS framework with a **Multiscale Stochastic Particle method based on the Fokker–Planck equation (MSP)**, specifically designed to simulate nonequilibrium flows of **diatomic gases** involving internal energy exchange.

## Highlights
- MSPD incorporates **Langevin-based integration schemes** to model internal mode relaxation (rotational and vibrational).
- The method introduces a **modified collision operator** to achieve **second-order temporal accuracy** even in near-continuum regimes.
- It is implemented on top of SPARTACUS, maintaining full compatibility with SPARTA infrastructure.

## Implementation Notes
The MSPD method is implemented by extending the original SPARTACUS framework with the following components:

- **Computation of required macroscopic quantities on each cell**: At every time step, cell-averaged quantities are computed and used to determine the coefficients needed for the Langevin update.

- **Langevin-based particle updates**: Each particle's velocity and internal energies (rotational and vibrational) are updated according to a Langevin-based scheme driven by the cell-local parameters.

- The main modifications are located in `collide_bgk.cpp`, where the MSPD-specific updates are implemented. All usage and configuration are consistent with the original SPARTACUS workflow via the `style` keyword in the input script.

- Details of how to configure the MSPD model can be found in the input script comments for each benchmark case.

## MSPD Input Configuration

The MSPD model is integrated as a **style within the `collide bgk` framework**. To activate the MSPD-based collision model, use the following syntax in the input script:

```text
collide              bgk air mspd air.bgk
```

This enables the MSPD collision operator for the `air` species.

Additional configuration options specific to MSPD can be set via the `collide_bgk_modify` command:

```text
collide_bgk_modify   time_ave 0.99 reset_wmax 0.98 relax_mod 0 0 vib_energy smooth
```

### Options explained:

- `relax_mod <rot_mode> <vib_mode>`  
  Controls how the **rotational and vibrational relaxation numbers** are specified:
  - `0`: constant relaxation number;
  - `1`: temperature-dependent relaxation number.

- `vib_energy <mode>`  
  Specifies the treatment of **vibrational energy**:
  - `none`: vibration not activated (e.g., low-temperature flows);
  - `smooth`: continuous vibrational energy (harmonic oscillator model);
  - `discrete`: discrete vibrational energy levels (quantum energy model).

These options allow flexible modeling of internal energy relaxation in diatomic gases. Examples of full configurations are provided in each benchmark case in the `examples_MSPD/` directory.

## Reproducibility
All benchmark cases reported in our paper are included in the repository:
- Homogeneous relaxation of internal modes
- Normal shock structure in diatomic gases
- Hypersonic flow over a cylinder
- 70-degree blunted cone test case

Please refer to the `examples_MSPD/` directory and follow the detailed usage instructions in each case subfolder.

----------------------------------------------------------------------
----------------------------------------------------------------------
*The following sections contain the original README files from **SPARTACUS** and **SPARTA** for reference.*

**SPARTACUS** is a Unified Stochastic Particle (*USP*) solver based on SPARTA.
**SPARTACUS** stands for **SPARTA** **C**ombined with **US**P.
The *USP* method is proposed by *Fei et al.* to simulate multi-scale gas flows.

> [examples_USP](./examples_USP/) contains Several benchmark cases and tutorial examples.

> [doc/collide.html](./doc/collide.html) & [doc/collide_bgk_modify.html](./doc/collide_bgk_modify.html)  contains essential documentation about this solver.

Downloading, installing and using SPARTACUS is exactly the same as SPARTA, which is fully described within [SPARTA Users Manual](https://sparta.github.io/doc/Manual.html), unless specified in the above files and directories

The current maintainer of this solver is [**Feng Kai**](https://scholar.google.com/citations?user=dk0a3ysAAAAJ) from [**Jun Zhang**](https://scholar.google.com/citations?user=6vjJtPsAAAAJ)'s Lab at Beihang University. For issue reporting or suggestions, please contact kfeng@buaa.edu.cn or jun.zhang@buaa.edu.cn.

*Coupled time step and grid adaptation based on local characteristic scale is available now!*

<img src=./doc/JPG/adaptation.gif width=50% />

**SPARTACUS Papers**

The following papar introduces the implementation of the USP method in the current version of SPARTACUS.

> Feng, K., Tian, P., **Zhang, J.***, Fei, F., Wen, D., 2022. [**SPARTACUS: An open-source unified stochastic particle solver for the simulation of multiscale nonequilibrium gas flows.**](https://doi.org/10.1016/j.cpc.2022.108607) Computer Physics Communications 284 (2023) 108607.

> Feng, K., Cui, Z., Tian, P., **Zhang, J.*** [**A Unified Stochastic Particle Method with Spatiotemporal Adaptation for Simulating Multiscale Gas Flows.**](https://doi.org/10.1016/j.jcp.2024.112915) , Journal of Computational Physics 505 (2024) 112915.

The following paper covers the theoretical foundation of SPARTACUS, namely the development of USP method.

> Fei, F., **Zhang, J.***, Li, J., **Liu, Z.***, 2020. [**A unified stochastic particle Bhatnagar-Gross-Krook method for multiscale gas flows.**](https://doi.org/10.1016/j.jcp.2019.108972) Journal of Computational Physics 400, 108972. 

> **Fei, F.***, Ma, Y., **Wu, J.***, Zhang, J., 2021. [**An efficient algorithm of the unified stochastic particle Bhatnagar-Gross-Krook method for the simulation of multi-scale gas flows.**](https://doi.org/10.1186/s42774-021-00069-8) Adv. Aerodyn. 3, 18. 

> **Fei, F.***, **Hu, Y.***, Jenny, P., 2022. [**A unified stochastic particle method based on the Bhatnagar-Gross-Krook model for polyatomic gases and its combination with DSMC.**](https://doi.org/10.1016/j.jcp.2022.111640) Journal of Computational Physics 471, 111640. 

The following are some partial papers that discuss the application of SPARTACUS or the USP method

> Tian, P., Feng, K., Ma, Q., Li, Z., **Zhang, J.***, 2023.[ **Unified stochastic particle simulation of polyatomic gas flows using SPARTACUS**](https://doi.org/10.1016/j.compfluid.2023.105987). Computers & Fluids 265, 105987. 

> **Zhang, J.***, Yao, S., **Fei, F.***, Ghalambaz, M., Wen, D., 2020. [**Competition of natural convection and thermal creep in a square enclosure.**](https://doi.org/10.1063/5.0022260) Physics of Fluids 32, 102001. 

> Ma, Q., Yang, C., Chen, S., Feng, K., Cui, Z., **Zhang, J.*** (2023).[**Effect of thermal fluctuations on spectra and predictability in compressible decaying isotropic turbulence.**](
> https://doi.org/10.1017/jfm.2024.342)  Journal of Fluid Mechanics 987 (2024) A29.

[junzhangmail]: mailto:jun.zhang@buaa.edu.cn 'junzhangmail'
[feifeimail]: mailto:ffei@hust.edu.cn 'feifeimail'

----------------------------------------------------------------------
----------------------------------------------------------------------

SPARTA on GitHub: https://github.com/sparta/sparta 

***Original SPARTA README:***



This is the SPARTA software package.

SPARTA stands for Stochastic PArallel Rarefied-gas Time-accurate
Analyzer.

Copyright (2014) Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software.  This software is distributed under
the GNU General Public License.

----------------------------------------------------------------------

SPARTA is a Direct Simulation Monte Carlo (DSMC) code designed to run
efficiently on parallel computers.  It was developed at Sandia
National Laboratories, a US Department of Energy facility, with
funding from the DOE.  It is an open-source code, distributed freely
under the terms of the GNU Public License (GPL).

The primary authors of the code are Steve Plimpton and Michael Gallis,
who can be emailed at sjplimp@sandia.gov and magalli@sandia.gov.  The
SPARTA web site at http://sparta.sandia.gov has more information about
the code and its uses.

The SPARTA distribution includes the following files and directories:

README			   this file
LICENSE			   the GNU General Public License (GPL)
bench                      benchmark problems
data                       files with species/reaction params, surface files
doc                        documentation
examples                   simple test problems
lib                        additional library files
python                     Python wrapper on SPARTA as a library
src                        source files
tools                      pre- and post-processing tools

Point your browser at any of these files to get started:

doc/Manual.html	           the SPARTA manual
doc/Section_intro.html	   hi-level introduction to SPARTA
doc/Section_start.html	   how to build and use SPARTA
doc/Developer.pdf          SPARTA developer guide
