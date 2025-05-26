
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
