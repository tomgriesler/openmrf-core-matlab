## Subfolders

- **src_main/**

Core utilities for sequence compilation, metadata backup and reading pulseq meas data 

- **src_mics/**

Auxiliary utilities and modules:

	- AFP: adiabatic full passage pulses

	- AHP: adiabatic half passage pulses

	- FOV: fov setup

    - GE: functions for TRID labels and receive gain adjustments

	- PNS: peripheral nerve stimulation

	- WASABI: simultaneous mapping of water shift and B1 (DOI: 10.1002/mrm.26133)

	- SIGPY: interface for SLR or adiabatic pulses

	- TRAJ: trajectory calibration tools

- **src_mrf/**

Research toolbox for Magnetic Resonance Fingerprinting:

    - src_sequences: functions for modular generation of MRF sequences

    - src_simulations: functions for the calculation of dictionaries based on .seq files

  

- **src_preparations/**

Preparation modules:

    - ADIASL: Adiabatic Spin-Locking
    - CEST: Chemical Exchange Saturation Transfer
    - FAT: Fat Suppression
    - INV: Inversion
    - MLEV: Malcolm-Levitt T2 or T2rho
    - SAT: Saturation
    - SL: Spin-Locking
    - T2
the spin-lock preparation uses by default the balanced SL module described in: https://doi.org/10.1002/mrm.28585


- **src_readouts/**

Imaging readouts:

    - EPI: Echo Planar Imaging
    - GRE: Gradient Echo
    - PRESS: Point RESolved Spectroscopy
    - RAD: Radial
    - SPI: Spiral
    - SPITSE: Spiral Turbo-Spin-Echo
    - TSE: Turbo-Spin-Echo
    - UTE: Ultra short TE

The following repositories were used for the basic EPI, RAD, SPITSE, TSE and UTE readout modules:

    - https://pulseq.github.io/writeEpi.html
    - https://pulseq.github.io/writeFastRadialGradientEcho.html
    - https://github.com/HennigJue/single-shot-spiral-TSE
    - https://pulseq.github.io/writeTSE.html, https://pulseq.github.io/writeUTE_rs.html

  
_Maximilian Gram: 09.03.2026_