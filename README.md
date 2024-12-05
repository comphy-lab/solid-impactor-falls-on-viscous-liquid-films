# Elastic Impactor Hitting a Viscous Liquid Film

This repository contains code for simulating the dynamics of an elastic spherical impactor falling onto a viscous liquid film. The simulation is implemented in 2D+axi (axisymmetric) configuration using the Basilisk framework.

## Background

The simulation models a viscoelastic scalar implementation of a spherical impactor falling onto a viscous liquid film. The physics includes:
- Two-phase flow with surface tension
- Viscoelastic effects
- Axisymmetric geometry
- Fluid-structure interaction

## Implementation

### Key Files
- `testCases/sphereFilm_catch.c`: Main simulation file implementing the physics
- `testCases/VideoAxi.py`: Post-processing script for visualization
- `src-local/log-conform-viscoelastic-scalar-2D.h`: Viscoelastic solver implementation
- `src-local/two-phase-TF-VE.h`: Two-phase flow implementation

### Key Parameters
- `We`: Weber number (ratio of inertial to surface tension forces)
- `Ohd`: Ohnesorge number for the droplet (ratio of viscous to inertial and surface tension forces)
- `Ec`: Elastic parameter
- `De`: Deborah number (ratio of relaxation time to observation time)
- `Ohf`: Ohnesorge number for the film
- `hf`: Film thickness
- `Bo`: Bond number (ratio of gravitational to surface tension forces)

### Numerical Methods
- Grid: Adaptive mesh refinement with levels from `MINlevel` (4) to `MAXlevel` (10)
- Time integration: Centered Navier-Stokes solver
- Interface tracking: Volume-of-Fluid method with tension
- Error tolerances:
  - Fraction error: 1e-3
  - Velocity error: 1e-3
  - Kinetic energy error: 1e-3
  - Area error: 1e-3

## Usage

### Compilation
The code requires the Basilisk framework. Compile using:
```bash
qcc -Wall -O2 -o sphereFilm_catch sphereFilm_catch.c -lm
```

### Running Simulations
Execute the compiled binary with appropriate parameters:
```bash
./sphereFilm_catch [parameters]
```

### Post-processing
Use the `VideoAxi.py` script for visualization:
```bash
python VideoAxi.py [options]
```

## Post-processing Scripts

### VideoAxi.py
A Python script for creating visualizations of the simulation results. Features:
- Custom color mapping for different phases
- Facet extraction for interface visualization
- Field data extraction and plotting
- Multi-processing capabilities for faster rendering

## Dependencies
- Basilisk framework
- Python libraries:
  - NumPy
  - Matplotlib
  - Subprocess
  - Multiprocessing
  - Argparse

## Author
Vatsal Sanjay (vatsalsanjay@gmail.com)
Physics of Fluids

## Version History
- Version 1.0 (December 5, 2024): Initial release
