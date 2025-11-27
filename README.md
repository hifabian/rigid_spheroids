# A (Galerkin) Spectral Method for Rigid Spheroids under Flow

This project provides a numerical simulation tool based on the (Galerkin) spectral method using real spherical harmonics to solve the Fokker-Planck equation describing rigid spheroids in 3D that follow Brownian motion and are under influence of a flow field.
Specifically, this simulation tool solves for the orientation distribution of rigid spheroids in a dilute solution, currently assuming the rotational diffusion coefficient for rods.
The system may either be a monodisperse or polydisperse system.
The order parameter and extinction angle are directly evaluated.
Currently, all types of planar flows are supported, which include:
- Simple shear,
- planar extensional flow,
- rotational flow, and
- any combination of these in the same plane.

The plane can be either set to `xy` or `xz` and the flow may be steady or unsteady, i.e. transient, both in terms of the flow field and the spheroid's motion.

For more details on the mathematics and references, see the documentation file in `doc/`.
