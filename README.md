# (Galerkin) Spectral Method for 3D Jeffrey Orbits

This project provides a numerical simulation tool based on the (Galerkin) spectral method using real spherical harmonics to solve the Fokker-Planck equation relating to the Jeffrey orbits in 3D.
Specifically, this simulation tool solves for the orientation distribution of rigid ellispoids in a dilute solution.
It may either be a monodisperse or polydisperse system.
The order parameter and extinction angle are directly evaluated.
Currently two types of flow are supported:
- Simple shear, and
- planar extension.
  
These flow types may be combined but must lie in the same plane.
The flow may be steady or transient.

For more details and references, see the documentation file in `doc/`.
