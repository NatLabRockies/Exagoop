Introduction
============

About ExaGOOP
-------------

ExaGOOP is a Material Point Method (MPM) based solver developed by the Scalable Algorithms,
Modeling and Simulation (SAMS) team at the National Laboratory of the Rockies,
Colorado. The solver can be used to simulate physical problems involving gaseous, liquid, and
solid phases, including scenarios with large deformations, history-dependent material
behaviour, and multi-phase interactions. ExaGOOP is built on the
`AMReX <https://amrex-codes.github.io/amrex/>`_ framework and has demonstrated excellent
scalability on both CPU and GPU architectures.

Background and Motivation
--------------------------

Numerical simulations have become an essential tool in engineering analysis today. They
provide an effective alternative to expensive experiments and are frequently used to
selectively identify optimal designs for experimental testing. The traditional numerical
methods utilized in continuum mechanics include Finite Volume Methods (FVMs), Finite Element
Methods (FEMs), and Finite Difference Methods (FDMs). The first two methods rely on the
integral form of the governing equations, while the last method is based on the differential
form, approximating derivatives using algebraic expressions. A shared characteristic of
these three methods is the requirement of a computational grid to solve the governing
equations. This grid consists of a collection of numbered points (or nodes) interconnected
by edges (in 2D) or faces (in 3D) to form grid cells. These grid cells together encompass
the entire computational domain. In solid mechanics, these methods are usually applied in
conjunction with the Lagrangian form of the governing equations. In this context, the grid
nodes are attached to the material and move as the material deforms. However, when the
material experiences severe deformation, the movement of the mesh nodes can lead to grid
entanglement, resulting in numerical instabilities that may cause the solver to halt
abruptly. In contrast, Eulerian-based methods do not require grid movement; instead, the
material deforms across fixed grid cells, decoupling it from the grid nodes. While this
formulation is numerically stable, it is less suited for problems involving moving
interfaces and materials with history-dependent properties.

MPM is a mesh-free technique that has garnered significant
attention for its ability to simulate problems with considerable deformations. MPM is a
Lagrangian, particle-based numerical method inspired by the Particle-In-Cell (PIC)
:cite:`osti_4769185` and Fluid-Implicit Particle (FLIP) methods :cite:`BRACKBILL1986314`.
Although MPM was initially developed for solid mechanics problems, it has since been
extended for use in fluid simulations as well. In MPM, the entire material domain is
discretized using particles or material points, where all material properties — such as
velocity, density, strain rates, and stresses — are stored. This aspect of MPM makes it
especially useful for modeling history-dependent constitutive models. The lack of a grid
connecting these material points allows MPM to effectively simulate scenarios with
significant material deformations, such as solid fractures, foam deformations, and granular
flows. Although a background grid is necessary for MPM, it is primarily used as a scratch
pad for gradient calculations and time integration. These features make MPM a promising
approach for studying materials with large deformations and history-dependent properties.

Key Features
------------

Key features of the ExaGOOP MPM solver include:

- Single-level, Cartesian background grid generation and manipulation using the AMReX library
- Material point functionalities implemented using the AMReX particle library
- CPU, GPU, and CPU+GPU implementation on multiple heterogeneous architectures
- Explicit time integration scheme
- Bilinear, Quadratic B-Spline, and Cubic B-Spline shape functions
- Particle in Cell (PIC) and Fluid Implicit Particle (FLIP) methods
- Update Stress Last (USL) and Modified Update Stress Last (MUSL) stress update schemes
- Multiple constitutive models for solids and fluids
- Optional coupled heat-transfer module for thermal problems
- Embedded boundary support via signed-distance level sets
- User-defined function (UDF) interface for arbitrary moving wall boundary conditions
