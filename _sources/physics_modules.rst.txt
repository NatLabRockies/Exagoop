.. _physics_modules:

Physics Modules
===============

ExaGOOP's core solver handles momentum and mass conservation for isothermal,
single-phase problems. Additional physics are compiled in through optional
modules, each controlled by a build-time flag. This page gives an overview of
the three available modules and guidance on when to use them. Full details —
including input file parameters, boundary condition syntax, and worked examples
— are in the individual module pages linked below.

.. list-table:: Module summary
   :header-rows: 1
   :widths: 30 20 20 30

   * - Module
     - CMake flag
     - GNUmake flag
     - Purpose
   * - :ref:`temperature_module`
     - ``EXAGOOP_USE_TEMP=ON``
     - ``USE_TEMP=TRUE``
     - Coupled heat transfer and thermal BCs
   * - :ref:`levelsets`
     - *(always ON for DIM ≥ 2)*
     - ``USE_EB=TRUE``
     - Solid obstacles via signed-distance level sets
   * - :ref:`udf_moving_wall`
     - *(runtime, no extra flag)*
     - *(runtime, no extra flag)*
     - Spatially and temporally varying velocity BCs


Temperature Module
------------------

The temperature module adds a Lagrangian energy equation solved alongside the
momentum equation. When enabled, each material point carries thermal state
(temperature, specific heat, thermal conductivity, heat flux, and internal heat
source), and the background grid carries five additional nodal fields that are
deposited and interpolated every time step.

**When to use it.** Enable this module for any problem where thermal effects
are significant: conjugate heat transfer, thermally driven flows, or
constitutive models whose parameters depend on temperature.

**Build flag.** ``EXAGOOP_USE_TEMP=ON`` (CMake) or ``USE_TEMP=TRUE`` (GNUmake).
The module is compiled *out* by default; simulations that do not need it incur
no runtime overhead.

:doc:`Full temperature module documentation <temperature_module>`


Embedded Boundaries / Level Sets
---------------------------------

The embedded boundary (EB) module allows solid obstacles to be placed inside
the computational domain using signed-distance level sets built on the AMReX
EB2 library. Particles that fall inside a body are removed at initialisation;
during time advancement the level-set field enforces momentum and thermal
boundary conditions on grid nodes near each body surface, and reflects or
absorbs particles that attempt to cross it. Multiple bodies can be active
simultaneously.

**When to use it.** Enable for any problem with a non-trivial internal
geometry — cylinders, spheres, arbitrary bodies — that cannot be handled by
simple domain-face boundary conditions alone.

**Build flag.** ``USE_EB=TRUE`` (GNUmake); EB is always compiled in by the
CMake build for ``DIM ≥ 2``. EB is not available in 1-D builds.

:doc:`Full embedded boundary documentation <levelsets>`


UDF Moving Wall Boundary Conditions
-------------------------------------

The UDF (user-defined function) interface allows arbitrary velocity boundary
conditions to be specified at any domain face at runtime, without recompiling
the solver. The user compiles a small shared library exporting a single
C-linkage function; ExaGOOP loads it at start-up and calls it every time step
to obtain the velocity at each boundary node.

**When to use it.** Use this for wall motions that cannot be expressed as a
simple constant: rotating walls, oscillating pistons, prescribed shear profiles,
or any spatially and temporally varying velocity field.

**Build flag.** No additional compile-time flag is required. The UDF is
activated at runtime through two extra parameters in the input file pointing
to the shared library path and function name.

:doc:`Full UDF moving wall documentation <udf_moving_wall>`


Combining Modules
-----------------

The modules are independent and can be combined freely. For example, a
simulation of heat transfer around a moving cylinder would require both
``USE_TEMP`` (for the energy equation) and ``USE_EB`` (for the cylinder
geometry). The tutorials section includes examples that exercise these
combinations.
