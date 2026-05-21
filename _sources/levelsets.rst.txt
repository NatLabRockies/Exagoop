.. highlight:: rst

.. _levelsets:

Embedded Boundaries via Level Sets (``USE_EB``)
================================================

ExaGOOP can represent solid obstacles inside the computational domain using a
*signed-distance level-set* approach built on top of the AMReX Embedded
Boundary (EB2) library.  Particles that fall inside a body are removed at
initialisation; during time advancement the level-set field is used to enforce
momentum and thermal boundary conditions on grid nodes that lie on or near each
body surface, and to reflect or absorb particles that attempt to cross it.


Build requirements
------------------

Level-set geometry is a compile-time feature controlled by the ``USE_EB`` flag.

**GNUmake**

.. code-block:: makefile

   USE_EB = TRUE    # enables embedded-boundary support
   USE_EB = FALSE   # disables (default)

**CMake**

.. code-block:: bash

   cmake -DEXAGOOP_USE_EB=ON  ..   # enables EB
   cmake -DEXAGOOP_USE_EB=OFF ..   # disables EB (default)

When ``USE_EB=TRUE``, the preprocessor symbol ``USE_EB=1`` is injected and
all ``#if USE_EB`` guarded blocks in the source are compiled in.  The AMReX
EB2 library is linked automatically via the AMReX submodule.


Multi-body implementation
--------------------------

ExaGOOP supports up to ``EXAGOOP_MAX_LS_BODIES = 8`` independent level-set
bodies simultaneously.  Each body is described by a *named block* in the
input file.  The set of active body names is declared first:

.. code-block:: bash

   eb2.body_names = cylinder plate

ExaGOOP reads this list and then processes a ParmParse namespace for every
name in it (e.g. ``cylinder.*``, ``plate.*``).  Each body owns an
independent signed-distance ``MultiFab`` (``lsphi``) built at a refinement
factor relative to the background grid, and an independent set of momentum
and temperature boundary condition parameters.

The ``LevelSetBody`` struct (``Source/mpm_eb.H``) aggregates all per-body
data:

- ``name`` — string identifier matching the entry in ``eb2.body_names``
- ``lsphi`` — pointer to the per-body nodal signed-distance ``MultiFab``
- ``ls_refinement`` — refinement factor for the level-set grid (1 = same
  resolution as the background grid)
- Momentum BC parameters: ``mom_bc_type``, ``wall_mu``, ``wall_vel``
- Temperature BC parameters: ``temp_bc_type``, ``T_wall``, ``heat_flux``,
  ``h_conv``, ``T_inf``

At startup, ``mpm_ebtools::init_eb()`` loops over ``eb2.body_names``, builds
the AMReX EB2 implicit function for each body, fills its ``lsphi`` ``MultiFab``,
and appends the populated ``LevelSetBody`` to ``mpm_ebtools::ls_bodies``.
Subsequent operations (particle removal, nodal BC enforcement, particle
reflection) iterate over this vector.

.. note::

   The AMReX ``EBFArrayBoxFactory`` — which provides cell-type flags used by
   some AMReX routines — is currently built from the *last* body processed.
   For single-body simulations this is exact.  Multi-body support for the
   factory (union of all body shops) is a planned enhancement; the
   per-body ``lsphi`` fields are already fully independent.


Supported geometry types
-------------------------

Each body must specify a ``<name>.geom_type`` entry.  The three analytic
types currently supported are:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - ``geom_type``
     - Description
   * - ``sphere``
     - Sphere (circle in 2-D) defined by a centre point and radius.
   * - ``plane``
     - Half-space defined by a point on the plane and an outward normal.
   * - ``cylinder``
     - Finite cylinder defined by a centre, direction axis, radius, and
       half-length.

A fourth type ``udf`` allows an arbitrary signed-distance function to be
supplied as a runtime-loaded shared library; see :ref:`udf_moving_wall` for
the shared-library mechanism (an analogous ``geom_type = udf`` path exists
in ``Source/mpm_eb_udf_build.cpp``).


Sphere (circle in 2-D)
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   <name>.geom_type                = sphere
   <name>.sphere_radius            = 0.15
   <name>.sphere_center            = 0.5 0.5 0.0   # x y z (z ignored in 2-D)
   <name>.sphere_has_fluid_inside  = false          # false => fluid outside
   <name>.ls_refinement            = 2              # level-set grid refinement


Plane (half-space)
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   <name>.geom_type        = plane
   <name>.plane_point      = 0.5 0.0 0.0   # a point on the plane
   <name>.plane_normal     = 1.0 0.0 0.0   # outward normal (need not be unit)
   <name>.ls_refinement    = 1


Cylinder
~~~~~~~~~

.. code-block:: bash

   <name>.geom_type             = cylinder
   <name>.cylinder_radius       = 0.1
   <name>.cylinder_center       = 0.5 0.5 0.5   # axis centre
   <name>.cylinder_direction    = 2              # axis direction: 0=x, 1=y, 2=z
   <name>.cylinder_half_height  = 0.3
   <name>.ls_refinement         = 2


Momentum level-set boundary conditions
----------------------------------------

Each body's momentum BC type is set with ``<name>.levelset_mom``.  The three
recognised values map to the same integer BC codes used for domain-wall
momentum BCs:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - ``levelset_mom``
     - Behaviour
   * - ``noslipwall``
     - Full no-slip: at contact, the particle velocity is set to the body
       wall velocity (default ``0.0`` for all components; can be overridden
       with ``<name>.wall_vel``).  The same condition is applied to grid
       nodes whose level-set value is negative.
   * - ``slipwall``
     - Slip wall: the normal velocity component is reversed on contact; the
       tangential component is unconstrained.  Grid nodes at the surface have
       only their normal velocity component overwritten.
   * - ``partialslip``
     - Partial-slip Coulomb friction: contact is treated as no-slip but a
       Coulomb slip correction is applied when the tangential stress exceeds
       ``<name>.wall_mu`` times the normal stress.

Optional per-body sub-parameters for momentum BCs:

.. list-table::
   :header-rows: 1
   :widths: 38 12 50

   * - Parameter
     - Type
     - Description
   * - ``<name>.levelset_mom``
     - ``string``
     - Momentum BC type: ``noslipwall``, ``slipwall``, or ``partialslip``.
       Default: ``noslipwall``.
   * - ``<name>.wall_mu``
     - ``Real``
     - Coulomb friction coefficient.  Only used when
       ``levelset_mom = partialslip``.  Default: ``0.0``.
   * - ``<name>.wall_vel``
     - ``Real[SPACEDIM]``
     - Constant wall velocity vector (m/s).  Applied when
       ``levelset_mom = noslipwall``.  Default: ``0.0 0.0 0.0``.

**Example** — stationary no-slip cylinder:

.. code-block:: bash

   cylinder.levelset_mom  = noslipwall
   cylinder.wall_vel      = 0.0 0.0 0.0

**Example** — slip sphere:

.. code-block:: bash

   sphere.levelset_mom    = slipwall


Temperature level-set boundary conditions
------------------------------------------

Thermal BCs on level-set bodies are specified under the same per-body
namespace and require ``USE_TEMP=TRUE`` at build time.  Four BC types are
supported:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - ``temp_bc_type``
     - Description
   * - ``adiabatic``
     - No thermal constraint (natural boundary condition, zero heat flux).
       Grid nodes inside the body are left at whatever temperature the
       particle-to-grid projection deposited.  This is the **default** when
       ``<name>.temp_bc_type`` is not specified.
   * - ``isothermal``
     - Dirichlet (fixed temperature).  Every grid node whose level-set value
       is negative is overwritten with ``<name>.lset_T_wall`` each time step.
   * - ``heatflux``
     - Neumann (prescribed flux).  The surface-node temperature is set by a
       one-sided finite-difference ghost-point formula:

       .. math::

          T_\text{surface} = T_\text{nb}
            + \frac{q_\text{prescribed}}{k_\text{node}} \, \Delta x_n

       where :math:`T_\text{nb}` is the first interior neighbour in the
       dominant surface-normal direction, :math:`k_\text{node}` is the
       nodal thermal conductivity interpolated from particles
       (``MASS_CONDUCTIVITY / MASS_INDEX``), and :math:`\Delta x_n` is
       the grid spacing in that direction.
   * - ``convective``
     - Robin (Newton cooling).  The surface-node temperature is set via a
       Biot-number weighting:

       .. math::

          T_\text{surface}
            = \frac{T_\text{nb} + \text{Bi}\, T_\infty}{1 + \text{Bi}},
          \qquad
          \text{Bi} = \frac{h\,\Delta x_n}{k_\text{node}}

       where :math:`h` is the convective heat-transfer coefficient and
       :math:`T_\infty` is the far-field fluid temperature.

The full set of temperature BC parameters for a level-set body is:

.. list-table::
   :header-rows: 1
   :widths: 38 12 50

   * - Parameter
     - Type
     - Description
   * - ``<name>.temp_bc_type``
     - ``string``
     - Thermal BC type: ``adiabatic``, ``isothermal``, ``heatflux``, or
       ``convective``.  Default: ``adiabatic``.
   * - ``<name>.lset_T_wall``
     - ``Real``
     - Surface temperature :math:`T_w`.  Used when
       ``temp_bc_type = isothermal``.  Default: ``0.0``.
   * - ``<name>.lset_heat_flux``
     - ``Real``
     - Prescribed surface heat flux :math:`q` (energy per unit area per
       unit time).  Used when ``temp_bc_type = heatflux``.  Default: ``0.0``.
   * - ``<name>.lset_h_conv``
     - ``Real``
     - Convective heat-transfer coefficient :math:`h`.  Used when
       ``temp_bc_type = convection``.  Default: ``0.0``.
   * - ``<name>.lset_T_inf``
     - ``Real``
     - Far-field fluid temperature :math:`T_\infty`.  Used when
       ``temp_bc_type = convection``.  Default: ``0.0``.

**Example** — isothermal cylinder at :math:`T = 1`:

.. code-block:: bash

   cylinder.temp_bc_type  = isothermal
   cylinder.lset_T_wall   = 1.0

**Example** — convectively cooled sphere:

.. code-block:: bash

   sphere.temp_bc_type    = convection
   sphere.lset_h_conv     = 5.0
   sphere.lset_T_inf      = 300.0


Implementation notes
---------------------

**Predictor–corrector consistency.**  When the MUSL stress-update scheme is
used, temperature BCs are applied twice per time step: once during the
predictor pass (Dirichlet-only, ``dirichlet_only = true``) and once during
the corrector pass (all BC types).  This matches the treatment of domain-wall
temperature BCs in ``Apply_Nodal_BCs_Temperature`` and ensures that the
isothermal constraint is active at both half-steps.

**Normal direction for flux and convection BCs.**  For ``heatflux`` and
``convective`` types, the interior-neighbour node is found by stepping one
cell in the direction of the largest-magnitude component of the surface
normal (the signed-distance gradient at the node).  This one-dimensional
stencil is consistent with the domain-wall flux and convection BCs and avoids
off-axis ambiguity near corners.

**Level-set grid refinement.**  The signed-distance field is built on a grid
refined by ``ls_refinement`` relative to the background MPM grid.  Higher
refinement gives a more accurate surface representation but increases memory.
A value of 1–2 is sufficient for smooth analytic bodies; sharp-edged bodies
may require 4 or higher.

**Particle removal.**  At initialisation, particles whose level-set value
(evaluated by multilinear interpolation of ``lsphi``) is negative are removed
from the simulation.  The removal loop iterates over all bodies so that
particles inside any body are eliminated.


Complete per-body input block example
--------------------------------------

The following block defines a single circular cylinder (2-D sphere) at domain
centre, with a slip momentum BC and an isothermal thermal BC.  This is the
configuration used in the ``2D_Heat_Conduction_Cylinder_Dirichlet`` tutorial
(see :ref:`tutorial-heat-cylinder`).

.. code-block:: bash

   eb2.body_names                      = cylinder

   cylinder.geom_type                  = sphere
   cylinder.sphere_radius              = 0.15
   cylinder.sphere_center              = 0.5 0.5 0.0
   cylinder.sphere_has_fluid_inside    = false
   cylinder.ls_refinement              = 2
   cylinder.levelset_mom               = slipwall
   cylinder.temp_bc_type               = isothermal
   cylinder.lset_T_wall                = 1.0
