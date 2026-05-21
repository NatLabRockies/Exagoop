.. highlight:: rst

.. _tutorial-heat-cylinder:

Tutorial: 2-D Steady Heat Conduction around a Cylinder (Dirichlet)
====================================================================

This tutorial walks through the ``Tests/2D_Heat_Conduction_Cylinder_Dirichlet``
test case.  It demonstrates how to use an embedded-boundary level-set body
with an isothermal (Dirichlet) thermal BC to compute steady-state heat
conduction in a 2-D domain containing a circular cylinder — the canonical
Laplace-equation problem with an analytic solution.

Required build flags: ``USE_EB=TRUE``, ``USE_TEMP=TRUE``, ``DIM=2``.


Physical problem
-----------------

A unit square domain :math:`[0,1]^2` contains a circular cylinder of radius
:math:`r = 0.15` centred at :math:`(0.5,\,0.5)`.  The four domain walls are
held at :math:`T = 0`, and the cylinder surface is isothermal at :math:`T = 1`.
There are no volumetric heat sources, and the initial elastic solid velocity is zero throughout.

The governing equation is the steady-state heat equation (Laplace equation):

.. math::

   \nabla^2 T = 0

with boundary conditions:

.. math::

   T &= 0 \quad \text{on all four domain walls} \\
   T &= 1 \quad \text{on the cylinder surface } (r = 0.15)

In the time-dependent MPM formulation, the simulation is advanced until
the transient decays and the solution converges to the steady state.  The
thermal diffusivity is :math:`\alpha = k / (\rho\,c_p) = 1`, so the
characteristic diffusion time across the domain is :math:`t_\alpha \sim 1`.
The simulation is run to :math:`t = 0.1`, by which point the solution has
relaxed to within plotting accuracy of the Laplace steady state.


Analytical solution
--------------------

The steady-state temperature in polar coordinates centred on the cylinder
axis is:

.. math::

   T(r, \theta) =
     \frac{\ln(R_\text{out}/r)}{\ln(R_\text{out}/r_0)}

where :math:`r_0 = 0.15` is the cylinder radius and :math:`R_\text{out}`
is an effective outer radius determined by the Dirichlet condition
:math:`T = 0` on the square boundary.  For a circular outer boundary of
radius :math:`R_\text{out}`, this reduces to:

.. math::

   T(r) = \frac{\ln(R_\text{out}/r)}{\ln(R_\text{out}/r_0)}

Because the outer boundary is square rather than circular, the exact solution
is a 2-D series; however, at large separation the logarithmic profile provides
an excellent approximation.  The post-processing script in
``PostProcess/validate.py`` computes the numerical angular average at each
radial band and compares it against a reference profile obtained by
high-resolution numerical integration of the Laplace equation on the same
geometry.


Material properties
--------------------

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - Property
     - Value
     - Notes
   * - Density :math:`\rho`
     - 1.0
     - uniform
   * - Specific heat :math:`c_p`
     - 1.0
     - uniform; gives :math:`\alpha = 1`
   * - Thermal conductivity :math:`k`
     - 1.0
     - uniform
   * - Volumetric heat source :math:`\dot{q}`
     - 0.0
     - no internal source
   * - Initial temperature :math:`T_0`
     - 0.0
     - uniform; below the cylinder BC


Domain and discretisation
--------------------------

The computational domain is the unit square :math:`[0,1]^2` discretised
with a :math:`128 \times 128` background grid (``mpm.ncells = 128 128 0``).
One material point per cell is placed throughout the domain excluding the
cylinder interior, giving approximately 15 000 active particles.  The time
step is fixed at :math:`\Delta t = 10^{-5}` and the simulation runs to
:math:`t_\text{final} = 0.1`.

The level-set grid is refined by a factor of 2 (``cylinder.ls_refinement = 2``),
giving an effective resolution of :math:`256 \times 256` for the signed-distance
field, which accurately represents the circular geometry.


Boundary conditions
--------------------

**Domain walls** — all four faces use Dirichlet temperature BCs at :math:`T = 0`
and slip momentum BCs (no momentum exchange with the walls):

.. code-block:: bash

   mpm.bc_xlo_mom = slip
   mpm.bc_xhi_mom = slip
   mpm.bc_ylo_mom = slip
   mpm.bc_yhi_mom = slip

   mpm.bc_xlo_temp          = dirichlet
   mpm.bc_xlo_temp.T_wall   = 0.0
   mpm.bc_xhi_temp          = dirichlet
   mpm.bc_xhi_temp.T_wall   = 0.0
   mpm.bc_ylo_temp          = dirichlet
   mpm.bc_ylo_temp.T_wall   = 0.0
   mpm.bc_yhi_temp          = dirichlet
   mpm.bc_yhi_temp.T_wall   = 0.0

**Cylinder surface** — isothermal at :math:`T = 1` with a slip momentum BC
(no drag on the particles from the cylinder):

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

For a full description of the level-set BC parameters see :ref:`levelsets`.


Setting up the test case
-------------------------

The particle file and input file are generated from the JSON configuration
in ``PreProcess/config.json``.  The key fields are the grid resolution, the
particles-per-cell count, and the thermal properties:

.. code-block:: json

   {
     "grid": { "nx": 128, "ny": 128 },
     "ppc": [1, 1],
     "bodies": [{
       "temperature": {
         "T": 0.0,
         "spheat": 1.0,
         "thermcond": 1.0,
         "heatsrc": 0.0
       }
     }]
   }

To (re)generate the particle and input files, run from the test directory:

.. code-block:: bash

   cd Tests/2D_Heat_Conduction_Cylinder_Dirichlet/   
   bash Generate_MPs_Inputfile_Generic.sh

This writes ``mpm_particles.dat`` (particle positions and thermal properties)
and ``Inputs_2DHeat_Conduction_Cylinder_Dirichlet.inp`` one level up.


Building
---------

Navigate to the test directory and build with the required flags:

.. code-block:: bash

   cd Tests/2D_Heat_Conduction_Cylinder_Dirichlet
   # create a temporary copy of the GNUmakefile here
   make USE_EB=TRUE USE_TEMP=TRUE DIM=2 -j4

The resulting executable is named
``ExaGOOP2d.<suffix>.ex`` (suffix depends on the build environment).

.. note::

   ``USE_EB=TRUE`` pulls in the AMReX EB2 library.  On some systems this
   significantly increases compile time; using ``-j4`` or higher
   parallelises the build across compilation units.


Running
--------

From the test directory:

.. code-block:: bash

   ./ExaGOOP2d.<suffix>.ex Inputs_2DHeat_Conduction_Cylinder_Dirichlet.inp

For a parallel run on four MPI ranks:

.. code-block:: bash

   mpirun -n 4 ./ExaGOOP2d.<suffix>.ex \
       Inputs_2DHeat_Conduction_Cylinder_Dirichlet.inp

ASCII particle snapshots are written at intervals of
``mpm.write_output_time = 0.01`` to the subdirectory
``2D_Heat_Conduction_Cylinder_Dirichlet/<solution folder name set in input file>/matpnt*``.  A checkpoint is written
to ``2D_Heat_Conduction_Cylinder_Dirichlet/solution folder name set in input file/chk*`` at the same frequency.


Post-processing and validation
--------------------------------

The validation script reads the final particle snapshot, bins particles by
radial distance from the cylinder axis, and compares the angularly-averaged
temperature profile against the analytical steady-state:

.. code-block:: bash

   cd Tests/2D_Heat_Conduction_Cylinder_Dirichlet/
   python3 PostProcess/validate.py --time 0.1
   

A passing run prints the RMS error and ``PASS``.  The acceptance criterion
is RMS :math:`< 5 \times 10^{-2}` (the finite-time run and square-vs-circle
outer boundary geometry contribute a background error of order
:math:`10^{-2}`).

The companion script ``Plot_Temperature.py`` generates a colour map of the
particle temperature field together with a contour overlay:

.. code-block:: bash

   python Plot_Temperature.py


Expected result
----------------

At :math:`t = 0.1` the temperature field has relaxed to a radially symmetric
profile that rises from :math:`T = 0` at the domain walls to :math:`T = 1`
at the cylinder surface.  The gradient is steepest in the gap between the
cylinder and the nearest wall.  Representative radially-averaged values are:

.. list-table::
   :header-rows: 1
   :widths: 20 30

   * - Radial distance from axis :math:`r`
     - :math:`T_\text{avg}(r)`
   * - 0.15 (cylinder surface)
     - 1.000
   * - 0.20
     - ≈ 0.82
   * - 0.30
     - ≈ 0.57
   * - 0.40
     - ≈ 0.32
   * - ≥ 0.50 (approaching walls)
     - → 0.0

The temperature profile is symmetric about the cylinder axis to within
numerical discretisation error.
