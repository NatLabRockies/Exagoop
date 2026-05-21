.. highlight:: rst

.. _udf_moving_wall:

UDF Moving Wall Boundary Conditions
=====================================

ExaGOOP supports spatially- and temporally-varying velocity boundary
conditions at any domain face through a *user-defined function* (UDF)
supplied as a runtime-loaded shared library.  This mechanism allows
arbitrary wall motions such as linearly translating walls,  rotating wall etc. to be specified without recompiling the solver.


Overview
---------

A UDF boundary condition is attached to a single named momentum BC face
(e.g. ``mpm.bc_zhi_mom``) by providing two additional sub-namespace
parameters in the input file: the path to the shared library and the name
of the function to call.  At every time step, ExaGOOP calls the UDF on
the CPU for every node on that face, stores the resulting velocity field
on a device-resident array, and uses that array when enforcing nodal and
particle boundary conditions on the GPU.


UDF interface
--------------

The shared library must export a single C-linkage function with the
following signature:

.. code-block:: c

   void func_name(double x, double y, double z, double t, double vel[3]);

where:

- ``x``, ``y``, ``z`` — physical coordinates of the boundary node (m)
- ``t`` — current simulation time (s)
- ``vel[3]`` — output velocity vector ``{vx, vy, vz}`` to be written by
  the function

The function must be declared ``extern "C"`` (or compiled as plain C) so
that the symbol name is not mangled.  A minimal example for a
solid-body rotation about the :math:`z`-axis:

.. code-block:: c

   #define OMEGA 0.5

   void wall_vel_twist(double x, double y, double z,
                       double t, double vel[3])
   {
       vel[0] = -OMEGA * y;
       vel[1] =  OMEGA * x;
       vel[2] =  0.0;
   }

The function is called once per boundary node per time step on a single
CPU thread; it must be thread-safe and must not call any AMReX or GPU
routines.


Building the shared library
----------------------------

The shared library must be compiled natively on the machine where the
solver will run.  The required compiler flags differ by platform:

.. list-table::
   :header-rows: 1
   :widths: 20 30 25 25

   * - Platform
     - Compiler flags
     - Output extension
     - Detected by
   * - macOS (Darwin)
     - ``-O2 -fPIC -dynamiclib``
     - ``.dylib``
     - ``$(shell uname) == Darwin``
   * - Linux
     - ``-O2 -fPIC -shared``
     - ``.so``
     - otherwise

A portable ``Makefile`` that detects the platform automatically:

.. code-block:: makefile

   CC  = gcc
   SRC = my_wall_vel.c

   UNAME := $(shell uname)
   ifeq ($(UNAME), Darwin)
     TARGET = libmy_wall_vel.dylib
     CFLAGS = -O2 -fPIC -dynamiclib
   else
     TARGET = libmy_wall_vel.so
     CFLAGS = -O2 -fPIC -shared
   endif

   $(TARGET): $(SRC)
   	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) -lm

   clean:
   	rm -f libmy_wall_vel.so libmy_wall_vel.dylib

.. warning::

   A library compiled on Linux (ELF format) cannot be loaded by macOS and
   will cause an ``Abort`` at startup with the message
   ``slice is not valid mach-o file``.  Always compile the UDF on the same
   machine (or the same OS) where the simulation will run.


Linking the solver against ``libdl``
--------------------------------------

The solver must be linked against the dynamic loader library.  Add the
following line to the test-specific ``GNUmakefile`` **before** the
``include Make.rules`` line:

.. code-block:: makefile

   LDFLAGS += -ldl


Input file parameters
----------------------

UDF wall velocity is specified per face.  The face namespace follows the
same pattern as the existing momentum BC keys (``mpm.bc_xlo_mom``,
``mpm.bc_xhi_mom``, etc.).  Two sub-keys are added under the face
namespace:

.. list-table::
   :header-rows: 1
   :widths: 38 12 50

   * - Parameter
     - Type
     - Description
   * - ``mpm.bc_<face>_mom.udf_lib``
     - ``string``
     - Path to the shared library (relative to the run directory or
       absolute).  If omitted, no UDF is loaded for this face.
   * - ``mpm.bc_<face>_mom.udf_func``
     - ``string``
     - Name of the exported C function in the library.  Must match
       exactly (case-sensitive); no default.

Both keys must be present together — specifying only one is an error.
The UDF is only activated for the face where both keys appear; all other
faces continue to use their static ``mpm.wall_vel_lo`` /
``mpm.wall_vel_hi`` values.

**Example** — no-slip rotating top wall in the :math:`z`-direction:

.. code-block:: bash

   mpm.bc_zhi_mom            = noslip
   mpm.bc_zhi_mom.udf_lib    = "./UDF/libwall_twist.dylib"
   mpm.bc_zhi_mom.udf_func   = "wall_vel_twist"

The face BC type (``noslip``, ``slip``, ``partialslip``) continues to
control the treatment of the normal and tangential velocity components;
the UDF replaces the constant wall velocity that would otherwise be read
from ``mpm.wall_vel_hi``.


Implementation details
-----------------------

The CPU–GPU split is necessary because UDF function pointers cannot be
called inside GPU kernels.  ExaGOOP resolves this as follows.

**Per-timestep CPU pre-computation** (``compute_udf_wall_vel_at_nodes``):
For each face that has a UDF loaded, the function is called once per
boundary node with the node's physical coordinates and the current time.
The results are stored in a host-side ``amrex::Vector<Real>`` and then
copied to a ``Gpu::DeviceVector<Real>`` via
``Gpu::copy(Gpu::hostToDevice, ...)``.  The device array is indexed as:

- 2D: ``[j * SPACEDIM + c]`` where ``j`` is the node index along the
  single perpendicular direction
- 3D: ``[(j * (ncells[p1]+1) + k) * SPACEDIM + c]`` where ``j`` and
  ``k`` are node indices in the two perpendicular directions
  (``p0 = (dir==0) ? 1 : 0``, ``p1 = (dir==2) ? 1 : 2``)

**Nodal velocity enforcement** (``apply_udf_nodal_bcs``): A GPU
``ParallelFor`` over the nodal ``MultiFab`` overwrites the velocity at
each boundary node using a direct lookup into the device array.  For a
``noslip`` face, all velocity components are set to the UDF value; for a
``slip`` face, only the normal component is set.

**Particle boundary conditions** (``moveParticles``): When a particle
crosses a UDF face, its wall velocity is computed by bilinear
interpolation (3D) or linear interpolation (2D) of the pre-computed
device array at the particle's transverse position.  This interpolated
value replaces the constant ``wall_vel`` used by the non-UDF path and is
then passed to the standard ``applybc`` kernel.


Supported face labels
----------------------

Any of the six domain-face labels may carry a UDF:
``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, ``zhi``.  Multiple faces
can each have an independent UDF.  The same library file can export
multiple functions, each assigned to a different face.


Limitations
------------

- The UDF is called on a single CPU thread.  For large 3D grids with many
  boundary nodes the call overhead is proportional to the number of nodes
  on the face; the cost is typically negligible compared to the particle
  and grid operations.
- The UDF must be **stateless** (no global mutable state) because the
  order in which nodes are visited is not guaranteed between time steps.
  Time dependence must be expressed solely through the ``t`` argument.
- Shared libraries must be in a format compatible with the OS on which the
  solver runs (Mach-O ``.dylib`` on macOS, ELF ``.so`` on Linux).
  Cross-compiled libraries will fail to load.
