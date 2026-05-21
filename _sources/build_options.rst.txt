.. highlight:: rst

.. _build_options:

Build Options Reference
=======================

ExaGOOP is compiled from source using either CMake or GNUmake.  Most
compile-time options are set as variables in the relevant build file; the
table below summarises all user-visible flags.  Detailed descriptions
follow.

.. list-table:: Build option summary
   :header-rows: 1
   :widths: 20 22 22 12 12 12

   * - Feature
     - GNUmake variable
     - CMake variable
     - GNUmake default
     - CMake default
     - Values
   * - Temperature module
     - ``USE_TEMP``
     - ``EXAGOOP_USE_TEMP``
     - ``TRUE``
     - ``ON``
     - ``TRUE``/``FALSE``
   * - Embedded boundary
     - ``USE_EB``
     - *(hardcoded ON)*
     - ``FALSE``
     - ``ON`` (DIM ≥ 2)
     - ``TRUE``/``FALSE``
   * - HDF5 I/O
     - ``AMREX_USE_HDF5``
     - ``EXAGOOP_USE_HDF5``
     - ``FALSE``
     - ``OFF``
     - ``TRUE``/``FALSE``
   * - Floating-point precision
     - ``PRECISION``
     - ``EXAGOOP_PRECISION``
     - ``DOUBLE``
     - ``DOUBLE``
     - ``SINGLE``/``DOUBLE``
   * - Spatial dimensions
     - ``DIM``
     - ``EXAGOOP_DIM``
     - ``2``
     - ``3``
     - ``1``, ``2``, ``3``
   * - Full profiling
     - ``PROFILE``
     - *(not exposed)*
     - ``FALSE``
     - —
     - ``TRUE``/``FALSE``
   * - Lightweight profiling
     - ``TINY_PROFILE``
     - *(not exposed)*
     - ``FALSE``
     - —
     - ``TRUE``/``FALSE``
   * - Debug build
     - ``DEBUG``
     - ``CMAKE_BUILD_TYPE``
     - ``FALSE``
     - *(unset)*
     - ``TRUE``/``FALSE``
   * - Address/UB sanitizers
     - ``FSANITIZER``
     - *(not exposed)*
     - ``FALSE``
     - —
     - ``TRUE``/``FALSE``


.. note::

   GNUmake ``TRUE``/``FALSE`` values are translated to ``1``/``0`` by
   the line ``DEFINES += -DTRUE=1 -DFALSE=0`` in
   ``Build_Gnumake/GNUmakefile``, so ``#if USE_TEMP`` guards work
   correctly for both build paths.


----

USE_TEMP — Temperature Module
------------------------------

Compiles in the optional coupled heat-transfer solver.  See
:ref:`temperature_module` for the full description of the physics and
input parameters.

**GNUmake**

.. code-block:: makefile

   USE_TEMP = TRUE    # enable  (default)
   USE_TEMP = FALSE   # disable

Set this in the relevant ``GNUmakefile`` (``Build_Gnumake/GNUmakefile``
or a test-specific file) **before** the ``include Make.defs`` line.
The build system injects ``-DUSE_TEMP=1`` (enabled) or ``-DUSE_TEMP=0``
(disabled) automatically.

**CMake**

.. code-block:: bash

   cmake -DEXAGOOP_USE_TEMP=ON  ..    # enable  (default)
   cmake -DEXAGOOP_USE_TEMP=OFF ..    # disable

**What it enables**

When ``USE_TEMP=1``:

- Five additional nodal fields are added to the ``MultiFab``
  (``NUM_STATES`` grows from 20 to 25).
- Each particle carries five thermal real-data slots
  (temperature, specific heat, thermal conductivity, heat flux,
  heat source).
- The thermal solve sequence (P2G deposition, nodal advance,
  G2P interpolation) runs every time step after the momentum update.


----

USE_EB — Embedded Boundaries
------------------------------

Enables AMReX Embedded Boundary (EB / cut-cell) support for domains
with complex geometric boundaries.

**GNUmake**

.. code-block:: makefile

   USE_EB = FALSE   # disable (default)
   USE_EB = TRUE    # enable

When ``TRUE``, ``Build_Gnumake/GNUmakefile`` additionally includes
``$(AMREX_HOME)/Src/EB/Make.package``.

**CMake**

EB is enabled by default for DIM ≥ 2 builds via a hardcoded ``set``
in ``CMakeLists.txt``; there is no command-line ``-D`` switch to
disable it without editing the file.  The effective preprocessor
symbol is set via:

.. code-block:: cmake

   # Inside CMakeLists.txt — not overridable from the command line
   set(EXAGOOP_ENABLE_EB ON)

**Enforcement rule**

EB is silently forced to ``FALSE``/``0`` in 1D builds regardless of
what was requested:

.. code-block:: makefile

   # GNUmake — automatic when DIM=1
   ifeq ($(DIM),1)
     USE_EB := FALSE
   endif

A warning is printed if ``USE_EB=TRUE`` was explicitly set for a 1D
build.

**What it enables**

- Registers the AMReX EB factory so cut-cell normals, apertures, and
  volume fractions are available on the background grid.
- Activates the ``mpm_eb.cpp`` routines for EB-aware boundary
  velocity enforcement.
- Requires AMReX to be compiled with ``AMREX_USE_EB`` (handled
  automatically by the build system when ``USE_EB=TRUE``).


----

AMREX_USE_HDF5 — HDF5 Checkpoint/Restart
------------------------------------------

Links against an HDF5 library and enables HDF5-format checkpoint and
restart files through AMReX's native HDF5 I/O layer.

**GNUmake**

HDF5 is configured in ``Build_Gnumake/Make.local`` so that machine- or
environment-specific library paths stay out of the main makefile.
A typical ``Make.local`` block looks like:

.. code-block:: makefile

   # --- Make.local (machine-specific) ---
   HDF5_HOME = /path/to/hdf5          # root of HDF5 installation
   HDF5_INC  = $(HDF5_HOME)/include
   HDF5_LIB  = $(HDF5_HOME)/lib

   USE_HDF5            = TRUE         # signal AMReX to enable HDF5
   AMREX_USE_HDF5      = TRUE         # define the preprocessor symbol
   AMREX_USE_HDF5_PARALLEL = TRUE     # enable parallel (collective) HDF5

   CXXFLAGS += -I$(HDF5_INC)
   LDFLAGS  += -L$(HDF5_LIB) -lhdf5 -lhdf5_hl

On systems using Cray modules (e.g. Kestrel at the National Laboratory of the Rockies), the module
environment automatically provides the correct paths, so only the
``USE_HDF5``/``AMREX_USE_HDF5`` flags are needed.

The default in a clean repository is ``FALSE`` (HDF5 disabled).

**CMake**

.. code-block:: bash

   cmake -DEXAGOOP_USE_HDF5=OFF ..           # disable (default)
   cmake -DEXAGOOP_USE_HDF5=ON  ..           # enable serial HDF5
   cmake -DEXAGOOP_USE_HDF5_PARALLEL=ON ..   # enable parallel HDF5

These options are forwarded to AMReX as ``AMReX_HDF5`` and
``AMReX_HDF5_PARALLEL``; CMake's ``find_package(HDF5)`` locates the
library automatically if ``HDF5_ROOT`` or the module environment is
set.

**What it enables**

- AMReX defines the preprocessor symbol ``AMREX_USE_HDF5``.
- ``mpm_init.cpp`` uses ``#ifdef AMREX_USE_HDF5`` guards to select
  HDF5 checkpoint reading instead of the native AMReX plot-file
  format.
- Parallel HDF5 (``AMREX_USE_HDF5_PARALLEL``) allows collective I/O
  across all MPI ranks for large-scale restart files.


----

PRECISION — Floating-Point Width
----------------------------------

Selects whether ``amrex::Real`` resolves to ``float`` (32-bit) or
``double`` (64-bit) throughout the solver.

**GNUmake**

.. code-block:: makefile

   PRECISION = DOUBLE   # 64-bit (default)
   PRECISION = SINGLE   # 32-bit

**CMake**

.. code-block:: bash

   cmake -DEXAGOOP_PRECISION=DOUBLE ..   # 64-bit (default)
   cmake -DEXAGOOP_PRECISION=SINGLE ..   # 32-bit

**Preprocessor symbols injected**

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Setting
     - Symbol defined
     - String constant
   * - ``DOUBLE``
     - ``AMREX_USE_DOUBLE``
     - ``EXAGOOP_PRECISION_STR="double"``
   * - ``SINGLE``
     - ``AMREX_USE_FLOAT``
     - ``EXAGOOP_PRECISION_STR="single"``

The string constant is printed in the startup banner.  Any value other
than ``SINGLE`` or ``DOUBLE`` is a fatal CMake error; the GNUmake
path silently applies no precision flags if the value is unrecognised.


----

DIM — Spatial Dimensions
--------------------------

Selects the number of spatial dimensions compiled into the executable.
This is a **compile-time constant**; a single binary cannot switch
between 1D, 2D, and 3D at run time.

**GNUmake**

.. code-block:: makefile

   DIM = 2   # 2-D build (default in Build_Gnumake/GNUmakefile)
   DIM = 3   # 3-D build
   DIM = 1   # 1-D build

Test-specific makefiles (e.g. ``Tests/1D_Heat_Conduction/GNUmakefile``)
override this to match the problem dimension.

**CMake**

.. code-block:: bash

   cmake -DEXAGOOP_DIM=3 ..   # 3-D (default)
   cmake -DEXAGOOP_DIM=2 ..   # 2-D
   cmake -DEXAGOOP_DIM=1 ..   # 1-D

**What it controls**

- AMReX defines ``AMREX_SPACEDIM`` to the chosen value; all
  dimensional arrays (``GpuArray``, loop bounds, ``IntVect``) are
  sized at compile time.
- The executable name includes the dimension suffix, e.g.
  ``ExaGOOP2d.gnu.MPI.ex`` or ``ExaGOOP3d.gnu.MPI.ex``.
- Setting ``DIM=1`` forces ``USE_EB=FALSE`` (AMReX EB is not
  supported in 1D).
- ``mpm.ppc``, ``mpm.prob_lo``, and all other spatially-indexed
  parameters expect exactly ``SPACEDIM`` entries in the input file.
  Per-face momentum BC keys (``mpm.bc_xlo_mom``, ``mpm.bc_ylo_mom``,
  etc.) exist for each of the ``2 * SPACEDIM`` domain faces.


----

PROFILE — Full AMReX Profiling
--------------------------------

Enables AMReX's detailed timer-based profiler, which wraps every
``BL_PROFILE`` region and produces a ``bl_prof/`` output directory
containing call-graph timing data.

**GNUmake**

.. code-block:: makefile

   PROFILE = FALSE   # disable (default)
   PROFILE = TRUE    # enable

Related flags (all default ``FALSE``):

.. code-block:: makefile

   COMM_PROFILE  = TRUE   # also profile MPI communication
   TRACE_PROFILE = TRUE   # record full call traces
   MEM_PROFILE   = TRUE   # track memory usage

**CMake**

Full profiling is an AMReX-layer flag and is **not** exposed as an
``-DEXAGOOP_…`` option in ``CMakeLists.txt``.  Pass it directly to
the AMReX sub-build or enable it inside ``CMakeLists.txt`` before the
``add_subdirectory(amrex)`` call.

**What it enables**

- AMReX defines ``AMREX_PROFILING`` (and optionally
  ``AMREX_COMM_PROFILING``, ``AMREX_TRACE_PROFILING``).
- A ``.PROF`` suffix is appended to the GNUmake executable name.
- ``PROFILE=TRUE`` and ``TINY_PROFILE=TRUE`` are **mutually
  exclusive**; AMReX will error if both are set.


----

TINY_PROFILE — Lightweight Profiling
--------------------------------------

Enables AMReX's low-overhead profiler, which accumulates per-region
wall-clock totals and prints a summary table at the end of the run.
Use this as a first-pass performance diagnostic when full ``PROFILE``
overhead is undesirable.

**GNUmake**

.. code-block:: makefile

   TINY_PROFILE = FALSE   # disable (default)
   TINY_PROFILE = TRUE    # enable (only when PROFILE = FALSE)

**CMake**

Not exposed as an ExaGOOP CMake option; configure at the AMReX level
if needed.

**What it enables**

- AMReX defines ``AMREX_TINY_PROFILING``.
- A ``.TPROF`` suffix is appended to the GNUmake executable name.
- Mutually exclusive with ``PROFILE=TRUE``.


----

DEBUG — Debug Build
---------------------

Switches from an optimised release build to a debug build: disables
most optimisations, enables full debug symbols, activates AMReX
runtime assertions, and adds array-bounds checking.

**GNUmake**

.. code-block:: makefile

   DEBUG = FALSE   # optimised release (default)
   DEBUG = TRUE    # debug

When ``TRUE``, AMReX's ``gnu.mak`` (or equivalent compiler file)
substitutes ``-g -O0 -ggdb -ftrapv`` for the release flags
``-g1 -O3``, and defines ``AMREX_DEBUG``.

**CMake**

.. code-block:: bash

   cmake -DCMAKE_BUILD_TYPE=Release ..   # optimised (recommended default)
   cmake -DCMAKE_BUILD_TYPE=Debug   ..   # debug

The CMake path defines ``EXAGOOP_BUILD_TYPE="Debug"`` or
``"Release"`` in the startup banner via ``add_compile_definitions``.
Standard ``-DNDEBUG`` suppression applies for Release builds.

**What it enables**

- AMReX array-bounds and precondition assertions are active.
- A ``.DEBUG`` suffix is appended to the GNUmake executable name.
- Incompatible with GPU profiling runs due to severe performance
  degradation.


----

FSANITIZER — Address and Undefined-Behaviour Sanitizers
---------------------------------------------------------

Instruments the binary with LLVM/GCC runtime sanitizers to detect
memory errors, undefined behaviour, and pointer misuse at run time.

**GNUmake**

.. code-block:: makefile

   FSANITIZER = FALSE   # disable (default)
   FSANITIZER = TRUE    # enable AddressSanitizer + UBSan

A related flag enables thread-safety checking independently:

.. code-block:: makefile

   THREAD_SANITIZER = FALSE   # disable (default)
   THREAD_SANITIZER = TRUE    # enable ThreadSanitizer

``FSANITIZER`` and ``THREAD_SANITIZER`` should not be set ``TRUE``
simultaneously — the two sanitizer runtimes conflict.

**CMake**

Not exposed as an ExaGOOP CMake option.  Pass the sanitizer flags
manually via ``CMAKE_CXX_FLAGS`` if needed:

.. code-block:: bash

   cmake -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined" ..

**Compiler flags injected (via AMReX's** ``gnu.mak`` **)**

.. code-block:: text

   -fsanitize=address
   -fsanitize=undefined
   -fsanitize=pointer-compare
   -fsanitize=pointer-subtract
   -fsanitize=builtin
   -fsanitize=pointer-overflow

**What it enables**

- AddressSanitizer catches heap/stack buffer overflows, use-after-free,
  and use-after-return.
- UndefinedBehaviorSanitizer detects signed-integer overflow, null
  pointer dereferences, misaligned accesses, and invalid enum values.
- Pointer sanitizers flag illegal pointer arithmetic.
- Requires a runtime library (``libasan``, ``libubsan``); the
  executable is significantly slower and uses more memory than a
  release build.

.. note::

   Sanitizers are incompatible with GPU backends (``USE_CUDA``,
   ``USE_HIP``).  Use them on CPU-only debug builds only.
