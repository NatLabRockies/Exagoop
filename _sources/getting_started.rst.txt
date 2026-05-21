.. _getting_started:

Getting Started
===============

This page covers everything needed to go from a fresh checkout to a running
ExaGOOP simulation: system requirements, obtaining the source, building with
CMake or GNUmake, and verifying the installation.


System Requirements
-------------------

ExaGOOP runs on macOS and Linux. Windows is not officially supported; Windows
users should use the Windows Subsystem for Linux (WSL).

**Compiler.** A C++17-capable compiler is required:

- GCC 8 or newer
- Clang 3.6 or newer

Microsoft Visual C++ (MSVC) is not offcially supported. However, the developers are not aware of any potential issues preventing Windows 11
builds on non-GPU architectures.

**Build system.** Either of:

- CMake 3.20 or newer (recommended)
- GNU Make 3.81 or newer

**MPI.** An MPI-2 implementation (e.g. OpenMPI, MPICH) is required for
multi-process runs. Single-process CPU builds can be compiled without MPI by
setting ``EXAGOOP_ENABLE_MPI=OFF`` (CMake) or ``USE_MPI=FALSE`` (GNUmake).

**GPU (optional).**

- NVIDIA GPUs: CUDA 11 or newer
- AMD GPUs: ROCm 5.2 or newer

The AMReX framework (bundled as a submodule) handles the low-level GPU
portability layer.


Obtaining the Source
--------------------

Clone the repository with its AMReX submodule in a single step:

.. code-block:: bash

   git clone --recurse-submodules https://github.com/NatLabRockies/Exagoop.git

This creates an ``Exagoop/`` directory. The AMReX sources reside under
``Exagoop/Submodules/amrex``.

Next, set the two environment variables that the build system and test scripts
rely on. Add these to your ``~/.bashrc`` or ``~/.zshrc`` to make them
persistent:

.. code-block:: bash

   export MPM_HOME=/path/to/Exagoop
   export AMREX_HOME=${MPM_HOME}/Submodules/amrex


Building with CMake
-------------------

The ``Build_Cmake/cmake.sh`` script
contains a template ``cmake`` invocation; edit it to match your environment
before running it.

.. code-block:: bash

   cd $MPM_HOME/Build_Cmake
   # Edit cmake.sh as described below, then:
   sh cmake.sh

Key options in ``cmake.sh``:

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Option
     - Default
     - Description
   * - ``-DEXAGOOP_ENABLE_MPI``
     - ``OFF``
     - Enable MPI parallelism
   * - ``-DEXAGOOP_ENABLE_CUDA``
     - ``OFF``
     - Enable NVIDIA GPU support (requires CUDA ≥ 11)
   * - ``-DEXAGOOP_ENABLE_HIP``
     - ``OFF``
     - Enable AMD GPU support (requires ROCm ≥ 5.2)
   * - ``-DAMReX_CUDA_ARCH``
     - ``Auto``
     - CUDA compute capability, e.g. ``80`` for A100
   * - ``-DAMReX_AMD_ARCH``
     - ``gfx90a``
     - ROCm target architecture
   * - ``-DEXAGOOP_USE_TEMP``
     - ``OFF``
     - Compile in the heat-transfer module
   * - ``-DEXAGOOP_PRECISION``
     - ``DOUBLE``
     - Floating-point precision: ``SINGLE`` or ``DOUBLE``
   * - ``-DEXAGOOP_USE_HDF5``
     - ``OFF``
     - Enable HDF5 output

``EXAGOOP_DIM`` is fixed at ``2`` in the CMake build. Users should set these variables corresponding to the dimensionality 
   of the problem they intend to solve.

On a successful build the executable ``ExaGOOP<dim>d.*.exe`` is placed in
``$MPM_HOME/Build_Cmake``. The executable name is determined by the dimensionality used as well other build variables. 


Building with GNUmake
---------------------

Open ``$MPM_HOME/Build_Gnumake/GNUmakefile`` and set the options near the
top of the file:

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Variable
     - Default
     - Description
   * - ``DIM``
     - ``1``
     - Spatial dimension: ``1``, ``2``, or ``3``
   * - ``COMP``
     - ``gnu``
     - Compiler: ``gnu`` or ``clang``
   * - ``USE_MPI``
     - ``TRUE``
     - Enable MPI
   * - ``USE_OMP``
     - ``FALSE``
     - Enable OpenMP threading
   * - ``USE_CUDA``
     - ``FALSE``
     - Enable NVIDIA GPU (CUDA)
   * - ``USE_HIP``
     - ``FALSE``
     - Enable AMD GPU (ROCm/HIP)
   * - ``USE_TEMP``
     - ``FALSE``
     - Compile in the heat-transfer module
   * - ``USE_EB``
     - ``FALSE``
     - Compile in embedded boundary support (requires ``DIM`` ≥ 2)
   * - ``DEBUG``
     - ``FALSE``
     - Build in debug mode
   * - ``PRECISION``
     - ``DOUBLE``
     - Floating-point precision: ``SINGLE`` or ``DOUBLE``

Then build:

.. code-block:: bash

   cd $MPM_HOME/Build_Gnumake
   make -j$(nproc)

The executable is named ``ExaGOOP<dim>d.<comp>[.<parallel>].ex`` and is
placed in the same directory.


Verifying the Build
-------------------

The quickest verification is to run the 1-D axial bar vibration test, which
has no physics-module dependencies and produces a known analytical solution:

.. code-block:: bash

   cd $MPM_HOME/Tests/1D_Axial_Bar_Vibration
   sh Generate_MPs_and_InputFiles.sh
   # This creates an initial material point file (mpm_particles.dat or mpm_particles.h5) and an input file- Inputs_1DAxialBarVibration.inp
   # Now copy the ExaGOOP executable from the build folder and run the solver
   ./ExaGOOP1d.*.ex Inputs_1DAxialBarVibration.ex

If the run completes without errors, the build is working correctly.  For a
richer check, the tutorials section walks through two multi-physics test cases
end-to-end.


Visualizing Output
------------------

ExaGOOP writes output in the AMReX plotfile format. Two types of output files
are produced:

- ``plt*`` files — particle (material point) data
- ``nplt*`` files — nodal grid data

Both can be loaded directly in `ParaView <https://www.paraview.org/>`_ using
the built-in AMReX/BoxLib reader. Open the top-level ``plt*/`` or ``nplt*/``
directory as a dataset.


Getting Help
------------

- For questions and discussion: `GitHub Discussions <https://github.com/NatLabRockies/Exagoop/discussions>`_
- To report a bug: `GitHub Issues <https://github.com/NatLabRockies/Exagoop/issues>`_
  (include your compiler version, MPI implementation, and the build flags used)
