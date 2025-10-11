1. Copy the ExaGOOP executable to the current folder
   ```bash
   cp $MPM_HOME/Build_Gnumake/ExaGOOP3d.<comp>.<MPI>.ex .
   ```

2. Generate the input file and initial material point file
   ```bash
   python3 generate_particle_and_inputfiles.py <no of cells> <buffery> <periodic> <number of particle per cell> <scheme order> <alpha_pic_flip> <stress update scheme> <CFL> <Solution Folder Name>
   ```
   For example, 

   ```bash
   python3 generate_particle_and_inputfiles.py 25 3 1 1 3 1 1 0.1 25Nx_3Buf_1Per_1PPC_3Order_0Alpha_1StressScheme_0.1CFL   
   ```
   This should generate an input file (inputs_axialbar.dat) and initial material point file (mpm_particles.dat). 

3. Run ExaaGOOP

   ```bash
   ExaGOOP3d.<comp>.<MPI>.ex inputs_axialbar.dat  
   ```

   At the end of the simulation, two output files are written- AxialBarVel.out.0 and AxialBarEnergy.out.0, which are ascii files showing time evolution of total energies and center of mass velocities
   
4. To plot the energy and center of mass velocity, run the following python scripts:
  ```bash
   python plot_energy.py <energy_pic_filename.png>
   python plot__vel.py <velocity_pic_filename.png>
   ```
