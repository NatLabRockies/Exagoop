1. Copy the ExaGOOP executable to the current folder.
   ```bash
   cp ../../Build_Gnumake/ExaGOOP3d.gnu.MPI.ex .
   ```

2. Run the python script generate_particle_and_inputfiles.py without any arguments:   
   
   ```bash
   python generate_particle_and_inputfiles.py
   ```
   This should generate one file: the initial material point file 'mpm_particles.dat'
   
3. Run the ExaGOOP executable

   ```bash
   ./ExaGOOP3d.gnu.MPI.ex inputs_elasticdisk.in
   ```
   
   At the end of the simulation, one should see the output diagnostic file: 'ElasticDiskCollisionEnergy.out.0'. This is an ASCII file providing details of the time evolution of the total energy of the colliding elastic disks.
   
4. To view these output, the user can use the python script to generate image file.
   ```bash
   python plot_energy.py   
   ```
   This produces an image file 'Energy_vs_time.png' that shows the energy evolution over time.
   
   
5. The user can also view the solution files in './particle_files' and './grid_files'. These files can be opened and viewed in Paraview using the AMREX particle and grid files respectively. The restart files are output in './checkpoint_files'. Th


