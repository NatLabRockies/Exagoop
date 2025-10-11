1. Copy the ExaGOOP executable to the current folder.
   ```bash
   cp ../../Build_Gnumake/ExaGOOP3d.gnu.MPI.ex .
   ```

2. Run the python script generate_particle_and_inputfiles.py using the arguments as given below:

   ```bash
   python generate_particle_and_inputfiles.py <no_of_cells_in_x> <buffery> <periodicity> <np_per_cell_x> <order_scheme> <alpha_pic_flip> <stress_update_scheme> <CFL> <Output folder name>   
   ```
   
   For example, the user can use the following statement,
   
   ```bash
   python generate_particle_and_inputfiles.py 25 3 1 1 1 1 1 0.1 test   
   ```
   This should generate two files: the input file 'inputs_axialbar.in' and initial material point file 'mpm_particles.dat'
   
3. Run the ExaGOOP executable

   ```bash
   ./ExaGOOP3d.gnu.MPI.ex inputs_axialbar.in
   ```
   
   At the end of the simulation, one should see the output files: 'AxialBarVel.out.0' and 'AxialBarEnergy.out.0'. These are ASCII files providing details of the time evolution of the total energies and center of mass velocity.
   
4. To view these outputs, the user can use the python scripts to generate image files.
   ```bash
   python plot_energy.py <energy_picture_filename>
   python plot_vel.py <vel_picture_filename>
   ```
   
5. The user can also view the solution files in './particle_files' and './grid_files'. The restart files are output in './checkpoint_files'
       