find . -mindepth 1 ! -name README.md ! -name 'clean_folder.sh' ! -name 'plot_energy.py' ! -name 'plot_vel.py' ! -name 'generate_particle_and_inputfiles.py' -exec rm -rf {} +
