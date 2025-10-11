find . -mindepth 1 ! -name 'clean_folder.sh' ! -name 'inputs_elasticdisk.in' ! -name 'plot_energy.py' ! -name 'generate_particle_and_inputfiles.py' ! -name 'README.md' -exec rm -rf {} +
