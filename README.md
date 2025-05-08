# AM5

# Description of the model:
TODO

# Getting the INPUT data:
TODO

# Compiling the model:
The script build_am5.sh can be run to compile the model.

The script exec/env.sh should be updated to load the environment in your system. <br>
Libraries required:
- NetCDF C and Fortran (77/90) headers and libraries
- Fortran 2003 standard compiler
- Fortran compiler that supports Cray Pointer
- MPI C and Fortran headers and librarie

The file exec/intel-classic.mk should be updated with the flags, based on your compile

# Running the model:
The script run_c96L65_am5f5b7r0_pdclim2010F_gwtest_cqa15.sh can be used to run the model. It is currently set up to run with 576 nodes and 2 openmp threads, but the script has guidelines on how to update this.
