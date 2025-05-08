#!/bin/bash

cd run
mkdir -p RESTART

# Copy over the executable
cp ../exec/am5f5b7r0_beres_cqa_compile.x .

# Copy or symlink the data

# Run the model
# In this case, the model is using 576 cores and 2 threads
# This can be updated by changing the input.nml
# coupler_nml::atmos_npes = 576
# coupler_nml::atmos_nthreads = 2

# layout(1) * layout(2) * 6 in fv_core_nml and land_model_nml must equal the number of atmos_npes
# fv_core_nml::layout   = 4,24
# land_model_nml::layout   = 4,24

# layout(1) * layout(2) in ice_model_nml must equal atmos_npes
# ice_model_nml::layout = 96,6

# NOTE: layout must be divisible by io_layout

srun --ntasks=576 --cpus-per-task=2 --export=ALL,OMP_NUM_THREADS=2 ./am5f5b7r0_beres_cqa_compile.x | & tee fms.out

echo "----------------- END OF MODEL RUN -----------------"
