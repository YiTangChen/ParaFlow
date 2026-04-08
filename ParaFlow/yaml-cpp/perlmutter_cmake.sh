# cmake -DCMAKE_C_COMPILER=/pscratch/sd/t/tpeterka/software/mpich-4.3.0/install/bin/mpicc -DCMAKE_CXX_COMPILER=/pscratch/sd/t/tpeterka/software/mpich-4.3.0/install/bin/mpicxx .
# cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx .

cmake -DCMAKE_C_COMPILER=/opt/cray/pe/mpich/9.0.1/ofi/gnu/12.3/bin/mpicc -DCMAKE_CXX_COMPILER=/opt/cray/pe/mpich/9.0.1/ofi/gnu/12.3/bin/mpicxx .
cmake -DCMAKE_INSTALL_PREFIX=~/third_bin -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=~/third_lib ..


# export PATH=~/third_bin/bin:$PATH
# export LD_LIBRARY_PATH=~/third_lib:$LD_LIBRARY_PATH