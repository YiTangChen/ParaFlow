# cmake -DCMAKE_C_COMPILER=/pscratch/sd/t/tpeterka/software/mpich-4.3.0/install/bin/mpicc -DCMAKE_CXX_COMPILER=/pscratch/sd/t/tpeterka/software/mpich-4.3.0/install/bin/mpicxx .
# cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx .


# cmake . -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx \
# -DnetCDF_DIR=$HOME/pstrach_folder/MOPS/third_lib/lib64/cmake/netCDF \
# -Dyaml-cpp_DIR=$HOME/pstrach_folder/MOPS/third_lib/lib64/cmake/yaml-cpp \
# -Dndarray_DIR=$HOME/pstrach_folder/MOPS/third_lib/lib/cmake/ndarray \
# -Dpybind11_DIR=$HOME/pstrach_folder/MOPS/third_lib/share/cmake/pybind11 \
# -DVTK_DIR=$HOME/pstrach_folder/MOPS/third_lib/lib64/cmake/vtk-9.2 \
# -DMOPS_DIR=$HOME/pstrach_folder/MOPS/third_lib/lib/cmake/MOPS

cmake -DCMAKE_C_COMPILER=/opt/cray/pe/mpich/9.0.1/ofi/gnu/12.3/bin/mpicc -DCMAKE_CXX_COMPILER=/opt/cray/pe/mpich/9.0.1/ofi/gnu/12.3/bin/mpicxx .