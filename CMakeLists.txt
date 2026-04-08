project                 (MPAS_OCEAN_Load)

cmake_minimum_required  (VERSION 3.16)

include                 (CMakePrintHelpers)

set                     (CMAKE_CXX_STANDARD 17)
add_compile_options     (-Wno-deprecated-declarations)
add_compile_options     (-Wno-macro-redefined)

find_package            (yaml-cpp REQUIRED)
include_directories     (${YAML_CPP_INCLUDE_DIR})

link_directories        (./ParaFlow)
link_directories        (./ParaFlow/OSUFlow)


# DIY
find_path                   (DIY_INCLUDE_DIRS diy/types.h)
include_directories         (${DIY_INCLUDE_DIRS}
                             SYSTEM ${MPI_INCLUDE_PATH})

# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
if(NOT NETCDF_DIR)
    message(FATAL_ERROR "NETCDF_DIR is not set. Pass -DNETCDF_DIR=<path> to cmake.")
endif()
include_directories         (${NETCDF_DIR}/include)
###### NERSC
# set                         (NetCDF_LIBRARIES "${NETCDF_DIR}/lib/libnetcdf.so" "${NETCDF_DIR}/lib/libnetcdf_c++4.so")
###### OSC
set                         (NetCDF_LIBRARIES "${NETCDF_DIR}/lib/libnetcdf.so" "${NETCDF_CXX_DIR}/lib64/libnetcdf-cxx4.so")


# ############################################################################## 

add_executable(DrawSubdomain)
target_sources(DrawSubdomain
    PRIVATE
        "DrawSubdomain.cpp"
)

target_include_directories(DrawSubdomain
    PRIVATE
        "./ParaFlow"
        "./ParaFlow/diy/include"
        "./ParaFlow/OSUFlow"
        ${YAML_CPP_INCLUDE_DIR}
)

target_link_libraries(DrawSubdomain
    PRIVATE
        pthread
        ${NetCDF_LIBRARIES}
        ParaFlow
        OSUFlow
        yaml-cpp
)

target_compile_options(DrawSubdomain PRIVATE 
  -pthread -Wdeprecated-declarations
)

target_link_options(DrawSubdomain PRIVATE
  -pthread -Wdeprecated-declarations
)

# # ############################################################################## 
add_executable(ParaFlow_streamline)
target_sources(ParaFlow_streamline
    PRIVATE
        "ParaFlow_streamline.cpp"
)

target_include_directories(ParaFlow_streamline
    PRIVATE
        "./ParaFlow"
        "./ParaFlow/diy/include"
        "./ParaFlow/OSUFlow"
        ${YAML_CPP_INCLUDE_DIR}
)

target_link_libraries(ParaFlow_streamline
    PRIVATE
        pthread
        ${NetCDF_LIBRARIES}
        ParaFlow
        OSUFlow
        yaml-cpp
)

target_compile_options(ParaFlow_streamline PRIVATE -g -rdynamic)
target_link_options(ParaFlow_streamline PRIVATE -rdynamic)

# # ############################################################################## 
add_executable(ParaFlow_pathline)
target_sources(ParaFlow_pathline
    PRIVATE
        "ParaFlow_pathline.cpp"
)

target_include_directories(ParaFlow_pathline
    PRIVATE
        "./ParaFlow"
        "./ParaFlow/diy/include"
        "./ParaFlow/OSUFlow"
        ${YAML_CPP_INCLUDE_DIR}
)

target_link_libraries(ParaFlow_pathline
    PRIVATE
        pthread
        ${NetCDF_LIBRARIES}
        ParaFlow
        OSUFlow
        yaml-cpp
)

target_compile_options(ParaFlow_pathline PRIVATE -g -rdynamic)
target_link_options(ParaFlow_pathline PRIVATE -rdynamic)
