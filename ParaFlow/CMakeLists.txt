project                     (ParaFlow)
cmake_minimum_required      (VERSION 3.13)
set                         (CMAKE_CXX_STANDARD 17)
add_compile_options         (-Wno-deprecated-declarations)
add_compile_options         (-Wno-macro-redefined)
set                         (OSUFLOW_SRC "./OSUFlow")
set                         (libraries pthread)

find_package                (MPI REQUIRED)
find_package                (yaml-cpp REQUIRED)
# find_package                (MOPS REQUIRED)
link_directories            (./OSUFlow)

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
set                         (NetCDF_LIBRARIES "${NETCDF_DIR}/lib/libnetcdf.so" "${NETCDF_CXX_DIR}/lib64/libnetcdf-c++4.so")
###### OSC
# set                         (NetCDF_LIBRARIES "${NETCDF_DIR}/lib/libnetcdf.so" "${NETCDF_CXX_DIR}/lib64/libnetcdf-cxx4.so")

# ##############################################################################  


add_library(ParaFlow "")
target_sources(ParaFlow
  PRIVATE
    ParaFlow.cpp
  PUBLIC
    ParaFlow.hpp
    block.hpp
    utils.hpp
)

target_include_directories(ParaFlow
    PUBLIC
        "./OSUFlow"
    PRIVATE
        "./diy/include"
        ${MPI_CXX_INCLUDE_DIRS}
        ${YAML_CPP_INCLUDE_DIR}
)

target_link_directories(ParaFlow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${OSUFLOW_SRC})

target_link_libraries(ParaFlow
    PUBLIC
        yaml-cpp
        ${NetCDF_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        ${OSUFLOW_SRC}/libOSUFlow.a
)


target_compile_options(ParaFlow PRIVATE
  -pthread -Wdeprecated-declarations
)

target_link_options(ParaFlow PRIVATE
  -pthread -Wdeprecated-declarations
)
