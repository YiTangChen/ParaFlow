# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.18)  # 3.18+ needed for first-class CUDA support
project (OSUFlow)
set (CMAKE_CXX_STANDARD 17)

# **************** optional CUDA backend ****************
# Enable with: -DOSUFLOW_ENABLE_CUDA=ON
# Host compiler (mpicxx) is inherited via CMAKE_CUDA_HOST_COMPILER so MPI
# headers remain visible inside .cu files.
option(OSUFLOW_ENABLE_CUDA "Build the CUDA particle-tracing backend" OFF)
if(OSUFLOW_ENABLE_CUDA)
    # Force nvcc to use the same host compiler (mpicxx) we're already using
    # for .C/.cpp. Must be set BEFORE enable_language(CUDA).
    if(NOT DEFINED CMAKE_CUDA_HOST_COMPILER)
        set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}" CACHE STRING
            "Host compiler used by nvcc" FORCE)
    endif()
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        # Perlmutter A100 = 80; H100/GH200 = 90. Override with -DCMAKE_CUDA_ARCHITECTURES=...
        set(CMAKE_CUDA_ARCHITECTURES 80)
    endif()
    message(STATUS "OSUFlow: CUDA backend ENABLED (arch=${CMAKE_CUDA_ARCHITECTURES}, host=${CMAKE_CUDA_HOST_COMPILER})")
endif()

# **************** link netcdf library ****************
if(NOT NETCDF_DIR)
    message(FATAL_ERROR "NETCDF_DIR is not set. Pass -DNETCDF_DIR=<path> to cmake.")
endif()
include_directories         (${NETCDF_DIR}/include)
###### NERSC
set                         (NetCDF_LIBRARIES "${NETCDF_DIR}/lib/libnetcdf.so" "${NETCDF_DIR}/lib/libnetcdf_c++4.so")
###### OSC
# set                         (NetCDF_LIBRARIES "${NETCDF_DIR}/lib/libnetcdf.so" "${NETCDF_CXX_DIR}/lib64/libnetcdf-cxx4.so")

# grab all the header files
file(GLOB HEADER_FILES "*.h" )

list(APPEND SRC_FILES
	Grid.C       
	polynomials.C  
	TimeVaryingFieldLine.C 
	eigenvals.C  
	Interpolator.C  
	Rake.C	    
	Topology.C 
	eigenvecs.C  
	IsoSurf.C
	Solution.C
	triangulator.C 
	Element.C
	StreakLine.C
	VectorMatrix.C
	Field.C 
	PathLine.C
	Streamline.C
	FieldLine.C
	Plot3DReader.C 
	TimeLine.C
	OSUFlow.C
	FileReader.C
	Plot3DReader.C
	CurvilinearGrid.C
	MPASOReader.C
	)

if(OSUFLOW_ENABLE_CUDA)
    list(APPEND SRC_FILES
        GPU/CUDA/kernels/MPASODeviceField.cu
        GPU/CUDA/kernels/Kernel/MPASOTracerKernels.cu
        GPU/CUDA/MPASOGPUTracer.cpp
        )
endif()

# build the project as a library
add_library(${PROJECT_NAME} ${LIB_TYPE} ${SRC_FILES} ${HEADER_FILES})

target_link_libraries       (${PROJECT_NAME} PRIVATE ${NetCDF_LIBRARIES})

if(OSUFLOW_ENABLE_CUDA)
    target_compile_definitions(${PROJECT_NAME} PUBLIC OSUFLOW_ENABLE_CUDA)
    # Source-root on include path so files under GPU/CUDA can #include
    # top-level headers like VectorMatrix.h / Grid.h / Solution.h.
    # GPU/CUDA/kernels is the shared kernel-core submodule; MPASOGPUTracer.cpp
    # #includes "MPASODeviceField.h" / "Kernel/MPASOTracerKernels.h" from it.
    target_include_directories(${PROJECT_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/GPU/CUDA
        ${CMAKE_CURRENT_SOURCE_DIR}/GPU/CUDA/kernels)
    # Link CUDA runtime. Static lib still records this as a usage requirement,
    # but downstream executables must also link cudart explicitly because they
    # consume libOSUFlow.a by raw path, not as a CMake target.
    target_link_libraries(${PROJECT_NAME} PUBLIC CUDA::cudart)
    # Leave separable compilation OFF: all __device__ symbols live in one .cu
    # file, and enabling it requires a device-link step that isn't compatible
    # with how ParaFlow links libOSUFlow.a as a raw archive.
endif()

# ADD-BY-LEETEN 09/01/2011-BEGIN
set_target_properties(${PROJECT_NAME} PROPERTIES 
	DEBUG_OUTPUT_NAME "${PROJECT_NAME}_d"
	RELEASE_OUTPUT_NAME "${PROJECT_NAME}"
	)

install(
	FILES ${HEADER_FILES}
	DESTINATION include/osuflow)
	
install(
	TARGETS ${PROJECT_NAME}
	DESTINATION lib/osuflow)
