# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.2)
project (OSUFlow) 
set (CMAKE_CXX_STANDARD 17)

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

# build the project as a library
add_library(${PROJECT_NAME} ${LIB_TYPE} ${SRC_FILES} ${HEADER_FILES})

target_link_libraries       (${PROJECT_NAME} PRIVATE ${NetCDF_LIBRARIES})

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
