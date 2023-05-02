set(CMAKE_VERBOSE_MAKEFILE TRUE)

add_executable              (Streamline Streamline.cpp)
target_link_libraries       (Streamline ParaFlow ${libraries})

add_library                 (ParaFlow "")
target_sources(ParaFlow
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/ParaFlow.cpp
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/ParaFlow.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/block.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/utils.hpp
  )
target_link_libraries       (ParaFlow OSUFlow)
target_link_directories     (ParaFlow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${OSUFLOW_SRC})

install(TARGETS Streamline
        DESTINATION ${CMAKE_INSTALL_PREFIX}/ParaFlow
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_WRITE WORLD_EXECUTE)