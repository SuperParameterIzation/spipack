cmake_minimum_required (VERSION 3.0)
set(CMAKE_MACOSX_RPATH 1)
project(SPIPACK)

# opening message
message(STATUS "")
message(STATUS "Welcome to the SuperParameterIzed software PACKage (SPIPACK)!")
message(STATUS "")

###############################################
# Check flags
###############################################

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wall --std=c++17 -pthread")

###############################################
# Find dependencies
###############################################

# tell cmake where to find the cmake source code
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)

# find dependencies
include(FindDependencies)

# print external libraries
message(STATUS "")
message(STATUS "External libraries: " ${SPIPACK_EXTERNAL_LIBRARIES})
message(STATUS "External includes: " ${SPIPACK_EXTERNAL_INCLUDE_DIRS})
message(STATUS "")

include_directories(${SPIPACK_EXTERNAL_INCLUDE_DIRS})

###############################################
# Build SPIPACK
###############################################

add_subdirectory(modules)

add_executable(run-tests ${GTEST_SOURCE})
target_link_libraries(run-tests ${SPIPACK_EXTERNAL_LIBRARIES})
add_dependencies(run-tests GTEST)

###############################################
# Build the documentation (doxygen)
###############################################

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
	message(STATUS "To build documentation, run 'make doc'")

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/documentation/spipack.doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/spipack.doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/spipack.doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)
