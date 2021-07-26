cmake_minimum_required (VERSION 3.4)
set(CMAKE_MACOSX_RPATH 1)
project(SPIPACK)

# opening message
message(STATUS "")
message(STATUS "Welcome to the SuperParameterIzed software PACKage (SPIPACK)!")
message(STATUS "")

###############################################
# Check flags
###############################################

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-w -g -pthread -fopenmp -O3 -ldl -fPIC")

###############################################
# Find dependencies
###############################################

# tell cmake where to find the cmake source code
list(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
	${CMAKE_SOURCE_DIR}/cmake/find/
	${CMAKE_SOURCE_DIR}/cmake/check/
	${CMAKE_SOURCE_DIR}/cmake/build/
	${CMAKE_SOURCE_DIR}/cmake/append/
	${CMAKE_SOURCE_DIR}/cmake/exahype/
	)
  
set(SPIPACK_LINK_DIRS )

# find dependencies
include(FindDependencies)

# build ExaHyPE, which we must build as a library ourselves
include(BuildExaHyPE)

# print external libraries
message(STATUS "")
message(STATUS "External libraries: " ${SPIPACK_EXTERNAL_LIBRARIES})
message(STATUS "External includes: " ${SPIPACK_EXTERNAL_INCLUDE_DIRS})
message(STATUS "")

include_directories(${CMAKE_SOURCE_DIR} ${SPIPACK_EXTERNAL_INCLUDE_DIRS})

###############################################
# Build SPIPACK
###############################################

# Create the SPIPACKDirectories.hpp file
configure_file(cmake/SPIPACKDirectories.hpp.in "${CMAKE_CURRENT_SOURCE_DIR}/spipack/SPIPACKDirectories.hpp" @ONLY)

add_subdirectory(modules)

include_directories(${GTEST_INCLUDE_DIR})

add_executable(run-tests ${GTEST_SOURCE})
target_link_libraries(run-tests ${SPIPACK_EXTERNAL_LIBRARIES} ${SPIPACK_LIBRARIES} ${GTEST_LIBRARIES})
add_dependencies(run-tests ${SPIPACK_LIBRARIES})
if( SPIPACK_BUILT_GTEST )
	add_dependencies(run-tests GTEST)
endif()

###############################################
# Install SPIPACK
###############################################

set(SPIPACK_EXPORT_LIBRARIES "")
FOREACH(BUILT_LIB ${SPIPACK_LIBRARIES})
	list(APPEND SPIPACK_EXPORT_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/${library_prefix}${BUILT_LIB}${shared_library_suffix})
ENDFOREACH()

# Create the SPIPACKConfig.cmake file
configure_file(cmake/SPIPACKConfig.cmake.in "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Config.cmake" @ONLY)
configure_file(cmake/SPIPACKConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${CMAKE_PROJECT_NAME}Config.cmake" @ONLY)

# Install the SPIPACKConfig.cmake file
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${CMAKE_PROJECT_NAME}Config.cmake"
        DESTINATION "${CMAKE_INSTALL_PREFIX}"
	    COMPONENT dev)

# install all the SPIPACK headers
install(DIRECTORY spipack
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
	      FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY spipack
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
	      FILES_MATCHING PATTERN "*.h")

# install all the ExaHyPE external files
install(DIRECTORY spipack/ConservationEquations/ExaHyPE
        DESTINATION "${CMAKE_INSTALL_PREFIX}"
	      FILES_MATCHING PATTERN "*.exahype")
install(DIRECTORY spipack/ConservationEquations/ExaHyPE
        DESTINATION "${CMAKE_INSTALL_PREFIX}"
	      FILES_MATCHING PATTERN "exahype.log-filter")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${INSTALL_LIB_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${INSTALL_LIB_DIR})

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

###############################################
# Build the examples
###############################################

#add_subdirectory(examples EXCLUDE_FROM_ALL)
