cmake_minimum_required (VERSION 2.8)
set(CMAKE_MACOSX_RPATH 1)
project(SPIPACK)

# opening message
message(STATUS "")
message(STATUS "Welcome to the SuperParameterIzed software PACKage (SPIPACK)!")
message(STATUS "")

# tell cmake where to find the cmake source code
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)

# find dependencies
include(FindDependencies)

# print external libraries
message(STATUS "")
message(STATUS "External libraries: " ${SPIPACK_EXTERNAL_LIBRARIES})
message(STATUS "External includes: " ${SPIPACK_EXTERNAL_INCLUDE_DIRS})
message(STATUS "")
