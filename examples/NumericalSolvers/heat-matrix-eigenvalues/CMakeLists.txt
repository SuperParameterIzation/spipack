cmake_minimum_required (VERSION 3.4)
set(CMAKE_MACOSX_RPATH 1)
project(HEAT_MATRIX_EIGENVALUES_EXAMPLE)

find_package(SPIPACK REQUIRED HINTS ${SPIPACK_DIR})

set(CMAKE_CXX_FLAGS ${SPIPACK_CXX_FLAGS})
set(CMAKE_CXX_COMPILER ${SPIPACK_CXX_COMPILER})

include_directories(${SPIPACK_INCLUDE_DIRS})

set(SOURCE
  _heat-matrix-eigenvalues.cpp
)

add_executable(run-heat-matrix-eigenvalues ${SOURCE})
target_link_libraries(run-heat-matrix-eigenvalues ${SPIPACK_EXTERNAL_LIBRARIES} ${SPIPACK_LIBRARIES})
