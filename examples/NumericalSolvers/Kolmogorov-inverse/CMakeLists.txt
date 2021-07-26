cmake_minimum_required (VERSION 3.4)
set(CMAKE_MACOSX_RPATH 1)
project(KOLMOGOROV_OPERATOR_EIGENVALUES_EXAMPLE)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # Only look for SPIPACK when building example as standalone project
  find_package(SPIPACK REQUIRED HINTS ${SPIPACK_DIR})

  set(CMAKE_CXX_FLAGS ${SPIPACK_CXX_FLAGS})
  set(CMAKE_CXX_COMPILER ${SPIPACK_CXX_COMPILER})

  include_directories(${SPIPACK_INCLUDE_DIRS} ${SPIPACK_EXTERNAL_INCLUDE_DIRS})
  link_directories(${SPIPACK_LINK_DIRS})
endif()

set(SOURCE
  _Kolmogorov-inverse.cpp
)

message(STATUS ${SPIPACK_LIBRARIES})

add_executable(run-Kolmogorov-inverse ${SOURCE})
target_link_libraries(run-Kolmogorov-inverse ${SPIPACK_EXTERNAL_LIBRARIES} ${SPIPACK_LIBRARIES})
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # add dependencies if this is not a stand alone project
  add_dependencies(run-Kolmogorov-inverse spiNS spiTS)
endif()
