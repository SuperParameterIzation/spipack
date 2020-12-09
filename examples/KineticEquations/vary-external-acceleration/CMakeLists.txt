cmake_minimum_required (VERSION 3.4)
set(CMAKE_MACOSX_RPATH 1)
project(RUN_MICROSCALE_MODEL)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # Only look for SPIPACK when building example as standalone project
  find_package(SPIPACK REQUIRED HINTS ${SPIPACK_DIR})

  set(CMAKE_CXX_FLAGS ${SPIPACK_CXX_FLAGS})
  set(CMAKE_CXX_COMPILER ${SPIPACK_CXX_COMPILER})

  include_directories(${SPIPACK_INCLUDE_DIRS} ${SPIPACK_EXTERNAL_INCLUDE_DIRS})
endif()

set(LINEAR_SOURCE
  LinearExternalAcceleration.cpp

  _run-linear-external-acceleration.cpp
)

add_executable(run-linear-external-acceleration ${LINEAR_SOURCE})
target_link_libraries(run-linear-external-acceleration ${SPIPACK_EXTERNAL_LIBRARIES} ${SPIPACK_LIBRARIES})
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # add dependencies if this is not a stand alone project
  add_dependencies(run-linear-external-acceleration spiNS spiTS spiKE)
endif()

set(QUADRATIC_SOURCE
  QuadraticExternalAcceleration.cpp

  _run-quadratic-external-acceleration.cpp
)

add_executable(run-quadratic-external-acceleration ${QUADRATIC_SOURCE})
target_link_libraries(run-quadratic-external-acceleration ${SPIPACK_EXTERNAL_LIBRARIES} ${SPIPACK_LIBRARIES})
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # add dependencies if this is not a stand alone project
  add_dependencies(run-quadratic-external-acceleration spiNS spiTS spiKE)
endif()

set(COMPONENTWISE_QUADRATIC_SOURCE
  ComponentwiseQuadraticExternalAcceleration.cpp

  _run-componentwise-quadratic-external-acceleration.cpp
)

add_executable(run-componentwise-quadratic-external-acceleration ${COMPONENTWISE_QUADRATIC_SOURCE})
target_link_libraries(run-componentwise-quadratic-external-acceleration ${SPIPACK_EXTERNAL_LIBRARIES} ${SPIPACK_LIBRARIES})
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # add dependencies if this is not a stand alone project
  add_dependencies(run-componentwise-quadratic-external-acceleration spiNS spiTS spiKE)
endif()
