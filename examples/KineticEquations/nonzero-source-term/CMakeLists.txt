cmake_minimum_required (VERSION 3.4)
set(CMAKE_MACOSX_RPATH 1)
project(RUN_MICROSCALE_MODEL)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # Only look for SPIPACK when building example as standalone project
  find_package(SPIPACK REQUIRED HINTS ${SPIPACK_DIR})

  set(CMAKE_CXX_FLAGS ${SPIPACK_CXX_FLAGS})
  set(CMAKE_CXX_COMPILER ${SPIPACK_CXX_COMPILER})

  include_directories(${SPIPACK_INCLUDE_DIRS} ${SPIPACK_EXTERNAL_INCLUDE_DIRS})
endif()

set(SOURCE
  _run-nonzero-source-term.cpp
)

add_executable(run-nonzero-source-term ${SOURCE})
target_link_libraries(run-nonzero-source-term ${SPIPACK_EXTERNAL_LIBRARIES} ${SPIPACK_LIBRARIES})
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # add dependencies if this is not a stand alone project
  add_dependencies(run-nonzero-source-term spiNS spiTS spiKE)
endif()
