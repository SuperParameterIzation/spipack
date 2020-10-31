cmake_minimum_required (VERSION 2.8)

project(GaussianProcess_CO2)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # Only look for MUQ when building example as standalone project
  find_package(MUQ REQUIRED)

  include_directories(${MUQ_INCLUDE_DIRS})

  set(CMAKE_CXX_FLAGS ${MUQ_CXX_FLAGS})
  set(CMAKE_CXX_COMPILER ${MUQ_CXX_COMPILER})
endif()


add_executable(GaussianProcess_CO2_exe GaussianProcess_CO2.cpp)

target_link_libraries(GaussianProcess_CO2_exe ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})
