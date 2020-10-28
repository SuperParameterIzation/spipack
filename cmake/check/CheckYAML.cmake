# make sure that YAML works
set(CMAKE_REQUIRED_INCLUDES ${YamlCpp_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
include(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <iostream>
  int main() {
    return 0;
  }
  "
  YAML_CODE_COMPILES)

set(YAML_COMPILES 1)
if( NOT YAML_CODE_COMPILES )
  set(YAML_COMPILES 0)
endif()
