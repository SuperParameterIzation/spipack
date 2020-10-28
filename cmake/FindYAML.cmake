find_package(PkgConfig)

if(NOT DEFINED ${CMAKE_PROJECT_NAME}_YamlCpp_DIR)
	pkg_check_modules(PC_YamlCpp QUIET libyaml-cpp)
	set(YamlCpp_DEFINITIONS ${PC_YamlCpp_CFLAGS_OTHER})

	find_path(YamlCpp_INCLUDE_DIR NAMES yaml-cpp/yaml.h
          HINTS ${PC_YamlCpp_INCLUDEDIR} ${PC_YamlCpp_INCLUDE_DIRS} /usr/local/include
          PATH_SUFFIXES yaml-cpp YamlCpp )

	find_library(YamlCpp_LIBRARY NAMES yaml-cpp
             HINTS ${PC_YamlCpp_LIBDIR} ${PC_YamlCpp_LIBRARY_DIRS} )

else()
	find_path(YamlCpp_INCLUDE_DIR NAMES yaml-cpp/yaml.h
	          HINTS ${${CMAKE_PROJECT_NAME}_YamlCpp_DIR}
		  PATH_SUFFIXES include NO_DEFAULT_PATH)

	find_library(YamlCpp_LIBRARY NAMES yaml-cpp
	             HINTS ${${CMAKE_PROJECT_NAME}_YamlCpp_DIR}
		     PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
endif()

set(SPIPACK_BUILT_YAML OFF)

if( NOT YamlCpp_LIBRARY )
  include(ExternalProject)
  set(SPIPACK_BUILT_YAML ON)

  ExternalProject_Add(YAML
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    STEP_TARGETS install
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/yaml
    -DCMAKE_INSTALL_LIBDIR=${CMAKE_BINARY_DIR}/external/yaml/lib
    -DCMAKE_CXX_FLAGS=-w --std=c++11
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DBUILD_SHARED_LIBS=ON
    -DYAML_CPP_BUILD_TESTS=OFF
  )

  if(APPLE)
    set(suffix ".dylib")
  else()
    set(suffix ".so")
  endif()

  set(SPIPACK_EXTERNAL_INCLUDE_DIRS
    ${CMAKE_BINARY_DIR}/external/yaml/include)

  set(SPIPACK_EXTERNAL_LIBRARIES
    ${CMAKE_BINARY_DIR}/external/yaml/lib/libyaml-cpp${suffix})
else()
  list(APPEND SPIPACK_EXTERNAL_INCLUDE_DIRS
    ${YamlCpp_INCLUDE_DIR}
    )

  list(APPEND SPIPACK_EXTERNAL_LIBRARIES
    ${YamlCpp_LIBRARY}
    )
endif()
