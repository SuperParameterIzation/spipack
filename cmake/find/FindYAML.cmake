find_package(PkgConfig)

if(NOT DEFINED SPIPACK_YAML_DIR)
	pkg_check_modules(PC_YamlCpp QUIET libyaml-cpp)
	set(YamlCpp_DEFINITIONS ${PC_YamlCpp_CFLAGS_OTHER})

	find_path(YamlCpp_INCLUDE_DIR NAMES yaml-cpp/yaml.h
          HINTS ${PC_YamlCpp_INCLUDEDIR} ${PC_YamlCpp_INCLUDE_DIRS} /usr/local/include
          PATH_SUFFIXES yaml-cpp YamlCpp )

	find_library(YamlCpp_LIBRARY NAMES yaml-cpp
             HINTS ${PC_YamlCpp_LIBDIR} ${PC_YamlCpp_LIBRARY_DIRS} )

else()
	find_path(YamlCpp_INCLUDE_DIR NAMES yaml-cpp/yaml.h
	          HINTS ${SPIPACK_YAML_DIR}
		  PATH_SUFFIXES include NO_DEFAULT_PATH)

	find_library(YamlCpp_LIBRARY NAMES yaml-cpp
	             HINTS ${SPIPACK_YAML_DIR}
		     PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
endif()

if( YamlCpp_LIBRARY )
	set(YAML_FOUND 1)
else()
	set(YAML_FOUND 0)
endif()
