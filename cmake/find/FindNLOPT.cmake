find_package(PkgConfig)

if(NOT DEFINED SPIPACK_NLOPT_DIR)
	pkg_check_modules(PC_NLOPT QUIET libnlopt)
	set(NLOPT_DEFINITIONS ${PC_NLOPT_CFLAGS_OTHER})

	find_path(NLOPT_INCLUDE_DIR nlopt.h
          HINTS ${PC_NLOPT_INCLUDEDIR} ${PC_NLOPT_INCLUDE_DIRS}
          PATH_SUFFIXES nlopt )

	find_library(NLOPT_LIBRARY NAMES nlopt nlopt_cxx
             HINTS ${PC_NLOPT_LIBDIR} ${PC_NLOPT_LIBRARY_DIRS} )
else()

  message(STATUS "HAS HINTS for NLOPT")
  message(STATUS ${SPIPACK_NLOPT_DIR})

	find_path(NLOPT_INCLUDE_DIR nlopt.h
	          HINTS ${SPIPACK_NLOPT_DIR}/include)

	find_library(NLOPT_LIBRARY NAMES nlopt nlopt_cxx
	             HINTS ${SPIPACK_NLOPT_DIR}/lib)
endif()

set(NLOPT_FOUND 0)
if( NLOPT_LIBRARY AND NLOPT_INCLUDE_DIR )
  set(NLOPT_FOUND 1)
endif()
