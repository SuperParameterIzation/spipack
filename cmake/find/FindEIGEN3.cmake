
find_package(PkgConfig)
if(NOT DEFINED SPIPACK_EIGEN3_DIR)
	pkg_check_modules(PC_EIGEN3 QUIET EIGEN3)
	set(EIGEN3_DEFINITIONS ${PC_EIGEN3_CFLAGS_OTHER})

	find_path(EIGEN3_INCLUDE_DIR Eigen/Core
    		  HINTS ${PC_EIGEN3_INCLUDEDIR} ${PC_EIGEN3_INCLUDE_DIRS} PATH_SUFFIXES eigen3)
else()
	find_path(EIGEN3_INCLUDE_DIR Eigen/Core
	          HINTS ${SPIPACK_EIGEN3_DIR} PATH_SUFFIXES eigen3 NO_DEFAULT_PATH)
endif()

set(EIGEN3_FOUND 0)
if( EIGEN3_INCLUDE_DIR )
	set(EIGEN3_FOUND 1)
endif()
