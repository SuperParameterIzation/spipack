find_package(PkgConfig)

if(NOT DEFINED SPIPACK_SPECTRA_DIR)
	find_path(SPECTRA_INCLUDE_DIR Spectra)
else()
	find_path(SPECTRA_INCLUDE_DIR Spectra
	          HINTS ${SPIPACK_SPECTRA_DIR})
endif()

set(SPECTRA_FOUND 1)
if( NOT SPECTRA_INCLUDE_DIR OR NOT EIGEN3_FOUND )
	set(SPECTRA_FOUND 0)
endif()
