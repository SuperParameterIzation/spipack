macro(ExaHyPE_System equation dim)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/exahype/${equation}/buildinfo.h.in "${CMAKE_CURRENT_SOURCE_DIR}/spipack/ConservationEquations/ExaHyPE/${equation}/buildinfo.h" @ONLY)

  # load the exahype file list
  include(${equation}/ExaHyPEFiles)

  add_library(spiEX_${equation} SHARED ${EXAHYPE_${equation}_SOURCE})

  target_include_directories(spiEX_${equation} PRIVATE
  ${PEANO_DIRS}
  spipack/ConservationEquations/ExaHyPE/${equation})

  set_target_properties(spiEX_${equation}
  PROPERTIES COMPILE_FLAGS "-DDim${dim} -DTrackGridStatistics"
  )

  install(TARGETS spiEX_${equation}
        EXPORT ${CMAKE_PROJECT_NAME}Depends
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

  list(APPEND SPIEX_LIBRARIES spiEX_${equation})
endmacro(ExaHyPE_System)

message(STATUS "Building the ExaHyPE library: locally named spiEX")

set(EXAHYPE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/ExaHyPE-Engine)
set(PEANO_DIRS
  ${EXAHYPE_DIR}/Peano/ ${EXAHYPE_DIR}/Peano/peano ${EXAHYPE_DIR}/Peano/tarch ${EXAHYPE_DIR}/Peano/multiscalelinkedcell ${EXAHYPE_DIR}/Peano/sharedmemoryoracles ${EXAHYPE_DIR}/Peano/mpibalancing
)

list(APPEND SPIPACK_EXTERNAL_INCLUDE_DIRS ${EXAHYPE_DIR}/ExaHyPE/)

set(SPIEX_LIBRARIES )

#ExaHyPE_System(CompressibleEuler 2)
ExaHyPE_System(Boltzmann 2)
