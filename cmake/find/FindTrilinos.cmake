find_package(Trilinos HINTS ${SPIPACK_TRILINOS_DIR})

list(APPEND SPIPACK_LINK_DIRS ${Trilinos_INSTALL_DIR}/lib)
