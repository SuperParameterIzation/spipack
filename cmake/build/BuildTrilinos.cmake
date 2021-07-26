include(ExternalProject)

ExternalProject_Add(
  TRILINOS
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/trilinos/
    GIT_REPOSITORY git@github.com:trilinos/Trilinos.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/external/trilinos
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -DCMAKE_CXX_STANDARD:STRING=17
    -DTrilinos_ENABLE_Sacado=ON
    -DTPL_ENABLE_BLAS=OFF
    -DBUILD_SHARED_LIBS=ON
  	BUILD_COMMAND make -j5
    INSTALL_COMMAND make -j5 install
)

set(Trilinos_DIR ${CMAKE_BINARY_DIR}/external/trilinos/)

list(APPEND Trilinos_INCLUDE_DIRS
  "${CMAKE_BINARY_DIR}/external/trilinos/include"
)

list(APPEND Trilinos_LIBRARIES
"${CMAKE_BINARY_DIR}/external/trilinos/lib/${library_prefix}teuchoscomm${shared_library_suffix}"
"${CMAKE_BINARY_DIR}/external/trilinos/lib/${library_prefix}teuchoskokkoscompat${shared_library_suffix}"
"${CMAKE_BINARY_DIR}/external/trilinos/lib/${library_prefix}kokkoscontainers${shared_library_suffix}"
"${CMAKE_BINARY_DIR}/external/trilinos/lib/${library_prefix}teuchosparameterlist${shared_library_suffix}"
"${CMAKE_BINARY_DIR}/external/trilinos/lib/${library_prefix}teuchoscore${shared_library_suffix}"
"${CMAKE_BINARY_DIR}/external/trilinos/lib/${library_prefix}kokkoscore${shared_library_suffix}"
"${CMAKE_BINARY_DIR}/external/trilinos/lib/${library_prefix}teuchosparser${shared_library_suffix}"
"${CMAKE_BINARY_DIR}/external/trilinos/lib/${library_prefix}teuchoskokkoscomm${shared_library_suffix}"
"${CMAKE_BINARY_DIR}/external/trilinos/lib/${library_prefix}teuchoskokkoscompat${shared_library_suffix}"
"${CMAKE_BINARY_DIR}/external/trilinos/lib/${library_prefix}sacado${shared_library_suffix}"
)
