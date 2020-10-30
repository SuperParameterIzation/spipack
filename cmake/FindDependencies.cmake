macro(Dependency name)
  # try to find the package
  find_package(${name})

  # assume that we did not build this package
  set(SPIPACK_BUILT_${name} OFF)

  if( ${${name}_FOUND} )
    # check that the package works
    include(Check${name})
    if( NOT ${${name}_COMPILES})
      # build and re-check the build
      include(Build${name})
      set(SPIPACK_BUILT_${name} ON)
    endif()
  else()
    # build and check the build
    include(Build${name})
    set(SPIPACK_BUILT_${name} ON)
  endif()

  # fail if we could not find this dependency
  if( NOT ${${name}_COMPILES} )
    message(FATAL_ERROR "\nSPI-PACK FAILED TO FIND OR BUILD ${name}\n")
  endif()

  # append the dependency information into the include directories and external libraries
  include(Append${name})
endmacro(Dependency)

set(library_prefix "lib")
set(static_library_suffix ".a")
if(APPLE)
  set(shared_library_suffix ".dylib")
else()
  set(shared_library_suffix ".so")
endif()

Dependency(YAML)
Dependency(EIGEN3)
Dependency(BOOST)
Dependency(HDF5)
Dependency(MUQ)
Dependency(GTEST)

# add the header only submodules
list(APPEND SPIPACK_EXTERNAL_INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/external/cereal/include
  ${CMAKE_SOURCE_DIR}/external/nanoflann/include/
  )
