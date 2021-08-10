# create the library
add_subdirectory(src)

# add the tests to the list of source files
add_subdirectory(tests)

set(CONSERVATION_EQUATIONS_TESTS
  ${CONSERVATION_EQUATIONS_TESTS}

  PARENT_SCOPE
)

set(PARALLEL_CONSERVATION_EQUATIONS_TESTS
  ${PARALLEL_CONSERVATION_EQUATIONS_TESTS}

  PARENT_SCOPE
)
