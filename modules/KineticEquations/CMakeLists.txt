# create the library
add_subdirectory(src)

# add the tests to the list of source files
add_subdirectory(tests)

set(KINETIC_EQUATIONS_TESTS
  ${KINETIC_EQUATIONS_TESTS}

  PARENT_SCOPE
)
