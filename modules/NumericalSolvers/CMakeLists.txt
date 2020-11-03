# create the library
add_subdirectory(src)

# add the tests to the list of source files
add_subdirectory(tests)

set(NUMERICAL_SOLVERS_TESTS
  ${NUMERICAL_SOLVERS_TESTS}

  PARENT_SCOPE
)

# add the examples to the list of source files
add_subdirectory(examples)

set(NUMERICAL_SOLVERS_EXAMPLES
  ${NUMERICAL_SOLVERS_EXAMPLES}

  PARENT_SCOPE
)
