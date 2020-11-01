# create the library
add_subdirectory(src)

# add the tests to the list of source files
add_subdirectory(tests)

set(TOOLS_TESTS
  ${TOOLS_TESTS}

  PARENT_SCOPE
)
