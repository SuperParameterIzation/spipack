set(SPIPACK_LIBRARIES )

add_subdirectory(Tools)
list(APPEND SPIPACK_LIBRARIES spiTS)

add_subdirectory(NumericalSolvers)
list(APPEND SPIPACK_LIBRARIES spiNS)

set(GTEST_SOURCE
      modules/_run-tests.cpp
      ${TOOLS_TESTS}
      ${NUMERICAL_SOLVERS_TESTS}
    PARENT_SCOPE)

set(EXAMPLES_SOURCE
      modules/_run-examples.cpp
      ${NUMERICAL_SOLVERS_EXAMPLES}
    PARENT_SCOPE)

set(SPIPACK_LIBRARIES ${SPIPACK_LIBRARIES} PARENT_SCOPE)
