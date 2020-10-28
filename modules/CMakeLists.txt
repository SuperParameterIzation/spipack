add_subdirectory(NumericalSolvers)

set(GTEST_SOURCE
      modules/_run-tests.cpp
      ${NUMERICAL_SOLVERS_TESTS}
    PARENT_SCOPE)
