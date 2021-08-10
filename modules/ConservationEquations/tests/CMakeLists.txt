set(CONSERVATION_EQUATIONS_TESTS
  #modules/ConservationEquations/tests/CompressibleEulerTests.cpp

  modules/ConservationEquations/tests/BoltzmannTests.cpp

  PARENT_SCOPE
)

set(PARALLEL_CONSERVATION_EQUATIONS_TESTS
  modules/ConservationEquations/tests/ParallelBoltzmannTests.cpp

  PARENT_SCOPE
)
