add_subdirectory(NumericalSolvers)
add_subdirectory(KineticEquations)

add_custom_target(examples DEPENDS
                    run-tests
                    run-density-estimation
                    run-Kolmogorov-eigen-decomposition
                    run-Kolmogorov-inverse
                    run-microscale-model
                  )
