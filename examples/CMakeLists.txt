add_subdirectory(NumericalSolvers)
add_subdirectory(KineticEquations)

add_custom_target(examples DEPENDS
                    run-tests
                    run-density-estimation
                    run-Kolmogorov-eigen-decomposition
                    run-Kolmogorov-inverse

                    run-elastic-collisions
                    run-inelastic-collisions
                    run-hyperelastic-collisions
                    run-linear-external-acceleration
                    run-quadratic-external-acceleration
                    run-componentwise-quadratic-external-acceleration
                    run-nonzero-source-term
                  )
