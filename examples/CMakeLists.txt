add_subdirectory(NumericalSolvers)

add_custom_target(examples DEPENDS
                    run-tests
                    run-heat-matrix-eigenvalues
                  )
