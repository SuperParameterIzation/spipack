add_subdirectory(NumericalSolvers)

add_custom_target(examples DEPENDS
                    run-heat-matrix-eigenvalues
                  )
