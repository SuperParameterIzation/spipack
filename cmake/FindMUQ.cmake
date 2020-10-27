find_package(MUQ REQUIRED HINTS ${SPIPACK_MUQ_DIR})

if( MUQ_FOUND )
  list(APPEND SPIPACK_EXTERNAL_INCLUDE_DIRS
    ${MUQ_INCLUDE_DIRS}
  )

  list(APPEND SPIPACK_EXTERNAL_LIBRARIES
    ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES}
  )
else()
  message(STATUS "NO MUQ FOUND")
endif()
