set(CMAKE_REQUIRED_INCLUDES ${PARCER_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${PARCER_LIBRARIES} ${MPI_LIBRARIES})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
include(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <mpi.h>
  #include <vector>
  #include <memory>
  #include <cereal/archives/binary.hpp>
  #include <cereal/types/vector.hpp>
  #include <cereal/types/memory.hpp>
  #include <parcer/Communicator.h>
  int main() {
  }
  "
  PARCER_CODE_COMPILES)

set(PARCER_COMPILES 1)
if( NOT PARCER_CODE_COMPILES )
  set(PARCER_COMPILES 0)
endif()
