# make sure that the HDF5 library is available
set(CMAKE_REQUIRED_LIBRARIES ${SUNDIALS_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${SUNDIALS_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")

# Try to compile without any additional flags
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <idas/idas.h>

  #include <nvector/nvector_serial.h>

  #include <sundials/sundials_dense.h>
  #include <sundials/sundials_types.h>
  #include <sundials/sundials_math.h>

  int main()
  {
    N_Vector state;
    state = N_VNew_Serial(10);
    N_VDestroy(state);
    return 0;
  }


  "
  SUNDIALS_IDA_COMPILES)

# If it failed, try compiling with an additional lapack flag
if(NOT SUNDIALS_IDA_COMPILES)
  set(CMAKE_REQUIRED_LIBRARIES ${SUNDIALS_LIBRARIES} -llapack)
  set(CMAKE_REQUIRED_INCLUDES ${SUNDIALS_INCLUDE_DIR})
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <idas/idas.h>

    #include <nvector/nvector_serial.h>

    #include <sundials/sundials_dense.h>
    #include <sundials/sundials_types.h>
    #include <sundials/sundials_math.h>

    int main()
    {
      N_Vector state;
      state = N_VNew_Serial(10);
      N_VDestroy(state);
      return 0;
    }


    "
    SUNDIALS_IDA_COMPILES_WITH_LAPACK)

endif()

if((NOT SUNDIALS_IDA_COMPILES) AND (NOT SUNDIALS_IDA_COMPILES_WITH_LAPACK))
  	set(SUNDIALS_TEST_FAIL 1)
else()
	  set(SUNDIALS_TEST_FAIL 0)
endif()
