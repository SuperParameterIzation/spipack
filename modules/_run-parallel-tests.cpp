// include the google testing header
#include <gtest/gtest.h>
#include <mpi.h>

using namespace std;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  int ierr;
  ierr = MPI_Init(NULL, NULL);

  int res = RUN_ALL_TESTS();

  MPI_Finalize();

  return res;
}
