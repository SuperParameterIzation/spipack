#include <gtest/gtest.h>

#include "spipack/Tools/SortVector.hpp"

using namespace spi::Tools;

TEST(SortVectorTests, Descending) {
  const std::size_t dim = 100;

  // create a vector of random integeters
  std::vector<std::size_t> rands(dim);
  for( std::size_t i=0; i<dim; ++i ) { rands[i] = rand(); }

  // get the indices of the elements in descending order
  const std::vector<std::size_t> inds = SortVector<std::vector<std::size_t> >::Descending(rands);

  for( std::size_t i=0; i<dim-1; ++i ) { EXPECT_TRUE(rands[inds[i]]>=rands[inds[i+1]]); }
}

TEST(SortVectorTests, Ascending) {
  const std::size_t dim = 100;

  // create a vector of random integeters
  std::vector<std::size_t> rands(dim);
  for( std::size_t i=0; i<dim; ++i ) { rands[i] = rand(); }

  // get the indices of the elements in descending order
  const std::vector<std::size_t> inds = SortVector<std::vector<std::size_t> >::Ascending(rands);

  //for( std::size_t i=0; i<dim; ++i ) { std::cout << rands[i] << std::endl; }
  for( std::size_t i=0; i<dim-1; ++i ) { EXPECT_TRUE(rands[inds[i]]<=rands[inds[i+1]]); }
}
