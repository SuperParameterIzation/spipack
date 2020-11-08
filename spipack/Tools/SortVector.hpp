#ifndef SORTVECTOR_HPP_
#define SORTVECTOR_HPP_

#include <numeric>

namespace spi {
namespace Tools {

/// Get the indices of a vector that order the vector in either ascending or descending order
/**
We can sort any time of vector with the accessor operator[].
*/
template<typename VEC>
class SortVector {
public:
  /// Sort the vector from smallest to largest
  /**
  @param[in] in The vector we wish to sort
  \return The indices that order the vector from smallest to largest
  */
  inline static std::vector<std::size_t> Ascending(VEC const& in) {
    // initialize original index locations
    std::vector<std::size_t> idx(in.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values
    std::sort(idx.begin(), idx.end(),
       [&in](std::size_t i1, std::size_t i2) { return in[i1]<in[i2]; });

    return idx;
  }

  /// Sort the vector from largest to smallest
  /**
  @param[in] in The vector we wish to sort
  \return The indices that order the vector from largest to smallest
  */
  inline static std::vector<std::size_t> Descending(VEC const& in) {
    // initialize original index locations
    std::vector<std::size_t> idx(in.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values
    std::sort(idx.begin(), idx.end(),
       [&in](std::size_t i1, std::size_t i2) { return in[i1]>in[i2]; });

    return idx;
  }
};

} // namespace tools
} // namespace spi

#endif
