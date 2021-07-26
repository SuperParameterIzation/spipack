set(CMAKE_REQUIRED_INCLUDES ${MUQ_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
include(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <MUQ/SamplingAlgorithms/SampleCollection.h>
  int main() {
    muq::SamplingAlgorithms::SampleCollection collection;
    Eigen::VectorXd temp = Eigen::VectorXd::Random(1024);
    collection.Add(std::make_shared<muq::SamplingAlgorithms::SamplingState>(temp));
    return 0;
  }
  "
  MUQ_CODE_COMPILES)

set(MUQ_COMPILES 1)
if( NOT MUQ_CODE_COMPILES )
  set(MUQ_COMPILES 0)
endif()
