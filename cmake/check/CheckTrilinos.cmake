set(CMAKE_REQUIRED_INCLUDES ${Trilinos_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${Trilinos_LIBRARIES})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
include(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <Sacado.hpp>
  template <typename ScalarT>
  ScalarT func(const ScalarT& a, const ScalarT& b, const ScalarT& c) {
  ScalarT r = c*std::log(b+1.)/std::sin(a);
  return r;
  }
  int main() {
  double a = std::atan(1.0);
  double b = 2.0;
  double c = 3.0;
  int num_deriv = 2;
  Sacado::Fad::DFad<double> afad(num_deriv, 0, a);
  Sacado::Fad::DFad<double> bfad(num_deriv, 1, b);
  Sacado::Fad::DFad<double> cfad(c);
  Sacado::Fad::DFad<double> rfad;
  double r = func(a, b, c);
  rfad = func(afad, bfad, cfad);
  double r_ad = rfad.val();
  double drda_ad = rfad.dx(0);
  double drdb_ad = rfad.dx(1);
  return 0;
  }
  "
  SACADO_CODE_COMPILES)

set(Trilinos_COMPILES 1)
if( NOT SACADO_CODE_COMPILES )
  set(Trilinos_COMPILES 0)
endif()
