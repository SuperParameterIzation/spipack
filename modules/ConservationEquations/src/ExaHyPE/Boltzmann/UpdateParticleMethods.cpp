// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "UpdateParticleMethods.h"

#include "kernels/GaussLegendreBasis.h"

Boltzmann::UpdateParticleMethods::UpdateParticleMethods() : exahype::plotters::LimitingADERDG2UserDefined::LimitingADERDG2UserDefined(){
  // @TODO Please insert your code here.
}

double interpolate_derivative(
    const double* offsetOfPatch,
    const double* sizeOfPatch,
    const double* x,
    int           numberOfUnknowns,
    int           unknown,
    int           order,
    const double* u
) {

  double result = 0.0;

  double xRef[DIMENSIONS];
  xRef[0] =  (x[0] - offsetOfPatch[0]) / sizeOfPatch[0];
  xRef[1] =  (x[1] - offsetOfPatch[1]) / sizeOfPatch[1];
  #if DIMENSIONS==3
  xRef[2] =  (x[2] - offsetOfPatch[2]) / sizeOfPatch[2];
  #endif

  double jac[DIMENSIONS];
  jac[0] = 1.0/sizeOfPatch[0];
  jac[1] = 1.0/sizeOfPatch[1];
  #if DIMENSIONS==3
  jac[2] = 1.0/sizeOfPatch[2];
  #endif

  std::vector<double> deriv(DIMENSIONS, 0.0);

  // The code below evaluates the basis functions at the reference coordinates
  // and multiplies them with their respective coefficient.
  dfor(ii,order+1) { // Gauss-Legendre node indices
    int iGauss = peano::utils::dLinearisedWithoutLookup(ii,order + 1);
    /*result += kernels::lobatto::basisFunction[order][ii(0)](xRef[0]) *
              kernels::lobatto::basisFunction[order][ii(1)](xRef[1]) *
              #if DIMENSIONS==3
              kernels::lobatto::basisFunction[order][ii(2)](xRef[2]) *
              #endif
               u[iGauss * numberOfUnknowns + unknown];
    assertion6(std::isfinite(result), result, unknown, iGauss, numberOfUnknowns, offsetOfPatch[0], sizeOfPatch[0]);*/

    deriv[0] += kernels::legendre::basisFunctionFirstDerivative[order][ii(0)] (xRef[0]) * jac[0] *
     kernels::legendre::basisFunction[order][ii(1)] (xRef[1]) *
     #if DIMENSIONS==3
              kernels::lobatto::basisFunction[order][ii(2)](xRef[2]) *
     #endif
     u[iGauss * numberOfUnknowns + unknown];
    ;

    deriv[1] += kernels::legendre::basisFunction[order][ii(0)] (xRef[0]) *
     kernels::legendre::basisFunctionFirstDerivative[order][ii(1)] (xRef[1]) * jac[1] *
     #if DIMENSIONS==3
              kernels::lobatto::basisFunction[order][ii(2)](xRef[2]) *
     #endif
     u[iGauss * numberOfUnknowns + unknown];
    ;

    #if DIMENSIONS==3
    deriv[2] += kernels::legendre::basisFunction[order][ii(0)] (xRef[0]) *
     kernels::legendre::basisFunction[order][ii(1)] (xRef[1]) *
              kernels::lobatto::basisFunctionFirstDerivative[order][ii(2)](xRef[2]) * jac[2] *
     u[iGauss * numberOfUnknowns + unknown];
    ;
    #endif
  }

  std::cout << "dervi: " << deriv[0] << " " << deriv[1] << std::endl;
  std::cout << "expected dervi: " << x[1] << " " << x[0] << std::endl;

  return result;

}

void Boltzmann::UpdateParticleMethods::plotADERDGPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* const u,
    double timeStamp) {
  // @TODO Please insert your code here.

  std::vector<double> x(2);
  x[0] = offsetOfPatch[0]+sizeOfPatch[0]/2.0;
  x[1] = offsetOfPatch[1]+sizeOfPatch[1]/2.0;

  const double f = kernels::legendre::interpolate(offsetOfPatch.data(), sizeOfPatch.data(), x.data(), 5, 0, 3, u);

  interpolate_derivative(offsetOfPatch.data(), sizeOfPatch.data(), x.data(), 5, 0, 3, u);

  assert(DIMENSIONS==2);

  std::cout << "offsetOfPatch: " << offsetOfPatch[0] << " " << offsetOfPatch[1] << std::endl;
  std::cout << "sizeOfPatch: " << sizeOfPatch[0] << " " << sizeOfPatch[1] << std::endl;
  std::cout << "timestamp: " << timeStamp << " coord: " << x[0] << " " << x[1] << " (ADERDG) f: " << f << " expected: " << 1.0+x[0]*x[1] << std::endl;
  std::cout << std::endl;
}

void Boltzmann::UpdateParticleMethods::plotFiniteVolumesPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* const u,
    double timeStamp) {
  // @TODO Please insert your code here.
/*
  std::vector<double> x(2);
  x[0] = 0.1; x[1] = 0.1;

  const double f = kernels::legendre::interpolate(offsetOfPatch.data(), sizeOfPatch.data(), x.data(), 5, 1, 3, u);

  std::cout << "(FV) f: " << f << std::endl;*/
}


void Boltzmann::UpdateParticleMethods::startPlotting( double time) {
  // @TODO Please insert your code here.
  std::cout << "start plotting: " << time << std::endl;
}


void Boltzmann::UpdateParticleMethods::finishPlotting() {
  // @TODO Please insert your code here.
  std::cout << "finish plotting" << std::endl;
}
