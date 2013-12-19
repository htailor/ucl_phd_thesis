#ifndef INCLUDE_INTEGRAL_CALCULATIONS_HH
#define INCLUDE_INTEGRAL_CALCULATIONS_HH

#include <complex>

using namespace std;

double Integrate_LAMBDA_ALPHA(int _LAMBDA, int _LAMBDA2);
double Integrate_ROE_ALPHA(int _ROE, int _ROE2);
complex<double> Integrate_R_ALPHA(int _R, int _R2);

double Integrate_ROE_N(int p2);
complex<double> Integrate_LAMBDA_1(int _R, int _LAMBDA);
complex<double> Integrate_ROE_1(int _R, int _ROE);
complex<double> Integrate_R_N(int _R, int _LAMBDA);


#endif
