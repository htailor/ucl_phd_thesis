#ifndef INCLUDE_INTEGRAL_CALCULATIONS_HH
#define INCLUDE_INTEGRAL_CALCULATIONS_HH

#include <complex>

using namespace std;

double Integrate_LAMBDA_ALPHA(int _t, int _t2);
double Integrate_ROE_ALPHA(int _t, int _t2);
complex<double> Integrate_R_ALPHA(int _s, int _s2);

complex<double> Integrate_R_N_ROE_N(int _s, int _t);
//complex<double> CIntegrate_R_N_ROE_N(int _s, int _t);
complex<double> Integrate_ROE_1_LAMBDA_1(int _s, int _t);

#endif
