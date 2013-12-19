#ifndef INCLUDE_EIGENFUNCTIONS_HH
#define INCLUDE_EIGENFUNCTIONS_HH


#include <complex>

using namespace std;

complex<double> PSI_R(int s, double eta);
complex<double> g(int s, double eta);
complex<double> dPSI_R(int ds, double eta);
double PSI_ROE(int roe, double eta);
double PSI_LAMBDA(int lambda, double eta);

#endif
