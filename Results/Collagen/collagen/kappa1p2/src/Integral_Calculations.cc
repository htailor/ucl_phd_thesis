#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include "Functions.hh"
#include "Potential.hh"
#include "Matrix.hh"
#include "Vector.hh"
#include "Eigenfunctions.hh"
#include "Transfer_Matrix_Definitions.hh"
#include "Eigenvalues.hh"
#include "Integral_Calculations.hh"

#define PI 3.1415926535897932384626433832795

using namespace std;

extern double L;
extern int m;
extern double Delta;
extern double kappa;

double Integrate_LAMBDA_ALPHA(int _LAMBDA, int _LAMBDA2)
{
	double _LAMBDA_ALPHA = 0;
	// i sums over LAMBDA_alpha integral (MXD)
	for(int i=-m;i<=m;i++){
		_LAMBDA_ALPHA = _LAMBDA_ALPHA
					+ Delta
					*PSI_LAMBDA(_LAMBDA,i*Delta)
					*i*Delta
					*PSI_LAMBDA(_LAMBDA2,i*Delta);
	}
	return _LAMBDA_ALPHA;
}

double Integrate_ROE_ALPHA(int _ROE, int _ROE2)
{
	double _ROE_ALPHA = 0;
	// i sums over ROE_alpha integral (MXD)
	for(int i=-m;i<=m;i++){
		_ROE_ALPHA = _ROE_ALPHA
			  + Delta
			  *PSI_ROE(_ROE,Delta*i)
		      	  *i*Delta
			  *PSI_ROE(_ROE2,Delta*i);
	}
	return _ROE_ALPHA;
}

complex<double> Integrate_R_ALPHA(int _R, int _R2)
{
  complex<double> Integrate_R_ALPHA (0,0);
  // i sums over R_alpha integral
  for(int i=-m;i<=m;i++){
	 Integrate_R_ALPHA = Integrate_R_ALPHA
				+Delta
				*PSI_R(_R,((double)i*Delta))
				*(double)i*Delta
				*conj(PSI_R(_R2,((double)i*Delta)));
  }
  return Integrate_R_ALPHA;
}


complex<double> Integrate_LAMBDA_1(int _R, int _LAMBDA)
{
	complex<double> _Integrate_LAMBDA_1 = 0;
	// i sums over LAMBDA_1 integral
	for(int i=-m;i<=m;i++){
		_Integrate_LAMBDA_1 = _Integrate_LAMBDA_1 + Delta*polar((double)1,(_R*Delta*i*2*PI)/(SquareRoot(2.)*L))*exp(-3*Potential(Delta*i)/kappa)*PSI_LAMBDA(_LAMBDA,i*Delta);
	}
	return _Integrate_LAMBDA_1;
}

complex<double> Integrate_ROE_1(int _R, int _ROE)
{
	complex<double> _Integrate_ROE_1 = 0;
	// i sums over ROE_1 integral
	for(int i=-m;i<=m;i++){
		_Integrate_ROE_1 = _Integrate_ROE_1 + Delta*PSI_R(_R,Delta*i*SquareRoot(1.5))*exp(-3*Potential(Delta*i)/kappa)*PSI_ROE(_ROE,i*Delta);
	}
	return _Integrate_ROE_1;
}

complex<double> Integrate_R_N(int _R, int _LAMBDA)
{
	complex<double> _Integrate_R_N = 0;
	// i sums over R'_N integral
	for(int i=-m;i<=m;i++){
		_Integrate_R_N = _Integrate_R_N	+ Delta*PSI_R(_R,SquareRoot(2.)*i*Delta)*exp(-3*Potential(i*Delta)/kappa)*PSI_LAMBDA(_LAMBDA,i*Delta);
	}
	return _Integrate_R_N;
}

double Integrate_ROE_N(int _ROE)
{
	double _ROE_N = 0;
	// i sums over ROE'_N integral
	for(int i=-m;i<=m;i++){
		_ROE_N = _ROE_N + Delta*PSI_ROE(_ROE,Delta*i)*exp(-(3*Potential(Delta*i))/kappa);
	}
	return _ROE_N;
}
