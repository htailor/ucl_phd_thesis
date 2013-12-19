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

#define PI 3.1415926535897932384626433832795

using namespace std;

extern double L;
extern int m;
extern double Delta;

extern Matrix TROE_LAMBDA_evec;

complex<double> Integrate_ROE_1_LAMBDA_1(int _s, int _t)
{
  complex<double> _Integrate_ROE_1_LAMBDA_1 = 0;
  // i sums over ROE_1 integral
  // j sums over LAMBDA_1 integral
  for(int i=-m;i<=m;i++){
	 for(int j=-m;j<=m;j++){
		_Integrate_ROE_1_LAMBDA_1 =
				_Integrate_ROE_1_LAMBDA_1
				+Squared(Delta)
				*TROE_LAMBDA_evec.get(_t,r(i+m+1,j+m+1)-1)
				*exp(-T_V(i*Delta,j*Delta)/2.)
				*PSI_R(_s,((SquareRoot(3.)*i*Delta+j*Delta)/SquareRoot(2.)));
	 }
  }
  return _Integrate_ROE_1_LAMBDA_1;
}

complex<double> Integrate_R_N_ROE_N(int _s, int _t)
{
  complex<double> _Integrate_R_N_ROE_N = 0;
  // i sums over R_N integral
  // j sums over ROE_N integral
  for(int i=-m;i<=m;i++){
	 for(int j=-m;j<=m;j++){
		 _Integrate_R_N_ROE_N =
				_Integrate_R_N_ROE_N
				+Squared(Delta)
				*PSI_R(_s,(SquareRoot(2.)*i*Delta))
				*exp(-T_V(j*Delta,i*Delta)/2.)
				*TROE_LAMBDA_evec.get(_t,r(j+m+1,i+m+1)-1);
	 }
  }
  return _Integrate_R_N_ROE_N;
}


double Integrate_LAMBDA_ALPHA(int _t, int _t2)
{
  double Integrate_LAMBDA_ALPHA = 0;
  // i sums over ROE_alpha integral
  // j sums over LAMBDA_alpha integral
  for(int i=-m;i<=m;i++){
	 for(int j=-m;j<=m;j++){
		 Integrate_LAMBDA_ALPHA = 
				Integrate_LAMBDA_ALPHA
				+Squared(Delta)
				*TROE_LAMBDA_evec.get(_t2,r(i+m+1,j+m+1)-1)
				*j*Delta //lambda
				*TROE_LAMBDA_evec.get(_t,r(i+m+1,j+m+1)-1);
	 }
  }
  return Integrate_LAMBDA_ALPHA;
}


double Integrate_ROE_ALPHA(int _t, int _t2)
{
  double Integrate_ROE_ALPHA = 0;
  // i sums over ROE_alpha integral
  // j sums over LAMBDA_alpha integral
  for(int i=-m;i<=m;i++){
	 for(int j=-m;j<=m;j++){
		 Integrate_ROE_ALPHA =
				 Integrate_ROE_ALPHA
				+Squared(Delta)
				*TROE_LAMBDA_evec.get(_t2,r(i+m+1,j+m+1)-1)
				*i*Delta //roe
				*TROE_LAMBDA_evec.get(_t,r(i+m+1,j+m+1)-1);
	 }
  }
  return Integrate_ROE_ALPHA;
}

complex<double> Integrate_R_ALPHA(int _s, int _s2)
{
  complex<double> Integrate_R_ALPHA (0,0);
  // i sums over R_alpha integral
  for(int i=-m;i<=m;i++){
	 Integrate_R_ALPHA = Integrate_R_ALPHA
				+Delta
				*PSI_R(_s,((double)i*Delta))
				*(double)i*Delta
				*conj(PSI_R(_s2,((double)i*Delta)));
  }
  return Integrate_R_ALPHA;
}

