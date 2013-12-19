#include "Eigenfunctions.hh"
#include "Functions.hh"
#include <cmath>
#include <complex>

#define PI 3.1415926535897932384626433832795

extern int m;
extern double T_ROE_LAMBDA_BETA;
extern double T_ROE_LAMBDA_MU;
extern double T_ROE_LAMBDA_C;
extern double T_ROE_LAMBDA_DELTA;
extern double T_ROE_LAMBDA_B;
extern double T_ROE_LAMBDA_K;
extern double T_ROE_LAMBDA_EV0;
extern double T_ROE_LAMBDA_C0;
extern double Delta;
extern double L;

using namespace std;

complex<double> PSI_R(int r, double eta)
{
   return pow(L,-0.5)*polar((double)1,(r*eta*2*PI)/L);
}

complex<double> g(int r, double eta)
{
    return polar((double)1,(r*eta*2*PI)/L);
}


double PSI_ROE(int roe, double eta)
{
   double nC = pow(1/(pow((double)2,roe)*Factorial(roe)),0.5);
   return 2*T_ROE_LAMBDA_C0*nC*Hermite(roe,T_ROE_LAMBDA_B*eta)*exp(-0.5*T_ROE_LAMBDA_C*Squared(eta));
}

double PSI_LAMBDA(int lambda, double eta)
{
   double nC = pow(1/(pow((double)2,lambda)*Factorial(lambda)),0.5);
   return 2*T_ROE_LAMBDA_C0*nC*Hermite(lambda,T_ROE_LAMBDA_B*eta)*exp(-0.5*T_ROE_LAMBDA_C*Squared(eta));
}
