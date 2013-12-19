#include <cmath>
#include "Eigenvalues.hh"
#include <stdlib.h>

#define PI 3.1415926535897932384626433832795

extern double T_ROE_LAMBDA_EV0;
extern double T_ROE_LAMBDA_K;
extern double L;

double LAMBDA_R(int x)
{
   return pow(PI,0.5)*exp(-pow((PI*x/L),2));
}

double LAMBDA_ROE(int x)
{
   return T_ROE_LAMBDA_EV0*exp(-T_ROE_LAMBDA_K*x);
}

double LAMBDA_LAMBDA(int x)
{
   return T_ROE_LAMBDA_EV0*exp(-T_ROE_LAMBDA_K*x);
}
