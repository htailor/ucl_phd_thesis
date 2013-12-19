#include "Transfer_Matrix_Definitions.hh"
#include "Functions.hh"
#include <cmath>

extern double eta_b;

using namespace std;

double T(double ea, double eb)
{
   return exp(-Squared(ea-eb));
}

double T_hat(double na, double nb)
{
   return exp(-Squared(na-nb))*exp(-(Potential(na) + Potential(nb))/2); 
}

double T_hat00(double na, double nb)
{
   return T_hat(na,nb);
}

