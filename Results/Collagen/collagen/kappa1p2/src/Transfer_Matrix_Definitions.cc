#include "Transfer_Matrix_Definitions.hh"
#include "Potential.hh"
#include "Functions.hh"
#include <cmath>

using namespace std;

extern double beta;
extern double kappa;

double T_R(double ea, double eb)
{
   return exp(-Squared(ea-eb));
}

double T_ROE(double na, double nb)
{
   return exp(-Squared(na-nb))*exp(-3*(Potential(na) + Potential(nb))/kappa);
}

double T_LAMBDA(double na, double nb)
{
   return exp(-Squared(na-nb))*exp(-3*(Potential(na) + Potential(nb))/kappa);
}
