#include "Transfer_Matrix_Definitions.hh"
#include "Functions.hh"
#include <cmath>

using namespace std;

extern double beta;
extern double kappa;

double T_R(double ea, double eb)
{
   return exp(-Squared(ea-eb));
}

double T_ROE_LAMBDA(double a1, double b1, double a2, double b2)
{
   return exp(-Squared(a1-a2))*exp(-Squared(b1-b2))*exp(-(T_V(a1,b1)+T_V(a2,b2))/2);
}

double T_V(double a, double b)
{
	return beta*(
			Potential((2*a)/SquareRoot(beta*kappa))
			+Potential((a+b*SquareRoot(3))/SquareRoot(beta*kappa))
			+Potential((SquareRoot(3)*b-a)/SquareRoot(beta*kappa))
			);
}
