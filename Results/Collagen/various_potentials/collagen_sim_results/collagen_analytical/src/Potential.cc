#include "Potential.hh"
#include "Functions.hh"

using namespace std;

extern double kappa;
extern double sigma;
extern double beta;

double Potential(double x){

	return 0.5*sigma*Squared(x);
  
}

