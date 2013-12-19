#include "Potential.hh"
#include "Functions.hh"

using namespace std;

extern double kappa;
extern double sigma;
extern double beta;


double Potential(double x){

	double etab = 0.5;
	double PotentialValueConst = 16*(0.5*sigma*Squared(etab)) - 0.0625*(0.5*sigma*Squared(etab));
	double PotentialValue;

	if(abs(x) <= etab){
	   PotentialValue = 16*(0.5*sigma*Squared(x));
	}
	else{
	   PotentialValue = 0.0625*(0.5*sigma*Squared(x)) + PotentialValueConst;
	}
	return PotentialValue;

}

