#include "Potential.hh"
#include "Functions.hh"

using namespace std;

extern double kappa;
extern double sigma;
extern double beta;


double Potential(double x){

	double etab = 1.75;
	double HardStiffConst = 8;
	double SoftStiffConst = (1/16);

	double PotentialValueConst = SoftStiffConst*(0.5*sigma*Squared(etab)) - HardStiffConst*(0.5*sigma*Squared(etab));
	double PotentialValue;

	if(abs(x) <= etab){
	   PotentialValue = SoftStiffConst*(0.5*sigma*Squared(x));
	}
	else{
	   PotentialValue = HardStiffConst*(0.5*sigma*Squared(x)) + PotentialValueConst;
	}
	return PotentialValue;
}
