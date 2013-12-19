#include "Potential.hh"
#include "Functions.hh"

using namespace std;

extern double kappa;
extern double sigma;

extern double kappa_sigma_r;

extern double e0;
extern double eta_d;
extern double beta;

double Potential(double x){

    return ((6*beta*e0)/Squared(eta_d))*Squared(x);

}

