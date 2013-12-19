#include "Potential.hh"
#include "Functions.hh"
#include <cmath>
#include <iostream>

using namespace std;

extern double eta_b;
extern double kappa_sigma_r;

double Potential(double _n){
  
  double PotentialValue;
  if(abs(_n) < eta_b){
      PotentialValue = 2*Squared(_n)/kappa_sigma_r;
  }
  else{
      PotentialValue = 2*Squared(eta_b)/kappa_sigma_r;
  }
  return PotentialValue;
}

