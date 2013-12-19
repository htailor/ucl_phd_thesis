#include "Eigenfunctions.hh"
#include "Functions.hh"
#include <cmath>
#include <complex>

#define PI 3.1415926535897932384626433832795

extern double Delta;
extern double L;

using namespace std;

complex<double> PSI_R(int r, double eta)
{
   return pow(L,-0.5)*polar((double)1,(r*eta*2.*PI)/L);
}

complex<double> g(int r, double eta)
{
    return polar((double)1,((double)r*eta*2.*PI)/L);
}

