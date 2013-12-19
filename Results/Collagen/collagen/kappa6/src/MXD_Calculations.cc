#include <cmath>
#include <complex>
#include <vector>
#include "Functions.hh"
#include "Eigenfunctions.hh"
#include "Eigenvalues.hh"
#include "MXD_Calculations.hh"

#define PI 3.1415926535897932384626433832795

using namespace std;

extern double L;
extern int m;
extern double Delta;
extern int N;

extern int rmax;
extern int pmax;
extern int vmax;

extern vector<vector<double > > MXDIntegral_LAMBDA_ALPHA;
extern vector<vector<double > > MXDIntegral_ROE_ALPHA;
extern vector<vector<complex<double> > > MXDIntegral_R_ALPHA;

extern vector<vector<complex<double> > > Integral_LAMBDA_1;
extern vector<vector<complex<double> > > Integral_R_N;
extern vector<vector<complex<double> > > Integral_ROE_1;
extern vector<double> Integral_ROE_N;

double MXD_LAMBDA_ALPHA_PF(int Extension,int alpha){

    double Sum_PF = 0;
    #pragma omp parallel for reduction(+:Sum_PF)
    for(int r=-rmax;r<=rmax;r++){
       for(int p=0;p<pmax;p++){
          for(int v=0;v<vmax;v++){
             for(int v2=0;v2<vmax;v2++){
                Sum_PF = Sum_PF
                		+pow(LAMBDA_R(r),(N-1))
                		*pow(LAMBDA_ROE(p),(N-1))
                		*pow(LAMBDA_LAMBDA(v),(N-alpha))
                		*pow(LAMBDA_LAMBDA(v2),(alpha-1))
                		*real(
                			polar((double)1,(r*Delta*2*SquareRoot(3)*Extension*PI)/L)
                			*Integral_R_N[r+rmax][v]
                			*Integral_ROE_1[r+rmax][p]
                			*Integral_LAMBDA_1[r+rmax][v2] //GAMMA
   		                      *Integral_ROE_N[p]
   		                      *MXDIntegral_LAMBDA_ALPHA[v][v2]
                		);
             }
	  }
       }
    }
    return Sum_PF;
}

double MXD_ROE_ALPHA_PF(int Extension,int alpha){

    double Sum_PF = 0;
    #pragma omp parallel for reduction(+:Sum_PF)
    for(int r=-rmax;r<=rmax;r++){
       for(int v=0;v<vmax;v++){
          for(int p=0;p<pmax;p++){
             for(int p2=0;p2<pmax;p2++){
                Sum_PF = Sum_PF
                		+pow(LAMBDA_R(r),(N-1))
                		*pow(LAMBDA_LAMBDA(v),(N-1))
                		*pow(LAMBDA_ROE(p),(N-alpha))
                		*pow(LAMBDA_ROE(p2),(alpha-1))
                		*real(
                			polar((double)1,(r*Delta*2*SquareRoot(3)*Extension*PI)/L)
                			*Integral_R_N[r+rmax][v]
                			*Integral_ROE_1[r+rmax][p2] //GAMMA
                			*Integral_LAMBDA_1[r+rmax][v]
   		                	*Integral_ROE_N[p]
   		                	*MXDIntegral_ROE_ALPHA[p][p2]
                		);
             }
	  }
       }
    }
    return Sum_PF;
}

double MXD_R_ALPHA_PF(int Extension,int alpha){

    double Sum_PF = 0;
    #pragma omp parallel for reduction(+:Sum_PF)
    for(int r=-rmax;r<=rmax;r++){
       for(int v=0;v<vmax;v++){
          for(int p=0;p<pmax;p++){
             for(int r2=-rmax;r2<=rmax;r2++){
                Sum_PF = Sum_PF
                		+pow(LAMBDA_R(r),(alpha-1))
                		*pow(LAMBDA_LAMBDA(v),(N-1))
                		*pow(LAMBDA_ROE(p),(N-1))
                		*pow(LAMBDA_R(r2),(N-alpha))
                		*real(
                			polar((double)1,((r+r2)*SquareRoot(3)*Delta*Extension*PI)/L)
                			*Integral_R_N[r2+rmax][v] //GAMMA
                			*Integral_ROE_1[r+rmax][p] 
                			*Integral_LAMBDA_1[r+rmax][v]
   		                	*Integral_ROE_N[p]
   		                	*MXDIntegral_R_ALPHA[r+rmax][r2+rmax]
                		);
             }
	  }
       }
    }
    return Sum_PF;
}

