#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include "Functions.hh"
#include "Potential.hh"
#include "Matrix.hh"
#include "Vector.hh"
#include "Eigenfunctions.hh"
#include "Transfer_Matrix_Definitions.hh"
#include "Eigenvalues.hh"

#define PI 3.1415926535897932384626433832795

using namespace std;

extern double L;
extern int m;
extern double Delta;
extern int N;

extern int smax;
extern int tmax;

extern double pfConstant;

extern Matrix TROE_LAMBDA_evec;
extern Vector TROE_LAMBDA_eval;

extern vector<vector<complex<double> > > Integral_ROE_1_LAMBDA_1;
extern vector<vector<complex<double> > > Integral_R_N_ROE_N;
extern vector<vector<complex<double> > > CIntegral_R_N_ROE_N;

extern vector<vector<double > > MXDIntegral_LAMBDA_ALPHA;
extern vector<vector<double > > MXDIntegral_ROE_ALPHA;
extern vector<vector<complex<double> > > MXDIntegral_R_ALPHA;

double MXD_LAMBDA_ALPHA_PF(int Extension,int alpha){

	double Sum_PF = 0;
	#pragma omp parallel for reduction(+:Sum_PF)
	for(int s=-smax;s<=smax;s++){
		for(int t=0;t<tmax;t++){
			for(int t2=0;t2<tmax;t2++){
			// i sums over R_N integral
			// j sums over ROE_N integral
			Sum_PF = Sum_PF	+ pow(LAMBDA_R(s),(N-1))
					*pow(TROE_LAMBDA_eval.get(t),(N-alpha))
					*pow(TROE_LAMBDA_eval.get(t2),(alpha-1))
					*real(
						g(s,Delta*(double)Extension*SquareRoot(3.))
						*Integral_R_N_ROE_N[s+smax][t]
						*Integral_ROE_1_LAMBDA_1[s+smax][t2] //GAMMA
						*MXDIntegral_LAMBDA_ALPHA[t][t2]
					);
			}
		}
	}
	//printf("***Calculation for extension %3i complete.***\n",Extension);
	return Sum_PF;
}


double MXD_ROE_ALPHA_PF(int Extension,int alpha){

	double Sum_PF = 0;
	#pragma omp parallel for reduction(+:Sum_PF)
	for(int s=-smax;s<=smax;s++){
		for(int t=0;t<tmax;t++){
			for(int t2=0;t2<tmax;t2++){
			// i sums over R_N integral
			// j sums over ROE_N integral
			Sum_PF = Sum_PF
					+ pow(LAMBDA_R(s),(N-1))
					*pow(TROE_LAMBDA_eval.get(t),(N-alpha))
					*pow(TROE_LAMBDA_eval.get(t2),(alpha-1))
					*real(
						g(s,Delta*(double)Extension*SquareRoot(3.))
						*Integral_R_N_ROE_N[s+smax][t]
						*Integral_ROE_1_LAMBDA_1[s+smax][t2] //GAMMA
						*MXDIntegral_ROE_ALPHA[t][t2]	
					);
			}
		}
	}
	//printf("***Calculation for extension %3i complete.***\n",Extension);
	return Sum_PF;
}

double MXD_R_ALPHA_PF(int Extension,int alpha){

	double Sum_PF = 0;
	#pragma omp parallel for reduction(+:Sum_PF)
	for(int s=-smax;s<=smax;s++){
		for(int t=0;t<tmax;t++){
			for(int s2=-smax;s2<=smax;s2++){
			Sum_PF = Sum_PF + pow(LAMBDA_R(s),(alpha-1))
					*pow(LAMBDA_R(s2),(N-alpha))
					*pow(TROE_LAMBDA_eval.get(t),(N-1))
					*real(
						g((s+s2),(Delta*(double)Extension*SquareRoot(3.))/2.)
						*Integral_ROE_1_LAMBDA_1[s+smax][t]
						*Integral_R_N_ROE_N[s2+smax][t] //THETA
						*MXDIntegral_R_ALPHA[s+smax][s2+smax]
					);
			}
		}
	}
	//printf("***Calculation for extension %3i complete.***\n",Extension);
	return Sum_PF;
}
