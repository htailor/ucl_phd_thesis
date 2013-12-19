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

extern double pfConstant;

extern int smax;
extern int tmax;

extern Matrix TROE_LAMBDA_evec;
extern Vector TROE_LAMBDA_eval;

extern vector<vector<complex<double> > > Integral_ROE_1_LAMBDA_1;
extern vector<vector<complex<double> > > Integral_R_N_ROE_N;

extern vector<vector<double > > MXDIntegral_LAMBDA_ALPHA;
extern vector<vector<double > > MXDIntegral_ROE_ALPHA;

double Partition_Function(int Extension){
	double Sum_PF = 0;
	#pragma omp parallel for reduction(+:Sum_PF)
	for(int s=-smax;s<=smax;s++){
		for(int t=0;t<tmax;t++){
			// i sums over R_N integral
			// j sums over ROE_N integral
			Sum_PF = Sum_PF
					+pow(LAMBDA_R(s),(N-1))
					*pow(TROE_LAMBDA_eval.get(t),(N-1))
					*real(
							g(s,Delta*Extension*SquareRoot(3.))
							*Integral_ROE_1_LAMBDA_1[s+smax][t]
							*Integral_R_N_ROE_N[s+smax][t]
					);
		}
	}
	//printf("***Calculation for extension %3i complete.***\n",Extension);
	return Sum_PF;
}



double dPartition_Function(int Extension){

	double Sum_PF = 0;
	#pragma omp parallel for reduction(+:Sum_PF)
	for(int s=-smax;s<=smax;s++){
		double s_variable = s;
		for(int t=0;t<tmax;t++){
			// i sums over R_N integral
			// j sums over ROE_N integral
			Sum_PF = Sum_PF
					+pow(LAMBDA_R(s),(N-1))
					*pow(TROE_LAMBDA_eval.get(t),(N-1))
					*real(
							g(s,Delta*Extension*SquareRoot(3.))
							*Integral_ROE_1_LAMBDA_1[s+smax][t]
							*Integral_R_N_ROE_N[s+smax][t]
							*polar((double)1,PI/2)*(s_variable*2*SquareRoot(3.)*PI)/L
					);
		}
	}
	//printf("***Calculating differential for extension %3i complete.***\n",Extension);
	return Sum_PF;
}
