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

extern int rmax;
extern int p2max;
extern int pmax;
extern int vmax;

extern vector<vector<complex<double> > > Integral_LAMBDA_1;
extern vector<vector<complex<double> > > Integral_R_N;
extern vector<vector<complex<double> > > Integral_ROE_1;
extern vector<double> Integral_ROE_N;

double Partition_Function(int Extension){

	double Sum_PF = 0;
#pragma omp parallel for reduction(+:Sum_PF)
    for(int r=-rmax;r<=rmax;r++){
       for(int p=0;p<pmax;p++){
          for(int v=0;v<vmax;v++){
             Sum_PF = Sum_PF
               		+pow(LAMBDA_R(r),(N-1))
               		*pow(LAMBDA_ROE(p),(N-1))
               		*pow(LAMBDA_LAMBDA(v),(N-1))
               		*real(
               				polar((double)1,(r*Delta*Extension*2*SquareRoot(3)*PI)/L)
               				*Integral_R_N[r+rmax][v]
               				*Integral_LAMBDA_1[r+rmax][v]
               				*Integral_ROE_1[r+rmax][p]
               				*Integral_ROE_N[p]
               		);
          }
       }
    }
    printf("***Calculation for extension %3i complete.***\n",Extension);
    return Sum_PF;

}

double dPartition_Function(int Extension){

    double Sum_PF = 0;
    #pragma omp parallel for reduction(+:Sum_PF)
    for(int s=-rmax;s<=rmax;s++){
       for(int p=0;p<pmax;p++){
             for(int v=0;v<vmax;v++){
            	double s_varaible = s;
                Sum_PF = Sum_PF
                		+pow(LAMBDA_R(s),(N-1))
                		*pow(LAMBDA_ROE(p),(N-1))
                		*pow(LAMBDA_LAMBDA(v),(N-1))
                		*real(
                				polar((double)1,(s*Delta*2*SquareRoot(3)*Extension*PI)/L)
                				*Integral_R_N[s+rmax][v]
                				*Integral_LAMBDA_1[s+rmax][v]
                				*Integral_ROE_1[s+rmax][p]
                				*Integral_ROE_N[p]
                				*polar((double)1,PI/2)*(2*SquareRoot(3)*s_varaible*PI)/L
                		);
             }


       }
    }
    printf("***Calculating differential for extension %3i complete.***\n",Extension);
    return Sum_PF;
}
