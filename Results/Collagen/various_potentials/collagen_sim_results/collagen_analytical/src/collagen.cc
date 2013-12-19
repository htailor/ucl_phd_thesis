#include "Matrix.hh"
#include "Vector.hh"
#include "Potential.hh"
#include "Matrix_Functions.hh"
#include "Calculations.hh"
#include "MXD_Calculations.hh"
#include "Integral_Calculations.hh"
#include "Transfer_Matrix_Definitions.hh"
#include "Functions.hh"
#include "Eigenfunctions.hh"
#include "Eigenvalues.hh"
#include "Output.hh"
#include "Menu.hh"
#include <cmath>
#include <getopt.h>
#include <ctime>
#include <complex>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <sstream>
#include <complex>

using namespace std;

#define PI 3.1415926535897932384626433832795

// Global Parameters//

int umax=0;
int umin;
int N=0; 
int m=0;
int MATRIX_SIZE;

int nIntact;
int nBroken;

int rmax=0;
int vmax=0;
int pmax=0;
int p2max=0;

double L;
double kappa;
double sigma;
double TOLERANCE;
double Delta;
double beta;

double T_ROE_LAMBDA_BETA;
double T_ROE_LAMBDA_MU;
double T_ROE_LAMBDA_C;
double T_ROE_LAMBDA_DELTA;
double T_ROE_LAMBDA_B;
double T_ROE_LAMBDA_K;
double T_ROE_LAMBDA_EV0;
double T_ROE_LAMBDA_C0;

Matrix TR;
Matrix TR_evec;
Vector TR_eval;

Matrix TROE;
Matrix TROE_evec;
Vector TROE_eval;

Matrix TLAMBDA;
Matrix TLAMBDA_evec;
Vector TLAMBDA_eval;

Matrix PSI_R_evec;
Matrix PSI_ROE_evec;
Matrix PSI_LAMBDA_evec;

Vector LAMBDA_R_eval;
Vector LAMBDA_ROE_eval;
Vector LAMBDA_LAMBDA_eval;

vector<vector<complex<double> > > Integral_LAMBDA_1;
vector<vector<complex<double> > > Integral_R_N;
vector<vector<complex<double> > > Integral_ROE_1;

vector<vector<double > > MXDIntegral_LAMBDA_ALPHA;
vector<vector<double > > MXDIntegral_ROE_ALPHA;
vector<vector<complex<double> > > MXDIntegral_R_ALPHA;

vector<double> Integral_ROE_N;

void NucleationRun(int intact,int broken){


	Matrix FE = Matrix(umax-umin+1,2);
	Matrix PF = Matrix(umax-umin+1,2);

	Matrix dFE = Matrix(umax-umin+1,2);
	Matrix dPF = Matrix(umax-umin+1,2);
    
    	Matrix pfMXD_r = Matrix(umax-umin+1,N);
	Matrix MXD_r = Matrix(umax-umin+1,N);

	Matrix pfMXD_roe = Matrix(umax-umin+1,N);
	Matrix MXD_roe = Matrix(umax-umin+1,N);

	Matrix pfMXD_lambda = Matrix(umax-umin+1,N);
	Matrix MXD_lambda = Matrix(umax-umin+1,N);

	Matrix AverageX = Matrix(umax-umin+1,N);
	Matrix AverageY = Matrix(umax-umin+1,N);
	Matrix AverageZ = Matrix(umax-umin+1,N);

	for(int i=umin;i<=umax;i++){

		double PARTITION_FUNCTION_VALUE = Partition_Function(i);
		double dPARTITION_FUNCTION_VALUE = dPartition_Function(i);
		double dFreeEnergy = (-1/PARTITION_FUNCTION_VALUE)*dPARTITION_FUNCTION_VALUE;

		PF.set(i-umin,1,PARTITION_FUNCTION_VALUE);
		PF.set(i-umin,0,Delta*i);

		FE.set(i-umin,1,Free_Energy(PARTITION_FUNCTION_VALUE));
		FE.set(i-umin,0,Delta*i);

		dPF.set(i-umin,1,dPARTITION_FUNCTION_VALUE);
		dPF.set(i-umin,0,Delta*i);

		dFE.set(i-umin,1,dFreeEnergy);
		dFE.set(i-umin,0,Delta*i);

		for(int alpha=1;alpha<=N;alpha++){
            
			double PARTITION_FUNCTION_VALUE_mxd_r = MXD_R_ALPHA_PF(i,alpha);
			double PARTITION_FUNCTION_VALUE_mxd_lambda = MXD_LAMBDA_ALPHA_PF(i,alpha);
			double PARTITION_FUNCTION_VALUE_mxd_roe = MXD_ROE_ALPHA_PF(i,alpha);

            		double mxd_r = PARTITION_FUNCTION_VALUE_mxd_r/PARTITION_FUNCTION_VALUE;
			double mxd_lambda = PARTITION_FUNCTION_VALUE_mxd_lambda/PARTITION_FUNCTION_VALUE;
			double mxd_roe = PARTITION_FUNCTION_VALUE_mxd_roe/PARTITION_FUNCTION_VALUE;

			double AVERAGE_X_VALUE = mxd_r*SquareRoot(1./3.) + mxd_lambda*SquareRoot(1./6.) + mxd_roe*SquareRoot(1./2.);
			double AVERAGE_Y_VALUE = mxd_r*SquareRoot(1./3.) + mxd_lambda*SquareRoot(1./6.) - mxd_roe*SquareRoot(1./2.);
			double AVERAGE_Z_VALUE = mxd_r*SquareRoot(1./3.) - mxd_lambda*SquareRoot(2./3.);
            
           		pfMXD_r.set(i,alpha-1,PARTITION_FUNCTION_VALUE_mxd_r);
			MXD_r.set(i,alpha-1,mxd_r);
            
			pfMXD_lambda.set(i,alpha-1,PARTITION_FUNCTION_VALUE_mxd_lambda);
			MXD_lambda.set(i,alpha-1,mxd_lambda);

			pfMXD_roe.set(i,alpha-1,PARTITION_FUNCTION_VALUE_mxd_roe);
			MXD_roe.set(i,alpha-1,mxd_roe);

			AverageX.set(i,alpha-1,AVERAGE_X_VALUE);
			AverageY.set(i,alpha-1,AVERAGE_Y_VALUE);
			AverageZ.set(i,alpha-1,AVERAGE_Z_VALUE);
            
		}

	}

	FE.print_to_file("FreeEnergy_0_0.out");
	PF.print_to_file("PartitionFunction_0_0.out");

	dFE.print_to_file("dFreeEnergy_0_0.out");
	dPF.print_to_file("dPartitionFunction_0_0.out");
    
  	pfMXD_r.print_to_file("pfmxd_r.out");
	MXD_r.print_to_file("mxd_r.out");

	pfMXD_lambda.print_to_file("pfmxd_lambda.out");
	MXD_lambda.print_to_file("mxd_lambda.out");

	pfMXD_roe.print_to_file("pfmxd_roe.out");
	MXD_roe.print_to_file("mxd_roe.out");

	AverageX.print_to_file("Average_X.out");
	AverageY.print_to_file("Average_Y.out");
	AverageZ.print_to_file("Average_Z.out");

	system("mv *.out ./results");

	printf("*** Intact State Complete ***\n");

}




int main(int argc,char *argv[]) 
{
   omp_set_num_threads(4);
   Menu(argc,argv);    

   time_t NUCLEATION_START,NUCLEATION_FINISH;
   time (&NUCLEATION_START);
   printf("*******************************************************************************************\n");
   printf("Starting Nucleation Simulation.....\n");
   printf("*******************************************************************************************\n\n");

   printf("Creating data for base-pair potential...\n");

   Matrix PotentialData = Matrix(MATRIX_SIZE,2);

   for(int i=-m;i<=m;i++){
  	  double POTENTIAL_VALUE = Potential(i*Delta);
   	  PotentialData.set(i+m,1,POTENTIAL_VALUE);
   	  PotentialData.set(i+m,0,Delta*i);
   }

   PotentialData.print_to_file("POTENTIAL_DATA");
   system("mv POTENTIAL_DATA ./results");
   printf("Complete.\n\n");

///////////////////////////////////////////////////////////////////////////////////////

   ///STAGE 1 - CREATE TRANSFER MATRICES T,T00,T11///

///////////////////////////////////////////////////////////////////////////////////////


   /*
   Create the transfer matrices T,T00,T11 using the Transfer Matrix Functions (These functions
   contain the mathematical definitions and hence calculate the values for each element in the
   matrix
   */

   printf("Initialising Transfer Matrices TR, TROE, TLAMBDA.....");

   Matrix TR (MATRIX_SIZE,MATRIX_SIZE);     
   Matrix TROE (MATRIX_SIZE,MATRIX_SIZE);    
   Matrix TLAMBDA (MATRIX_SIZE,MATRIX_SIZE);

   for(int i=0;i<=2*m;i++){
      for(int j=0;j<=2*m;j++){
         TR.set(i,j,Delta*T_R(Delta*(i-m),Delta*(j-m)));          // TR definition is in Transfer_Matrix_Definitions.cc
         TROE.set(i,j,Delta*T_ROE(Delta*(i-m),Delta*(j-m)));      // TROE definition is in Transfer_Matrix_Definitions.cc
         TLAMBDA.set(i,j,Delta*T_LAMBDA(Delta*(i-m),Delta*(j-m)));// TLAMBDA definition is in Transfer_Matrix_Definitions.cc
      }
   }

   printf("COMPLETE.\n\n");  
   
   

///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 2 - CALCULATE EIGENVALUES AND EIGENVECTORS OF MATRICES T,T00,T11///

///////////////////////////////////////////////////////////////////////////////////////

   /*
   STAGE 2 solves the eigensystem for matrices T,T00 and T11. The Eigenvalues are stored in
   a vector and the eigenvectors are stored in a matrix. 
   */

   printf("Calculating Eigenvalues and Eigenvectors of TR, TROE, TLAMBDA.....");

   time_t Eigenvalue_Problem_START,Eigenvalue_Problem_FINISH;
   time (&Eigenvalue_Problem_START);

   TROE_eval = Eigenvalues(TROE);
   TLAMBDA_eval = Eigenvalues(TLAMBDA);
    
   //finds the maximum number of elements to sum over and ignores
   //the elements that have zero eigenvalues.

   rmax = m;//-TR_eval.Zero_Elements();
   vmax = MATRIX_SIZE-TLAMBDA_eval.Zero_Elements();
   pmax = MATRIX_SIZE-TROE_eval.Zero_Elements();

   printf("\n\n------------------------------------------------\n");
   printf("rmax: %4i vmax: %4i pmax: %4i",rmax,vmax,pmax);
   printf("\n------------------------------------------------\n\n");

   time (&Eigenvalue_Problem_FINISH);
   printf("COMPLETED in %6.2f seconds.\n\n", difftime(Eigenvalue_Problem_FINISH,Eigenvalue_Problem_START));

 
///////////////////////////////////////////////////////////////////////////////////////

   ///STAGE 4 - PRE-CALCULATING INTEGRALS FOR THE PARTITION FUNCTION///

///////////////////////////////////////////////////////////////////////////////////////

   Integral_LAMBDA_1		= vector<vector<complex<double> > > (2*rmax+1, vector<complex<double> > (vmax) );
   Integral_R_N			= vector<vector<complex<double> > > (2*rmax+1, vector<complex<double> > (vmax) );
   Integral_ROE_1		= vector<vector<complex<double> > > (2*rmax+1, vector<complex<double> > (pmax) );
   Integral_ROE_N		= vector<double> (pmax);

   MXDIntegral_LAMBDA_ALPHA	= vector<vector<double> > (vmax,vector<double> (vmax));
   MXDIntegral_ROE_ALPHA	= vector<vector<double> > (pmax,vector<double> (pmax));
   MXDIntegral_R_ALPHA		= vector<vector<complex<double> > > (2*rmax+1, vector<complex<double> > (2*rmax+1) );
   

   #pragma omp parallel for
   for(int v=0;v<vmax;v++){
      for(int v2=0;v2<vmax;v2++){
    	  MXDIntegral_LAMBDA_ALPHA[v][v2] = Integrate_LAMBDA_ALPHA(v,v2);
      }
   }

   #pragma omp parallel for
   for(int p=0;p<pmax;p++){
      for(int p2=0;p2<pmax;p2++){
    	  MXDIntegral_ROE_ALPHA[p][p2] = Integrate_ROE_ALPHA(p,p2);
      }
   }
   
   #pragma omp parallel for
   for(int r=-rmax;r<=rmax;r++){
      for(int r2=-rmax;r2<=rmax;r2++){
    	  MXDIntegral_R_ALPHA[r+rmax][r2+rmax] = Integrate_R_ALPHA(r,r2);
      }
   }

   #pragma omp parallel for
   for(int r=-rmax;r<=rmax;r++){
      for(int v=0;v<vmax;v++){
    	  Integral_LAMBDA_1[r+rmax][v] = Integrate_LAMBDA_1(r,v);
      }
   }

   #pragma omp parallel for
   for(int r=-rmax;r<=rmax;r++){
	   for(int v=0;v<vmax;v++){
		   Integral_R_N[r+rmax][v] = Integrate_R_N(r,v);
	   }
   }

   #pragma omp parallel for
   for(int r=-rmax;r<=rmax;r++){
      for(int p=0;p<pmax;p++){
    	  Integral_ROE_1[r+rmax][p] = Integrate_ROE_1(r,p);
      }
   }

   #pragma omp parallel for
   for(int p=0;p<pmax;p++){
      Integral_ROE_N[p] = Integrate_ROE_N(p);
   }

///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 4 - CALCULATE PARTITION FUNCTION///

///////////////////////////////////////////////////////////////////////////////////////

   printf("Calculating Partition Functions.....\n\n\n");

   NucleationRun(0,0);
   printf("\n");

///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 5 - CREATE OUTPUT FILES///

///////////////////////////////////////////////////////////////////////////////////////

   /*
   	   This section writes all the data and logs to files. The data output includes
   	   eigenvalue, eigenvectors, partition functions, free energies, Transfer matrices.
   */

   time (&NUCLEATION_FINISH);
   printf("---------------------------\n");
   printf("NUCLEATION RUNTIME: %6.2f\n",difftime(NUCLEATION_FINISH,NUCLEATION_START));
   printf("---------------------------\n\n");
   
   return 0;
}
