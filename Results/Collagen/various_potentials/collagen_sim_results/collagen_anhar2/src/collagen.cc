#include "Matrix.hh"
#include "Vector.hh"
#include "Potential.hh"
#include "Matrix_Functions.hh"
#include "Calculations.hh"
#include "Integral_Calculations.hh"
#include "MXD_Calculations.hh"
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

// Global Parameters//

int umax;
int umin=0;
int N=0; 
int m;
int MATRIX_SIZE;

int nIntact;
int nBroken;

int smax=0;
int tmax=0;

double L;
double eta_b=0;
double kappa;
double sigma;

double pfConstant;

double beta;
double TOLERANCE;
double Delta;

Matrix TR;
Matrix TR_evec;
Vector TR_eval;

Matrix TROE_LAMBDA;
Matrix TROE_LAMBDA_evec;
Vector TROE_LAMBDA_eval;

Matrix PSI_R_evec;
Vector LAMBDA_R_eval;

Vector LAMBDA_ROE_LAMBDA_eval;

vector<vector<double> > MXDIntegral_LAMBDA_ALPHA;
vector<vector<double> > MXDIntegral_ROE_ALPHA;
vector<vector<complex<double> > > MXDIntegral_R_ALPHA;

vector<vector<complex<double> > > Integral_ROE_1_LAMBDA_1;
vector<vector<complex<double> > > Integral_R_N_ROE_N;

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

   printf("Initialising Transfer Matrices TR, TROE_LAMBDA...\n");

   Matrix TR(MATRIX_SIZE,MATRIX_SIZE);
   Matrix TROE_LAMBDA((int)Squared(MATRIX_SIZE),(int)Squared(MATRIX_SIZE));

   #pragma omp parallel for
   for(int i=1;i<=MATRIX_SIZE;i++){
      for(int j=1;j<=MATRIX_SIZE;j++){
         for(int k=1;k<=MATRIX_SIZE;k++){
            for(int l=1;l<=MATRIX_SIZE;l++){
               TROE_LAMBDA.set(r(i,j)-1,r(k,l)-1,Squared(Delta)*T_ROE_LAMBDA(Delta*(i-m-1),Delta*(j-m-1),Delta*(k-m-1),Delta*(l-m-1)));
            }
         }
      }
   }

   printf("Complete.\n\n");


///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 2 - CALCULATE EIGENVALUES AND EIGENVECTORS OF MATRICES T,T00,T11///

///////////////////////////////////////////////////////////////////////////////////////

   /*
   STAGE 2 solves the eigensystem for matrices T,T00 and T11. The Eigenvalues are stored in
   a vector and the eigenvectors are stored in a matrix. 
   */

   printf("Calculating Eigenvalues and Eigenvectors of TR, TROE_LAMBDA...\n");

   time_t Eigenvalue_Problem_START,Eigenvalue_Problem_FINISH;
   time (&Eigenvalue_Problem_START);

   TROE_LAMBDA_eval = Eigenvalues(TROE_LAMBDA);
   TROE_LAMBDA_evec = Eigenvectors(TROE_LAMBDA);

   for(int i=0;i<Squared(MATRIX_SIZE);i++){
      for(int j=0;j<Squared(MATRIX_SIZE);j++){
         double evec_element = TROE_LAMBDA_evec.get(i,j)*pow(Delta,-1);
         TROE_LAMBDA_evec.set(i,j,evec_element);
      }
   }

   //smax = MATRIX_SIZE-TR_eval.Zero_Elements();
   smax = m;
   tmax = (int)Squared(MATRIX_SIZE)-TROE_LAMBDA_eval.Zero_Elements();

   printf("\n\n------------------------------------------------\n");
   printf("smax: %4i tmax: %4i",smax,tmax);
   printf("\n------------------------------------------------\n\n");

   time (&Eigenvalue_Problem_FINISH);
   printf("COMPLETED in %6.2f seconds.\n\n", difftime(Eigenvalue_Problem_FINISH,Eigenvalue_Problem_START));


///////////////////////////////////////////////////////////////////////////////////////

   ///STAGE 4 - PRE-CALCULATING INTEGRALS FOR THE PARTITION FUNCTION///

///////////////////////////////////////////////////////////////////////////////////////

   printf("Pre-calculating Integrals for the Partition Function...\n");

   MXDIntegral_LAMBDA_ALPHA = vector<vector<double> > (tmax, vector<double> (tmax) );
   MXDIntegral_ROE_ALPHA = vector<vector<double> > (tmax, vector<double> (tmax) );
   MXDIntegral_R_ALPHA = vector<vector<complex<double> > > (2*smax+1, vector<complex<double> > (2*smax+1) );

   Integral_ROE_1_LAMBDA_1 = vector<vector<complex<double> > > (2*smax+1, vector<complex<double> > (tmax) );
   Integral_R_N_ROE_N = vector<vector<complex<double> > > (2*smax+1, vector<complex<double> > (tmax) );


   cout << "** Calculating double integral MXDIntegral_LAMBDA_ALPHA...\n";
#pragma omp parallel for
   for(int t=0;t<tmax;t++){
	   for(int t2=0;t2<tmax;t2++){
		   MXDIntegral_LAMBDA_ALPHA[t][t2] = Integrate_LAMBDA_ALPHA(t,t2);
	   }
   }

   cout << "** Calculating double integral MXDIntegral_ROE_ALPHA...\n";
#pragma omp parallel for
   for(int t=0;t<tmax;t++){
	   for(int t2=0;t2<tmax;t2++){
		   MXDIntegral_ROE_ALPHA[t][t2] = Integrate_ROE_ALPHA(t,t2);
	   }
   }

   cout << "** Calculating double integral MXDIntegral_R_ALPHA...\n";
#pragma omp parallel for
   for(int s=-smax;s<=smax;s++){
	   for(int s2=-smax;s2<=smax;s2++){
		   MXDIntegral_R_ALPHA[s+smax][s2+smax] = Integrate_R_ALPHA(s,s2);
	   }
   }

   cout << "** Calculating double integral Integral_ROE_1_LAMBDA_1...\n";
#pragma omp parallel for
   for(int s=-smax;s<=smax;s++){
	   for(int t=0;t<tmax;t++){
		   Integral_ROE_1_LAMBDA_1[s+smax][t] = Integrate_ROE_1_LAMBDA_1(s,t);
	   }
   }

   cout << "** Calculating double integral Integral_R_N_ROE_N...\n";
#pragma omp parallel for
   for(int s=-smax;s<=smax;s++){
	   for(int t=0;t<tmax;t++){
		   Integral_R_N_ROE_N[s+smax][t] = Integrate_R_N_ROE_N(s,t);
	   }
   }
    
   printf("Complete.\n\n");

///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 4 - CALCULATE PARTITION FUNCTION///

///////////////////////////////////////////////////////////////////////////////////////

   printf("Calculating Partition Functions.....\n\n\n");

   NucleationRun(0,0);

   printf("Complete.\n\n");

///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 5 - CREATE OUTPUT FILES///

///////////////////////////////////////////////////////////////////////////////////////

   /*
   	   This section writes all the data and logs to files. The data output includes
   	   eigenvalue, eigenvectors, partition functions, free energies, Transfer matrices.
   */

   system("mv *.log ./logs 2> /dev/null");

   time (&NUCLEATION_FINISH);
   printf("---------------------------\n");
   printf("NUCLEATION RUNTIME: %6.2f\n",difftime(NUCLEATION_FINISH,NUCLEATION_START));
   printf("---------------------------\n\n");
   
   return 0;
}
