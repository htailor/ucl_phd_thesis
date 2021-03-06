#include "Matrix.hh"
#include "Vector.hh"
#include "Potential.hh"
#include "Matrix_Functions.hh"
#include "Partition_Function.hh"
#include "Transfer_Matrix_Definitions.hh"
#include "Functions.hh"
#include "Eigenfunctions.hh"
#include "Eigenvalues.hh"
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

using namespace std;

// Global Parameters//

double Beta = 38.6817; // eV^-1
double eta_d = 20.; // diameter of dna AA
double e0;

int umax=0;
int umin;
int N=0; 
int m=0;
int MATRIX_SIZE;

int nIntact;
int nBroken;

int smax=0;
int tmax=0;
int pmax=0;


double L=0;
double eta_b=0;
double kappa=0;
double sigma=0;
double kappa_sigma_r=0;
double TOLERANCE;
double Delta;

double beta;
double mu;
double c;
double delta;
double b;
double k;
double ev0;
double C0;


Matrix t;
Matrix t_evec;
Vector t_eval;

Matrix t00;
Matrix t00_evec;
Vector t00_eval;

vector<vector<complex<double> > > sVt00;
vector<vector<complex<double> > > CsVt00;
vector<vector<double> > t_eta_t2;
vector<vector<complex<double> > > s1_xi_s2;

void NucleationRun(int intact, int broken){
    
	string FE_filename  = "FreeEnergy_" + intToString(intact) + "_" + intToString(broken) + ".out";
	string PF_filename  = "PartitionFunction_" + intToString(intact) + "_" + intToString(broken) + ".out";

	Matrix FE = Matrix(umax-umin+1,2);
	Matrix PF = Matrix(umax-umin+1,2);

	string dFE_filename = "dFreeEnergy_" + intToString(intact) + "_" + intToString(broken) + ".out";
	string dPF_filename = "dPartitionFunction_" + intToString(intact) + "_" + intToString(broken) + ".out";
  
	Matrix dFE = Matrix(umax-umin+1,2);
	Matrix dPF = Matrix(umax-umin+1,2);

	string ZNJ_xi_filename = "znj_xi.out";
	string ZNJ_eta_filename = "znj_eta.out";

	Matrix ZNJ_xi = Matrix(umax-umin+1,N);
	Matrix ZNJ_eta = Matrix(umax-umin+1,N);

	string MXD_xi_filename = "MeanAxialDisp_xi.out";
	string MXD_eta_filename = "MeanAxialDisp_eta.out";

	Matrix MXD_xi = Matrix(umax-umin+1,N);
	Matrix MXD_eta = Matrix(umax-umin+1,N);

	string Average_X_filename = "Average_X.out";
	string Average_Y_filename = "Average_Y.out";

	Matrix AverageX = Matrix(umax-umin+1,N);
	Matrix AverageY = Matrix(umax-umin+1,N);
	
	for(int i=umin;i<=umax;i++){
		
		string ext_status_create = "touch RUNNING_EXT_" + intToString(i) + ".status";
		string ext_status_destroy = "rm -f RUNNING_EXT_" + intToString(i) + ".status";

		system(ext_status_create.c_str());
		vector<double> PARTITION_FUNCTION_DATA(Partition_Function(i));
		double PARTITION_FUNCTION_VALUE = PARTITION_FUNCTION_DATA[0];
		double dPARTITION_FUNCTION_VALUE = PARTITION_FUNCTION_DATA[1];
		double dFreeEnergy = -(1/PARTITION_FUNCTION_VALUE)*dPARTITION_FUNCTION_VALUE;

		PF.set(i-umin,1,PARTITION_FUNCTION_VALUE);
		PF.set(i-umin,0,Delta*i);

		FE.set(i-umin,1,Free_Energy(PARTITION_FUNCTION_VALUE));
		FE.set(i-umin,0,Delta*i);
			
		dPF.set(i-umin,1,dPARTITION_FUNCTION_VALUE);
		dPF.set(i-umin,0,Delta*i);	

		dFE.set(i-umin,1,dFreeEnergy);
		dFE.set(i-umin,0,Delta*i);

    	PF.print_to_file(PF_filename.c_str());
    	dPF.print_to_file(dPF_filename.c_str());

    	FE.print_to_file(FE_filename.c_str());
		dFE.print_to_file(dFE_filename.c_str());

 // Intact State Only - Calculate MXD
	for(int j=1;j<=N;j++){ // j represents the base-pair, i refers the extension applied.

		string mxd_status_create = "touch RUNNING_MXD_" + intToString(j) + ".status";
		string mxd_status_destroy = "rm -f RUNNING_MXD_" + intToString(j) + ".status";

		system(mxd_status_create.c_str());

            double ZNJ_XI_VALUE = ZNJ_XI(i,j);
            double ZNJ_ETA_VALUE = ZNJ_ETA(i,j); // refer to BasePairCalculations.cc for definition.

            double MXD_XI_VALUE = ZNJ_XI_VALUE/PARTITION_FUNCTION_VALUE;
            double MXD_ETA_VALUE = ZNJ_ETA_VALUE/PARTITION_FUNCTION_VALUE;

            double AVERAGE_X_VALUE = 0.5*(MXD_XI_VALUE+MXD_ETA_VALUE);
            double AVERAGE_Y_VALUE = 0.5*(MXD_XI_VALUE-MXD_ETA_VALUE);

            ZNJ_xi.set(i,j-1,ZNJ_XI_VALUE);
            MXD_xi.set(i,j-1,ZNJ_XI_VALUE/PARTITION_FUNCTION_VALUE);

            ZNJ_eta.set(i,j-1,ZNJ_ETA_VALUE);
            MXD_eta.set(i,j-1,ZNJ_ETA_VALUE/PARTITION_FUNCTION_VALUE);

            AverageX.set(i,j-1,AVERAGE_X_VALUE);
            AverageY.set(i,j-1,AVERAGE_Y_VALUE);

				system(mxd_status_destroy.c_str());		

        }

			system(ext_status_destroy.c_str());
	}

	ZNJ_xi.print_to_file(ZNJ_xi_filename.c_str());
	ZNJ_eta.print_to_file(ZNJ_eta_filename.c_str());

	MXD_eta.print_to_file(MXD_eta_filename.c_str());
    MXD_xi.print_to_file(MXD_xi_filename.c_str());

    AverageX.print_to_file(Average_X_filename.c_str());
    AverageY.print_to_file(Average_Y_filename.c_str());

    printf("\n***INTACT STATE COMPLETE***\n\n");
    system("mv *.out ./results/Intact");
	

}

int main(int argc,char *argv[]) 
{
   //omp_set_num_threads(4);
   omp_set_num_threads(8);
   Menu(argc,argv);    

   system("rm -f *.status");
   system("touch SIM_STARTED.status");
   system("touch RUNNING.status");

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

   PotentialData.print_to_file("./results/Potential_Data.out");
   printf("Complete.\n\n");

///////////////////////////////////////////////////////////////////////////////////////

   ///STAGE 1 - CREATE TRANSFER MATRICES T,T00,T11///

///////////////////////////////////////////////////////////////////////////////////////

   /*
   Create the transfer matrices T,T00,T11 using the Transfer Matrix Functions (These functions
   contain the mathematical definitions and hence calculate the values for each element in the
   matrix
   */

   printf("Initialising Transfer Matrices T, T00.....");

   Matrix t (MATRIX_SIZE,MATRIX_SIZE);     
   Matrix t00 (MATRIX_SIZE,MATRIX_SIZE);

   for(int i=0;i<=2*m;i++){
      for(int j=0;j<=2*m;j++){
         t00.set(i,j,Delta*T_hat00(Delta*(i-m),Delta*(j-m)));  // T_hat00 definition is in Transfer_Matrix_Definitions.cc
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

   printf("Calculating Eigenvalues and Eigenvectors of T, T00 .....");  

   time_t Eigenvalue_Problem_START,Eigenvalue_Problem_FINISH;
   time (&Eigenvalue_Problem_START);

   t00_eval = Eigenvalues(t00);
   t00_evec = Eigenvectors(t00);

   for(int i=0;i<MATRIX_SIZE;i++){
      for(int j=0;j<MATRIX_SIZE;j++){
         double NewElement_T00 = pow(Delta,-0.5)*t00_evec.get(i,j);
         t00_evec.set(i,j,NewElement_T00);
      }
   }

   smax = m;
   tmax = MATRIX_SIZE-t00_eval.Zero_Elements();
   
   printf("\n\n------------------------------------------------\n");
   printf("smax: %4i tmax: %4i",smax,tmax);
   printf("\n------------------------------------------------\n\n");

   time (&Eigenvalue_Problem_FINISH);
   printf("COMPLETED in %6.2f seconds.\n\n", difftime(Eigenvalue_Problem_FINISH,Eigenvalue_Problem_START));

///////////////////////////////////////////////////////////////////////////////////////

   ///STAGE 3 - CALCULATE DOUBLE INTEGRALS FOR T10 and T01 Matrices///

///////////////////////////////////////////////////////////////////////////////////////

   sVt00 = vector<vector<complex<double> > > (2*smax+1, vector<complex<double> > (tmax) );
   CsVt00 = vector<vector<complex<double> > > (2*smax+1, vector<complex<double> > (tmax) );
   t_eta_t2 = vector<vector<double > > (tmax, vector<double> (tmax) );
   s1_xi_s2 = vector<vector<complex<double> > > (2*smax+1, vector<complex<double> > (2*smax+1) );

   time_t DI_START, DI_FINISH;
   time (&DI_START);

   printf("Pre-calculating Standard Integrals.....");

/// Pre Calculating the Double Intergrals ///

   #pragma omp parallel for
   for(int t=0;t<tmax;t++){
      for(int t2=0;t2<tmax;t2++){
    	 double Integral_1=0;
    	 for(int i=-m;i<=m;i++){
    	    Integral_1 = Integral_1 + Delta*t00_evec.get(t,i+m)*Delta*i*t00_evec.get(t2,i+m);
    	 }
    	 t_eta_t2[t][t2] = Integral_1;
      }
   }

   #pragma omp parallel for
   for(int s1=-smax;s1<=smax;s1++){
      for(int s2=-smax;s2<=smax;s2++){
    	 complex<double> Integral_2 = (0,0);
    	 for(int i=-m;i<=m;i++){
    	    Integral_2 = Integral_2 + Delta*PSI_S(s1,i*Delta)*(Delta*i)*PSI_S(s2,i*Delta);
    	 }
    	 s1_xi_s2[s1+smax][s2+smax] = Integral_2;
      }
   }

   printf("COMPLETE.\n\nPre-calculating Complex Integrals.....COMPLETE\n\n");

   #pragma omp parallel for
   for(int s=-smax;s<=smax;s++){
      for(int t=0;t<tmax;t++){
    	 sVt00[s+smax][t] = EXP_sVt00(s,t);
      }
   }

   #pragma omp parallel for
   for(int s=-smax;s<=smax;s++){
      for(int t=0;t<tmax;t++){
    	 CsVt00[s+smax][t] = EXP_CsVt00(s,t);
      }
   }

   time (&DI_FINISH);
   printf("COMPLETED in %6.2f seconds.\n\n",difftime(DI_FINISH,DI_START));

///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 4 - CALCULATE PARTITION FUNCTION///

///////////////////////////////////////////////////////////////////////////////////////

   printf("Calculating Partition Functions.....\n\n");   

   time_t PF_TOTAL_START, PF_TOTAL_FINISH;
   time (&PF_TOTAL_START);

   //This calculates the Intact State
   system("touch RUNNING_INTACT_STATE.status");
   NucleationRun(0,0);
   system("rm -f RUNNING_INTACT_STATE.status");
   system("touch INTACT_STATE_COMPLETED.status");
   

   time(&PF_TOTAL_FINISH);
   printf("COMPLETED in %6.2f seconds\n\n",difftime(PF_TOTAL_FINISH,PF_TOTAL_START));


///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 5 - CREATE OUTPUT FILES///

///////////////////////////////////////////////////////////////////////////////////////

  
   //This section writes all the data and logs to files. The data output includes
   //eigenvalue, eigenvectors, partition functions, free energies, Transfer matrices.
   
   cout << "Writing output data to files....\n";
   OutputDataToFiles();

   time (&NUCLEATION_FINISH);
   printf("---------------------------\n");
   printf("NUCLEATION RUNTIME: %6.2f\n",difftime(NUCLEATION_FINISH,NUCLEATION_START));
   printf("---------------------------\n\n");

   system("rm -f RUNNING.status");
   system("touch SIM_FINISHED.status");
   
   system("find . -name 'core.*' | xargs rm -f");
    
   return 0;
}
