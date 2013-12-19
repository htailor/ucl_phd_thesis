#include "Matrix.hh"
#include "Vector.hh"
#include "Table.hh"
#include "Potential.hh"
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
#include <stdlib.h>
#include <getopt.h>
#include <ctime>
#include <complex>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <fstream>
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

double kappa_sigma_r;

double beta;
double eta_d; // diameter of dna AA
double e0;

double T_ROE_LAMBDA_BETA;
double T_ROE_LAMBDA_MU;
double T_ROE_LAMBDA_C;
double T_ROE_LAMBDA_DELTA;
double T_ROE_LAMBDA_B;
double T_ROE_LAMBDA_K;
double T_ROE_LAMBDA_EV0;
double T_ROE_LAMBDA_C0;



Table PSI_ROE_evec;
Table PSI_LAMBDA_evec;

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
    
	Table FE = Table(2,umax-umin+1);
	Table PF = Table(2,umax-umin+1);
    
	Table dFE = Table(2,umax-umin+1);
	Table dPF = Table(2,umax-umin+1);
    
    Table pfMXD_r = Table(N,umax-umin+1);
	Table MXD_r = Table(N,umax-umin+1);
    
	Table pfMXD_roe = Table(N,umax-umin+1);
	Table MXD_roe = Table(N,umax-umin+1);
    
	Table pfMXD_lambda = Table(N,umax-umin+1);
	Table MXD_lambda = Table(N,umax-umin+1);
    
	Table AverageX = Table(N,umax-umin+1);
	Table AverageY = Table(N,umax-umin+1);
	Table AverageZ = Table(N,umax-umin+1);
    
	for(int i=umin;i<=umax;i++){
        
		double PARTITION_FUNCTION_VALUE = Partition_Function(i);
		double dPARTITION_FUNCTION_VALUE = dPartition_Function(i);
		double dFreeEnergy = (-1/PARTITION_FUNCTION_VALUE)*dPARTITION_FUNCTION_VALUE;
        
        double eta = Delta*i;
        
		PF.set(0,i-umin,eta);
		PF.set(1,i-umin,PARTITION_FUNCTION_VALUE);
        
		FE.set(0,i-umin,eta);
		FE.set(1,i-umin,Free_Energy(PARTITION_FUNCTION_VALUE));
        
		dPF.set(0,i-umin,eta);
		dPF.set(1,i-umin,dPARTITION_FUNCTION_VALUE);
        
		dFE.set(0,i-umin,eta);
		dFE.set(1,i-umin,dFreeEnergy);
        
        
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
            
            pfMXD_r.set(alpha-1,i,PARTITION_FUNCTION_VALUE_mxd_r);
			MXD_r.set(alpha-1,i,mxd_r);
            
			pfMXD_lambda.set(alpha-1,i,PARTITION_FUNCTION_VALUE_mxd_lambda);
			MXD_lambda.set(alpha-1,i,mxd_lambda);
            
			pfMXD_roe.set(alpha-1,i,PARTITION_FUNCTION_VALUE_mxd_roe);
			MXD_roe.set(alpha-1,i,mxd_roe);
            
			AverageX.set(alpha-1,i,AVERAGE_X_VALUE);
			AverageY.set(alpha-1,i,AVERAGE_Y_VALUE);
			AverageZ.set(alpha-1,i,AVERAGE_Z_VALUE);
            
		}
        printf("***Calculation for extension %3i complete.***\n",i);
        
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
    
	printf("*** INTACT STATE COMPLETE ***\n");
    
    
}




int main(int argc,char *argv[])
{
   
   Menu(argc,argv);

    system("rm -f *.status");
    system("touch SIM_STARTED.status");
    system("touch RUNNING.status");


    char * num_threads;
    num_threads = getenv ("OMP_NUM_THREADS");
    if (num_threads!=NULL){
        printf ("*Number of processing threads: %s\n",num_threads);
    }

    int eigen_limit_tmax = (int)MATRIX_SIZE;
    char * eigen_limit_tmax_string;
    eigen_limit_tmax_string = getenv("EIGEN_LIMIT_TMAX");

    if (eigen_limit_tmax_string != NULL){
        eigen_limit_tmax = atoi(eigen_limit_tmax_string);
        if(eigen_limit_tmax==0){
           eigen_limit_tmax = (int)MATRIX_SIZE;
        }
    }
    printf ("*Setting EIGEN_LIMIT_TMAX = %i\n",eigen_limit_tmax);


    int eigen_limit_smax = m;
    char * eigen_limit_smax_string;
    eigen_limit_smax_string = getenv("EIGEN_LIMIT_SMAX");

    if (eigen_limit_smax_string != NULL){
        eigen_limit_smax = atoi(eigen_limit_smax_string);
        if(eigen_limit_smax==0){
           eigen_limit_smax = m;
        }
    }
    printf ("*Setting EIGEN_LIMIT_SMAX = %i\n",eigen_limit_smax);


   time_t NUCLEATION_START,NUCLEATION_FINISH;
   time (&NUCLEATION_START);
   printf("*******************************************************************************************\n");
   printf("Starting Nucleation Simulation.....\n");
   printf("*******************************************************************************************\n\n");

   printf("Creating data for base-pair potential...\n");

    Table PotentialData = Table(2,MATRIX_SIZE);
    
    for(int i=-m;i<=m;i++){
        double POTENTIAL_VALUE = Potential(i*Delta);
        PotentialData.set(0,i+m,Delta*i);
        PotentialData.set(1,i+m,POTENTIAL_VALUE);
        
    }
    
    PotentialData.print_to_file("./results/potential_data.out");
    printf("Complete.\n\n");

///////////////////////////////////////////////////////////////////////////////////////

   ///STAGE 1 - CREATE TRANSFER MATRICES T,T00,T11///

///////////////////////////////////////////////////////////////////////////////////////


   /*
   Create the transfer matrices T,T00,T11 using the Transfer Matrix Functions (These functions
   contain the mathematical definitions and hence calculate the values for each element in the
   matrix
   */
   
    
    
    

///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 2 - CALCULATE EIGENVALUES AND EIGENVECTORS OF MATRICES T,T00,T11///

///////////////////////////////////////////////////////////////////////////////////////

   /*
   STAGE 2 solves the eigensystem for matrices T,T00 and T11. The Eigenvalues are stored in
   a vector and the eigenvectors are stored in a matrix. 
   */
    
   rmax = eigen_limit_smax;
    
   vmax = eigen_limit_tmax;
   pmax = eigen_limit_tmax;

   printf("\n\n------------------------------------------------\n");
   printf("rmax: %4i vmax: %4i pmax: %4i",rmax,vmax,pmax);
   printf("\n------------------------------------------------\n\n");

    
   PSI_ROE_evec = Table(eigen_limit_tmax, MATRIX_SIZE);
   PSI_LAMBDA_evec = Table(eigen_limit_tmax, MATRIX_SIZE);

   for(int i=0;i < PSI_ROE_evec.n_column();i++){
      for(int j=0; j< PSI_ROE_evec.n_row();j++){
          PSI_ROE_evec.set(i,j,PSI_ROE(i, Delta*(j-m)));
          PSI_LAMBDA_evec.set(i,j,PSI_LAMBDA(i, Delta*(j-m)));
          
      }
   }
    
   LAMBDA_ROE_eval = Vector(eigen_limit_tmax);
   LAMBDA_LAMBDA_eval = Vector(eigen_limit_tmax);
   
   for (int i=0; i < LAMBDA_ROE_eval.dimension(); i++) {
        LAMBDA_ROE_eval.set(i,LAMBDA_ROE(i));
        LAMBDA_LAMBDA_eval.set(i,LAMBDA_LAMBDA(i));
   }


   OutputAnalyticalValues();

 
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

    string ndata;
    ifstream nfile_list("N.list");
    
    
    if(nfile_list.is_open()){
        
        getline (nfile_list,ndata);
        stringstream stringStream(ndata);
        string line;
        while(getline(stringStream, line))
        {
            vector<string> wordVector;
            size_t prev = 0, pos;
            while ((pos = line.find_first_of("=,", prev)) != string::npos){
                if(pos > prev){
                    wordVector.push_back(line.substr(prev, pos-prev));
                }
                prev = pos+1;
            }
            if(prev < line.length()){
                wordVector.push_back(line.substr(prev, string::npos));
            }
            
            for(int i=1;i < (int)wordVector.size();i++){
                
                string dir_name = "N" + (string)wordVector[i].c_str();
                string mkdir_command = "mkdir -p " + dir_name;
                
                system(mkdir_command.c_str());
                
 	            N = atoi(wordVector[i].c_str());
                
                printf("\nRunning N=%i\n",N);
                NucleationRun(0,0);
                
                string cp_command = "cp -r results/* ./" + dir_name;
                
                system(cp_command.c_str());
                
            }
        }
        nfile_list.close();
    }
    else{
        printf("*File N.list could not be opened.\n");
        
    }

   printf("\n");

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
    
    system("rm -f RUNNING.status");
    system("touch SIM_FINISHED.status");
   
   return 0;
}
