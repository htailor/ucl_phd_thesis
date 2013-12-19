#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include <getopt.h>
#include "Functions.hh"
#include "Matrix.hh"
#include "Vector.hh"

#define PI 3.1415926535897932384626433832795

using namespace std;

extern double L;
extern int m;
extern int N;
extern double Delta;
extern double TOLERANCE;
extern double eta_b;
extern double kappa;
extern double beta;
extern double sigma;

extern double kappa_sigma_r;

extern int umin;
extern int umax;
extern int MATRIX_SIZE;

extern double eta_d;
extern double e0;
extern double beta;


extern double T_ROE_LAMBDA_BETA;
extern double T_ROE_LAMBDA_MU;
extern double T_ROE_LAMBDA_C;
extern double T_ROE_LAMBDA_DELTA;
extern double T_ROE_LAMBDA_B;
extern double T_ROE_LAMBDA_K;
extern double T_ROE_LAMBDA_EV0;
extern double T_ROE_LAMBDA_C0;


void Menu(int argc, char *argv[])
{
    char z;
    int long_opt_index = 0;
    int longval;
    struct option long_options[] = {        /* long options array. Items are all caSe SensiTivE! */
        { "L", 1, NULL, 'L'  },             /* --add or -a  */
        { "m", 1, NULL, 'm'  },             /* --back or -b */
        { "kappa", 1, &longval, 'k' },/* return 'c', or return 0 and set longval to 'c' if "check" is parsed */
        { "sigma", 1, &longval, 's' },/* return 'c', or return 0 and set longval to 'c' if "check" is parsed */
        { "umin", 1, &longval, 'u' },
        { "umax", 1, &longval, 'U' },
        { 0,    0,    0,    0   }           /* terminating -0 item */
    };


    while ((z = getopt_long(argc, argv, "L:m:k:s:u:U:h", long_options, &long_opt_index)) != -1) {
       switch (z) {
           case 'L':   /* long_opt_index does not make sense for these */
               L = atof(optarg);
               break;
           case 'm':
               m = atoi(optarg);
               break;
           case 0:     /* this is returned for long options with option[i].flag set (not NULL). */
                       /* the flag itself will point out the option recognized, and long_opt_index is now relevant */
               switch (longval) {
                   case 'k':
                       kappa = atof(optarg);
                       break;
                   case 's':
                       sigma = atof(optarg);
                       break;
                   case 'u':
                       umin = atoi(optarg);
                       break;
                   case 'U':
                       umax = atoi(optarg);
                       break;
                   /* there's no default here */
               }
               break;
           case 'h':   /* mind that h is not described in the long option list */
               printf("Usage: Nucleation -L [] -m [] --kappa [] --sigma [] --umin [] --umax []\n");
               exit(1);
               break;
           default:
               printf("ERROR!\n");
       }
    }


   if(L==0 || m==0 || kappa==0 || sigma==0 || umax==0){

   printf("===========================================================================\n");
   printf("Nucleation Multi-Stranded Polymer (Collagen) Simulation\n");
   printf("===========================================================================\n\n");

   cout << "Please enter simulation parameters:\n\n";
   cout << "L: ";
   cin >> L;
   cout << "m: ";
   cin >> m;

   cout << "\nPlease enter potential parameters:\n\n";

   cout << "sigma: ";
   cin >> sigma;

   cout << "kappa: ";
   cin >> kappa;

   cout << "\nPlease enter Partition Function parameters:\n\n";
   cout << "Minimum Extension: ";
   cin >> umin;
   cout << "Maximum Extension: ";
   cin >> umax;

}

   kappa_sigma_r = kappa/sigma;


   MATRIX_SIZE = 2*m+1;

   Delta = L/(MATRIX_SIZE);

   beta = 38.6817; // eV^-1
   eta_d = 15.; // diameter of dna AA
   
   //beta = 1.;
   e0 = (sigma*Squared(eta_d))/12.; 
   
   TOLERANCE = 1e-5;



   T_ROE_LAMBDA_BETA = 4;
   T_ROE_LAMBDA_MU = 4*(3/kappa)*((6*beta*e0)/Squared(eta_d));
   T_ROE_LAMBDA_C = 0.5*pow((Squared(T_ROE_LAMBDA_MU)+2*T_ROE_LAMBDA_BETA*T_ROE_LAMBDA_MU),0.5);
   T_ROE_LAMBDA_DELTA = 2*T_ROE_LAMBDA_C + T_ROE_LAMBDA_MU + T_ROE_LAMBDA_BETA;
   T_ROE_LAMBDA_B = 0.5*pow(((Squared(T_ROE_LAMBDA_DELTA)-Squared(T_ROE_LAMBDA_BETA))/T_ROE_LAMBDA_DELTA),0.5);
   T_ROE_LAMBDA_K = log(T_ROE_LAMBDA_DELTA/T_ROE_LAMBDA_BETA);
   T_ROE_LAMBDA_EV0 = pow((4*PI/T_ROE_LAMBDA_DELTA),0.5);
   T_ROE_LAMBDA_C0 = 0.5*pow(PI/T_ROE_LAMBDA_C,-0.25);
     
   system("rm ./results/Intact/*.out ./results/Frayed/*.out ./results/Bubble/*.out ./results/*.out 2> ./logs/DUMP.log");
   system("rm ./logs/T_R/*.R ./logs/T_ROE/*.ROE ./logs/T_LAMBDA/*.LAMBDA ./logs/*.log 2>> ./logs/DUMP.log");

   printf("-----------------\n");
   printf("Global Parameters\n");
   printf("-----------------\n");

   FILE * Simulation_Parameters = fopen("Parameters","w");

   fprintf(Simulation_Parameters,"L:\t%3.6f\nm:\t%i\nkappa:\t%3.6f\nsigma:\t%3.6f\nkappa_sigma_r:\t%3.6f\nDelta:\t%3.6f\nExtension Minimum:\t%i\nExtension Maximum:\t%i\nbeta:\t%-5.2f\nmu:\t%-5.2f\nc:\t%-5.2f\ndelta:\t%-5.2f\nb:\t%-5.2f\nk:\t%-5.2f\nev0:\t%-5.2f\nC0:\t%-5.2f\ne0:\t%-8.6f\nbeta:\t%-8.6f",L,m,kappa,sigma,kappa_sigma_r,Delta,umin,umax,T_ROE_LAMBDA_BETA,T_ROE_LAMBDA_MU,T_ROE_LAMBDA_C,T_ROE_LAMBDA_DELTA,T_ROE_LAMBDA_B,T_ROE_LAMBDA_K,T_ROE_LAMBDA_EV0,T_ROE_LAMBDA_C0,e0,beta);
   fclose(Simulation_Parameters);

   printf("\nL:\t%f\nm:\t%i\nkappa:\t%f\nsigma:\t%f\nkappa_sigma_r:\t%f\nDelta:\t%f\n\nExtension Minimum:\t%i\nExtension Maximum:\t%i\n\n\nbeta:\t%-5.2f\nmu:\t%-5.2f\nc:\t%-5.2f\ndelta:\t%-5.2f\nb:\t%-5.2f\nk:\t%-5.2f\nev0:\t%-5.2f\nC0:\t%-5.2f\ne0:\t%-8.6f\nbeta:\t%-8.6f\n\n",L,m,kappa,sigma,kappa_sigma_r,Delta,umin,umax,T_ROE_LAMBDA_BETA,T_ROE_LAMBDA_MU,T_ROE_LAMBDA_C,T_ROE_LAMBDA_DELTA,T_ROE_LAMBDA_B,T_ROE_LAMBDA_K,T_ROE_LAMBDA_EV0,T_ROE_LAMBDA_C0,e0,beta);
  

   system("mv Parameters ./results");

}



