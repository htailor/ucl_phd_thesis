#include <stdio.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include <getopt.h>
#include "Matrix.hh"
#include "Vector.hh"
#include "Functions.hh"


#define PI 3.1415926535897932384626433832795

using namespace std;

extern double eta_d;
extern double e0;
extern double Beta;

extern double L;
extern int m;
extern int N;
extern double Delta;
extern double TOLERANCE;

extern double sigma;
extern double kappa;
extern double kappa_sigma_r;
extern int umin;
extern int umax;
extern int MATRIX_SIZE;

extern Matrix t;
extern Matrix t_evec;
extern Vector t_eval;

extern Matrix t00;
extern Matrix t00_evec;
extern Vector t00_eval;

void Menu(int argc, char *argv[])
{
    char s;
    int long_opt_index = 0;
    int longval;
    struct option long_options[] = {        /* long options array. Items are all caSe SensiTivE! */
        { "L", 1, NULL, 'L'  },             /* --add or -a  */
        { "m", 1, NULL, 'm'  },             /* --back or -b */
        { "N", 1, NULL, 'N'  },             /* --back or -b */
        { "kappa", 1, &longval, 'K' },/* return 'c', or return 0 and set longval to 'c' if "check" is parsed */
        { "sigma", 1, &longval, 'S' },/* return 'c', or return 0 and set longval to 'c' if "check" is parsed */
        { "umin", 1, &longval, 'u' },
        { "umax", 1, &longval, 'U' },
        { 0,    0,    0,    0   }           /* terminating -0 item */
    };


    while ((s = getopt_long(argc, argv, "L:m:N:K:S:u:U:h", long_options, &long_opt_index)) != -1) {
       switch (s) {
           case 'L':   /* long_opt_index does not make sense for these */
               L = atof(optarg);
               break;
           case 'm':
               m = atoi(optarg);
               break;
           case 'N':
               N = atoi(optarg);
               break;
           case 0:     /* this is returned for long options with option[i].flag set (not NULL). */
                       /* the flag itself will point out the option recognized, and long_opt_index is now relevant */
               switch (longval) {
                   case 'K':
                       kappa = atof(optarg);
                       break;
                   case 'S':
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
               printf("Usage: Nucleation -L [] -m [] -N [] --kappa [] --sigma [] --umin [] --umax []\n");
               exit(1);
               break;
           default:
               printf("ERROR!\n");
       }
    }


   if(L==0 || m==0 || N==0 || kappa==0 || sigma==0 || umax==0){

   system("clear");

   printf("===========================================================================\n");
   printf("Nucleation DNA Simulation\n");
   printf("===========================================================================\n\n");

   cout << "Please enter simulation parameters:\n";
   cout << "\nL: ";
   cin >> L;
   cout << "m: ";
   cin >> m;
   cout << "N: ";
   cin >> N;
   cout << "\nPlease enter potential parameters:\n";
       
   cout << "\nkappa: ";
   cin >> kappa;
   cout << "sigma: ";
   cin >> sigma;

   cout << "\nPlease enter Partition Function parameters:\n";
   cout << "\nMinimum Extension: ";
   cin >> umin;
   cout << "Maximum Extension: ";
   cin >> umax;

}
    
   MATRIX_SIZE = 2*m+1;

   Delta = L/(2*m+1);
   TOLERANCE = 1e-5;
   kappa_sigma_r = kappa/sigma;

   e0 = (sigma*Squared(eta_d))/12.;
   
   system("rm ./results/Intact/*.out ./results/Frayed/*.out ./results/Bubble/*.out ./results/*.out 2> ./logs/DUMP.log");
   system("rm ./logs/T/*.T ./logs/T00/*.T00 ./logs/T11/*.T11 ./logs/*.log 2>> ./logs/DUMP.log");

   printf("\n\n-----------------\n");
   printf("Global Parameters\n");
   printf("-----------------\n");

   FILE * Simulation_Parameters = fopen("Parameters","w");

   fprintf(Simulation_Parameters,"L:\t%-5.2f\nm:\t%-3i\nN:\t%-3i\nkappa:\t%-8.6f\nsigma:\t%-8.6f\nkappa_sigma_r:\t%-8.6f\nDelta:\t%-5.4f\nExtension Minimum:\t%-3i\nExtension Maximum:\t%-3i\ne0:\t%-8.6f\nBeta:\t%-8.6f\n",L,m,N,kappa,sigma,kappa_sigma_r,Delta,umin,umax,e0,Beta);
   fclose(Simulation_Parameters);

   printf("L:\t%-5.2f\nm:\t%-3i\nN:\t%-3i\nkappa:\t%-8.6f\nsigma:\t%-8.6f\nkappa_sigma_r:\t%-8.6f\nDelta:\t%-5.4f\nExtension Minimum:\t%-3i\nExtension Maximum:\t%-3i\ne0:\t%-8.6f\n",L,m,N,kappa,sigma,kappa_sigma_r,Delta,umin,umax,e0);

   printf("\nTOLERANCE LEVEL:\t%-5.2e\n\n\n",TOLERANCE);

}


void OutputDataToFiles(){

   printf("Creating Output Files.....");

   t00.print_to_file("t00_Matrix.log");

   for(int i=0;i<MATRIX_SIZE;i++){
      string filename;
      stringstream out;
      out << "Phi_t_" << i << ".T00";
      filename = out.str();

      t00_evec.print_row(filename.c_str(),i);
   }

   t00_eval.print_to_file("t00_eval.log");
   //t00_evec.print_to_file("t00_evec.log");

   system("mv *.T00 ./logs/T00");
  
   system("mv *.log ./logs");
   system("mv Parameters ./results");

   printf("COMPLETE.\n\n\n");

}

