/*
 * Output_Data.cc
 *
 *  Created on: 26 Oct 2010
 *      Author: ht
 */
#include "Output.hh"
#include "Matrix.hh"
#include "Vector.hh"
#include "Eigenfunctions.hh"
#include "Eigenvalues.hh"
#include <complex>
#include <iostream>
#include <sstream>
#include <stdio.h>

using namespace std;

extern int MATRIX_SIZE;
extern int smax;
extern int vmax;
extern int pmax;
extern int p2max;

extern int m;
extern double Delta;


extern Matrix TR;
extern Matrix TR_evec;
extern Vector TR_eval;

extern Matrix TROE_LAMBDA;
extern Matrix TROE_LAMBDA_evec;
extern Vector TROE_LAMBDA_eval;

extern Matrix PSI_R_evec;
extern Vector LAMBDA_R_eval;

void OutputTransferMatrices(){

	TROE_LAMBDA.print_to_file("TROE_LAMBDA_Matrix.log");
	system("mv *.log ./logs");

}

void OutputEigensystemValues(){

	for(int i=0;i<MATRIX_SIZE;i++){
		string filename;
		stringstream out;
		out << "Phi_roe_lambda" << i << ".TROE_LAMBDA";
		filename = out.str();
		TROE_LAMBDA_evec.print_row(filename.c_str(),i);
	}

	TROE_LAMBDA_eval.print_to_file("TROE_LAMBDA_eval.log");
	TROE_LAMBDA_evec.print_to_file("TROE_LAMBDA_evec.log");

	system("mv *.log ./logs");
	system("mv *.TROE_LAMBDA ./logs/TROE_LAMBDA");

}

void OutputAnalyticalValues(){

	for(int r=-smax;r<smax;r++){
	   for(int i=-m;i<=m;i++){
	      PSI_R_evec.set(r+smax,i+m,real(PSI_R(r,i*Delta)));
	   }
	}

	for(int i=0;i<smax;i++){
		string filename;
 		stringstream out;
   		out << "PSI_R_" << i << ".RA";
   		filename = out.str();
   		PSI_R_evec.print_row(filename.c_str(),i);
	}

	for(int s=1;s<smax;s++){
	   LAMBDA_R_eval.set(s-1,LAMBDA_R(s));
	}

	PSI_R_evec.print_to_file("PSI_R.log");
	LAMBDA_R_eval.print_to_file("LAMBDA_R.log");
	system("mv *.log ./logs");
	system("mv *.RA ./logs/R_Analytical");

}
