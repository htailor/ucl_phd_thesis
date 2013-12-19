/*
 * Output_Data.cc
 *
 *  Created on: 26 Oct 2010
 *      Author: ht
 */
#include "Output.hh"
#include "Matrix.hh"
#include "Vector.hh"
#include "Functions.hh"
#include <iostream>
#include <sstream>
#include <stdio.h>

using namespace std;

extern int MATRIX_SIZE;
extern int rmax;
extern int vmax;
extern int pmax;
extern int p2max;

extern double Delta;

extern Matrix TR;
extern Matrix TR_evec;
extern Vector TR_eval;

extern Matrix TROE;
extern Matrix TROE_evec;
extern Vector TROE_eval;

extern Matrix TLAMBDA;
extern Matrix TLAMBDA_evec;
extern Vector TLAMBDA_eval;

extern Matrix PSI_R_evec;
extern Matrix PSI_ROE_evec;
extern Matrix PSI_LAMBDA_evec;

extern Vector LAMBDA_R_eval;
extern Vector LAMBDA_ROE_eval;
extern Vector LAMBDA_LAMBDA_eval;

void OutputTransferMatrices(){

	TR.print_to_file("TR_Matrix.log");
	TROE.print_to_file("TROE_Matrix.log");
	TLAMBDA.print_to_file("TLAMBDA_Matrix.log");
	system("mv *.log ./logs");

}

void OutputEigensystemValues(){

	for(int i=0;i<MATRIX_SIZE;i++){
		string filename;
		stringstream out;
		out << "Phi_roe_" << i << ".TROE";
		filename = out.str();
		TROE_evec.print_row(filename.c_str(),i);
	}

	for(int i=0;i<MATRIX_SIZE;i++){
		string filename;
		stringstream out;
		out << "Phi_lambda_" << i << ".TLAMBDA";
		filename = out.str();
		TLAMBDA_evec.print_row(filename.c_str(),i);
	}

	TR_eval.print_to_file("TR_eval.log");
	TR_evec.print_to_file("TR_evec.log");

	TROE_eval.print_to_file("TROE_eval.log");
	TROE_evec.print_to_file("TROE_evec.log");

	TLAMBDA_eval.print_to_file("TLAMBDA_eval.log");
	TLAMBDA_evec.print_to_file("TLAMBDA_evec.log");

	system("mv *.TROE ./logs/TROE");
	system("mv *.TLAMBDA ./logs/TLAMBDA");

}

void OutputAnalyticalValues(){

	for(int i=0;i<rmax;i++){
		string filename;
		stringstream out;
		out << "PSI_R_" << i << ".RA";
		filename = out.str();
		PSI_R_evec.print_row(filename.c_str(),i);
	}

	for(int i=0;i<pmax;i++){
		string filename;
		stringstream out;
		out << "PSI_ROE_" << i << ".ROEA";
		filename = out.str();
		PSI_ROE_evec.print_row(filename.c_str(),i);
	}

	for(int i=0;i<vmax;i++){
		string filename;
		stringstream out;
		out << "PSI_LAMBDA_" << i << ".LAMBDAA";
		filename = out.str();
		PSI_LAMBDA_evec.print_row(filename.c_str(),i);
	}

	PSI_R_evec.print_to_file("PSI_R.log");
	PSI_ROE_evec.print_to_file("PSI_ROE.log");
	PSI_LAMBDA_evec.print_to_file("PSI_LAMBDA.log");

	LAMBDA_R_eval.print_to_file("LAMBDA_R.log");
	LAMBDA_ROE_eval.print_to_file("LAMBDA_ROE.log");
	LAMBDA_LAMBDA_eval.print_to_file("LAMBDA_LAMBDA.log");

	system("mv *.RA ./logs/R_Analytical");
	system("mv *.ROEA ./logs/ROE_Analytical");
	system("mv *.LAMBDAA ./logs/LAMBDA_Analytical");

}

void CheckNormalisation(){

	   cout << "\n";
	   for(int i=0;i<MATRIX_SIZE;i++){
	   	  double intNum=0;
	   	  for(int j=0;j<MATRIX_SIZE;j++){
	         intNum = intNum + Squared(TROE_evec.get(i,j))*Delta;
	      }
	   	  cout << "Normalisation Integral for TROE Eigenvector (" << i << ") = " << intNum << "\n";
	   }

	   cout << "\n";
	   for(int i=0;i< MATRIX_SIZE;i++){
	      double intNum=0;
	      for(int j=0;j<MATRIX_SIZE;j++){
	         intNum = intNum + Squared(TR_evec.get(i,j))*Delta;
	      }
	      cout << "Normalisation Integral for TR Eigenvector (" << i << ") = " << intNum << "\n";
	   }


}
