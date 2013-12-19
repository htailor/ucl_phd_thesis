/*
 * Output_Data.cc
 *
 *  Created on: 26 Oct 2010
 *      Author: ht
 */
#include "Output.hh"
#include "Table.hh"
#include "Vector.hh"
#include "Functions.hh"
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdio.h>

using namespace std;

extern int MATRIX_SIZE;
extern int rmax;
extern int vmax;
extern int pmax;

extern double Delta;

extern Table PSI_ROE_evec;
extern Table PSI_LAMBDA_evec;

extern Vector LAMBDA_ROE_eval;
extern Vector LAMBDA_LAMBDA_eval;


void OutputAnalyticalValues(){

	for(int i=0;i<pmax;i++){
		string filename;
		stringstream out;
		out << "./logs/ROE_Analytical/PSI_ROE_" << i << ".ROEA";
		filename = out.str();
		PSI_ROE_evec.print_column(filename.c_str(),i);
	}

	for(int i=0;i<vmax;i++){
		string filename;
		stringstream out;
		out << "./logs/LAMBDA_Analytical/PSI_LAMBDA_" << i << ".LAMBDAA";
		filename = out.str();
		PSI_LAMBDA_evec.print_column(filename.c_str(),i);
	}

	PSI_ROE_evec.print_to_file("./logs/PSI_ROE.log");
	PSI_LAMBDA_evec.print_to_file("./logs/PSI_LAMBDA.log");

	LAMBDA_ROE_eval.print_to_file("./logs/LAMBDA_ROE.log");
	LAMBDA_LAMBDA_eval.print_to_file("./logs/LAMBDA_LAMBDA.log");

}
