#include "Table.hh"
#include <stdio.h>
#include <stdlib.h>

extern double Delta;

using namespace std;

Table::Table(){}

Table::Table(int num_column, int num_row){
   
   _N_ROW = num_row;

   _N_COLUMN = num_column;

   _ELEMENT_ARRAY_SIZE = num_row * num_column;
    
   _TABLE = new double [_ELEMENT_ARRAY_SIZE];	

}



int Table::n_row(){

   return _N_ROW;

}

int Table::n_column(){

   return _N_COLUMN;

}



void Table::set(int _column, int _row, double _value){

	_TABLE[_column * _N_ROW + _row] = _value ; 

}

double Table::get(int _column, int _row){

	return _TABLE[_column * _N_ROW + _row];
 
}


void Table::print_to_file(const char* filename){
    
    FILE * table_file = fopen(filename,"w");

    for(int _i=0;_i<_N_ROW;_i++){    
        for (int _j=0;_j<_N_COLUMN;_j++){
            fprintf(table_file,"%e\t", _TABLE[_j*_N_ROW+_i]);
        }
        fprintf(table_file,"\n");
    }
    fclose(table_file);
}



//the ith row of the table
void Table::print_row(const char* filename, int _row){

    FILE * row_file = fopen(filename,"w");
    
    for (int _j=0;_j<_N_COLUMN;_j++){
            fprintf(row_file,"%e ", _TABLE[_j*_N_ROW + _row ]);
        }
    fclose(row_file);

}



// the ith column of the table
void Table::print_column(const char* filename, int _column){
    
    FILE * column_file = fopen(filename,"w");
    
    for (int _i=0;_i<_N_ROW;_i++){
            fprintf(column_file,"%e\n", _TABLE[_column*_N_ROW + _i]);
        }
    fclose(column_file);
}



double* Table::return_array(){

   return _TABLE;

}



