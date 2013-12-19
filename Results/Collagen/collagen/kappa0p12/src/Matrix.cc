#include "Matrix.hh"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

extern "C" void dsyevx_ ( const char* jobz, const char* range, const char* uplo, int* n, double* a, int* lda, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, double* w, double* z, int* ldz, double* work, int* lwork, int* iwork, int* ifail, int* info );

extern double Delta;

using namespace std;

Matrix::Matrix(){}

Matrix::Matrix(int _n){
    
    _N = _n;
    
    _ELEMENT_ARRAY_SIZE = _N * _N;
    
    _MATRIX = new double [_ELEMENT_ARRAY_SIZE];

    _EIGEN_LIMIT = _N;
    
}



int Matrix::dimension(){
    
    return _N;
    
}

int Matrix::eigen_limit(){
    
    return _EIGEN_LIMIT;
    
}



double Matrix::get(int _i,int _j){
    
    return _MATRIX[_i*_N+_j];
    
}

void Matrix::set(int _i,int _j, double _value){
    
    _MATRIX[_i*_N+_j] = _value;
    
}



void Matrix::print_to_file(const char* filename){
    
    FILE * matrix_file = fopen(filename,"w");
    
    for(int _i=0;_i<_N;_i++){
        for (int _j=0;_j<_N;_j++){
            fprintf(matrix_file,"%e\t", _MATRIX[_i*_N+_j]);
        }
        fprintf(matrix_file,"\n");
    }
    fclose(matrix_file);
}

void Matrix::print_row(const char* filename, int _row){
    
    FILE * row_file = fopen(filename,"w");
    
    for (int _j=0;_j<_N;_j++){
            fprintf(row_file,"%e ", _MATRIX[_row*_N+_j]);
        }
    fclose(row_file);
}

void Matrix::print_column(const char* filename, int _column){
    
    FILE * row_file = fopen(filename,"w");
    
    for (int _i=0;_i<_N;_i++){
            fprintf(row_file,"%e\n", _MATRIX[_i*_N+_column]);
        }
    fclose(row_file);
}



void Matrix::calculate_eigensystem(int _eigen_limit){
    
    _EIGEN_LIMIT = _eigen_limit;
    
    printf("Overriding eigen limit parameter: %i\n",_EIGEN_LIMIT);
    
    int LDA=_N;
    int LDZ=_N;
    
    int il, iu, m, info, lwork;
    double abstol, vl, vu;
    double wkopt;
    double* work;
    
    /* Local arrays */
    /* iwork dimension should be at least 5*n */
    int iwork[5*_N];
    int ifail[_N];

    double* w = new double[_EIGEN_LIMIT];
    double* z = new double[LDZ*_EIGEN_LIMIT];
   
    _EVAL = new double [_EIGEN_LIMIT];
    _EVEC = new double [_EIGEN_LIMIT*LDZ];
    
    /* Executable statements */
    printf("Lapack DSYEVX Results\n");
    
    /* Negative abstol means using the default value */
    abstol = -1.0;
    
    /* Set il, iu to compute NSELECT smallest eigenvalues */
    il = _N - _EIGEN_LIMIT + 1;
    iu = _N;
    
    /* Query and allocate the optimal workspace */
    lwork = -1;
    
    dsyevx_("Vectors", "Indices", "Upper", &_N, _MATRIX, &LDA, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &LDZ, &wkopt, &lwork, iwork, ifail, &info);
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );
    
    /* Solve eigenproblem */
    dsyevx_("Vectors", "Indices", "Upper", &_N, _MATRIX, &LDA, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &LDZ, work, &lwork, iwork, ifail, &info);
    
    /* Check for convergence */
    if( info > 0 ) {
        printf("The algorithm failed to compute eigenvalues.\n");
        exit( 1 );
    }
    
    /* Print the number of eigenvalues found */
    printf("\nThe total number of eigenvalues found:%2i\n", m );
    /* Free workspace */

    free((void*)work);

    for (int _eval = 0;_eval < _EIGEN_LIMIT; _eval++){
    	_EVAL[_eval] = w[_EIGEN_LIMIT-_eval-1];
    }
 
    for( int x = 0; x < _N; x++ ){
	for( int y = (_EIGEN_LIMIT-1); y > -1; y-- ){
        	_EVEC[((_EIGEN_LIMIT-1)-y)*_N+x]= z[y*_N+x];
        }
    }

}

void Matrix::calculate_eigensystem(){
    
  
    printf("Overriding eigen limit parameter: %i\n",_EIGEN_LIMIT);
    
    int LDA=_N;
    int LDZ=_N;
    
    int il, iu, m, info, lwork;
    double abstol, vl, vu;
    double wkopt;
    double* work;
    
    /* Local arrays */
    /* iwork dimension should be at least 5*n */
    int iwork[5*_N];
    int ifail[_N];

    double* w = new double[_EIGEN_LIMIT];
    double* z = new double[LDZ*_EIGEN_LIMIT];
   
    _EVAL = new double [_EIGEN_LIMIT];
    _EVEC = new double [_EIGEN_LIMIT*LDZ];
    
    /* Executable statements */
    printf("Lapack DSYEVX Results\n");
    
    /* Negative abstol means using the default value */
    abstol = -1.0;
    
    /* Set il, iu to compute NSELECT smallest eigenvalues */
    il = _N - _EIGEN_LIMIT + 1;
    iu = _N;
    
    /* Query and allocate the optimal workspace */
    lwork = -1;
    
    dsyevx_("Vectors", "Indices", "Upper", &_N, _MATRIX, &LDA, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &LDZ, &wkopt, &lwork, iwork, ifail, &info);
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );
    
    /* Solve eigenproblem */
    dsyevx_("Vectors", "Indices", "Upper", &_N, _MATRIX, &LDA, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &LDZ, work, &lwork, iwork, ifail, &info);
    
    /* Check for convergence */
    if( info > 0 ) {
        printf("The algorithm failed to compute eigenvalues.\n");
        exit( 1 );
    }
    
    /* Print the number of eigenvalues found */
    printf("\nThe total number of eigenvalues found:%2i\n", m );
    
    /* Free workspace */
    free((void*)work);

    for (int _eval = 0;_eval < _EIGEN_LIMIT; _eval++){
    	_EVAL[_eval] = w[_EIGEN_LIMIT-_eval-1];
    }
   
    for( int x = 0; x < _N; x++ ){
	for( int y = (_EIGEN_LIMIT-1); y > -1; y-- ){
        	_EVEC[((_EIGEN_LIMIT-1)-y)*_N+x]= z[y*_N+x];
        }
    }
    
}



double Matrix::eigenvalue(int _i){
    
    return _EVAL[_i];
    
}

double Matrix::eigenvector(int _evec, int _element){
    
    return _EVEC[_evec*_N+_element];
    
}



void Matrix::print_eigenvalues(const char* filename){

    FILE * file = fopen(filename,"w");
    
    for (int _i=0;_i<_EIGEN_LIMIT;_i++){
           fprintf(file,"%f\n",_EVAL[_i]);
        }

    fclose(file);
}

void Matrix::print_eigenvectors(const char* filename){

	//Eigenvectors printed in columns 

	FILE * file = fopen(filename,"w");
    
	for( int _element = 0; _element < _N; _element++ ){
		for( int _evec = 0; _evec < _EIGEN_LIMIT; _evec++ ){
			fprintf(file,"%f ", _EVEC[_evec*_N+_element] );
		}
                fprintf(file,"\n");
        }

	fclose(file);


}

void Matrix::screen_print_eigenvalue(int _i){

	if(_i<_EIGEN_LIMIT){
		printf("eval %i: %f\n",_i,_EVAL[_i]);
	}
	else{
		printf("\nERROR: index out of _EIGEN_LIMIT range\n");
	}
}

void Matrix::screen_print_eigenvector(int _evec){

	if(_evec<_EIGEN_LIMIT){
		printf("evec %i: ",_evec);
		for( int _element = 0; _element < _N; _element++ ){
			printf("%f\t", _EVEC[_evec*_N+_element] );
		}
		printf("\n");
	}
	else{
		printf("\nERROR: index out of _EIGEN_LIMIT range\n");
	}

}



double* Matrix::return_array(){
    
    return _MATRIX;
    
}



