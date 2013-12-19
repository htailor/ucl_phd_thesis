#include "Vector.hh"
#include <stdio.h>
#include <iostream>

using namespace std;

Vector::Vector(){}

Vector::Vector(int _n){
    
    _N = _n;
    _VECTOR = new double[_N];
    
}

int Vector::dimension(){
    
    return _N;
    
}

void Vector::set(int _i, double _value){
    
    _VECTOR[_i] = _value;
    
}

double Vector::get(int _i){
    
    return _VECTOR[_i];
    
}

void Vector::print_to_file(const char* filename){
    
    FILE * vector_file = fopen(filename,"w");
    
    for(int i=0;i<_N;i++){
        fprintf(vector_file,"%e \n", _VECTOR[i]);
    }
    
    fclose(vector_file);
    
}

int Vector::num_zero_elements(){
    
    int count=0;
    for(int i=0 ; i<_N ; i++){
        if( _VECTOR[i] == 0 ){
            count = count + 1;
        }
    }
    return count;
}

void Vector::print_to_screen(){
    
    for(int i=0;i<_N;i++){
        printf("%e \n", _VECTOR[i]);
    }
    
}


