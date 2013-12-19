#ifndef INCLUDE_VECTOR_HH
#define INCLUDE_VECTOR_HH


class Vector{
    
public:
    
    Vector();
    Vector(int _n);
    
    int dimension();
    int num_zero_elements();
    
    double get(int _i);
    
    void set(int _i, double _value);
    void print_to_file(const char* filename);
    void print_to_screen();
    
    
private:
    
    int _N;
    double* _VECTOR;
    
};

#endif
