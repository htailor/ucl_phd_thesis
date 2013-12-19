#ifndef INCLUDE_MATRIX_HH
#define INCLUDE_MATRIX_HH

class Matrix{
    
public:
    
    Matrix();
    
    Matrix(int _n);

    

    int dimension();
    
    int eigen_limit();

    
    
    double get(int _i, int _j);

    void set(int _i, int _j, double _value);



    void print_to_file(const char* filename);

    void print_row(const char* filename, int _row);

    void print_column(const char* filename, int _column);


    
    void calculate_eigensystem(int _eigen_limit);
    
    void calculate_eigensystem();



    double eigenvalue(int _i);
    
    double eigenvector(int _i, int _j);
    
    
    
    void print_eigenvalues(const char* filename);

    void print_eigenvectors(const char* filename);
  
    void screen_print_eigenvalue(int _i);

    void screen_print_eigenvector(int _i);

    
    double* return_array();
    
    
private:
    
    double* _MATRIX;
    
    double* _EVEC;
    
    double* _EVAL;
    
    int _N, _EIGEN_LIMIT, _ELEMENT_ARRAY_SIZE;
    
};

#endif
