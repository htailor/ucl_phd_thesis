#ifndef INCLUDE_TABLE_HH
#define INCLUDE_TABLE_HH


class Table{
 
public:

 Table();

 Table(int num_row, int num_column);



 int n_row();

 int n_column();



 double get(int i, int j);

 void set(int i, int j, double value);



 void print_to_file(const char* filename);



 void print_row(const char* filename, int _row);

 void print_column(const char* filename, int _column);



 double* return_array();



private:

 int _N, _EIGEN_LIMIT, _ELEMENT_ARRAY_SIZE , _N_COLUMN, _N_ROW;

 double* _TABLE;

};

#endif
