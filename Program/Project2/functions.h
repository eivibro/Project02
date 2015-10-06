#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <iostream>
#include <armadillo>

using namespace arma;

mat transformation_matrix(int row, int column, int dimensions);
mat main_matrix(int dimensions);
void find_largest_element(mat A, int *row, int *column);

#endif // FUNCTIONS_H
