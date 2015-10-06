#include "functions.h"

mat transformation_matrix(int row, int column, int dimensions)
{
    mat A(dimensions,dimensions, fill::zeros);
    A(row,column) = 1;
    return A;
}


mat main_matrix(int dimensions)
{
    mat A(dimensions, dimensions, fill::zeros);
    A(0,0) = 2;
    A(0,1) = -1;
    A(dimensions-1,dimensions-1) = 2;
    A(dimensions-1, dimensions-2) = -1;
    for(int i = 1; i < dimensions-1; i++){
        A(i,i) = 2;
        A(i,i+1) = -1;
        A(i,i-1) = -1;
    }
    return A;
}

//Finding the largest off-diagonal element in matrix A
//Matrices in Armadillo is stored in a column-major order, and
//thus we loop through the columns in the outer for-loop
void find_largest_element(mat A, int *row, int *column)
{
    double max = 0;
    int n = A.n_rows;
    for(int j = 0; j < n; j++){
        for(int i = 0; i < n; i++){
            if(j != i){
                double k = A(i,j)*A(i,j);
                if(k > max){
                    max = k; *row = i; *column = j;
                }
            }
        }
    }
}
