#include "functions.h"
#include "lib.h"
#include <iostream>
#include <fstream>
#include <iomanip>

//Creating matrix for calculations with Armadillo
mat mainMatrix(int dimensions, double step_length, double(*potential)(double,int,double), double omega_r)
{
    double diagonal_fixed = 2./(step_length*step_length);
    double off_diagonal = -1./(step_length*step_length);
    mat A(dimensions, dimensions, fill::zeros);
    A(0,0) = diagonal_fixed+potential(step_length, 1, omega_r);
    A(0,1) = off_diagonal;
    A(dimensions-1,dimensions-1) = diagonal_fixed + potential(step_length,dimensions,omega_r);
    A(dimensions-1, dimensions-2) = off_diagonal;
    for(int i = 1; i < dimensions-1; i++){
        A(i,i) = diagonal_fixed + potential(step_length,i+1,omega_r);
        A(i,i+1) = off_diagonal;
        A(i,i-1) = off_diagonal;
    }
    return A;
}

//Potential for two electrons repulsed by Coulomb potential
double potentialRepulsiveColumbInteractions(double step_length, int i, double omega_r)
{
    return ((step_length*i)*(step_length*i)*omega_r*omega_r+1./(step_length*i));
}

//Finding the largest off-diagonal element in matrix A
//Matrices in Armadillo is stored in a column-major order, and
//thus we loop through the columns in the outer for-loop.
//We assume that we have a symmetric matrix so that we only need
//to loop throug the elements on either side of the diagonal
//Returns the maximum element squared
double findLargestOffDiagonalElement(mat A, int *row, int *column)
{
    double max = 0;
    int n = A.n_rows;
    for(int j = 0; j < n; j++){
        for(int i = j+1; i < n; i++){
            double k = A(i,j)*A(i,j);
            if(k > max){
                max = k; *row = i; *column = j;
            }
        }
    }
    return max;
}

//Jacobi Rotation with armadillo matrix multiplication
//Few lines of code, but really slow.
mat jacobiRotationArmadillo(mat A, int *row, int *column){
    int dim = A.n_rows;
    double epsilon = 1e-8;
    double max_squared = 100;
    while(max_squared > epsilon){
        max_squared = findLargestOffDiagonalElement(A, row, column);
        double tau = (A(*column, *column)-A(*row,*row))/(2*A(*row,*column));
        double t = -tau-sqrt(1+tau*tau);
        double c = 1/sqrt(1+t*t);
        double s = t*c;
        mat S = transformationMatrix(*row,*column,dim, c, s);
        A = trans(S)*A*S;
    }
    return A;
}

//The dimensionless quantum harmonical oscillator
double harmonicOscillatorPotential(double step_length, int i, double omega)
{
    return (step_length*i)*(step_length*i);
}

//Single rotation for the Jacobi rotation method
//More lines than Armadillo version, but faster
mat jacobiRotation(mat A, int k, int l){
    int dim = A.n_rows;
    //Calculating sin and cos of theta to cancel the
    //largest off diagonal elements A(k,l) and A(l,k)
    double s, c;
    if(A(k, l)!=0){
        double t1, t2, t, tau;
        tau = (A(l,l)-A(k,k))/(2*A(k,l));
        double sqrt_root = sqrt(1+tau*tau);
        t1 = -tau-sqrt_root;
        t2 = -tau+sqrt_root;
        if(t2>t1){t = t1;}
        else{t = t2;}
        c = 1/sqrt(1+t*t);
        s = t*c;
    }else{c = 1.; s=0.;}
    //Calculating the new matrix elements given by the rotation
    double Akk, All, Akl, Ail, Aik;
    Akk = A(k,k);
    All = A(l,l);
    Akl = A(k,l);

    A(k,k) = Akk*c*c - 2.*c*s*Akl + All*s*s;
    A(l,l) = All*c*c + 2.*c*s*Akl + Akk*s*s;
    //From the definition we get
    A(l,k) = 0;  A(k,l) = 0;
    for(int i = 0; i < dim; i++){
        if(i != l && i != k){
            //A(i,i) = A(i,i) but obviously already true
            Aik = A(i,k);
            Ail = A(i,l);
            A(i,k) = Aik*c-Ail*s;
            A(i,l) = Ail*c+Aik*s;
            A(l,i) = A(i,l);
            A(k,i) = A(i,k);
        }
    }return A;
}

//The Jacobi method utilizing jacobiRotation(
mat jacobiRotationMethod(mat A)
{
    int dim = A.n_rows;
    int row = 0;
    int column = 0;
    double epsilon = 1e-8;
    double max_iterations = dim*dim*dim;
    double max_squared = 10;
    int iterations = 0;
    while(max_squared>epsilon && iterations < max_iterations){
        max_squared = findLargestOffDiagonalElement(A, &row, &column);
        A = jacobiRotation(A, row, column);
        iterations++;
    }
    cout << "Number of Jacobi rotations: " << iterations << endl;
    return A;
}

//Creates transformation matrix for the armadillo method
mat transformationMatrix(int row, int column, int dimensions, double c, double s){
    mat A(dimensions,dimensions, fill::eye);
    A(row,row) = c;
    A(column, column) = c;
    A(row, column) = s;
    A(column, row) = -s;
    return A;
}

//Function utilizing the function tqli() in lib.cpp to calculate eigenvalues and eigenvectors
void tqliEigensolver(double *off_diagonal, double *diagonal, double **output_matrix, int n, double step_length,
                     double omega_r,double(*potential)(double,int,double))
{
    double constant_diagonal = 2./(step_length*step_length);
    double constant_off_diagonal = -1./(step_length*step_length);
    for(int i = 0; i < n; i++){
        diagonal[i] = constant_diagonal+potential(step_length,i+1,omega_r);
        off_diagonal[i] = constant_off_diagonal;
        output_matrix[i][i] = 1;
        for(int j = i + 1; j < n; j++){
            output_matrix[i][j] = 0;
        }
    }
    tqli(diagonal, off_diagonal, n, output_matrix);
}



void writeResultsToFile(string filename, vec eigen_values, double *diagonal, double *rho,
                        double **output_matrix)
{
    int dim = eigen_values.n_elem + 2;
    double **eigen_vector_matrix;
    eigen_vector_matrix = (double **) matrix(dim,3,sizeof(double));
    for(int j = 0; j < 3; j++){
        int index = -1;
        for(int i = 0; i < dim-2; i++){
            if(eigen_values(j) == diagonal[i]){index = i;}
        }
        if(index == -1){
            cout << "Could not find eigenvalue, something went wrong" << endl;
        };
        eigen_vector_matrix[0][j] = 0;
        eigen_vector_matrix[dim-1][j] = 0;
        for(int i = 1; i < dim-1; i++){
            eigen_vector_matrix[i][j] = output_matrix[i-1][index];
        }
    }

    ofstream outfile;
    outfile.open(filename);
    outfile << dim << endl;
    outfile << eigen_values[0] << endl;
    outfile << eigen_values[1] << endl;
    outfile << eigen_values[2] << endl;
    for(int i = 0; i < dim; i++){
        outfile << left << setw(15) << setprecision(4) << setw(15) <<
                   eigen_vector_matrix[i][0]*eigen_vector_matrix[i][0]
                << setw(15) << eigen_vector_matrix[i][1]*eigen_vector_matrix[i][1]
                << setw(15) << eigen_vector_matrix[i][2]*eigen_vector_matrix[i][2]
                << "\t" << rho[i]<< endl;
    }
    free_matrix((void**) eigen_vector_matrix);
}
