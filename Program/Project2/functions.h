#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <fstream>

using namespace std;
using namespace arma;

mat transformationMatrix(int row, int column, int dimensions, double c, double s);
mat mainMatrix(int dimensions, double step_length, double (*potential)(double, int,double), double omega_r);
double findLargestOffDiagonalElement(mat A, int *row, int *column);
mat jacobiRotation(mat A, int row, int column);
mat jacobiRotationMethod(mat A);
mat jacobiRotationArmadillo(mat A, int *row, int *column);
double harmonicOscillatorPotential(double step_length, int i, double omega);
double potentialRepulsiveColumbInteractions(double step_length, int i, double omega_r);
vec simple_triDiag_solver(mat A, double eigen_value);
void tqliEigensolver(double *off_diagonal,double *diagonal,double **output_matrix,int n, double step_length,
                     double omega_r, double(*potential)(double,int,double));
void writeResultsToFile(string s, vec eigen_values, double *diagonal, double *rho, double **output_matrix);

#endif // FUNCTIONS_H
