#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <iostream>
#include <armadillo>

using namespace arma;

mat transformationMatrix(int row, int column, int dimensions, double c, double s);
mat mainMatrix(int dimensions, double step_length, double (*potential)(double, int));
double findLargestOffDiagonalElement(mat A, int *row, int *column);
mat jacobiRotation(mat A, int row, int column);
mat jacobiRotationMethod(mat A);
mat jacobiRotationArmadillo(mat A, int *row, int *column);
double harmonicOscillatorPotential(double step_length, int i);

#endif // FUNCTIONS_H
