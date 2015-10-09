#include <iostream>
#include "functions.h"
#include "test_functions.h"

using namespace std;

int main()
{
    int N = 400;
    double rho_min = 0;
    double rho_max = 7;
    double step_length = (rho_min-rho_max)/(N-1);
    //testArmadilloBasedJacobiRotation();
    vec eigen_values = zeros<vec>(N-2);
    mat B = mainMatrix(N-2, step_length, &harmonicOscillatorPotential);
    mat A = jacobiRotationMethod(B);
    for(int i = 0; i < N-2;i++){
        eigen_values(i) = A(i,i);
    }
    eigen_values = sort(eigen_values);
    for(int i = 0; i<3;i++){
        cout << "Lambda" << i+1 << " = " << eigen_values(i) << endl;
    }
    return 0;
}

