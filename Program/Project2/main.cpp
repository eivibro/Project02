#include <iostream>
#include <cstdlib>
#include "functions.h"
#include "test_functions.h"
#include "lib.h"
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;

int main()
{
    //Running unit test before calculations
    testJacobiRotation();
    testfindLargestOffDiagonalElement();
    //Variables for measuring time
    clock_t start, finish;
    //initializing values
    int N;
    double rho_min, rho_max;
    ifstream infile("inValue.txt");;
    infile >> N >> rho_min >> rho_max;
    infile.close();
    double step_length = (rho_max-rho_min)/(N-1);
    //First case, single harmonic oscillator
    //testArmadilloBasedJacobiRotation();
    vec eigen_values1 = zeros<vec>(N-2);
    mat B = mainMatrix(N-2, step_length, &harmonicOscillatorPotential, 0);
    start = clock();
    mat A = jacobiRotationMethod(B);
    finish = clock();
    cout << "Time my algo\t" << setprecision(7) << (finish-start)/double(CLOCKS_PER_SEC) << endl;
    //Extraxting and sorting the eigenvalues
    for(int i = 0; i < N-2;i++){
        eigen_values1(i) = A(i,i);
    }
    eigen_values1 = sort(eigen_values1);
    for(int i = 0; i<3;i++){
        cout << "Lambda" << i+1 << " = " << eigen_values1(i) << endl;
    }
//    //Running tqli()
//    double omega_r = 5;
//    double *off_diagonal, *diagonal, **output_matrix;
//    diagonal = new double[N-2];
//    off_diagonal = new double[N-2];
//    //Initializing output matrix
//    output_matrix = (double **) matrix(N-2,N-2,sizeof(double));
//    //Solving
//    start = clock();
//    tqliEigensolver(off_diagonal,diagonal,output_matrix,N-2, step_length,
//                    omega_r, &potentialRepulsiveColumbInteractions);
//    finish = clock();
//    cout << "Time tqli()\t" << setprecision(7) << (finish-start)/double(CLOCKS_PER_SEC) << endl;

    //Second case, two electrons interacting with Columb potential
//    double omega_r = 5;
//    vec eigen_values2 = zeros<vec>(N-2);
//    mat C = mainMatrix(N-2, step_length, &potentialRepulsiveColumbInteractions, omega_r);
//    mat D = jacobiRotationMethod(C);
//    //Extraxting and sorting the eigenvalues
//    for(int i = 0; i < N-2;i++){
//        eigen_values2(i) = D(i,i);
//    }
//    eigen_values2 = sort(eigen_values2);
//    for(int i = 0; i<3;i++){
//        cout << "Lambda" << i+1 << " = " << eigen_values2(i) << endl;
//    }

    //Part d
    //Calculating the eigenvectors with the tqli() method from lib.cpp
    //Initializing the matrix elements
//    double omega_r = 5;
//    double *off_diagonal, *diagonal, **output_matrix;
//    diagonal = new double[N-2];
//    off_diagonal = new double[N-2];
//    //Initializing output matrix
//    output_matrix = (double **) matrix(N-2,N-2,sizeof(double));
//    //Solving
//    tqliEigensolver(off_diagonal,diagonal,output_matrix,N-2, step_length,
//                    omega_r, &potentialRepulsiveColumbInteractions);
//    vec eigen_values1 = zeros<vec>(N-2);
//    for(int i = 0; i < N-2; i++){
//        eigen_values1(i) = diagonal[i];
//    }
//    eigen_values1 = sort(eigen_values1);
//    //Finding the index of the lowest eigenvalue in the diagonal elements
//    //so that we can extract the eigenvector
//    int index = -1;
//    for(int i = 0; i < N-2; i++){
//        if(eigen_values1(0) == diagonal[i]){index = i;}
//    }
//    if(index == -1){
//        cout << "Could not find eigenvalue, something went wrong" << endl;
//    }
//    double *rho;
//    rho = new double[N];
//    rho[0] = 0;
//    rho[N-1] = rho_max;
//    for(int i = 1; i < N-1; i++){
//        rho[i] = step_length*i;
//    }
//    writeResultsToFile("results.txt", eigen_values1, diagonal, rho, output_matrix);
//    free_matrix((void**) output_matrix);
//    delete [] diagonal; delete [] off_diagonal; delete [] rho;
    return 0;
}

