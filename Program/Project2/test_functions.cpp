#include <test_functions.h>


void testJacobiRotation(){
    mat B(3,3,fill::zeros);
    B(0,0) = 2; B(0,1) = -1; B(0,2) = 3;
    B(1,0) = -1; B(1,1) = 1; B(1,2) = -1;
    B(2,0) = 3; B(2,2) = 3; B(2,1) = -1;
    mat A = jacobiRotationMethod(B);
    vec eigen_values = zeros<vec>(3);
    for(int i = 0; i < 3; i++){
        eigen_values(i) = A(i,i);
    }
    eigen_values = sort(eigen_values);
    //Eigen values found with numpy in python
    double *eigen_values_python = new double[3];
    eigen_values_python[0] = -0.55247461;
    eigen_values_python[1] = 0.6090937;
    eigen_values_python[2] = 5.94338091;
    double sum = 0;
    eigen_values.print();
    for(int i = 0; i < 3; i++){
        sum += (eigen_values_python[i]-eigen_values(i))*
                (eigen_values_python[i]-eigen_values(i));

    }
    if(sum > 1e-6){cout << "Jacobi rotation is not functioning correctly!" << endl;}
    else{cout << "Test for Jacobi rotation passed!" << endl;}
}

void testfindLargestOffDiagonalElement(){
    double max_squared = 10;
    int row = 0;
    int column = 0;
    mat B(20,20,fill::ones);
    B(13,17) = -4;
    B(17,13) = -4;
    max_squared = findLargestOffDiagonalElement(B, &row, &column);
    if(max_squared == 16 && row == 17 && column ==13){
        cout << "findLargestOffDiagonalElement passed!" << endl;
    }
}



