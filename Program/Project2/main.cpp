#include <iostream>
#include "functions.h"

using namespace std;

int main()
{
    int N = 10;
    int row = 0;
    int column = 0;
    mat B = main_matrix(N);
    B(3,7) = 9;
    B(4,4) = 100;
    B.print();
    find_largest_element(B, &row, &column);
    cout << row << ", " << column << endl;
    return 0;
}

