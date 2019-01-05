# Matrix

Matrix header-only library (C++) from numerical optimization and machine learning.   
Ease of using is the emphasis of this library. Supported features are as following:
* templated type
* copy-free matrix-view
* flexiable matrix element access/selection
* variety of operators, including matrix-matrix operators and matrix-number operators
* broadcasting
* type casting
* linear algebra functions

Author: Liu Xiaodong (liuxiaodong008008@gmail.com) QQ:2410018191

Examples are listed in main.cpp
```
#include <iostream>
#include "matrix.h"
#include "matrix_operation.h"
using namespace std;

int main() {

    // declaration
    Matrix<double> m(3,3,0.0);
    cout<<"m "<<m<<endl;
    // m Matrix: 3x3
    // 0,      0,      0
    // 0,      0,      0
    // 0,      0,      0


    // initialization
    m<<1,2,3,4,5,6,7,8,9;
    cout<<"m "<<m<<endl;
    // m Matrix: 3x3
    // 1,      2,      3
    // 4,      5,      6
    // 7,      8,      9


    // change value (1)
    m<<11,12;
    cout<<"m "<<m<<endl;
    // m Matrix: 3x3
    // 11,     12,     3
    // 4,      5,      6
    // 7,      8,      9


    // change value (2)
    m(1,1)=15;
    cout<<"m "<<m<<endl;
    // m Matrix: 3x3
    // 11,     12,     3
    // 4,      15,     6
    // 7,      8,      9


    // change value (3)
    m(2)=13;
    cout<<"m "<<m<<endl;
    // m Matrix: 3x3
    // 11,     12,     13
    // 4,      15,     6
    // 7,      8,      9

    // traverse (1)
    for (int r=0;r<m.rows();r++) {
        for (int c=0;c<m.cols();c++) {
            m(r,c)/=2;
        }
    }
    cout<<"m "<<m<<endl;
    // m Matrix: 3x3
    // 5.5,    6,      6.5
    // 2,      7.5,    3
    // 3.5,    4,      4.5


    // traverse (2)
    for (int idx=0;idx<m.count();idx++) {
        m(idx) = (int) m(idx);
    }
    cout<<"m "<<m<<endl;
    cout<<"m "<<m<<endl;
    // m Matrix: 3x3
    // 5,      6,      6
    // 2,      7,      3
    // 3,      4,      4


    // traverse (3)
    auto it = m.rowwise().iterator();
    for (;it.valid();it.forward()) {
        cout<<"row "<<it.get()<<endl;
    }
    // row MatrixView: 1x3
    // 5,      6,      6
    //
    // row MatrixView: 1x3
    // 2,      7,      3
    //
    // row MatrixView: 1x3
    // 3,      4,      4


    // selection view (1)
    // a matrixview is a view of a matrix, no element copy
    auto v1 = m(Range(Point(0,0),Size(2,2)));
    cout<<"v1 "<<v1<<endl;
    // v2 MatrixView: 2x2
    // 5,      6
    // 2,      7


    // selection view (2)
    // view from view()
    auto v2 = m.view(Range(Point(0,0),Size(2,2)));
    cout<<"v2 "<<v2<<endl;
    // v2 MatrixView: 2x2
    // 5,      6
    // 2,      7

    // selection view (3)
    // row/col selection
    auto rowv1 = m.row(1);
    cout<<"rowv1 "<<rowv1<<endl;
    // rowv1 MatrixView: 1x3
    // 2,      7,      3

    // selection view (4)
    // view of view
    auto v2c2 = v2.col(2);
    cout<<"v2c2 "<<v2c2<<endl;
    // v2c2 MatrixView: 2x1
    // 6
    // 3

    // selection view (5)
    // view to matrix
    auto m1 = v2c2.copy();
    cout<<"m1 "<<m1<<endl;
    // m1 Matrix: 2x1
    // 6
    // 3

    // selection view (6)
    // change to view is change to the viewed matrix
    auto v3 = m.col(1);
    cout<<"m "<<m<<endl;
    v3(1) = 100;
    cout<<"m "<<m<<endl;
    // m Matrix: 3x3
    // 5,      6,      6
    // 2,      7,      3
    // 3,      4,      4
    //
    // m Matrix: 3x3
    // 5,      6,      6
    // 2,      100,    3
    // 3,      4,      4


    // matrix operations (1)
    // operators + , - , * , / , % , & , | , ^ , > , >=, < , <=, ==, !=, &&,
    // +=, -=, *=, /=, %=, &=, |=, ^=, ... almost all operators are supported
    Matrix<double> n(3,3,6);
    cout<<(m+n)<<endl;
    cout<<(m-n)<<endl;
    cout<<(m*n)<<endl;
    cout<<(m/n)<<endl;
    cout<<boolalpha<<(m>n)<<endl;
    cout<<boolalpha<<(m<=n)<<endl;
    // Matrix: 3x3
    // 11,     12,     12
    // 8,      106,    9
    // 9,      10,     10
    //
    // Matrix: 3x3
    // -1,     0,      0
    // -4,     94,     -3
    // -3,     -2,     -2
    //
    // Matrix: 3x3
    // 30,     36,     36
    // 12,     600,    18
    // 18,     24,     24
    //
    // Matrix: 3x3
    // 0.833333,       1,              1
    // 0.333333,       16.6667,        0.5
    // 0.5,            0.666667,       0.666667
    //
    // Matrix: 3x3
    // false,  false,  false
    // false,  true,   false
    // false,  false,  false
    //
    // Matrix: 3x3
    // true,   true,   true
    // true,   false,  true
    // true,   true,   true

    // matrix operations (2)
    // with numbers
    cout<<"m+1 "<<(m+1)<<endl;
    // m+1 Matrix: 3x3
    // 6,      7,      7
    // 3,      101,    4
    // 4,      5,      5

    // matrix operations (3)
    // copy
    auto l = m.copy();
    cout<<"l "<<l<<endl;
    // l Matrix: 3x3
    // 5,      6,      6
    // 2,      100,    3
    // 3,      4,      4

    // matrix operations (4)
    // type casting
    auto lb = (m-6).cast<bool>();
    cout<<"lb "<<lb<<endl;
    // lb Matrix: 3x3
    // true,   false,  false
    // true,   true,   true
    // true,   true,   true

    // matrix operations (5)
    // transpose
    cout<<"mT "<<m.t()<<endl;
    // mT Matrix: 3x3
    // 5,      2,      3
    // 6,      100,    4
    // 6,      3,      4

    // matrix operations (6)
    // broadcasting
    n = m.colwise()+m.col(0).copy();
    cout<<"n "<<n<<endl;
    // n Matrix: 3x3
    // 10,     11,     11
    // 4,      102,    5
    // 6,      7,      7

    // matrix operations (7)
    // inverse, matmul, any, all, sum, reduce ...
    cout<<"inverse(m) "<<inverse(m)<<endl;
    cout<<"matmul(m,m) "<<matmul(m,m)<<endl;
    cout<<"sum(m) "<<sum(m,REDUCE_DIRECTION::VERTICAL)<<endl;
    // inverse(m) Matrix: 3x3
    // 2,      0,      -3
    // 0.00515464,     0.0103093,      -0.0154639
    // -1.50515,       -0.0103093,     2.51546
    //
    // matmul(m,m) Matrix: 3x3
    // 55,     654,    72
    // 219,    10024,  324
    // 35,     434,    46
    //
    // sum(m) Matrix: 1x3
    // 10,     110,    13

    return 0;
}
```
