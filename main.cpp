#include <iostream>
#include "ComplexNumber.h"
#include "SquareMatrix.h"
#include "ComplexVector.h"
#include <map>
#include "RationalNumber.h"


using namespace std;
const double pi=4*atan(1); // pi value

const complex i(0,1); // imaginary number i

int main() {
    vector<vector<complex>> vv= {
            {1+i,2-3*i},
            {3*i ,4  }
    };

    SquareMatrix A=SquareMatrix(vv);
    cout<<"A="<<endl<<A<<endl;
    cout<<"------"<<endl;
    cout<<"A^-1="<<endl <<A.inverse()<<endl;
    cout<<"det(A) = "<<endl;
    cout<<A.det()<<endl;

    complex a(1,3),b(2,5);
    cout<<"a = "<<a<<endl;
    cout<<"b = "<<b<<endl;
    cout<<"a+b = "<<a+b<<endl;
    cout<<"a*b = "<<a*b<<endl;
    cout<<"a/b = "<<a/b<<endl;
    cout<<"cis(pi/4)= "<<complex::cis(pi/4)<<endl;
    complex z(1,2);
    cout<<"Log(z) = "<<complex::Log(z)<<endl;
    cout<<"Exp(z) = "<<complex::Exp(z)<<endl;
    cout<<"Sin(z) = "<<complex::Sin(z)<<endl;
    cout<<"Arg(z) = "<<z.Arg()<<endl;
    complex one(1);
    cout<<"cube roots of unity : "<< ComplexVector(one.nthRoots(3))<<endl;
    cout<<"fourth roots of unity : "<< ComplexVector(one.nthRoots(4))<<endl;


    rational q(3,2);
    cout<<ComplexVector(complex::solve_qudratic_eqtn(1,1,1))<<endl;

    cout<<ComplexVector(complex::solve_cubic_eqtn(i,0,0,1))<<endl;
    cout<<SquareMatrix::fft(ComplexVector(vector<complex>({2,3,5,7})));
    ComplexVector v(vector<complex>({1,2,3}));
    cout<<v<<endl;


    return 0;
}

/*
Output::

A=
[ 1 + i,2-3i ]
[ 3i,4 ]

------
A^-1=
[ -0.689655 + 0.275862i,0.137931-0.655172i ]
[ 0.206897 + 0.517241i,-0.241379-0.103448i ]

det(A) =
-5-2i
a = 1 + 3i
b = 2 + 5i
a+b = 3 + 8i
a*b = -13 + 11i
a/b = 0.586207 + 0.034483i
cis(pi/4)= 0.707107 + 0.707107i
Log(z) = 0.804719 + 1.107149i
Exp(z) = -1.131204 + 2.471727i
Sin(z) = 2.032723-3.051898i
Arg(z) = 1.10715
cube roots of unity : [ 1,-0.5 + 0.866025i,-0.5-0.866025i ]
fourth roots of unity : [ 1,i,-1,-i ]
[ -0.5 + 0.866025i,-0.5-0.866025i ]
[ 0.866025 + 0.5i,-0.866025 + 0.5i,-i ]
[ 17,-3 + 4i,-3,-3-4i ][ 1,2,3 ]


*/





