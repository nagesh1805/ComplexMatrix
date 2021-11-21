


//
// Created by Nagesh Talagani on 21/11/21.
//

#ifndef COMPLEXMATRIX_COMPLEXNUMBER_H
#define COMPLEXMATRIX_COMPLEXNUMBER_H

#include <iostream>
#include <vector>
#include <string>
#include<cmath>
#include "RationalNumber.h"
const double pi_val = std::atan(1.0)*4;
const int Max_Iterations=255;
using namespace std;
class complex{


public:
    complex();
    complex(double r);
    complex(int r);
    complex(double r,double i);
    complex(int r,int i);
    ~complex();

    complex conjugate() const;
    complex rotate(double thrta) const;
    double get_real() const;
    double get_img() const;
    static double Abs(const complex &z);
    double Abs() const;
    double Arg() const;
    double get_polar_theta()const;
    bool equalsTo(const complex &z);
    bool isReal();
    bool isPureImaginary();
    vector<complex> nthRoots(int n) const;
    vector<complex> Pow(int m, int n) const;
    static complex cis(double t);
    static complex Log(const complex &z);
    static complex Exp(const complex &z);
    static complex Cos(const complex &z);
    static complex Sin(const complex &z);
    static complex Tan(const complex &z);
    static complex Cot(const complex &z);
    static complex Cosec(const complex &z);
    static complex Sec(const complex &z);
    static vector<complex> solve_qudratic_eqtn(const complex &a,const complex &b, const complex &c);
    static vector<complex> solve_cubic_eqtn(const complex &a,const complex &b, const complex &c,const complex &d);

    string toString() const;

    friend ostream& operator<<(ostream &os, const complex &z);
    friend complex operator-(const complex &z1);

    friend complex operator+(const complex &z1, const complex &z2);

    friend complex operator-(const complex &z1, const complex &z2);

    friend complex operator*(const complex &z1, const complex &z2);

    friend complex operator/(const complex &z1, const complex &z2);

    friend complex operator^(const complex &z1, const complex &z2);
    vector<complex> operator^(const rational &q) const;

    friend void operator+=(complex &z1, const complex &z2);

    friend void operator-=(complex &z1, const complex &z2);

    friend void operator*=(complex &z1, const complex &z2);

    friend void operator/=(complex &z1, const complex &z2);

    friend bool operator<(const complex z1,const complex z2);
    friend bool operator<=(const complex z1, const complex z2);
    friend bool operator>(const complex z1,const  complex z2);
    friend bool operator>=(const complex z1, const complex z2);

    friend bool operator==(const complex &z1, const complex &z2);

    friend bool operator!=(const complex &z1, const complex &z2);

private:
    double real_part, img_part;
    static string to_myString(double a);
    static vector<complex> cardonas_helper(const complex &p, const complex &q);
};





#endif //COMPLEXMATRIX_COMPLEXNUMBER_H
