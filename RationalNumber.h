//
// Created by Nagesh Talagani on 21/11/21.
//

#ifndef COMPLEXMATRIX_RATIONALNUMBER_H
#define COMPLEXMATRIX_RATIONALNUMBER_H

#include<iostream>
#include <string>
#include <cmath>
#include <string>
using namespace std;
class rational{

public:
    rational();
    //RationalNumber(int n);
    rational(int n, int d);
    rational(const rational &q);
    ~rational(){};
    int get_numerator() const;
    int get_denominator() const;
    double to_double() const;
    string toString() const;
    friend rational operator-(const rational &q);
    friend rational operator+(const rational &q1, const rational &q2);
    friend rational operator-(const rational &q1, const rational &q2);
    friend rational operator*(const rational &q1, const rational &q2);
    friend rational operator/(const rational &q1, const rational &q2);
    friend ostream& operator<<(ostream &os, rational &q);
private:
    int num,den;
};



#endif //COMPLEXMATRIX_RATIONALNUMBER_H
