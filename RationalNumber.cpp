//
// Created by Nagesh Talagani on 21/11/21.
//

#include "RationalNumber.h"

using namespace std;


//RationalNumber::RationalNumber(int n):RationalNumber(n,1){}

rational::rational(int n, int d){
    num=n,den=d;
}

int rational::get_numerator() const{
    return num;
}

int rational::get_denominator() const{
    return den;
}
double rational::to_double() const{
    return (double(num)/den);
}
rational::rational(const rational &q){
    num=q.num;
    den=q.den;
}
string rational::toString() const{
    return ""+to_string(num)+"/"+to_string(den);
}

ostream& operator<<(ostream &os, rational &q){
    return os<<q.toString();
}

rational operator-(const rational &q){
    return rational(-q.num,-q.den);
}
rational operator+(const rational &q1, const rational &q2){
    return rational(q1.num*q2.den + q1.den*q2.num,q1.den*q2.den);
}
rational operator-(const rational &q1, const rational &q2){
    return q1 + (-q2);
}
rational operator*(const rational &q1, const rational &q2){
    return rational(q1.num*q2.num,q1.den*q2.den);
}
rational operator/(const rational &q1, const rational &q2){
    return rational(q1.num*q2.den,q1.den*q2.num);
}
