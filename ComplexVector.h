//
// Created by Nagesh Talagani on 21/11/21.
//

#ifndef COMPLEXMATRIX_COMPLEXVECTOR_H
#define COMPLEXMATRIX_COMPLEXVECTOR_H

#include<cmath>
#include <vector>
#include "ComplexNumber.h"

class ComplexVector{


public:
    ComplexVector();
    ComplexVector(int size);
    ComplexVector(const vector<complex> &vec);
    ComplexVector(const vector<double> &vec);
    ComplexVector(const ComplexVector &mv);
    ComplexVector(const ComplexVector &v1, const ComplexVector &v2);
    ~ComplexVector();
    int size() const;
    void set_ith(int i, const complex &z);
    ComplexVector conjugate() const;
    complex get_ith(int i) const;
    ComplexVector point_wise_mul(const ComplexVector &v) const;
    static ComplexVector point_wise_mul(const ComplexVector &v1, const ComplexVector &v2);
    ComplexVector operator+(const ComplexVector &v) const;
    ComplexVector operator-(const ComplexVector &v) const;
    complex operator*(const ComplexVector &v) const;
    ComplexVector operator*(const complex &c) const;
    ComplexVector join(const ComplexVector &v) const;
    friend ComplexVector operator*(const complex &c, const ComplexVector &v);
    friend bool operator==(const ComplexVector &u, const ComplexVector &v);
    friend bool operator!=(const ComplexVector &u, const ComplexVector &v);
    ComplexVector operator-() const;
    complex& operator[](int k) const;
    complex inner_product(const ComplexVector &v) const;
    complex dot_product(const ComplexVector &v) const;
    static complex inner_product(const ComplexVector &v1, const ComplexVector &v2);
    static ComplexVector zeros(int n);
    static ComplexVector ones(int n);
    string toString() const;
    double norm() const;
    static double norm(const ComplexVector &v);
    vector<complex> to_vector() const;
    friend ostream& operator<<(ostream &os,const ComplexVector& v);
    bool isZero() const;
    bool isOrthogonalTo(const ComplexVector &v) const;
    ComplexVector operator=(const ComplexVector &v);


private:
    complex *ar;
    int len;

};



#endif //COMPLEXMATRIX_COMPLEXVECTOR_H
