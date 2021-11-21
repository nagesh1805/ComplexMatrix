//
// Created by Nagesh Talagani on 21/11/21.
//

#ifndef COMPLEXMATRIX_MATRIX_H
#define COMPLEXMATRIX_MATRIX_H

#include <iostream>
#include <vector>
#include "ComplexNumber.h"
#include "ComplexVector.h"
using namespace std;

class Matrix{

public:
    Matrix();
    Matrix(int nr,int nc);
    Matrix(const vector<vector<complex>> &vec);
    Matrix(const vector<vector<double>> &vec);
    Matrix(const vector<complex> &vec);
    Matrix(const vector<double> &vec);
    Matrix(const ComplexVector &V);
    Matrix(const vector<ComplexVector> & vc);
    Matrix(const Matrix &m);
    //  Matrix(const Matrix &A11,const Matrix &A12, const Matrix &A21, const Matrix &A22);
    ~Matrix(){};

    static Matrix In(int n);
    static Matrix Pn(vector<int> &p);
    Matrix Eij(int i, int j);
    Matrix Eic(int i, const complex &c);
    Matrix Eijc(int i, int j, const complex &c);
    Matrix Fij(int i, int j);
    Matrix Fic(int i, const complex &c);
    Matrix Fijc(int i, int j, const complex &c);
    int get_nrows() const;
    int get_ncols() const;
    vector<ComplexVector> get_array() const;
    bool isZero() const;
    static Matrix zeros(int m, int n);
    static Matrix ones(int m, int n);
    ComplexVector getRow(int i) const;
    ComplexVector getColumn(int i) const;
    void setRow(int i, const ComplexVector v);
    complex get_ijth(int i, int j) const;
    void set_aij(int i, int j, const complex z);

    Matrix rref();
    int rank();
    Matrix transpose() const;
    Matrix conjugate() const;
    Matrix ctranspose() const;

    string toString() const;
    friend ostream& operator<<(ostream &os,const Matrix& A);
    Matrix operator-() const;
    Matrix operator+(const Matrix &B) const;
    Matrix operator-(const Matrix &B) const;
    Matrix operator*(const Matrix &B) const;
    ComplexVector operator[](int i) const;
    Matrix operator*(const complex &c) const;
    friend Matrix operator*(const complex &c,const Matrix &A);
    //friend Matrix operator*(Matrix &A,ComplexNumber &c);
    friend Matrix operator*(double c,const Matrix &A);
    Matrix operator*(double c) const;
    bool operator==(const Matrix &A) const;
    bool operator!=(const Matrix &A) const;
    friend ComplexVector operator*(const Matrix &A, const ComplexVector &v);
    friend ComplexVector operator*(const ComplexVector &v,const Matrix &A);

protected:
    //ComplexNumber **arr;
    vector<ComplexVector> arr;
    int nrow,ncol;

};




#endif //COMPLEXMATRIX_MATRIX_H
