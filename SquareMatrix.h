//
// Created by Nagesh Talagani on 21/11/21.
//

#ifndef COMPLEXMATRIX_SQUAREMATRIX_H
#define COMPLEXMATRIX_SQUAREMATRIX_H

#include "Matrix.h"
#include "ComplexNumber.h"

class SquareMatrix:public Matrix{

public:
    SquareMatrix();
    SquareMatrix(int size);
    SquareMatrix(vector<vector<complex>> &vec);
    SquareMatrix(vector<vector<double>> &vec);
    SquareMatrix(const Matrix &m);
    SquareMatrix(const SquareMatrix &m);
    //  SquareMatrix(const SquareMatrix &A11,const SquareMatrix &A12, const SquareMatrix &A21, const SquareMatrix &A22);
    ~SquareMatrix(){};
    SquareMatrix operator^(int n) const;
    int size() const;
    bool isSymmetric() const;
    bool isSkew_Symmetric() const;
    bool isHermition() const;
    bool isSkew_Hermition() const;
    bool isOrthoganal() const;
    bool isUnitary() const;
    complex det() const;
    complex trace() const;
    SquareMatrix inverse() const;
    static SquareMatrix zeros(int n);
    static SquareMatrix ones(int n);
    static SquareMatrix D(const ComplexVector &v);
    static ComplexVector Diagonal( const SquareMatrix &A);
    static SquareMatrix D(int n);
    static SquareMatrix DFT(int n);
    //static SquareMatrix FFT(int n);
    static ComplexVector dft(const ComplexVector &f);
    static ComplexVector fft(const ComplexVector &f);
    static vector<complex> dft(vector<complex> &f);
    static vector<complex> fft(vector<complex> &f);
private:
    static ComplexVector fft_helper(const ComplexVector &f);
    int n;

};




#endif //COMPLEXMATRIX_SQUAREMATRIX_H
