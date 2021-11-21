//
// Created by Nagesh Talagani on 21/11/21.
//

#include "SquareMatrix.h"

SquareMatrix::SquareMatrix():Matrix(1,1){
    this->n=1;
}

SquareMatrix::SquareMatrix(int size):Matrix(size,size){
    this->n=size;
}

SquareMatrix::SquareMatrix(vector<vector<complex>> &vec):Matrix(vec){
    this->n=vec.size();
}
SquareMatrix::SquareMatrix(vector<vector<double>> &vec):Matrix(vec){
    this->n=vec.size();
}

SquareMatrix::SquareMatrix(const Matrix &m):Matrix(m){
    this->n=m.get_nrows();
}

SquareMatrix::SquareMatrix(const SquareMatrix &m){
    *this=m;
}

/*SquareMatrix::SquareMatrix(const SquareMatrix &A11,const SquareMatrix &A12, const SquareMatrix &A21, const SquareMatrix &A22){
    nrow=A11.get_nrows()+A21.get_nrows();
    ncol=A11.get_ncols()+A12.get_ncols();
    n=nrow;
    arr = new ComplexNumber* [n];
    for(int i=0;i<nrow;i++){
        arr[i]= new ComplexNumber[n]();
    }

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i<A11.get_nrows() and j<A11.get_ncols()){
                arr[i][j]=A11[i][j];
            }
            else if(i<A11.get_nrows() and j>=A11.get_ncols()){
                arr[i][j]=A12[i][j-A11.get_ncols()];
            }
            else if(i>=A11.get_nrows() and j<A11.get_ncols()){
                arr[i][j]=A21[i-A11.get_nrows()][j];
            }
            else{
                arr[i][j]=A22[i-A11.get_nrows()][j-A11.get_ncols()];
            }
        }
    }

}*/
/*SquareMatrix::~SquareMatrix(){
    for(int i=0;i<n;i++){
        delete[] arr[i];
    }
}*/
int SquareMatrix::size() const{
    return n;
}

bool SquareMatrix::isSymmetric() const{
    return (*this == ((*this).transpose()));
}

bool SquareMatrix::isSkew_Symmetric() const{
    return (*this == -((*this).transpose()));
}
bool SquareMatrix::isHermition() const{
    return (*this == ((*this).ctranspose()));
}
bool SquareMatrix::isSkew_Hermition() const{
    return (*this == -((*this).ctranspose()));
}

bool SquareMatrix::isOrthoganal() const{
    return ((*this)*(this->transpose()) == In(this->size()));
}
bool SquareMatrix::isUnitary() const{
    return ((*this)*(this->ctranspose()) == In(this->size())); // A*Abar^t=I
}

complex SquareMatrix::trace() const{
    complex s=0;
    for(int i=0;i<n;i++){
        s+=(*this)[i][i];
    }
    return s;
}


complex SquareMatrix::det() const{
    SquareMatrix A=*this;
    complex d=1;
    for(int j=0;j<n;j++){
        if(A[j][j]==0){
            bool colum_j_is_zero=true;
            for(int k=j+1;k<n;k++){
                if(A[k][j]!=0){
                    colum_j_is_zero=false;
                    A.Eij(k,j);
                    d=d*(-1);
                    break;
                }
            }
            if(colum_j_is_zero) return 0.0;
        }
        for(int i=j+1;i<n;i++){
            A.Eijc(i,j,-A[i][j]/A[j][j]);
        }
    }
    for (int i=0;i<n;i++){
        d=d*A[i][i];
    }
    return d;
}

SquareMatrix SquareMatrix::inverse() const{
    SquareMatrix A=*this;
    SquareMatrix A_inv=In(A.size());
    for(int j=0;j<this->size();j++){
        if(A[j][j]==0){
            for(int k=j+1;k<this->size();k++){
                if(A[k][j]!=0){
                    A.Eij(k,j);
                    A_inv.Eij(k,j);
                    break;
                }
            }

        }
        if(A[j][j] != 0){
            complex ajj=A[j][j];
            A.Eic(j,1/ajj);
            A_inv.Eic(j,1/ajj);
            for(int i=0;i<this->size();i++){
                if(i!=j){
                    complex aij=A[i][j];
                    A.Eijc(i,j,-aij);
                    A_inv.Eijc(i,j,-aij);
                }

            }
        } // else raise exception

    }
    return A_inv;
}

SquareMatrix SquareMatrix::operator^(int n) const{
    if(n==0){
        return In(n);
    }
    else if(n==1){
        return (*this);
    }
    else if (n>0){
        SquareMatrix A=(*this)^(n/2);
        return (n%2==0)?A*A:A*A*(*this);
    }
    else{
        return (this->inverse())^(-n);
    }
}

SquareMatrix SquareMatrix::zeros(int n){
    return SquareMatrix(n);
}
SquareMatrix SquareMatrix::ones(int n){
    return Matrix::ones(n,n);
}

SquareMatrix SquareMatrix::D(const  ComplexVector &v){
    SquareMatrix d(v.size());
    for(int i=0;i<v.size();i++){
        d.set_aij(i,i,v[i]);
    }
    return d;
}

SquareMatrix SquareMatrix::D(int n){
    ComplexVector vd(n);
    complex wn=complex::cis(-2*pi_val/(2*n));
    for(int i=0;i<n;i++){
        vd[i]=wn^i;
    }
    return D(vd);
}

ComplexVector SquareMatrix::Diagonal(const SquareMatrix &A){
    ComplexVector v(A.size());
    for(int i=0;i<A.size();i++){
        v.set_ith(i,A[i][i]);
    }
    return v;
}

SquareMatrix SquareMatrix::DFT(int n){
    SquareMatrix dft(n);
    complex wn=complex::cis(-2*pi_val/n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            dft.set_aij(i,j,wn^(i*j));
        }
    }
    return dft;
}



ComplexVector SquareMatrix::dft(const ComplexVector &f){
    return DFT(f.size())*f;
}

vector<complex> SquareMatrix::dft(vector<complex> &f){
    return dft(ComplexVector(f)).to_vector();
}

ComplexVector SquareMatrix::fft_helper(const ComplexVector &f){
    int n=f.size();
    int m=n/2;
    if(n==1){
        return DFT(1)*f;
    }
    else{
        ComplexVector f_e(m),f_o(m);
        for(int i=0;i<n;i++){
            if(i%2==0){
                f_e.set_ith(i/2,f[i]);
            }
            else{
                f_o.set_ith((i-1)/2,f[i]);
            }
        }

        ComplexVector Dm=Diagonal(D(m));   // m=n/2

        ComplexVector ff_0=fft_helper(f_e);
        ComplexVector ff_1=fft_helper(f_o);

        ComplexVector Df=Dm.point_wise_mul(ff_1);

        ComplexVector ff(ff_0+Df,ff_0-Df);

        return ff;
    }

}

ComplexVector SquareMatrix::fft(const ComplexVector &f){
    int n=pow(2,ceil(log2(f.size())));
    ComplexVector zero_part(n-f.size());
    ComplexVector ff=f.join(zero_part);
    return fft_helper(ff);
}

vector<complex> SquareMatrix::fft(vector<complex> &f){
    return fft(ComplexVector(f)).to_vector();
}




