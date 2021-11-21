//
// Created by Nagesh Talagani on 21/11/21.
//

#include "Matrix.h"

Matrix::Matrix(){
    Matrix(1,1);
}

Matrix::Matrix(int nr,int nc){
    nrow=nr;
    ncol=nc;
    for(int i=0;i<nrow;i++){
        ComplexVector v(ncol);
        arr.push_back(v);
    }
}

Matrix::Matrix(const vector<vector<complex>> &vec){
    nrow=vec.size();
    ncol=vec[0].size();
    for(int i=1;i<nrow;i++){
        if(vec[i].size()>ncol) ncol=vec[i].size();
    }
    for(int i=0;i<nrow;i++){
        ComplexVector v(vec[i]);
        arr.push_back(v);
    }

}

Matrix::Matrix(const vector<complex> &vec){
    nrow=1;
    ncol=vec.size();
    ComplexVector v(vec);
    arr.push_back(v);
}

Matrix::Matrix(const ComplexVector &V):Matrix(V.to_vector()){
}

Matrix::Matrix(const vector<vector<double>> &vec){
    nrow=vec.size();
    ncol=vec[0].size();
    for(int i=1;i<nrow;i++){
        if(vec[i].size()>ncol) ncol=vec[i].size();
    }
    for(int i=0;i<nrow;i++){
        ComplexVector v(vec[i]);
        arr.push_back(v);
    }
}

Matrix::Matrix(const vector<double> &vec){
    nrow=1;
    ncol=vec.size();
    ComplexVector v(vec);
    arr.push_back(v);
}
Matrix::Matrix(const vector<ComplexVector> & vc){
    nrow = vc.size();
    ncol = vc[0].size();
    arr=vc;
}

Matrix::Matrix(const Matrix &m){
    *this=m;
}

/*Matrix::Matrix(const Matrix &A11,const Matrix &A12, const Matrix &A21, const Matrix &A22){
    nrow=A11.get_nrows()+A21.get_nrows();
    ncol=A11.get_ncols()+A12.get_ncols();
    arr = new ComplexNumber* [nrow];
    for(int i=0;i<nrow;i++){
        arr[i]= new ComplexNumber[ncol]();
    }
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            if(i<A11.get_nrows() and j<A11.get_ncols()){
                arr[i][j]=A11.get_ijth(i,j);
            }
            else if(i<A11.get_nrows() and j>=A11.get_ncols()){
                arr[i][j]=A12.get_ijth(i,j-A11.get_ncols());
            }
            else if(i>=A11.get_nrows() and j<A11.get_ncols()){
                arr[i][j]=A12.get_ijth(i-A11.get_nrows(),A11.get_ncols());
            }
            else{
                arr[i][j]=A22.get_ijth(i-A11.get_nrows(),j-A11.get_ncols());
            }
        }
    }

}*/
/*Matrix::~Matrix(){
    for(int i=0;i<nrow;i++){
        delete[] arr[i];
    }
    delete[] arr;
}*/
Matrix Matrix::In(int n){
    Matrix I(n,n);
    for(int i=0;i<n;i++){
        I.set_aij(i,i,1);
    }
    return I;
}

Matrix Matrix::Pn(vector<int> &p){
    Matrix P_n(p.size()),I_n=In(p.size());
    for(int i=0;i<p.size();i++){
        P_n.setRow(i,I_n.getRow(p[i]));
    }
    return P_n;
}

string Matrix::toString() const{
    string s="";
    for(int i=0;i<nrow;i++){
        s += (arr[i]).toString()+"\n";
    }
    return  s;
}

int Matrix::get_nrows() const{
    return nrow;
}

int Matrix::get_ncols() const{
    return ncol;
}

vector<ComplexVector> Matrix::get_array() const{
    return arr;
}

complex Matrix::get_ijth(int i, int j) const{
    return arr[i][j];
}
void Matrix::set_aij(int i, int j, complex z){
    arr[i][j]=z;
}

Matrix Matrix::Eij(int i,int j){
    ComplexVector Row_i =this->getRow(i);
    arr[i]=arr[j];
    arr[j]=Row_i;
    return *this;
}

Matrix Matrix::Eic(int i,const complex &c){
    arr[i]=c*arr[i];
    return *this;
}
Matrix Matrix::Eijc(int i, int j, const complex &c){
    arr[i]=arr[i]+c*arr[j];
    return *this;
}

Matrix Matrix::Fij(int i,int j){
    *this=this->transpose().Eij(i,j).transpose();
    return *this;

}

Matrix Matrix::Fic(int i,const complex &c){
    *this=this->transpose().Eic(i,c).transpose();
    return *this;
}
Matrix Matrix::Fijc(int i, int j,const complex &c){
    *this=this->transpose().Eijc(i,j,c).transpose();
    return *this;
}

Matrix Matrix::transpose() const{
    Matrix At(this->get_ncols(),this->get_nrows());
    for(int i=0;i<At.get_nrows();i++){
        for(int j=0;j<At.get_ncols();j++){
            At.set_aij(i,j,arr[j][i]);
        }
    }
    return At;
}

Matrix Matrix::conjugate() const{
    Matrix Abar(this->get_nrows(),this->get_ncols());
    for(int i=0;i<Abar.get_nrows();i++){
        //   Abar[i]=arr[i].conjugate();
        Abar.setRow(i,arr[i].conjugate());
    }
    return Abar;
}
Matrix Matrix::ctranspose() const{
    return (this->conjugate()).transpose();
}
ComplexVector Matrix::getRow(int i) const{
    return arr[i];

}
ComplexVector Matrix::getColumn(int i) const{
    return (this->transpose()).getRow(i);
}
void Matrix::setRow(int i, const ComplexVector v){
    arr[i]=v;
}

bool Matrix::isZero() const{
    for(int i=0;i<this->get_nrows();i++){
        if(not (*this)[i].isZero()) return false;
    }
    return true;
}

int Matrix::rank(){
    Matrix R=this->rref();
    int r=0;
    for(int k=0;k<R.get_nrows();k++){
        if(not R.getRow(k).isZero()) r++;
    }
    return r;
}

Matrix Matrix::rref(){
    Matrix R=*this;
    for(int j=0;j<R.get_nrows();j++){
        if(R[j][j]==0){
            for(int k=j+1;k<R.get_nrows();k++){
                if(R[k][j]!=0){
                    R.Eij(k,j);
                    break;
                }
            }
        }
        if(R[j][j]!=0){
            R.Eic(j,1.0/(R[j][j]));
            for(int i=0;i<R.get_nrows();i++){
                if(i!=j){
                    R.Eijc(i,j,-R[i][j]);
                }
            }
        }
    }
    return R;
}


ostream& operator<<(ostream &os,const Matrix &A){
    return os<<A.toString();
}

Matrix Matrix::operator-() const{
    return ((-1)*(*this));
}
Matrix Matrix::operator+(const Matrix &B) const{
    Matrix C(B.get_nrows(),B.get_ncols());
    for(int i=0;i<B.get_nrows();i++){
        C.setRow(i,arr[i]+B[i]);
    }
    return C;
}
Matrix Matrix::operator-(const Matrix &B) const{
    return (*this)+(-B);
}
Matrix Matrix::operator*(const Matrix &B) const{
    Matrix C(this->get_nrows(),B.get_ncols());
    for(int i=0;i<this->get_nrows();i++){
        for(int j=0;j<B.get_ncols();j++){
            C.set_aij(i,j,arr[i]*B[j]);
        }
    }
    return C;

}

ComplexVector Matrix::operator[](int i) const{
    return arr[i];
}

Matrix operator*(const complex &c,const Matrix &A){
    Matrix C(A.get_nrows(),A.get_ncols());
    for(int i=0;i<A.get_nrows();i++){
        C.setRow(i,c*A[i]);
    }
    return C;
}
Matrix Matrix::operator*(const complex &c)  const{
    return c*(*this);
}

Matrix operator*(double c,const Matrix &A){
    return  complex(c)*A;
}
Matrix Matrix::operator*(double c) const{
    return complex(c)*(*this);
}
bool Matrix::operator==(const Matrix &A) const{
    if(this->get_nrows()!=A.get_nrows() or this->get_ncols() !=A.get_ncols()) return false;
    else{
        for(int i=0;i<this->get_nrows();i++){
            if (arr[i] != A[i]) return false;
        }
        return true;
    }
}

bool Matrix::operator!=(const Matrix &A) const{
    return not((*this)==A);
}

ComplexVector operator*(const Matrix &A, const ComplexVector &v){
    ComplexVector vv(v.size());
    for (int i=0;i<v.size();i++){
        vv.set_ith(i,A[i]*v);
    }
    return vv;
}
ComplexVector operator*(const ComplexVector &v,const Matrix &A){
    return A.transpose()*v;
}

Matrix Matrix::zeros(int m, int n){
    return Matrix(m,n);
}
Matrix Matrix::ones(int m, int n){
    vector<ComplexVector> vc(n,ComplexVector::ones(n));
    return Matrix(vc);
}

