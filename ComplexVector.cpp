
#include "ComplexVector.h"

using namespace std;

ComplexVector::ComplexVector(){
    ComplexVector(1);
}
ComplexVector::ComplexVector(int size){
    len=size;
    ar = new complex [len]();

}
ComplexVector::ComplexVector(const vector<complex> &vec){
    len=vec.size();
    ar = new complex [len]();
    for(int i=0;i<len;i++){
        ar[i]=vec[i];
    }
}
ComplexVector::ComplexVector(const vector<double> &vec){
    len=vec.size();
    ar = new complex [len]();
    for(int i=0;i<len;i++){
        ar[i]=vec[i];
    }
}

ComplexVector::ComplexVector(const ComplexVector &mv){
    len=mv.size();
    ar = new complex [len];

    for(int i=0;i<len;i++){
        ar[i]=mv[i];
    }
}

ComplexVector::ComplexVector(const ComplexVector &v1, const ComplexVector &v2){
    len=v1.size()+v2.size();

    ar = new complex [len];
    for(int i=0;i<len;i++){
        if(i<v1.size()){
            ar[i]=v1[i];
        }
        else{
            ar[i]=v2[i-v1.size()];
        }
    }
}
ComplexVector::~ComplexVector(){
    delete[] ar;
}

int ComplexVector::size() const{
    return len;
}

ComplexVector ComplexVector::conjugate() const{
    ComplexVector v(this->size());
    for(int i=0;i<v.size();i++){
        v[i]=(*this)[i].conjugate();
    }
    return v;
}
void ComplexVector::set_ith(int i, const complex &z){
    ar[i]=z;
}
complex ComplexVector::get_ith(int i) const{
    return ar[i];
}

double ComplexVector::norm() const{
    complex z;
    z=this->inner_product(*this);
    return sqrt(complex::Abs(z));
}
double ComplexVector::norm(const ComplexVector &v){
    return v.norm();
}
complex ComplexVector::inner_product(const ComplexVector &v) const{
    complex s=0;
    for(int i=0;i<v.size();i++){
        s=s+((*this)[i]*(v[i].conjugate()));
    }
    return s;
}

complex ComplexVector::dot_product(const ComplexVector &v) const{
    complex s=0;
    for(int i=0;i<v.size();i++){
        s=s+((*this)[i]*(v[i]));
    }
    return s;
}

complex inner_product(const ComplexVector &v1, const ComplexVector &v2){
    return v1.inner_product(v2);
}

ComplexVector ComplexVector::operator+(const ComplexVector &v) const{
    ComplexVector vv(v.size());
    for(int i=0;i<v.size();i++){
        vv[i]=(*this)[i]+v[i];
    }
    return vv;
}

ComplexVector ComplexVector::operator-(const ComplexVector &v) const{
    return (*this) + (-v);
}
ComplexVector ComplexVector::operator-() const{
    return (-1)*(*this);
}

complex ComplexVector::operator*(const ComplexVector &v) const{
    return dot_product(v);
}
string ComplexVector::toString() const{
    string s="[ ";
    for(int i=0;i<this->size();i++){
        if(i!=(this->size()-1)) s=s+(*this)[i].toString()+",";
        else s=s+(*this)[i].toString()+" ]";
    }
    return s;
}
ostream& operator<<(ostream &os,const ComplexVector& v){
    return os<<v.toString();
}

vector<complex> ComplexVector::to_vector() const{
    vector<complex> v;
    for(int i=0;i<this->len;i++){
        v.push_back((*this)[i]);
    }
    return v;
}


complex& ComplexVector::operator[](int k) const{
    return ar[k];
}

ComplexVector ComplexVector::point_wise_mul(const ComplexVector &v) const{
    ComplexVector vv(v.size());
    for(int i=0;i<v.size();i++){
        vv[i]=(*this)[i]*v[i];
    }
    return vv;
}

ComplexVector ComplexVector::join(const ComplexVector &v) const{
    return ComplexVector(*this,v);
}

ComplexVector ComplexVector::operator*(const complex &c) const{
    ComplexVector v(this->size());
    for (int i=0;i<this->size();i++) {
        v[i]=c*(*this)[i];
    }
    return v;
}

ComplexVector ComplexVector::point_wise_mul(const ComplexVector &v1, const ComplexVector &v2) {
    return v1.point_wise_mul(v2);
}

bool ComplexVector::isZero() const{
    for(int i=0;i<this->size();i++){
        complex z=(*this)[i];
        if(z != 0){
            return false;
        }
    }
    return true;
}

bool ComplexVector::isOrthogonalTo(const ComplexVector &v) const{
    return (this->inner_product(v) == 0);
}

ComplexVector ComplexVector::operator=(const ComplexVector &v){
    len=v.size();
    ar = new complex [len];

    for(int i=0;i<len;i++){
        ar[i]=v[i];
    }
    return *this;
}

ComplexVector operator*(const complex &c, const ComplexVector &v){
    ComplexVector vv(v.size());
    for (int i=0;i<vv.size();i++) {
        vv[i]=c*v[i];
    }
    return vv;
}

bool operator==(const ComplexVector &u, const ComplexVector &v){
    if (u.size() != v.size()) return  false;
    else{
        for(int i=0;i<u.size();i++){
            if (u[i]!=v[i]) return  false;
        }
        return true;
    }
}
bool operator!=(const ComplexVector &u, const ComplexVector &v){
    return not(u==v);
}

ComplexVector ComplexVector::zeros(int n){
    return ComplexVector(n);
}
ComplexVector ComplexVector::ones(int n){
    vector<complex> v(n,1);
    return ComplexVector(v);
}


