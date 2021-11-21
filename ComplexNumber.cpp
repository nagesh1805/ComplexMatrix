//
// Created by Nagesh Talagani on 21/11/21.
//


#include<iostream>
#include <string>
#include <cmath>
#include"ComplexNumber.h"
using namespace std;

complex::complex(double r,double i){
    real_part=r,img_part=i;
}
complex::complex(int r, int i){
    real_part =r, img_part =i;
}
complex::complex(int r):complex(r,0){}
complex::complex():complex(0,0){}
complex::complex(double r):complex(r,0.0){}
complex::~complex(){}

string complex::to_myString(double a){

    string as = to_string(a);
    int ix =as.find(".");
    string ip = as.substr(0,ix);
    string fp = as.substr((ix+1));
    while(fp.back()=='0'){
        fp=fp.substr(0,fp.length()-1);
    }
    return (fp=="")?ip:ip+"."+fp;

}

vector<complex> complex::cardonas_helper(const complex &p, const complex &q)
{
    vector<complex> roots;
    vector<complex> sds = ((q^2)/4 +(p^3)/27)^rational(1,2);
    vector<complex> a = ((-q/2) +sds[0])^rational(1,3);
    vector<complex> b = ((-q/2) +sds[1])^rational(1,3);
    for(int i=0;i<a.size();i++){
        for(int j=0;j<b.size();j++){
            if ((a[i]*b[j]) == (-p/3)) {
                roots.push_back(a[i]+b[j]);
            }
        }
    }

    return roots;
}
string complex::toString() const{
    double rp=this->get_real();
    double ip=this->get_img();
    if(abs(rp)>0.0000001 and abs(ip)>0.0000001 and  ip > 0) return to_myString(rp)+" + "+((to_myString(ip)=="1")?"i": to_myString(ip)+"i");
    else if(abs(rp)>0.0000001 and abs(ip)>0.0000001 and ip< 0) return to_myString(rp)+((to_myString(ip)=="-1")?"-i": to_myString(ip)+"i");
    else if((abs(rp)< 0.0000001 and abs(ip)>0.0000001 and ip >0)) return ((to_myString(ip)=="1")?"i": to_myString(ip)+"i");
    else if(abs(rp)< 0.0000001 and abs(ip)>0.0000001 and ip <0) return ((to_myString(ip)=="-1")?"-i": to_myString(ip)+"i");
    else if(abs(rp) >0.0000001 and abs(ip)<0.0000001 ) return to_myString(rp);
    else return "0";
}

vector<complex> complex::Pow(int m, int n) const{
    complex z=(*this)^m;
    double r=z.Abs();
    double theta=z.Arg();
    vector<complex> nth_roots;
    for(int k=0;k<n;k++){
        complex zz(cos((theta+2*k*pi_val)/n),sin((theta+2*k*pi_val)/n));
        nth_roots.push_back(pow(r,1.0/n)*zz);
    }
    return nth_roots;
}

vector<complex> complex::nthRoots(int n) const{
    return Pow(1,n);
}

double complex::get_polar_theta() const{
    return atan2(img_part,real_part);
}

bool complex::equalsTo(const complex &z)
{
    return ((real_part==z.get_real())&&(img_part==z.get_img()));
}

complex complex::rotate(double theta) const{
    return ((*this)*cis(theta));
}

complex complex::conjugate() const{
    complex z(real_part,-img_part);
    return z;
}
double complex::get_real() const{
    return real_part;
}
double complex::get_img() const{
    return img_part;
}

bool complex::isReal(){
    return (img_part == 0);
}
bool complex::isPureImaginary(){
    return (real_part == 0);
}




double complex::Arg() const{
    if(real_part >0) return atan(img_part/real_part);
    else if(real_part<0 and img_part>=0) return  atan(img_part/real_part)+4*atan(1);
    else if(real_part<0 and img_part<0) return atan(img_part/real_part) - 4*atan(1);
    else if(real_part==0 and img_part >0) return 2*atan(1);
    else return -2*atan(1);
}
complex complex::cis(double t){
    return complex(cos(t),sin(t));
}


complex complex::Log(const complex &z){
    return complex(log(Abs(z)),z.Arg());
}

complex complex::Exp(const complex &z){
    return exp(z.get_real())*complex(cos(z.get_img()),sin(z.get_img()));
}

complex complex::Cos(const complex &z){
    complex i(0,1);
    return (Exp(i*z) - Exp(-i*z))/2;
}
complex complex::Sin(const complex &z){
    complex i(0,1);
    return (Exp(i*z) + Exp(-i*z))/2;
}
complex complex::Tan(const complex &z){
    return Sin(z)/Cos(z);
}
complex complex::Cot(const complex &z){
    return Sin(z)/Cos(z);
}
complex complex::Cosec(const complex &z){
    return 1/Sin(z);
}
complex complex::Sec(const complex &z){
    return 1/Cos(z);
}

vector<complex> complex::solve_qudratic_eqtn(const complex &a, const complex &b, const complex &c)
{
    vector<complex> roots;
    complex d = (b^2) - 4*a*c;
    vector<complex> sds=d^rational(1,2);
    return vector<complex>({(-b+sds[0])/2,(-b+sds[1])/2});
}

vector<complex> complex::solve_cubic_eqtn(const complex &a, const complex &b, const complex &c, const complex &d)
{
    /*                   x=(y-b/3a)
     *  ax^3+bx^2+cx+d=0  ------------->  Ay^3+Cy+D=0   ----> y^3+py+q=0   *
     *
     * */
    vector<complex> t_roots=complex::cardonas_helper((1/a)*(c-((b^2)/(3*a))),(1/a)*(d + (2*(b^3))/(27*(a^2)) - (b*c)/(3*a)));
    vector<complex> roots;
    for(int i=0;i<3;i++){
        roots.push_back(t_roots[i] - (b/(3*a)));
    }
    return roots;

}
double complex::Abs(const complex &z){
    return sqrt(z.real_part*z.real_part + z.img_part*z.img_part);
}

double complex::Abs() const{
    return sqrt(real_part*real_part + img_part*img_part);
}

ostream& operator<<(ostream &os, const complex &z){
    return os<<z.toString();
}
complex operator+(const complex &z1, const complex &z2){
    return complex(z1.real_part+z2.real_part,z1.img_part+z2.img_part);
}

complex operator-(const complex &z1, const complex &z2){
    return z1+(-z2);
}

complex operator*(const complex &z1, const complex &z2){
    return complex(z1.real_part*z2.real_part - z1.img_part*z2.img_part,z1.real_part*z2.img_part + z1.img_part*z2.real_part);
}

complex operator/(const complex &z1, const complex &z2){
    return z1*(z2.conjugate())*(1/(z2.Abs()*z2.Abs()));
}

complex operator^(const complex &z1, const complex &z2){
    if (z1==0) return 0;
    return complex::Exp(z2*(complex::Log(z1)));
}

complex operator-(const complex &z1){
    return complex(-z1.get_real(),-z1.get_img());
}
bool operator==(const complex &z1, const complex &z2){
    return (complex::Abs(z1-z2) < 0.0000001);
}

bool operator!=(const complex &z1, const complex &z2){
    return not(z1==z2);
}

void operator+=(complex &z1, const complex &z2){
    z1=z1+z2;
}

void operator-=(complex &z1, const complex &z2){
    z1=z1-z2;
}

void operator*=(complex &z1, const complex &z2){
    z1=z1*z2;
}

void operator/=(complex &z1, const complex &z2){
    z1=z1/z2;
}
bool operator<(const complex z1, const complex z2){
    return z1.Abs()<z2.Abs();
}
bool operator>(const complex z1, const complex z2){
    return z1.Abs()>z2.Abs();
}
bool operator<=(const complex z1, const complex z2){
    return not(z1.Abs()>z2.Abs());
}
bool operator>=(const complex z1, const complex z2){
    return not(z1.Abs()<z2.Abs());
}

vector<complex> complex::operator^(const rational &q) const{
    return Pow(q.get_numerator(),q.get_denominator());
}


