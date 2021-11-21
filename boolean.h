//
// Created by Nagesh Talagani on 21/11/21.
//

#ifndef COMPLEXMATRIX_BOOLEAN_H
#define COMPLEXMATRIX_BOOLEAN_H

#include<iostream>
#include<string>

using namespace std;
struct boolean{
    bool flag;
    boolean(){
        flag=false;
    }
    boolean(const bool a){
        flag =a;
    }
    boolean(const boolean &a){
        *this=a;
    }
    boolean operator=(const boolean a){
        flag=a.flag;
        return *this;
    }
    string toString() const{
        return (flag)?"true":"false";
    }
    friend ostream& operator<<(ostream &os,const boolean &a) {
        return os<<a.toString();
    }

    friend boolean operator&&(const boolean &a, const boolean &b){
        return boolean(a.flag&&b.flag);
    }
    friend boolean operator||(const boolean &a, const boolean &b){
        return boolean(a.flag||b.flag);
    }
    friend boolean operator!(const boolean &a){
        return boolean(!a.flag);
    }


};


#endif //COMPLEXMATRIX_BOOLEAN_H
