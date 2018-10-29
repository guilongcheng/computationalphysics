#ifndef LINALG
#define LINALG
#include"../Eigen/Dense"
#include<iostream>
using namespace Eigen;
using namespace std;
template <typename Derived, typename Derived2, typename Derived3>
void gauss_jordan(const EigenBase<Derived>& A, const EigenBase<Derived2>& b, EigenBase<Derived3>& x)
{
    double m;
    double tmp;
    int tmp_idx;
    int n = b.size();
    for (int k = 0; k < n-1; ++k) {

//选主元
        tmp_idx=k;
        for (int i = k+1; i < n; ++i) {
            if ( A(i,k) > A(tmp_idx,k) ) tmp_idx=i;
        }
        if ( tmp_idx != k){
            for (int j = k; j < n; ++j) {
                tmp = A(k,j);
                A(k,j) = A(tmp_idx,j);
                A(tmp_idx,j) = tmp;
            }
        }
// 消元
        for (int i = k+1; i < n; ++i) {
            m = A(i,k) / A(k,k);
            A(i,k)= 0;
            for (int j = k+1; j < n; ++j) {
                A(i,j) = A(i,j) - A(k,j) * m;
            }
            b(i) = b(i) - b(k) * m;
        }
    }

    x(n-1) = b(n-1) / A(n-1,n-1);
    for (int i = n-2; i > -1; --i) {
        tmp = 0;
        for (int k = i+1; k < n; ++k) {
            tmp = tmp + x(k) * A(i,k);
        }
        x(i) = (b(i) - tmp) / A(i,i);
    }

}

template <typename Derived>
void showM(const EigenBase<Derived>& A)
{
    cout<<"Matrix A is : "<<endl;
     for (int i = 0; i < A.rows(); ++i) {
         for (int j = 0; j < A.cols(); ++j) {
             cout<<A(i,j)<<",";
         }
         cout<<endl;
     }
     cout<<endl;
}

template <typename Derived>
void showV(const EigenBase<Derived>& V)
{
    cout<<"Vector V is : "<<endl;
    for (int i = 0; i < V.size(); ++i) {
        cout<<V(i)<<",";
    }
    cout<<endl;
}

//template <typename Derived, typename Derived2, typename Derived3>
//void gauss_jordan(EigenBase<Derived>& A, EigenBase<Derived2>& b, EigenBase<Derived3>& x);
//
//template <typename Derived>
//void showM(const EigenBase<Derived>& A);
//
//template <typename Derived>
//void showV(const EigenBase<Derived>& V);

#endif
