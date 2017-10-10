#include<iostream>
#include<cmath>
using namespace std;



template<int N>
void lu_decom(double (&A)[N][N],double (&b)[N], double (&x)[N],double (&L)[N][N],double (&U)[N][N]){
    // LU decompositon of matrix A
    // and solve A*x =b
    // A :Matrix
    // b :vectro
    // return: x(solve),L(Matrix),U(Matrix)
    double y[N];

    int i,j,k;
    for (i = 0; i < N; ++i) {
        x[i] = 0.0;
        y[i] = 0.0;
        for (j = 0; j < N; ++j) {
            L[i][j] = 0.0;
            U[i][j] = 0.0;
        }
        L[i][i] = 1.0;
    }

    for (i = 0; i < N; ++i)U[0][i] = A[0][i];
    for (i=1;i<N;++i)L[i][0] = A[i][0] / U[0][0];

    for (i = 1; i < N; ++i) {
        for (j = i; j < N; ++j) {

            double tmp = 0.0;
            for (k = 0; k < i; ++k) {
                tmp = tmp + L[i][k]*U[k][j];
            }

            U[i][j] = A[i][j] - tmp;

            tmp = 0.0;
            for (k = 0; k < i; ++k) {
                tmp = tmp + L[j+1][k]*U[k][i];
            }
            if (j!=N-1) {
                L[j+1][i] = (A[j+1][i] - tmp)/U[i][i];
            }
        }
    }

    y[0] = b[0];

    for (i = 1; i < N; ++i) {

        double tmp = 0.0;
        for (k = 0; k < i; ++k) {
            tmp = tmp + L[i][k]*y[k];
        }

        y[i] = b[i] - tmp;
    }

    x[N-1] = y[N-1] /U[N-1][N-1];

    for (i = N-2; i > -1; --i) {

        double tmp = 0.0;
        for (k = i+1; k < N; ++k) {
            tmp = tmp + U[i][k]*x[k];
        }

        x[i] = (y[i] - tmp) / U[i][i];
    }
}

int main(){


    const int N=3;
    double A[N][N] = {{2,2,3},{4,7,7},{-2,4,5}};
    double b[N] = {3,1,-7};
    double x[N],L[N][N],U[N][N];

    lu_decom(A,b,x,L,U);

    cout<<"solve is "<<x[0]<<","<<x[1]<<","<<x[2]<<endl;
    cout<<"Matrix L is "<<endl;
    for (int i = 0; i < N; ++i) {
        cout<<L[i][0]<<","<<L[i][1]<<","<<L[i][2]<<endl;
    }
    cout<<"Matrix U is "<<endl;
    for (int i = 0; i < N; ++i) {
        cout<<U[i][0]<<","<<U[i][1]<<","<<U[i][2]<<endl;
    }

}
