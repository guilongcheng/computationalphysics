#include<stdio.h>
#include<stdlib.h>
// c++ 可以方便的利用函数模版重载各类不同类型的参数。
// c 可以利用宏来实现函数模版的功能，或者就对不同的数据类型写对应的函数
void gauss_jordan(double **A, double *b, double *x, int n){
// A 为 nxn 的二维数组， b为 n 的一维数组
// 注意会改变原来A矩阵的值。
    double m;
    double tmp;
    int tmp_idx;
    for (int k = 0; k < n-1; ++k) {

//选主元
        tmp_idx=k;
        for (int i = k+1; i < n; ++i) {
            if ( A[i][k] > A[tmp_idx][k] ) tmp_idx=i;
        }
        if ( tmp_idx != k){
            for (int j = k; j < n; ++j) {
                tmp = A[k][j];
                A[k][j] = A[tmp_idx][j];
                A[tmp_idx][j] = tmp;
            }
        }
// 消元
        for (int i = k+1; i < n; ++i) {
            m = A[i][k] / A[k][k];
            A[i][k]= 0.0;
            for (int j = k+1; j < n; ++j) {
                A[i][j] = A[i][j] - A[k][j] * m;
            }
            b[i] = b[i] - b[k] * m;
        }
    }

    x[n-1] = b[n-1] / A[n-1][n-1];
    for (int i = n-2; i > -1; --i) {
        tmp = 0.0;
        for (int k = i+1; k < n; ++k) {
            tmp = tmp + x[k] * A[i][k];
        }
        x[i] = (b[i] - tmp) / A[i][i];
    }
}

void showV(double *x, int n){
     printf("vector x is :\n");
     for (int i = 0; i < n; ++i) {
         printf("%lf,",x[i]);
     }
     printf("\n");
}

void showM(double **A, int n){
     printf("Matrix A is : \n");
     for (int i = 0; i < n; ++i) {
         for (int j = 0; j < n; ++j) {
             printf("%lf,",A[i][j]);
         }
         printf("\n");
     }
     printf("\n");
}

int main(){


    const int N=3;
    double **A, *b, *x;
    A = (double **) malloc(N * sizeof(double *));
    for (int i = 0; i < N; ++i) {
        A[i] = (double *) malloc(N * sizeof(double ));
    }
    A[0][0] = 1;A[0][1] = 2;A[0][2] = 1;
    A[1][0] = 2;A[1][1] = 2;A[1][2] = 3;
    A[2][0] =-1;A[2][1] = 0;A[2][2] =-3;
    b = (double *) malloc(N * sizeof(double));
    x = (double *) malloc(N * sizeof(double));
    b[0] = 0; b[1] = 3; b[2] = 0;

//    lu_decom(A,b,x,L,U);
    gauss_jordan(A, b, x, N);

    showV(x, N);
    showM(A, N);

    for (int i = 0; i < N; ++i) {
        free(A[i]);
    }
    free(b);
    free(x);
}

