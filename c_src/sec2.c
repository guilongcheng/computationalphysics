#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"libsec2.h"


int main(){
    printf("============ Interpolation example 1==================\n");
    double xi[] = {0, 0.5, 1, 1.5, 2};
    double fi[] = {1, 0.938470, 0.765198, 0.511828, 0.223891};

    double x = 0.9;

    double f = aitken(x,5,xi,fi);
    printf("Interpolated value: %lf \n",f);

    printf("============ Interpolation example 2==================\n");

    f = updown(x,5,xi,fi);
    printf("Interpolated value: %lf \n",f);

    printf("============ Interpolation example 3 Millikan==================\n");
    double k[15]={4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
    double q[15]={6.558, 8.206, 9.880, 11.50, 13.14, 
    14.81, 16.40, 18.04, 19.68, 21.32, 22.96, 24.60,
    26.24, 27.88, 29.52};
    double alpha[15],*pu;
    int n = 15;
    int m = 1;
    pu =(double *)malloc(n*n*sizeof(double));
    orthogonalPolynomialFit(m,n,k,q,pu,alpha);
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += k[i];
    }
    
    printf("基本电荷为 %lf, 估计误差为 %lf \n",alpha[1], alpha[0] - alpha[1]*sum/n);

    return 0;
}