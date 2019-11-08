#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "libsec5.h"
#include "GLC.h"

int main()
{
    printf("======================矩阵 例题1=================\n");
    {
        int n = 5;
        double *a;
        a = (double *)malloc(n * n * sizeof(double));
        double tmpa[] = {1, -3, 2, -1, -2,
                         -2, 2, -1, 2, 3,
                         3, -3, -2, 1, -1,
                         1, -2, 1, -3, 2,
                         -3, -1, 2, 1, -3};

        for (int i = 0; i < n * n; i++)
        {
            *(a + i) = tmpa[i];
        }

        printf("矩阵A的行列式为 %4.0f \n", DetGauss(a, n));
    }

    printf("======================矩阵 例题2=================\n");
    {
        int n = 3;
        double *a;
        a = (double *)malloc(n * n * sizeof(double));
        double tmpa[] = {100, 100, 100,
                         -100, 300, -100,
                         -100, -100, 300};
        double b[] = {200, 0, 0};
        double x[n];

        for (int i = 0; i < n * n; i++)
        {
            *(a + i) = tmpa[i];
        }
        LinearSolveGauss(a, b, x, n);

        printf("矩阵Ax=b的解为: \n");
        for (int i = 0; i < n; i++)
        {
            printf("%3.1f \n", x[i]);
        }
    }

    printf("======================矩阵 例题3=================\n");
    {
        int n = 3;
        double *a;
        a = (double *)malloc(n * n * sizeof(double));
        double tmpa[] = {100, 100, 100,
                         -100, 300, -100,
                         -100, -100, 300};
        double *x;
        x = (double *)malloc(n * n * sizeof(double));

        for (int i = 0; i < n * n; i++)
        {
            *(a + i) = tmpa[i];
        }
        InvertGauss(a, x, n);

        printf("矩阵A的逆为: \n");
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                printf("%5.4f ", *(x + i * n + j));
            }
            printf("\n");
        }
    }

    return 0;
}