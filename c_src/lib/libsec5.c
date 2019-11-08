#include "math.h"

void Gaussian(double *a, int index[], int n)
// 采用列主元的高斯消去法，index里存储着主元的顺序。
{
    double c[n];

    for (int i = 0; i < n; i++)
    {
        index[i] = i;
    }

    //寻找每一行的比例因子
    for (int i = 0; i < n; i++)
    {
        double c1 = 0;
        for (int j = 0; j < n; j++)
        {
            double c0 = fabs(*(a + n * i + j));
            if (c0 > c1)
                c1 = c0;
        }
        c[i] = c1;
    }

    //寻找每一列的主元
    int k = 0;
    for (int j = 0; j < n - 1; j++)
    {
        double pi1 = 0;
        for (int i = j; i < n; i++)
        {
            double pi0 = fabs(*(a + index[i] * n + j));
            pi0 = pi0 / c[index[i]];
            if (pi0 > pi1)
            {
                pi1 = pi0;
                k = i;
            }
        }
        // 根据主元位置交换行
        int itmp = index[j];
        index[j] = index[k];
        index[k] = itmp;

        for (int i = j + 1; i < n; i++)
        {
            double pj = *(a + n * index[i] + j) / *(a + n * index[j] + j);
            // 将主元的比值记录在A矩阵元为0的位置
            *(a + n * index[i] + j) = pj;
            // 修改相应的矩阵元
            for (int l = j + 1; l < n; l++)
            {
                *(a + n * index[i] + l) = *(a + n * index[i] + l) - pj * (*(a + n * index[j] + l));
            }
        }
    }
}

double DetGauss(double *a, int n)
// 计算矩阵的行列式，
{
    int index[n];

    Gaussian(a, index, n);

    double d = 1;

    for (int i = 0; i < n; i++)
    {
        d = d * (*(a + n * index[i] + i));
    }

    int sgn = 1;
    for (int i = 0; i < n; i++)
    {
        if (i != index[i])
        {
            sgn = -sgn;
            int j = index[i];
            index[i] = index[j];
            index[j] = j;
        }
    }

    return d * sgn;
}

void LinearSolveGauss(double *a, double b[], double x[], int n)
// 利用高斯消元法求解线性方程组
{
    int index[n];
    Gaussian(a, index, n);

    //利用存在A矩阵下三角力的系数，对b做相应的变换。
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            b[index[j]] = b[index[j]] - *(a + n * index[j] + i) * b[index[i]];
        }
    }

    //回代求解x
    x[n - 1] = b[index[n - 1]] / *(a + n * index[n - 1] + n - 1);
    for (int i = n - 2; i >= 0; --i)
    {
        x[i] = b[index[i]];
        for (int j = i + 1; j < n; ++j)
        {
            x[i] -= *(a + index[i] * n + j) * x[j];
        }
        x[i] /= *(a + index[i] * n + i);
    }

    // 将x的顺序调整回去
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            if (index[i] > index[j])
            {
                int itmp = index[i];
                index[i] = index[j];
                index[j] = itmp;
                double xtmp = x[i];
                x[i] = x[j];
                x[j] = xtmp;
            }
        }
    }
}

void InvertGauss(double *a, double *x, int n)
{
    int index[n];
    double b[n][n];

    for (int i = 0; i < n; ++i)
        b[i][i] = 1;

    Gaussian(a, index, n);

    for (int i = 0; i < n - 1; ++i)
        for (int j = i + 1; j < n; ++j)
            for (int k = 0; k < n; ++k)
                b[index[j]][k] -= *(a + index[j] * n + i) * b[index[i]][k];

    for (int i = 0; i < n; ++i)
    {
        *(x + n * (n - 1) + i) = b[index[n - 1]][i] / *(a + index[n - 1] * n + n - 1);
        for (int j = n - 2; j >= 0; --j)
        {
            *(x + j * n + i) = b[index[j]][i];
            for (int k = j + 1; k < n; ++k)
            {
                *(x + j * n + i) -= (*(a + index[j] * n + k)) * (*(x + k * n + i));
            }
            *(x + j * n + i) /= *(a + index[j] * n + j);
        }
    }
}


void LUDecomp(double **a, int index[], int n)
//该程序用来将一个矩阵分解成一个上三角和一个下三角矩阵 {
{

}