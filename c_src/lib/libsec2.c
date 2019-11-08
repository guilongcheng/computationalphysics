#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "GLC.h"


double aitken(double x, int n, double xi[], double fi[])
// 使用aitken算法来进行拉格朗日插值，输入数据点xi，fi，以及数据点个数，返回插值点x处的插值值。
{
    double ft[n];
    for (int i = 0; i < n; i++)
    {
        ft[i] = fi[i];
    }
    // 注意循环的次数,第一列的结果就是fi的值，所以循环计算时，第一次循环的次数为n-1;
    for (int i = 0; i < (n - 1); ++i)
    {
        for (int j = 0; j < (n - 1) - i; ++j)
        {
            ft[j] = (x - xi[j]) / (xi[i + j + 1] - xi[j]) * ft[j + 1] + (x - xi[i + j + 1]) / (xi[j] - xi[i + j + 1]) * ft[j];
        }
    }
    return ft[0];
}

double updown(double x, int n, double xi[], double fi[])
// 采用上下修正的方法来计算拉格朗日插值
{
    int i, i0, j, j0, k, it;
    double dx0, dx1, dt;
    double dp[NMAX][NMAX], dm[NMAX][NMAX];
    double f, df;

    if (n > NMAX)
    {
        printf("数据点的个数太多，超过了21 ！！！\n");
        exit(1);
    }
    // 初始化，dp,dm 分别为上修正和下修正。用fi初始化它们，并寻找离x最近的点xi[i0];
    i0 = 0;
    dx0 = fabs(xi[n - 1] - xi[0]);
    for (i = 0; i < n; ++i)
    {
        dp[i][i] = fi[i];
        dm[i][i] = fi[i];
        dx1 = fabs(x - xi[i]);
        if (dx1 < dx0)
        {
            i0 = i;
            dx0 = dx1;
        }
    }
    j0 = i0;

    // 计算修正矩阵
    // 见修正矩阵的计算公式
    for (i = 0; i < n - 1; ++i)
    {
        for (j = 0; j < n - i - 1; ++j)
        {
            k = j + i + 1;
            dt = (dp[j][k - 1] - dm[j + 1][k]) / (xi[k] - xi[j]);
            dp[j][k] = dt * (xi[k] - x);
            dm[j][k] = dt * (xi[j] - x);
        }
    }
    // 改进近似结果
    f = fi[i0];
    it = 0;
    if (x < xi[i0])
        it = 1;
    for (i = 0; i < n - 1; ++i)
    {
        if ((it == 1) || (j0 == n - 1))
        {
            i0 = i0 - 1;
            df = dp[i0][j0];
            f = f + df;
            it = 0;
            if (j0 == n - 1)
                it = 1;
        }
        else if ((it == 0) || (i0 == 0))
        {
            j0 = j0 + 1;
            df = dm[i0][j0];
            f = f + df;
            it = 1;
            if (i0 == 0)
                it = 0;
        }
    }
    df = fabs(df);
    return f;
}

void orthogonalPolynomialFit(int m, int n, double x[], double f[], double* pu, double alpha[])
// n个数据点用m次多项式进行拟合，返回正交多项式u_k(x_i),以及系数alpha_i
{
    double s[n], g[n], h[n], u[n][n];
    if (m > n - 1)
    {
        printf("多项式阶数大于数据点个数，将使用m=n -1");
            m = n - 1;
    }
    //计算0阶多项式
    for (int i = 0; i < n; i++)
    {
        s[i] = 0;
        g[i] = 0;
        h[i] = 0;
        alpha[i] = 0;
    }

    double stmp;
    for (int i = 0; i <= n - 1; i++)
    {
        u[0][i] = 1;
        stmp = u[0][i] * u[0][i];
        s[0] += stmp;
        g[0] += x[i] * stmp;
        alpha[0] += u[0][i] * f[i];
    }
    g[0] = g[0] / s[0];
    alpha[0] = alpha[0] / s[0];

    //计算1阶多项式
    for (int i = 0; i <= n - 1; i++)
    {
        u[1][i] = x[i] * u[0][i] - g[0] * u[0][i];
        s[1] += u[1][i] * u[1][i];
        g[1] += x[i] * u[1][i] * u[1][i];
        h[1] += x[i] * u[1][i] * u[0][i];
        alpha[1] += u[1][i] * f[i];
    }
    g[1] = g[1] / s[1];
    h[1] = h[1] / s[0];
    alpha[1] = alpha[1] / s[1];

    //计算其他阶多项式
    if (m >= 2)
    {
        for (int i = 1; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                u[i + 1][j] = x[j] * u[i][j] - g[i] * u[i][j] - h[i] * u[i - 1][j];
                s[i+1] += u[i+1][j] * u[i+1][j];
                g[i+1] += x[j] * u[i+1][j] * u[i+1][j];
                h[i+1] += x[j] * u[i+1][j] * u[i][j];
                alpha[i+1] += u[i+1][j] * f[j];
            }
            g[i+1] = g[i+1] / s[i+1];
            h[i+1] = h[i+1] / s[i];
            alpha[i+1] = alpha[i+1] / s[i+1];
        }
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            *(pu + n * i + j) = u[i][j];
        }
        
    }
    
}