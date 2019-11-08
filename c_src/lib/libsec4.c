#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "libsec3.h"
#include "GLC.h"

void Eular(void (*g)(double, double[MMAX], double[MMAX]), double y0[MMAX], double *y, int m, int n, double h)
{
    double ya[n][MMAX];
    double dy1[MMAX], dy2[MMAX];

    for (int i = 0; i < m; i++)
    {
        ya[0][i] = y0[i];
    }

    for (int i = 1; i < n; i++)
    {
        g((i - 1) * h, ya[i - 1], dy1);
        for (int j = 0; j < m; j++)
        {
            ya[i][j] = ya[i - 1][j] + h * dy1[j];
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            *(y + i * m + j) = ya[i][j];
        }
    }
}

void PreEular(void (*g)(double, double[MMAX], double[MMAX]), double y0[MMAX], double *y, int m, int n, double h)
{
    double ya[n][MMAX];
    double dy1[MMAX], dy2[MMAX];

    for (int i = 0; i < m; i++)
    {
        ya[0][i] = y0[i];
    }

    for (int i = 1; i < n; i++)
    {
        g((i - 1) * h, ya[i - 1], dy1);
        for (int j = 0; j < m; j++)
        {
            ya[i][j] = ya[i - 1][j] + h * dy1[j];
        }
        g(i * h, ya[i], dy2);
        for (int j = 0; j < m; j++)
        {
            ya[i][j] = ya[i - 1][j] + h / 2 * (dy1[j] + dy2[j]);
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            *(y + i * m + j) = ya[i][j];
        }
    }
}

void TwoPointPre(void (*g)(double, double[MMAX], double[MMAX]), double y0[MMAX], double y1[MMAX], double *y, int m, int n, double h)
{
    double ya[n][MMAX];
    double dy1[MMAX], dy2[MMAX], dy3[MMAX];

    for (int i = 0; i < m; i++)
    {
        ya[0][i] = y0[i];
        ya[1][i] = y1[i];
    }

    for (int i = 2; i < n; i++)
    {
        g((i - 1) * h, ya[i - 1], dy2);
        for (int j = 0; j < m; j++)
        {
            ya[i][j] = ya[i - 2][j] + 2 * h * dy2[j];
        }
        g((i - 2) * h, ya[i - 2], dy1);
        g(i * h, ya[i], dy3);
        for (int j = 0; j < m; j++)
        {
            ya[i][j] = ya[i - 2][j] + h / 3 * (dy3[j] + 4 * dy2[j] + dy1[j]);
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            *(y + i * m + j) = ya[i][j];
        }
    }
}

void RungeKutta(void (*g)(double, double[MMAX], double[MMAX]), double y0[MMAX], double *y, int m, int n, double h)
//RungeKutta 算法求解常微分方程的初值问题。g为每次迭代y的变换函数，y0为初始值，指针y存放结果,m为方程个数，n为数据点个数，h为求解步长。
{
    double ya[n][MMAX];
    double gy[MMAX];
    double ypc[MMAX] = {0};
    double c1[MMAX], c2[MMAX], c3[MMAX], c4[MMAX];

    for (int i = 0; i < m; i++)
    {
        ya[0][i] = y0[i];
    }

    for (int i = 1; i < n; i++)
    {
        g((i - 1) * h, ya[i - 1], gy);
        for (int j = 0; j < m; j++)
        {
            c1[j] = h * gy[j];
            ypc[j] = ya[i - 1][j] + c1[j] / 2;
        }
        g((i - 1) * h + h / 2, ypc, gy);
        for (int j = 0; j < m; j++)
        {
            c2[j] = h * gy[j];
            ypc[j] = ya[i - 1][j] + c2[j] / 2;
        }
        g((i - 1) * h + h / 2, ypc, gy);
        for (int j = 0; j < m; j++)
        {
            c3[j] = h * gy[j];
            ypc[j] = ya[i - 1][j] + c3[j];
        }
        g(i * h, ypc, gy);
        for (int j = 0; j < m; j++)
        {
            c4[j] = h * gy[j];
        }

        for (int j = 0; j < m; j++)
        {
            ya[i][j] = ya[i - 1][j] + 1.0 / 6.0 * (c1[j] + 2 * c2[j] + 2 * c3[j] + c4[j]);
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            *(y + i * m + j) = ya[i][j];
        }
    }
}

void Shooting(void (*g)(double, double[2], double[2]), double y0[2], double y1[2], double *y, int n, int k, double h)
//利用打靶法求解一维二阶微分方程的边值问题,这里的k取0或者1说明未知初值的维度，y0 是 边界a上值，y1 是边界b上的值.
{
    double guessy0 = 1;        //初始猜测y[0]
    double guessy1 = 1 + 0.01; //初始猜测
    double guessy2;            //迭代使用
    double del = 1e-8;         //控制精度
    double dy = guessy1 - guessy0;
    double f1, f2;
    int step = 0;

    y0[k] = guessy0;
    RungeKutta(g, y0, y, 2, n, h);
    f1 = *(y + (n - 1) * 2 + 1 - k) - y1[1 - k];

    while (fabs(f1) > del && step < NMAXITER)
    {
        y0[k] = guessy1;
        RungeKutta(g, y0, y, 2, n, h);
        f2 = *(y + (n - 1) * 2 + 1 - k) - y1[1 - k];

        guessy2 = guessy1 - f2 * dy / (f2 - f1);
        guessy0 = guessy1;
        guessy1 = guessy2;
        f1 = f2;
        dy = guessy1 - guessy0;
        step++;
    }
    printf("总迭代次数为 %d \n", step);
    printf("初始值应为 %15.7e , %15.7e \n", y0[0], y0[1]);
}

void LinearDeqSolve(void (*g)(double, double[2], double[2]), double y0[2], double y1[2], double *y, int n, int k, double h)
//利用打靶法求解一维二阶微分方程的边值问题,这里的k取0或者1说明未知初值的维度，y0 是 边界a上值，y1 是边界b上的值.
{
    double f1, f2, a, b;
    double *solvey1, *solvey2;
    solvey1 = (double *)malloc(n * 2 * sizeof(double));
    solvey2 = (double *)malloc(n * 2 * sizeof(double));

    y0[k] = 1;
    RungeKutta(g, y0, solvey1, 2, n, h);
    f1 = *(solvey1 + (n - 1) * 2 + 1 - k);
    printf("y0 is %lf %lf \n", y0[0], y0[1]);

    y0[k] = 2;
    RungeKutta(g, y0, solvey2, 2, n, h);
    f2 = *(solvey2 + (n - 1) * 2 + 1 - k);
    printf("y0 is %lf %lf \n", y0[0], y0[1]);

    a = (f2 - y1[1 - k]) / (f2 - f1);

    b = (y1[1 - k] - f1) / (f2 - f1);

    printf("a and b is %lf, %lf \n", a, b);
    for (int i = 0; i < n; i++)
    {
        for (int ij = 0; ij < 2; ij++)
        {
            *(y + i * 2 + ij) = a * (*(solvey1 + i * 2 + ij)) + b * (*(solvey2 + i * 2 + ij));
        }
    }

    free(solvey1);
    free(solvey2);
}

void sturmLiouville(double p[], double p1[], double q[], double s[], double u0[2], double *u, double l, int n, double h)
//求解sturm Liouville 方程, p p1 q s
{
    double h2 = h * h;
    *(u) = u0[0];
    *(u + 1) = u0[1];

    for (int i = 1; i < n - 1; ++i)
    {
        double c2 = 2 * p[i] + h * p1[i];
        double c1 = 4 * p[i] - 2 * h2 * q[i];
        double c0 = 2 * p[i] - h * p1[i];
        double d = 2 * h2 * s[i];
        *(u + i + 1) = (c1 * (*(u + i)) - c0 * (*(u + i - 1)) + d) / c2;
    }
}

void Numerov(double q[], double s[], double u0[2], double *u, int n, double h)
//求解sturm Liouville 方程, p p1 q s
{
    double h2 = h * h;
    u[0] = u0[0];
    u[1] = u0[1];
    double g = h * h / 12;
    for (int i = 1; i < n - 1; ++i)
    {
        double c0 = 1 + g * q[i - 1];
        double c1 = 2 - 10 * g * q[i];
        double c2 = 1 + g * q[i + 1];
        double d = g * (s[i + 1] + s[i - 1] + 10 * s[i]);
        u[i + 1] = (c1 * u[i] - c0 * u[i - 1] + d) / c2;
    }
}

static double *ur, *ul;
static int nr, nl;

void wave(double (*v)(double, double[NPARS]), double pars[NPARS], double u0[2], double u[], double x0[2], double en, int n)
{
    double y[n];
    double ql[n], qr[n], s[n];
    double h = (x0[1] - x0[0]) / (n - 1);

    for (int i = 0; i < n; i++)
    {
        double x = x0[0] + i * h;
        ql[i] = 2 * (en - v(x, pars));
        qr[n - i - 1] = ql[i];
        s[i] = 0;
        //        printf("ql,v = %lf, %lf \n",ql[i],v(x,pars));
    }

    int im = 0;
    for (int i = 0; i < n - 1; i++)
    {
        if (((ql[i] * ql[i + 1]) < 0) && (ql[i] > 0))
        {
            im = i;
        }
    }

    nl = im + 2;
    nr = n - im + 1;
    Numerov(ql, s, u0, ul, nl, h);
    Numerov(qr, s, u0, ur, nr, h);

    double ratio = ur[nr - 2] / ul[im];
    for (int i = 0; i < im + 1; i++)
    {
        u[i] = ratio * ul[i];
        y[i] = u[i] * u[i];
    }

    for (int i = 0; i < nr - 1; i++)
    {
        u[i + im] = ur[nr - i - 2];
        y[i + im] = u[i + im] * u[i + im];
    }

    double sum = SimpsonArray(y, n, h);
    sum = sqrt(sum);

    for (int i = 0; i < n; i++)
    {
        u[i] = u[i] / sum;
    }
}

void SolveSchroedingerEQ(double (*v)(double, double[NPARS]), double pars[NPARS], double u0[2], double u[], double e0[2], double en, int n, double x0[2])
//求解一维薛定谔方程
{
    double del = 1e-6;
    double h = (x0[1] - x0[0]) / (n - 1);

    ur = (double *)malloc(sizeof(double) * n);
    ul = (double *)malloc(sizeof(double) * n);
    double f1, f2;
    int step = 0, nmaxstep = 20;
    double en0 = e0[0], en1 = e0[1], deltaen;
    deltaen = en1 - en0;

#ifdef DEBUG
    printf("初始能量为 ：%lf,%lf \n", e0[0], e0[1]);
    printf("求解边界为 ：%lf,%lf \n", x0[0], x0[1]);
    printf("解的初始值为 ：%lf,%lf \n", u0[0], u0[1]);
    printf("求解点数及步长为 %d, %lf \n", n, h);
    printf("初始能差 %lf\n", deltaen);
#endif

    wave(v, pars, u0, u, x0, en0, n);
    f1 = *(ur + nr - 1) + *(ul + nl - 1) - *(ur + nr - 3) - *(ul + nl - 3);
    f1 = f1 / (2 * h * (*(ur + nr - 2)));

    printf("f1 的值为 %lf \n", f1);
    printf("nl , nr 的值为 %d  %d \n", nl, nr);
    while ((fabs(deltaen) > del) && (step < nmaxstep))
    //割线法求解
    {
        wave(v, pars, u0, u, x0, en1, n);
        f2 = *(ur + nr - 1) + *(ul + nl - 1) - *(ur + nr - 3) - *(ul + nl - 3);
        f2 = f2 / (2 * h * (*(ur + nr - 2)));

        en = en1 - f2 * deltaen / (f2 - f1);

#ifdef DEBUG
        printf("第%d次迭代：nl,nr,f1,f2,en0,en1,en is %d, %d, %lf, %lf, %lf, %lf,%lf \n", step, nl, nr, f1, f2, en0, en1, en);
#endif
        en0 = en1;
        en1 = en;
        f1 = f2;
        deltaen = en1 - en0;
        step++;
    } //while

    printf("迭代次数为 %d \n", step);
    printf("能级为 %lf \n", en);
    free(ur);
    free(ul);
}