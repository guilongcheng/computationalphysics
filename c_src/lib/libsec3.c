#include "stdio.h"
#include "math.h"
#include "libsec2.h"
#include "GLC.h"

void FirstOrderDerivative(double h, double f[], double df[], int n, int k)
// 利用离散数据f来计算各点的导数值。 h 为步长，df为计算得到的导数值，n为数据点个数，k为用来做边界插值点的个数。
{
    double xl[k], fl[k], fr[k];

    for (int i = 1; i < n - 1; i++)
    {
        df[i] = (f[i + 1] - f[i - 1]) / (2 * h);
    }

    for (int i = 0; i < k; i++)
    {
        xl[i] = h * (i + 1);
        fl[i] = df[i + 1];
        fr[i] = df[n - i - 1];
    }
    df[0] = aitken(0, k, xl, fl);
    df[n - 1] = aitken(0, k, xl, fr);
}

void SecondOrderDerivative(double h, double f[], double ddf[], int n, int k)
// 利用离散数据f来计算各点的导数值。 h 为步长，df为计算得到的导数值，n为数据点个数，k为用来做边界插值点的个数。
{
    double xl[k], fl[k], fr[k];

    for (int i = 1; i < n - 1; i++)
    {
        ddf[i] = (f[i + 1] + f[i - 1] - 2 * f[i]) / (h * h);
    }

    for (int i = 0; i < k; i++)
    {
        xl[i] = h * (i + 1);
        fl[i] = ddf[i + 1];
        fr[i] = ddf[n - i - 1];
    }
    ddf[0] = aitken(0, k, xl, fl);
    ddf[n - 1] = aitken(0, k, xl, fr);
}

double FODAdaptive(double x, double (*f)(double), double h, double del)
//使用自适应方法来求解1阶导数
{
    double d = (f(x + h) - f(x - h)) / (2 * h);
    h = h / 2;
    double d1 = (f(x + h) - f(x - h)) / (2 * h);
    int step = 1;
    while (h * h * fabs(d - d1) > del && step < NMAXITER)
    {
        d = d1;
        h = h / 2;
        d1 = (f(x + h) - f(x - h)) / (2 * h);
        step++;
    }
    if (step == NMAXITER)
    {
        printf("达到了设定的最大迭代次数\n");
    }
    printf("迭代次数为 %d \n", step);

    return (4 * d1 - d) / 3;
}

double Simpson(double (*f)(double), double a, double b, long int n)
// simpson 算法计算 函数 f(x) 在[a,b]区间上的积分。n为区间划分的个数。
{
    double h = (b - a) / n;
    double quad = 0.0;
    double s0 = 0.0, s1 = 0.0, s2 = 0.0;
    for (long int i = 1; i < n; i += 2)
    {
        s0 = s0 + f(a + i * h);
        s1 = s1 + f(a + (i - 1) * h);
        s2 = s2 + f(a + (i + 1) * h);
    }
    quad = (s1 + 4 * s0 + s2) / 3;
    if (n % 2 == 0)
    {
        return h * quad;
    }
    else
    {
        return h * (quad + (5 * f(b) + 8 * f(b - h) - f(b - 2 * h)) / 12);
    }
}

double SimpsonArray(double y[], long int n, double h)
{
    double s0= 0 , s1 =0, s2=0;
    for (int i = 1; i < n - 1; i+=2)
    {
        s0 += y[i];
        s1 += y[i-1];
        s2 +=y[i+1];
    }
    double s = (s1 + 4*s0 + s2) /3;

    if( n % 2 ==0)
    {
        return h * (s + (5 * y[n-1] + 8*y[n-2] - y[n-3]) / 12);
    }
    else
    {
        return h * s;
    }
}

double SimpsonAdaptive(double (*f)(double), double a, double b, double del)
{
    double h = b - a;
    long int n = 2;
    double s0 = Simpson(f, a, b, n);
    n = n * 2;
    double s1 = Simpson(f, a, b, n);
    int step = 1;

    while (fabs(s1 - s0) > 15 * del && step < NMAXITER)
    {
        s0 = s1;
        n = n * 2;
        s1 = Simpson(f, a, b, n);
        step++;
    }
    if (step == NMAXITER)
    {
        printf("达到了设定的最大迭代次数\n");
    }

    printf("迭代次数为 %d \n", step);

    return s1;
}

double SimpsonAdaptive2(double (*f)(double), double a, double b, double del, int step)
// 采用递归计算积分
{
    double c = (a + b) / 2;
    double s0 = Simpson(f, a, b, 2);
    double s1 = Simpson(f, a, b, 4);
    step++;
    //    printf("第%d次迭代，积分区间为(%e,%e),积分结果为%e, %e, \n", step, a, b, s0, s1);
    //    printf("s1 - s0 = %e, del = %e \n ", fabs(s1 - s0), 15 * del);

    if (step >= NMAXITER)
    {
        printf("达到最大迭代次数，不收敛");
        return s1;
    }
    else
    {
        if (fabs(s1 - s0) < 15 * del)
        {
            return s1;
        }
        else
        {
            return SimpsonAdaptive2(f, a, c, del / 2, step) + SimpsonAdaptive2(f, c, b, del / 2, step);
        }
    }
}

double Bisect(double (*f)(double), double a, double b, double del)
{
    double c;
    double dx = b - a;
    int step = 0;
    while (fabs(dx) > del && step < NMAXITER)
    {
        c = (a + b) / 2;
        if ((f(a) * f(c)) < 0)
        {
            b = c;
            dx = b - a;
        }
        else
        {
            a = c;
            dx = b - a;
        }
        step++;
    }
    printf("二分法:总迭代次数为 %d \n", step);
    printf("二分法:方程的根为 %15.7e \n", c);
    printf("二分法:误差为 %15.7e \n", dx);
    return c;
}

double Newton(double (*f)(double), double (*pf)(double), double x0, double del)
{
    double dx = f(x0) / pf(x0);
    int step = 0;
    double x = x0;
    while (fabs(dx) > del && step < NMAXITER)
    {
        x = x - dx;
        dx = f(x) / pf(x);
        step++;
    }
    printf("牛顿法:总迭代次数为 %d \n", step);
    printf("牛顿法:方程的根为 %15.7e \n", x);
    printf("牛顿法:误差为 %15.7e \n", dx);
    return x;
}

double Secant(double (*f)(double), double x0, double del)
{
    double x1 = x0 ;
    double x2 = x1 + 10*del;
    double dx = x2 - x1;
    int step = 0;
    double x;
    double f1,f2;
    f1 = f(x1); 
    while (fabs(dx) > del && step < NMAXITER)
    {
        f2 = f(x2);
        x = x2 - f(x2) * dx / (f2 - f1); 
        f1 = f2;
        x1 = x2;
        x2 = x;
        dx = x2 - x1;
        step++;
    }
    printf("割线法:总迭代次数为 %d \n", step);
    printf("割线法:方程的根为 %15.7e \n", x);
    printf("割线法:误差为 %15.7e \n", dx);
    return x;
}