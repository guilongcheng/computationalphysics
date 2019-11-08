#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <libsec3.h>
#include <GLC.h>
#include <time.h>

double integ_examp_f1(double x)
{
    double a0 = 5;
    double s = sin(x);
    double c = cos(x);
    double f = (1 + a0 * (1 - c) * (1 - c)) / ((1 + a0 * s * s) * sqrt(1 + 2 * a0 * (1 - c)));
    return f;
}

double findroot_examp_f1(double x)
{
    return exp(x) * log(x) - x * x;
}

double findroot_examp_df1(double x)
{
    return exp(x) * (log(x) + 1 / x) - 2 * x;
}
int main()
{
    printf("============ Derivative example 1==================\n");
    int n = 101, m = 5;
    double x[n], f[n], df[n], ddf[n];
    int k = 2;
    double h = PI / (2 * n);
    for (int i = 0; i < n; i++)
    {
        x[i] = h * i;
        f[i] = sin(x[i]);
    }
    FirstOrderDerivative(h, f, df, n, k);
    SecondOrderDerivative(h, f, ddf, n, k);

    printf("result as : \n");
    for (int i = 0; i < n; i += 20)
    {
        printf("%lf    %lf    %lf    %lf\n", x[i], f[i], df[i], ddf[i]);
    }

    printf("============ 积分 example 1==================\n");
    time_t start, end;
    double del = 1e-6;
    start = time(NULL);
    double s = SimpsonAdaptive(integ_examp_f1, 0, PI, del);
    //double s = SimpsonAdaptive2(integ_examp_f1,0,PI,del,0);
    end = time(NULL);
    printf("用时为 %.2f 秒 \n", difftime(end, start));
    printf("辛普森自适应积分结果为 %.9f \n", s);
    printf("误差为 %e \n", fabs(s - PI));

    printf("============ 方程求根 example 1==================\n");
    Bisect(findroot_examp_f1, 1, 2, 1e-6);
    Newton(findroot_examp_f1,findroot_examp_df1,1,1e-6);
    Secant(findroot_examp_f1,1,1e-6);

    return 0;
}