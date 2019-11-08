#include <stdio.h>
#include <math.h>
double db=0.5;
double f(double t, double b)
{
    double result; /* db=0.005, a=100 */
    result = 1 - (b * b / t) - (1 / t) * exp( -t / 100);
    return result;
}

double fx(double xr, double b)
{
    double r2;
    r2 = xr * xr;
    return b / (r2 * sqrt(1 - b * b / r2 - 1 / xr * exp( -xr / 100)));
} /*积分函数*/

double theta(double b)
{
    double x0 = b, x1 = b + 0.5, x2;
    int n = 0;
    while ((fabs(x1 - x0) >= 0.000001) && (n <= 100))
    {
        x2 = x1 - (x1 - x0) / (f(x1, b) - f(x0, b)) * f(x1, b);
        x0 = x1;
        x1 = x2;
        n++;
    };

    int N = 100;
    double h, t = 0, s = 0, y[2 * N + 1], x[2 * N + 1];
    x[2 * N] = 10000, x[0] = x2;
    y[2 * N] = fx(10000,b), y[0] = fx(x2,b);
    h = (10000 - x2) / (2 * N);
    for (int i = 1; i < 2 * N + 1; i++)
    {
        x[i] = x2 + i * h;
        y[i] = fx(x[i],b);
    }
    for (int i = 2; i <= 2 * N; i = i + 2)
    {
        s = s + h * 4.0 / 3.0 * y[i - 1];
        s = s + h * 2.0 / 3.0 * y[i];
    }
    s = s + h / 3.0 * (y[0] + y[2 * N]) - h * 2.0 / 3.0 * y[2 * N];
    return (3.1415926 - 2 * s);
} /*求theta*/

double sigma(double bn)
{
    return (theta(bn + db) - theta(bn - db) / (2 * db));
}

int main()
{
    float b[20] = {0.01}, th[20], si[20];
    for (int i = 1; i < 20; i++)
        b[i] = b[i - 1] + db;
    for (int i = 0; i < 20; i++)
    {
        th[i] = theta(b[i]);
        si[i] = sigma(b[i]);
        printf("b[%d]=%f,th[%d]=%f,si[%d]=%f\n", i, b[i], i, th[i], i, si[i]);
    }
}