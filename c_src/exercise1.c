#include "libsec3.h"
#include "stdio.h"
#include "math.h"
#define PI 3.1415926535898

static double theta0;

double intf(double theta)
{
    return 1.0 / sqrt(cos(theta) - cos(theta0));
}

int main()
{
    double coeff = 4 * sqrt(1.0 / 2.0 * 9.8);

    theta0 = PI / 2.0;
    int stepmax = 1000;
    double deltaf = 1;
    double result = 0;
    double int1 = 0, int2;
    int step = 4;
    double a = 0, b = (theta0 + a) / 2;
    double del = 1e-6;
    int Nmax = 100;

    while ((deltaf > 1e-6) && (step < stepmax))
    {
        //int1 = int2;
        int2 = SimpsonAdaptive(intf, a, b, del);
        //int2 = Simpson(intf,a,b,Nmax);
        step++;
        printf("第 %d 次计算积分,积分结果为%15.7e, 积分区间为(%f,%f) \n", step, int2, a, b);
        a = b;
        b = (theta0 + b) / 2.0;
        result = result + int2;
        deltaf = fabs(int2);
    }
    if (step == stepmax)
        printf("达到了最大迭代次数\n");
    printf("积分结果为：%15.7e \n", result);

    return 0;
}