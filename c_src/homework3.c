#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "GLC.h"
#include "libsec3.h"

double En, b, kappa, alpha;

double V(double r)
{
    return kappa / r * exp(-r / alpha);
}

double frm(double rm)
{
    return 1 - b * b / (rm * rm) - V(rm) / En;
}

double intf(double r)
{
    return 1.0 / (r * r * sqrt(1 - b * b / (r * r) - V(r) / En));
}

double GeneralIntInfty(double (*intf)(double), double a0)
{
    int stepmax = 1000;
    double result = 0;
    double int0 = 1, int1 = 1, int2 = 1;
    int step = 0;
    double down = a0 + 0.5;
    double b0 = 10, up;
    double del = 1e-6;
    int Nmax = 100;

    int0 = SimpsonAdaptive(intf, down, b0, del);
    result = result + int0;

    up = down;
    int i = 2;
    down = a0 + pow(0.5, i);
    while ((fabs(int1) > 1e-6) && (step < stepmax))
    {
        //int1 = int2;
        int1 = SimpsonAdaptive(intf, down, up, del);
        //int2 = Simpson(intf,a,b,Nmax);
        step++;
        //printf("第 %d 次计算积分,积分结果为%15.7e, 积分区间为(%f,%f) \n", step, int1, down, up);
        up = down;
        i++;
        down = a0 + pow(0.5, i);
        result = result + int1;
    }

    down = b0;
    up = down * 10;
    while ((fabs(int2) > 1e-6) && (step < stepmax))
    {
        int2 = SimpsonAdaptive(intf, down, up, del);
        step++;
        //printf("第 %d 次计算积分,积分结果为%15.7e, 积分区间为(%f,%f) \n", step, int2, down, up);
        down = up;
        up = down * 10;
        result = result + int2;
    }

    if (step == stepmax)
        printf("达到了最大迭代次数\n");
    printf("积分结果为：%15.7e \n", result);

    return result;
}
int main()
{
    int n = 20;
    double rm;
    double theta[n];
    double sigma[n];
    double b0 = 0.01, db = 0.5;
    double alphas[4] = {0.1, 1, 10, 100};
    double *plotdata;
    plotdata = (double *)malloc(2 * n * sizeof(double));

    En = 1, kappa = 1; //参数初始化

    for (int j = 0; j < 4; j++) // 计算不同的alpha下的结果
    {
        alpha = alphas[j];
        for (int i = 0; i < n; i++) // 对不同的b计算出相应的theta
        {
            b = b0 + i * db;

            rm = Secant(frm, b, 1e-8);

            theta[i] = PI - 2 * b * GeneralIntInfty(intf, rm);
            int k = (int)(theta[i] / PI);
            if (k >= 1)
            {
                theta[i] = theta[i] - k * PI;
            }
            if (k <= -1)
            {
                theta[i] = theta[i] + k * PI;
            }
        }

        double dtheta[n];
        FirstOrderDerivative(db, theta, dtheta, n, 2);

        for (int i = 0; i < n; i++)
        {
            b = b0 + i * db;
            sigma[i] = b / sin(theta[i]) * fabs(1.0 / dtheta[i]);
            printf("theta[i],dtheta[i],sigma[i] = %15.7f  %15.7f  %15.7f \n", theta[i], dtheta[i], sigma[i]);
        }

        for (int i = 0; i < n; i++)
        {
            *(plotdata + 2 * i) = theta[i];
            *(plotdata + 2 * i + 1) = log(sigma[i]);
        }
        char filename[256];
        sprintf(filename, "homework3_data_alpha_%04.1f", alpha);
        savef(filename, plotdata, n, 2, 2);
    }
}