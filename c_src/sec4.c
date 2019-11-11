#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "libsec3.h"
#include "libsec4.h"
#include "GLC.h"

void Twofuncg(double x, double y[MMAX], double dy[MMAX])
{
    dy[0] = y[1];
    dy[1] = -y[0];
}

static double k, g;
void jump(double x, double y[MMAX], double dy[MMAX])
{
    double v = sqrt(y[2] * y[2] + y[3] * y[3]);

    dy[0] = y[2];
    dy[1] = y[3];
    dy[2] = -k * v * y[2];
    dy[3] = -g - k * v * y[3];
}

static double b, omega0, q;
void pendulum(double x, double y[MMAX], double dy[MMAX])
{
    dy[0] = y[1];
    dy[1] = -sin(y[0]) + b * cos(omega0 * x);
    dy[1] = dy[1] - q * y[1];
}

void shootf(double x, double y[MMAX], double dy[MMAX])
{
    dy[0] = y[1];
    dy[1] = -PI * PI / 4.0 * (y[0] + 1);
}

inline double SL_p(double x, double l)
{
    return 1 - x * x;
}
inline double SL_p1(double x, double l)
{
    return -2 * x;
}
inline double SL_q(double x, double l)
{
    return l * (l + 1);
}
inline double SL_s(double x, double l)
{
    return 0;
}

double SL_f(double l)
{
    double *u;
    int n = 101;
    double h = 1.0 / (n - 1);
    double u0[] = {0, h};
    double p[n], p1[n], q[n], s[n];
    u = (double *)malloc(sizeof(double) * n);

    for (int i = 0; i < n; i++)
    {
        double x = i * h;
        p[i] = SL_p(x, l);
        p1[i] = SL_p1(x, l);
        q[i] = SL_q(x, l);
        s[i] = SL_s(x, l);
    }

    sturmLiouville(p, p1, q, s, u0, u, n, h);

    return *(u + n - 1) - 1;
}

double schroedingerV(double x, double pars[])
{
    double alpha = pars[0];
    double lambda = pars[1];
    return alpha * alpha * lambda * (lambda - 1) * (0.5 - 1.0 / pow(cosh(alpha * x), 2)) / 2;
}

double schroedingerV2(double x, double pars[])
{
    double k = pars[0];
    return 1.0 * x*x / 2 ;
}

int main()
{
    printf("================求解微分方程 例一==============\n");
    {
        double y0[2] = {0, 1};
        double dy[2];
        int n = 100;
        int m = 2;
        double h = 2 * PI / n;
        double *y1, *y2;
        y1 = (double *)malloc(n * m * sizeof(double));
        y2 = (double *)malloc(n * m * sizeof(double));
        Eular(Twofuncg, y0, y1, m, n, h);
        PreEular(Twofuncg, y0, y2, m, n, h);
#ifdef SHOW
        printf("欧拉法的结果: \n");
        for (int i = 0; i < n; i++)
        {
            printf("%lf    %lf    %lf \n", i * h, *(y1 + m * i), *(y1 + m * i + 1));
        }

        printf("预处理欧拉法的结果: \n");
        for (int i = 0; i < n; i++)
        {
            printf("%lf    %lf    %lf \n", i * h, *(y2 + m * i), *(y2 + m * i + 1));
        }
#endif

#ifdef SAVEFILE
        double *y3;
        y3 = (double *)malloc(n * 5 * sizeof(double));
        for (int i = 0; i < n; i++)
        {
            *(y3 + i * 5) = i * h;
            *(y3 + i * 5 + 1) = *(y1 + i * m);
            *(y3 + i * 5 + 2) = *(y1 + i * m + 1);
            *(y3 + i * 5 + 3) = *(y2 + i * m);
            *(y3 + i * 5 + 4) = *(y2 + i * m + 1);
        }
        char filename[256] = "result_tmp1.txt";
        printf("结果保存在文件%s中 \n", filename);
        savef(filename, y3, n, 5, 5);
        free(y3);
#endif
        free(y1);
        free(y2);
    } //例题一

    printf("================求解微分方程 例二==============\n");
    {
        int n = 100, m = 4;
        double area = 0.93, mass = 250, density = 1.2;
        k = area * density / (2 * mass), g = 9.8;
        double y0[4], y1[4];
        double angle = 42.5 * PI / 180;
        double v0 = 67;
        double h = 2 * v0 * sin(angle) / (9.8 * n);
        y0[0] = 0, y0[1] = 0, y0[2] = v0 * cos(angle), y0[3] = v0 * sin(angle);
        double v = sqrt(y0[2] * y0[2] + y0[3] * y0[3]);
        double ax = -k * v * y0[2];
        double ay = -g - k * v * y0[3];
        double hh2 = h * h / 2;
        double p = y0[2] * ax + y0[3] * ay;
        y1[0] = y0[0] + h * y0[2] + hh2 * ax;
        y1[1] = y0[1] + h * y0[3] + hh2 * ay;
        y1[2] = y0[2] + h * ax - hh2 * k * (v * ax + p * y0[2] / v);
        y1[3] = y0[3] + h * ay - hh2 * k * (v * ay + p * y0[3] / v);

        double *y;
        y = (double *)malloc(n * m * sizeof(double));
        TwoPointPre(jump, y0, y1, y, m, n, h);
        double *y2;
        y2 = (double *)malloc(n * m * sizeof(double));
        Eular(jump, y0, y2, m, n, h);
#ifdef SHOW
        printf("两点预处理算法的结果: \n");
        for (int i = 0; i < n; i++)
        {
            printf("%lf    %lf    %lf \n", i * h, *(y + m * i), *(y + m * i + 1));
        }
#endif
#ifdef SAVEFILE
        char filename[256] = "result_tmp3.txt";
        printf("两点预处理算法的结果保存在文件%s中 \n", filename);
        savef(filename, y, n, m, 2);
        char filename2[256] = "result_tmp2.txt";
        printf("欧拉算法的结果保存在文件%s中 \n", filename2);
        savef(filename2, y2, n, m, 2);
#endif
        free(y);
        free(y2);
    }

    printf("================求解微分方程 例三==============\n");
    {
        int n = 2000, m = 2, nt = 10;
        double h = 3 * PI / nt;
        double y0[] = {0, 2};
        double *y;
        y = (double *)malloc(n * m * sizeof(double));
        q = 0.5, b = 1.15, omega0 = 2.0 / 3;
        RungeKutta(pendulum, y0, y, m, n, h);
        int np;
        //重标度theta使得在[-\pi,\pi]之间
        for (int i = 0; i < n; i++)
        {
            np = (int)(*(y + i * m) / (2 * PI) + 0.5);
            *(y + i * m) = *(y + i * m) - 2 * PI * np;
        }

#ifdef SHOW
        printf("RungeKutta算法的结果: \n");
        for (int i = 0; i < n; i++)
        {
            printf("%lf    %lf    %lf \n", i * h, *(y + m * i), *(y + m * i + 1));
        }
#endif

#ifdef SAVEFILE
        char filename[256] = "result_tmp4.txt";
        printf("RungeKutta算法的结果保存在文件%s中 \n", filename);
        savef(filename, y, n, m, m);
#endif
        free(y);
    }

    printf("================求解微分方程 例四==============\n");
    {
        int n = 100, m = 2;
        double h = 1.0 / n;
        double y0[] = {0, 0};
        double y1[] = {1, 0};
        double *y;
        y = (double *)malloc(n * m * sizeof(double));
        Shooting(shootf, y0, y1, y, n, 1, h);
#ifdef SHOW
        printf("RungeKutta算法的结果: \n");
        for (int i = 0; i < n; i++)
        {
            printf("%lf    %lf    %lf \n", i * h, *(y + m * i), *(y + m * i + 1));
        }
#endif

#ifdef SAVEFILE
        double *plot;
        int m2 = m + 2;
        plot = (double *)malloc(n * m2 * sizeof(double));
        for (int i = 0; i < n; i++)
        {
            *(plot + i * m2) = i * h;
            *(plot + i * m2 + 1) = cos(PI / 2 * i * h) + 2 * sin(PI / 2 * i * h) - 1;
            for (int ij = 0; ij < m; ij++)
            {
                *(plot + i * m2 + ij + 2) = *(y + i * m + ij);
            }
        }
        char filename[256] = "result_tmp5.txt";
        printf("RungeKutta算法的结果保存在文件%s中 \n", filename);
        savef(filename, plot, n, m2, m + 1);
        free(plot);
#endif
        free(y);
    }

    printf("================求解微分方程 例五==============\n");
    {
        int n = 100, m = 2;
        double h = 1.0 / n;
        double y0[] = {0, 0};
        double y1[] = {1, 0};
        double *y;
        y = (double *)malloc(n * m * sizeof(double));
        LinearDeqSolve(shootf, y0, y1, y, n, 1, h);
#ifdef SHOW
        printf("RungeKutta算法的结果: \n");
        for (int i = 0; i < n; i++)
        {
            printf("%lf    %lf    %lf \n", i * h, *(y + m * i), *(y + m * i + 1));
        }
#endif

#ifdef SAVEFILE
        double *plot;
        int m2 = m + 2;
        plot = (double *)malloc(n * m2 * sizeof(double));
        for (int i = 0; i < n; i++)
        {
            *(plot + i * m2) = i * h;
            *(plot + i * m2 + 1) = cos(PI / 2 * i * h) + 2 * sin(PI / 2 * i * h) - 1;
            for (int ij = 0; ij < m; ij++)
            {
                *(plot + i * m2 + ij + 2) = *(y + i * m + ij);
            }
        }
        char filename[256] = "result_tmp6.txt";
        printf("线性常微分方程的解保存在%s中 \n", filename);
        savef(filename, plot, n, m2, m + 1);
        free(plot);
#endif
        free(y);
    }

    printf("================求解微分方程 例六==============\n");
    {
        double l;
        l = Secant(SL_f, 0.5, 1e-6);
        printf("求解勒让德方程的本征值为 %lf \n", l);
    }

    printf("================求解微分方程 例七==============\n");
    {
        double u0[2] = {0, 0.01};
        double x0[2] = {-10, 10};
        int n = 501;
        double pars[2] = {1, 4};
        double u[n];
        double x[n];
        double en0[2];
        double en;
        printf("请输入初始能量:\n");
        scanf("%lf,%lf",&en0[0],&en0[1]);
        printf("初始能量为:%lf, %lf \n",en0[0],en0[1]);
        double h = (x0[1] - x0[0]) / (n - 1);
        for (int i = 0; i < n; i++)
        {
            x[i] = x0[0] + i * h;
            //   printf("v = %lf \n",schroedingerV(x[i],pars));
        }

        SolveSchroedingerEQ(schroedingerV2, pars, u0, u, en0, en, n, x0);
        for (int i = 0; i < 10; i++)
        {
            printf("能级 %d 的解析结果为 %lf \n", i, 1.0 / 2.0 * pars[0] * pars[0] * (pars[1] * (pars[1] - 1) / 2 - pow(pars[1] - 1 - i, 2)));
        }

#ifdef SHOW
        printf("RungeKutta算法的结果: \n");
        for (int i = 0; i < n; i++)
        {
            printf("%lf    %lf    %lf \n", i * h, *(y + m * i), *(y + m * i + 1));
        }
#endif

#ifdef SAVEFILE

        double *plot;
        int m2 = 3;
        plot = (double *)malloc(n * m2 * sizeof(double));
        for (int i = 0; i < n; i++)
        {
            *(plot + i * m2) = x[i];
            *(plot + i * m2 + 1) = u[i];
            *(plot + i * m2 + 2) = schroedingerV(x[i], pars);
        }
        char filename[256] = "result_tmp7.txt";
        printf("求解薛定谔方程的结果保存在文件%s中 \n", filename);
        savef(filename, plot, n, m2, m2);
        free(plot);
#endif
    }

    return 0;
}
