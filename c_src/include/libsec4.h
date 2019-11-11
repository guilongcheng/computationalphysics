#ifndef GLCSEC4
#define GLCSEC4
#include "GLC.h"
void Eular(void (*g)(double, double[MMAX], double[MMAX]), double y0[MMAX], double *y, int m, int n, double h);
void PreEular(void (*g)(double, double[MMAX], double[MMAX]), double y0[MMAX], double *y, int m, int n, double h);
void TwoPointPre(void (*g)(double, double[MMAX], double[MMAX]), double y0[MMAX], double y1[MMAX], double *y, int m, int n, double h);
void RungeKutta(void (*g)(double, double[MMAX], double[MMAX]), double y0[MMAX], double *y, int m, int n, double h);
void Shooting(void (*g)(double, double[2], double[2]), double y0[2],double y1[2],double *y,int n, int k, double h);
void LinearDeqSolve(void (*g)(double, double[2], double[2]), double y0[2],double y1[2],double *y, int n, int k, double h);
void sturmLiouville(double p[], double p1[], double q[], double s[], double u0[2], double *u, int n, double h);
void wave(double (*v)(double,double[NPARS]),double pars[NPARS],double u0[2],double u[], double x0[2] ,double en,int n);
void SolveSchroedingerEQ(double (*v)(double,double[NPARS]),double pars[NPARS],double u0[2],double u[],double e0[2],double en,int n,double x0[2]);
#endif