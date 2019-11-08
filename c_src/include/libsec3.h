#ifndef GLCSEC3
#define GLCSEC3
void FirstOrderDerivative(double h, double f[], double df[], int n, int k);
void SecondOrderDerivative(double h, double f[], double ddf[], int n, int k);
double FODAdaptive(double x, double (*f)(double), double h, double del);
double Simpson(double (*f)(double), double a, double b, long int n);
double SimpsonArray(double y[], long int n, double h);
double SimpsonAdaptive(double (*f)(double), double a, double b, double del);
double SimpsonAdaptive2(double (*f)(double), double a, double b, double del,int step);
double Bisect(double (*f)(double),double a,double b,double del);
double Newton(double (*f)(double), double (*pf)(double), double x0, double del);
double Secant(double (*f)(double), double x0, double del);
#endif