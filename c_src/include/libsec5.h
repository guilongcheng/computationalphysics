#ifndef GLCSEC5
#define GLCSEC5
void Gaussian(double *a, int index[], int n);
double DetGauss(double *a, int n);
void LinearSolveGauss(double *a, double b[], double x[], int n);
void InvertGauss(double *a, double *x, int n);
#endif