#ifndef GLCSEC2
#define GLCSEC2
double aitken(double x, int n, double xi[], double fi[]);
double updown(double x, int n, double xi[], double fi[]);
void orthogonalPolynomialFit(int m, int n, double x[], double f[], double* pu, double alpha[]);
#endif
