#ifndef NM_GAUSSIAN_QUADRATURE_HH
#define NM_GAUSSIAN_QUADRATURE_HH

double gaussian_quadrature_point (const int n, const int i);
double gaussian_quadrature_weight(const int n, const int i);

double gauss_kronrod_quadrature_point (const int n, const int i);
double gauss_kronrod_quadrature_weight(const int n, const int i);

#endif