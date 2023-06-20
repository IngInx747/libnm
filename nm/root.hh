#ifndef NM_ROOT_SEARCH_HH
#define NM_ROOT_SEARCH_HH

int root_search_bisectional(double (*)(double), double, double, double, double&);

int root_search_bisectional(double (*)(double, void*), void*, double, double, double, double&);

int root_search_newton_raphson(
    double (*)(double),
    double (*)(double),
    double, double, double, double&);

int root_search_newton_raphson(
    double (*)(double, void*), void*,
    double (*)(double, void*), void*,
    double, double, double, double&);

#endif