#ifndef NUMERICAL_INTEGRAL_HH
#define NUMERICAL_INREGARL_HH

/// @brief Integral f(x) within [x0, x1]
/// @param  Function pointer of integrand
/// @param  Lower bound of intergration interval
/// @param  Upper bound of intergration interval
/// @param  Number of quadrature nodes
/// @return Integration result
double integrate_gaussian_quadrature(double (*)(double), double, double, int);

/// @brief Integral f(x, ...) within [x0, x1]
/// @param  Function pointer
/// @param  Extra data handle of integrand
/// @param  Lower bound of intergration interval
/// @param  Upper bound of intergration interval
/// @param  Number of quadrature nodes
/// @return Integration result
double integrate_gaussian_quadrature(double (*)(double, void*), void*, double, double, int);

enum INTERNAL_ERROR_CODE
{
    NO_ERROR = 0,
    STACK_OVERFLOW = 1,
    INCOMPLETE_INTEGRAL = 2,
};

/// @brief Integral f(x) within [x0, x1] with adaptive subdivided intervals
/// @param  Function pointer of integrand
/// @param  Lower bound of intergration interval
/// @param  Upper bound of intergration interval
/// @param  Precision tolerance
/// @param  Integration result
/// @return Error code
int integrate_binary_adaptive(double (*)(double), double, double, double, double&);

/// @brief Integral f(x, ...) within [x0, x1] with adaptive subdivided intervals
/// @param  Function pointer of integrand
/// @param  Extra data handle
/// @param  Lower bound of intergration interval
/// @param  Upper bound of intergration interval
/// @param  Precision tolerance
/// @param  Integration result
/// @return Error code
int integrate_binary_adaptive(double (*)(double, void*), void*, double, double, double, double&);

int integrate_binary_adaptive(double (*)(double), double, double, double, double&, double*, double*, int&);

int integrate_binary_adaptive(double (*)(double, void*), void*, double, double, double, double&, double*, double*, int&);

#endif