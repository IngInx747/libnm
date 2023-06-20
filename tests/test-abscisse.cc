#include "integral.hh"
#include "root.hh"
#include <stdio.h>
#include <cmath>
#include <vector>

static double normal_distribution(double x)
{
    static const double KNrm = std::sqrt(std::acos(-1));
    static const double kSgm = 100;
    x -= 0.5;
    return std::exp(-x*x*kSgm*kSgm) * kSgm / KNrm;
}

static int binary_search_sorted_data(const double *vals, const int n, const double target)
{
    int i, j, m;

    if (vals[0] > target)
        return -1;
    if (vals[n-1] <= target)
        return n;

    i = 0;
    j = n - 1;

    while (i < j-1)
    {
        m = i + (j - i) / 2;
        if (vals[m] > target)
            j = m;
        else
            i = m;
    }

    return i;
}

static double abscisse_interpole_gaussian_quadrature(double (*fp)(double), double x0, double x1, const int n, double st, double tol)
{
    double xm, xs = x0, ds { 1e10 };
    const int max_num_iter = 1024;
    int num_iter {};

    for (; ds > tol && num_iter < max_num_iter; ++num_iter)
    {
        xm = (x0 + x1) * .5;
        double s = integrate_gaussian_quadrature(fp, xs, xm, n);
        if (s < st)
            x0 = xm;
        else
            x1 = xm;
        ds = s - st; ds = ds>0 ? ds : -ds;
    }

    return xm;
}

struct integral_function_data
{
    double (*fp)(double); // integrand
    double xo; // x offset
    double yo; // y offset
    int n;
};

static double integral_function(double x, void *ext)
{
    const auto &I = *(integral_function_data*)(ext);
    double s = integrate_gaussian_quadrature(I.fp, I.xo, x, I.n);
    return s - I.yo;
}

static double derivative_function(double x, void *ext)
{
    const auto &I = *(integral_function_data*)(ext);
    return I.fp(x);
}

static void test_abscisse()
{
    auto fp = normal_distribution;
    double x0 = 0;
    double x1 = 1;
    double t = 0.75;
    const double tol = 1e-10;

    std::vector<double> xs(4096), ss(4096);
    double s {}; int nd = 4096; // number of intervals

    int err = integrate_adaptive_bisection(fp, x0, x1, tol, s, xs.data(), ss.data(), nd);
    if (err) { printf("error = %d\n", err); return; }

    xs[nd] = x1;
    ss[nd] = s;
    ++nd; // number of points

    for (int i = 0; i < nd; ++i)
    {
        printf("x = %.10lf, s = %.10lf\n", xs[i], ss[i]);
    }
    printf("number of points: %d\n", nd);

    int id = binary_search_sorted_data(ss.data(), nd, t*s);
    if (id < 0)
        printf("For t = %lf, x in (-inf, %lf)\n", t, xs[0]);
    else if (id >= nd)
        printf("For t = %lf, x in [%lf, +inf)\n", t, xs[nd-1]);
    else
        printf("For t = %lf, x in [%lf, %lf)\n", t, xs[id], xs[id+1]);

    if (id >= 0 && id < nd)
    {
        double x2 = xs[id], x3 = xs[id+1];
        double sg = (s*t - ss[id]); // target integral
#if 0
        double xt = abscisse_interpole_gaussian_quadrature(fp, x2, x3, 7, sg, tol);
#else
        integral_function_data I { fp, x2, sg, 7 };
    #if 0
        double xt; err = root_search_bisectional(integral_function, (void*)(&I), x2, x3, tol, xt);
    #else
        double xt;// = (x2+x3)*.5;
        err = root_search_newton_raphson(
            integral_function, (void*)(&I),
            derivative_function, (void*)(&I),
            x2, x3, tol, xt);
    #endif
#endif
        printf("For t = %lf, x = %lf, r(x) = %lf\n", t, xt, (xt-x0)/(x1-x0));
        double st; err = integrate_adaptive_bisection(fp, x0, xt, tol, st);
        printf("For x = %lf, s = %lf, r(s) = %lf\n", xt, st, st/s);
    }
}

int main(int argc, char **argv)
{
    test_abscisse();

    return 0;
}