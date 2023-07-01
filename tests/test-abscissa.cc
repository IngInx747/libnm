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

struct Integral
{
    inline double operator()(double x) const
    {
        double s = integrate_g7(func, xo, x);
        return s - yo;
    }

    double (*func)(double); // integrand
    double xo; // x offset
    double yo; // y offset
};

static void test_abscissa()
{
    auto fp = normal_distribution;
    double x0 = 0;
    double x1 = 1;
    double t = 0.75;
    const double tol = 1e-10;

    int nd = 4096; // cap
    std::vector<double> xs(nd), ys(nd);

    // Integrate and build precision-truct intervals
    double s; int err = integrate_adaptive_bisection(fp, x0, x1, tol, s, xs.data(), ys.data(), nd);
    if (err) { printf("error = %d\n", err); return; }

    for (int i = 0; i < nd; ++i)
    {
        printf("x = %.10lf, s = %.10lf\n", xs[i], ys[i]);
    }
    printf("number of points: %d\n", nd);

    int id = binary_search_sorted_data(ys.data(), nd, t*s);
    if (id < 0)
        printf("For t = %lf, x in (-inf, %lf)\n", t, xs[0]);
    else if (id >= nd)
        printf("For t = %lf, x in [%lf, +inf)\n", t, xs[nd-1]);
    else
        printf("For t = %lf, x in [%lf, %lf)\n", t, xs[id], xs[id+1]);

    if (id >= 0 && id < nd)
    {
        double x2 = xs[id], x3 = xs[id+1];
        double yg = (s*t - ys[id]); // target integral
        Integral I { fp, x2, yg };
        #if 0
        double xt; err = root_search_bisectional(I, x2, x3, tol, xt);
        #else
        double xt; err = root_search_newton_raphson(I, fp, x2, x3, tol, xt);
        #endif
        double st; err = integrate_adaptive_bisection(fp, x0, xt, tol, st);
        printf("xt = %lf, xt/(x1-x0) = %lf\n", xt, (xt-x0)/(x1-x0));
        printf("It = %lf, It/I = %lf\n", st, st/s);
    }
}

int main(int argc, char **argv)
{
    test_abscissa();

    return 0;
}