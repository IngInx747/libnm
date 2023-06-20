#include "root.hh"
#include <stdio.h>
#include <cmath>

static double p2(double x)
{
    x -= 0.5;
    return x*x;
}

static double dp2(double x)
{
    x -= 0.5;
    return x*2;
}

static double p5(double x)
{
    x -= 0.5;
    return x*x*x*x*x;
}

static double dp5(double x)
{
    x -= 0.5;
    return x*x*x*x*5;
}

static void test_root()
{
    auto fp = p2;
    auto dp = dp2;
    double x0 = -1;
    double x1 = 1;
    const double tol = 1e-10;
    int err;

#if 0
    double xt; err = root_search_bisectional(fp, x0, x1, tol, xt);
#else
    double xt; err = root_search_newton_raphson(fp, dp, x0, x1, tol, xt);
#endif
    printf("root(%d), x = %lf, f(x) = %lf\n", err, xt, fp(xt));
}

int main(int argc, char **argv)
{
    test_root();

    return 0;
}