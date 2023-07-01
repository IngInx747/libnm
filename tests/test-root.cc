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
    double x0 = -1;
    double x1 =  1;
    const double tol = 1e-10;
    int err;

    auto f0 = p2;
    auto f1 = dp2;

#if 0
    double xt; err = root_search_bisectional(f0, f1, x1, tol, xt);
#else
    double xt; err = root_search_newton_raphson(f0, f1, x0, x1, tol, xt);
#endif
    printf("err=%d, f(%lf) = %lf\n", err, xt, f0(xt));
}

int main(int argc, char **argv)
{
    test_root();

    return 0;
}