#include "integral.hh"
#include <stdio.h>
#include <cmath>

static double p2(double x)
{
    return x*x;
}

static double normal_distribution(double x)
{
    static const double KNrm = std::sqrt(std::acos(-1));
    static const double kSgm = 100;
    x -= 0.5;
    return std::exp(-x*x*kSgm*kSgm) * kSgm / KNrm;
}

// 0.2108027355 in [0,1]
static double sech6(double x)
{
    double csh1 = std::cosh((x-0.2)*10);
    double csh2 = std::cosh((x-0.4)*100);
    double csh3 = std::cosh((x-0.6)*1000);
    return std::pow(csh1, -2) + std::pow(csh2, -4) + std::pow(csh3, -6);
}

// 0.34740017 in [0,8]
static double siegfried_rump_osc_1(double x)
{
    return std::sin(x + std::exp(x));
}

// 0.0986517 in [0,8]
static double siegfried_rump_osc_2(double x)
{
    double ex = std::exp(x);
    return std::sin(x + ex) * (ex - std::floor(ex));
}

// 0.504067 in [0,1]
static double oscillation_1(double x)
{
    return std::sin(1. / x);
}

// 0.3785300171242 in [0,1]
static double oscillation_2(double x)
{
    return std::sin(1. / x) * x;
}

// 3.4365636569 in [0,1]
static double helfgott(double x)
{
    double x2 = x*x, x3 = x*x*x, x4 = x*x*x*x;
    double xt = x4 + x3*10 + x2*19 - x*6 - 6;
    return xt * std::exp(x);
}

// 5050 in [0,1]
static double gauss_sum(double x)
{
    return std::floor(x*100 + 1) * 100;
}

static void test_integral()
{
    printf("sqrt(pi): %.20lf\n", std::sqrt(std::acos(-1)) *.5);

    auto fp = normal_distribution;
    double x0 = 0;
    double x1 = 1;

    for (int i = 1; i <= 10; ++i)
    {
        double d = integrate_gaussian_quadrature(fp, x0, x1, i);
        printf("gq: %.10lf\n", d);
    }

    for (int i = 1; i <= 20; ++i)
    {
        double eps = std::pow(0.1, i);
        double d; int err = integrate_adaptive_bisection(fp, x0, x1, eps, d);
        printf("ba(%d): %.20lf\n", err, d);
    }
}

int main(int argc, char **argv)
{
    test_integral();

    return 0;
}