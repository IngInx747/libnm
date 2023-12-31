#ifndef NM_ROOT_SEARCH_HH
#define NM_ROOT_SEARCH_HH

#include "error.hh"

#ifndef ROOT_SEARCH_MAX_ITER_NUM
#define ROOT_SEARCH_MAX_ITER_NUM 64
#endif

template <class Func>
inline int root_search_bisectional(const Func &func, double x0, double x1, double tol, double &root)
{
    const int max_num_iter = ROOT_SEARCH_MAX_ITER_NUM;
    double xm, dy = 1e10;
    int num_iter {};

    for (; dy > tol && num_iter < max_num_iter; ++num_iter)
    {
        xm = (x0 + x1) * 0.5;
        double ym = func(xm);
        if (ym < 0)
            x0 = xm;
        else // ym >= 0
            x1 = xm;
        dy = ym>0 ? ym : -ym;
    }

    root = xm;

    int errcode = 0; // handle exceptions
    if (num_iter >= max_num_iter)
        errcode |= ITERATION_OVERLIMIT;

    return errcode;
}

template <class Func0, class Func1>
inline int root_search_newton_raphson(const Func0 &f0, const Func1 &f1, const double x0, const double x1, const double tol, double &root)
{
    const int max_num_iter = ROOT_SEARCH_MAX_ITER_NUM;
    int num_iter {};

    // initial guess
    double xr = (root<x0 || root>x1) ? (x0+x1)*.5 : root;
    double dy = 1e10;

    for (; num_iter < max_num_iter; ++num_iter)
    {
        double y = f0(xr);
        dy = y>0 ? y : -y;
        if (dy < tol) break;
        double y1 = f1(xr);
        xr = xr - y / y1;
        xr = xr<x0 ? x0 : xr>x1 ? x1 : xr;
    }

    root = xr;

    int errcode = 0; // handle exceptions
    if (num_iter >= max_num_iter)
        errcode |= ITERATION_OVERLIMIT;

    return errcode;
}

#endif