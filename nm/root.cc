#include "root.hh"
#include "error.hh"
#include "functor.hh"

template <class Func>
int root_search_bisectional(const Func &func, double x0, double x1, double tol, double &root)
{
    const int max_num_iter = 4096;
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

int root_search_bisectional(double (*fp)(double), double x0, double x1, double tol, double &root)
{
    SingleVariableFunction func { fp };
    return root_search_bisectional(func, x0, x1, tol, root);
}

int root_search_bisectional(double (*fp)(double, void*), void *ext, double x0, double x1, double tol, double &root)
{
    ExtendedSingleVariableFunction func { fp, ext };
    return root_search_bisectional(func, x0, x1, tol, root);
}

template <class Func0, class Func1>
int root_search_newton_raphson(const Func0 &f0, const Func1 &f1, double x0, double x1, double tol, double &root)
{
    const int max_num_iter = 4096;
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

int root_search_newton_raphson(
    double (*fp0)(double),
    double (*fp1)(double),
    double x0, double x1, double tol, double &root)
{
    SingleVariableFunction func { fp0 };
    SingleVariableFunction drvt { fp1 };
    return root_search_newton_raphson(func, drvt, x0, x1, tol, root);
}

int root_search_newton_raphson(
    double (*fp0)(double, void*), void *ext0,
    double (*fp1)(double, void*), void *ext1,
    double x0, double x1, double tol, double &root)
{
    ExtendedSingleVariableFunction func { fp0, ext0 };
    ExtendedSingleVariableFunction drvt { fp1, ext1 };
    return root_search_newton_raphson(func, drvt, x0, x1, tol, root);
}
