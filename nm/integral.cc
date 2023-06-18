#include "integral.hh"
#include "quadrature.hh"
//#include <stdio.h>

////////////////////////////////////////////////////////////////
/// Function decorator
////////////////////////////////////////////////////////////////

struct SingleVariableFunction
{
    inline double operator()(double x) const { return fp(x); }
    double (*fp)(double);
};

struct ExtendedSingleVariableFunction
{
    inline double operator()(double x) const { return fp(x, handle); }
    double (*fp)(double, void*);
    void *handle;
};

////////////////////////////////////////////////////////////////
/// Stack data tracker
////////////////////////////////////////////////////////////////

struct DummyTracker
{
    inline void operator() (double) {}
};

struct LeafNodeDataTracker
{
    inline void operator() (double val)
    { fp(val, handle); }

    void (*fp)(double, void*);
    void *handle;
};

////////////////////////////////////////////////////////////////
/// Gaussian quadrature
////////////////////////////////////////////////////////////////

template <class Func>
inline double integrate_gaussian_quadrature(const Func &func, double x0, double x1, const int n)
{
    double sum {};
    const double scale  = (x1 - x0) * .5;
    const double offset = (x0 + x1) * .5;

    for (int i = 0; i < n; ++i)
    {
        const double t = gaussian_quadrature_point (n, i);
        const double w = gaussian_quadrature_weight(n, i);
        const double x = t*scale + offset;
        const double y = func(x);
        sum += y*w;
    }

    return sum * scale;
}

double integrate_gaussian_quadrature(double (*fp)(double), double x0, double x1, const int n)
{
    SingleVariableFunction func { fp };
    return integrate_gaussian_quadrature(func, x0, x1, n);
}

double integrate_gaussian_quadrature(double (*fp)(double, void*), void *handle, double x0, double x1, const int n)
{
    ExtendedSingleVariableFunction func { fp, handle };
    return integrate_gaussian_quadrature(func, x0, x1, n);
}

////////////////////////////////////////////////////////////////
/// Evaluation decorator
////////////////////////////////////////////////////////////////

template <class Func>
inline double integrate_middle_point(const Func &func, double x0, double x1, double &err)
{
    double xm = (x0 + x1) * .5;
    double y0 = func(x0);
    double y1 = func(x1);
    double ym = func(xm);
    double e = (y1 + y0) * .5 - ym;
    err = e>0 ? e : -e;
    return ym * (x1 - x0);
}

template <class Func>
struct MidPointIntegralEvaluation
{
    inline double operator() (const Func &func, double x0, double x1, double &err) const
    { return integrate_middle_point(func, x0, x1, err); }
};

template <class Func>
inline double integrate_g7_k15(const Func &func, double x0, double x1, double &err)
{
    double sg {}, sk {};
    const double scale  = (x1 - x0) * .5;
    const double offset = (x0 + x1) * .5;
    double ys[15];

    for (int i = 0; i < 15; ++i)
    {
        const double t = gauss_kronrod_quadrature_point (15, i);
        const double w = gauss_kronrod_quadrature_weight(15, i);
        const double x = t*scale + offset;
        const double y = func(x);
        sk += y*w;
        ys[i] = y;
    }

    for (int i = 0; i < 7; ++i)
    {
        //const double t = gaussian_quadrature_point (7, i);
        const double w = gaussian_quadrature_weight(7, i);
        const double y = ys[i*2 + 1];
        sg += y*w;
    }

    double e = (sg - sk) * scale;
    err = e>0 ? e : -e;

    return sg * scale;
}

template <class Func>
struct G7K15IntegralEvaluation
{
    inline double operator() (const Func &func, double x0, double x1, double &err) const
    { return integrate_g7_k15(func, x0, x1, err); }
};

////////////////////////////////////////////////////////////////
/// Binary adaptive
////////////////////////////////////////////////////////////////

template <class Func, class Eval>
inline double integrate_binary_adaptive_recursive(const Func &func, const Eval &eval, double x0, double x1, double tol, int level)
{
    double sum, err;
    sum = eval(func, x0, x1, err);

    //printf("%d\n",level);

    if (level < 32 && err > tol)
    {
        double xm = (x0 + x1) * .5;
        sum =
            integrate_binary_adaptive_recursive(func, x0, xm, tol, level + 1)+
            integrate_binary_adaptive_recursive(func, xm, x1, tol, level + 1);
    }

    return sum;
}

template <class Func, class Eval, class Tracker>
inline int integrate_binary_adaptive(const Func &func, const Eval &eval, Tracker &track, double x0, double x1, double tol, double &result)
{
    const unsigned long long max_stack_depth = 64; // up to 2^-(n-1) subdivided intervals
    const unsigned long long max_num_iter = 65536; // Warning: cap being too small will cause incomplete result

    double xs0[max_stack_depth];
    double xs1[max_stack_depth];

    double sum {}, err {};
    double xl = x0;
    double xr = x1;

#if 0 // do-while works but looks ugly
    int stack_top = 0;
    int num_iter {};

    do
    {
        double s = eval(func, xl, xr, err);

        //printf("%d\n",stack_top);
        //printf("%d\n",num_iter);

        if (stack_top < max_stack_depth && err > tol) // push right-half and make left-half current
        {
            double xm = (xl + xr) * .5;
            xs0[stack_top] = xm;
            xs1[stack_top] = xr;
            ++stack_top;
            xr = xm;
        }
        else // accept current result and make top(popped) current
        {
            sum += s;
            --stack_top;
            xl = xs0[stack_top];
            xr = xs1[stack_top];
        }

        ++num_iter;
    }
    while (stack_top >= 0 && num_iter < max_num_iter);
#endif

    long long top = 0;
    long long max_top = 0;
    long long num_iter {};

    for (; top >= 0 && num_iter < max_num_iter; ++num_iter)
    {
        double s = eval(func, xl, xr, err);

        if (max_top < top)
            max_top = top;

        //printf("%zd\n", top);
        //printf("%zd\n", num_iter);
        //printf("%lf\n", err);

        if (top < max_stack_depth && err > tol) // error exceeds tolerance
        {
            double xm = (xl + xr) * .5;
            xs0[top] = xm;
            xs1[top] = xr;
            xr = xm;
            ++top; // push right-half and make left-half current
        }
        else // result is accepted or too many levels
        {
            track(xl); track(sum); // track leaf node data(optional)
            --top;
            sum += s;
            xl = xs0[top];
            xr = xs1[top]; // pop stack top and make it current
        }
    }

    result = sum; // write out

    int errcode = 0; // handle exceptions
    if (max_top >= max_stack_depth)
        errcode |= STACK_OVERFLOW;
    if (num_iter >= max_num_iter)
        errcode |= INCOMPLETE_INTEGRAL;

    return errcode;
}

int integrate_binary_adaptive(double (*fp)(double), double x0, double x1, double tol, double &result)
{
    DummyTracker tracker {};
    SingleVariableFunction func { fp };
    G7K15IntegralEvaluation<decltype(func)> eval {};
    return integrate_binary_adaptive(func, eval, tracker, x0, x1, tol, result);
}

int integrate_binary_adaptive(double (*fp)(double, void*), void *ext, double x0, double x1, double tol, double &result)
{
    DummyTracker tracker {};
    ExtendedSingleVariableFunction func { fp, ext };
    G7K15IntegralEvaluation<decltype(func)> eval {};
    return integrate_binary_adaptive(func, eval, tracker, x0, x1, tol, result);
}

struct BinaryAdaptiveMemory
{
    inline void track(double val)
    {
        if (n >= cap) return;
        (i%2 == 0) ? xs[n] = val : ss[n] = val;
        (i%2 != 0) ? ++n : 0;
        i = (i + 1) % 2;
    }

    double *xs; // xl
    double *ss; // sum
    int cap { 1024 }, n {}, i {};
};

static void binary_adaptive_track_callback(double val, void *mem)
{
    using T_Tr = BinaryAdaptiveMemory;
    auto &data = (*(static_cast<T_Tr*>(mem)));
    data.track(val);
}

int integrate_binary_adaptive(double (*fp)(double), double x0, double x1, double tol, double &result, double *xs, double *ss, int &nd)
{
    SingleVariableFunction func { fp };
    G7K15IntegralEvaluation<decltype(func)> eval {};
    BinaryAdaptiveMemory mem { xs, ss, nd };
    LeafNodeDataTracker tracker { binary_adaptive_track_callback, (void*)(&mem) };
    int err = integrate_binary_adaptive(func, eval, tracker, x0, x1, tol, result);
    nd = mem.n;
    return err;
}

int integrate_binary_adaptive(double (*fp)(double, void*), void *ext, double x0, double x1, double tol, double &result, double *xs, double *ss, int &nd)
{
    ExtendedSingleVariableFunction func { fp, ext };
    G7K15IntegralEvaluation<decltype(func)> eval {};
    BinaryAdaptiveMemory mem { xs, ss, nd };
    LeafNodeDataTracker tracker { binary_adaptive_track_callback, (void*)(&mem) };
    int err = integrate_binary_adaptive(func, eval, tracker, x0, x1, tol, result);
    nd = mem.n;
    return err;
}
