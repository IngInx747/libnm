#ifndef NM_INTEGRAL_HH
#define NM_INREGARL_HH

#include "error.hh"

#ifndef INTEGRAL_MAX_ITER_NUM
#define INTEGRAL_MAX_ITER_NUM 65536
#endif

#ifndef INTEGRAL_MAX_STACK_DEPTH
#define INTEGRAL_MAX_STACK_DEPTH 64
#endif

#ifndef INTEGRAL_TRACKER_BUFFER_SIZE
#define INTEGRAL_TRACKER_BUFFER_SIZE 1024
#endif

inline double gaussian_quadrature_point(const int i)
{
    const double __gq7_x[] = {
        -0.9491079123427585245262,
        -0.7415311855993944398639,
        -0.4058451513773971669066,
         0,
         0.4058451513773971669066,
         0.7415311855993944398639,
         0.9491079123427585245262,
    };

    return __gq7_x[i];
}

inline double gaussian_quadrature_weight(const int i)
{
    const double __gq7_w[] = {
        0.1294849661688696932706,
        0.2797053914892766679015,
        0.38183005050511894495,
        0.417959183673469387755,
        0.38183005050511894495,
        0.2797053914892766679015,
        0.1294849661688696932706,
    };

    return __gq7_w[i];
}

inline double gauss_kronrod_quadrature_point(const int i)
{
    const double __gkq15_x[] = {
        -0.9914553711208126392069,
        -0.9491079123427585245262,
        -0.8648644233597690727897,
        -0.7415311855993944398639,
        -0.5860872354676911302941,
        -0.4058451513773971669066,
        -0.2077849550078984676007,
         0,
         0.2077849550078984676007,
         0.4058451513773971669066,
         0.5860872354676911302941,
         0.7415311855993944398639,
         0.8648644233597690727897,
         0.9491079123427585245262,
         0.9914553711208126392069,
    };

    return __gkq15_x[i];
}

inline double gauss_kronrod_quadrature_weight(const int i)
{
    const double __gkq15_w[] = {
        0.0229353220105292249637,
        0.0630920926299785532907,
        0.1047900103222501838399,
        0.1406532597155259187452,
        0.1690047266392679028266,
        0.1903505780647854099133,
        0.2044329400752988924142,
        0.209482141084727828013,
        0.2044329400752988924142,
        0.1903505780647854099133,
        0.1690047266392679028266,
        0.1406532597155259187452,
        0.1047900103222501838399,
        0.0630920926299785532907,
        0.0229353220105292249637,
    };

    return __gkq15_w[i];
}

template <class Func>
inline double integrate_g7(const Func &func, double x0, double x1)
{
    double sum {};
    const double scale  = (x1 - x0) * .5;
    const double offset = (x0 + x1) * .5;

    for (int i = 0; i < 7; ++i)
    {
        const double t = gaussian_quadrature_point (i);
        const double w = gaussian_quadrature_weight(i);
        const double x = t*scale + offset;
        const double y = func(x);
        sum += y*w;
    }

    return sum * scale;
}

template <class Func>
inline double integrate_g7_k15(const Func &func, double x0, double x1, double &err)
{
    double sg {}, sk {};
    const double scale  = (x1 - x0) * .5;
    const double offset = (x0 + x1) * .5;
    double ys[15];

    for (int i = 0; i < 15; ++i)
    {
        const double t = gauss_kronrod_quadrature_point (i);
        const double w = gauss_kronrod_quadrature_weight(i);
        const double x = t*scale + offset;
        const double y = func(x);
        sk += y*w;
        ys[i] = y;
    }

    for (int i = 0; i < 7; ++i)
    {
        const double w = gaussian_quadrature_weight(i);
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

template <class Func, class Eval, class Tracker>
inline int integrate_adaptive_bisection(const Func &func, const Eval &eval, Tracker &track, double x0, double x1, double tol, double &result)
{
    const long long max_stack_depth = INTEGRAL_MAX_STACK_DEPTH; // up to 2^-(n-1) subdivided intervals
    const long long max_num_iter = INTEGRAL_MAX_ITER_NUM; // Warning: cap being too small will cause incomplete result

    double xs0[max_stack_depth];
    double xs1[max_stack_depth];

    double sum {}, err {};
    double xl = x0;
    double xr = x1;

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
        errcode |= ITERATION_OVERLIMIT;

    return errcode;
}

struct DummyTracker
{
    inline void operator() (double) {}
};

struct AdaptiveBisectionTracker
{
    inline void operator()(double val)
    {
        if (n >= cap - 1) return;
        (i%2 == 0) ? xs[n] = val : ss[n] = val;
        (i%2 != 0) ? ++n : 0;
        i = (i + 1) % 2;
    }

    double *xs; // xl
    double *ss; // sum
    int cap { INTEGRAL_TRACKER_BUFFER_SIZE };
    int n {}, i {};
};

template <class Func, class Tracker>
inline int integrate_adaptive_bisection(const Func &func, Tracker &track, double x0, double x1, double tol, double &result)
{
    G7K15IntegralEvaluation<Func> eval {};
    return integrate_adaptive_bisection(func, eval, track, x0, x1, tol, result);
}

/// @brief Integrate the function
/// @tparam Func the definition of the integrand(function pointer or structure)
/// @param func the integrand
/// @param x0 domain lower bound
/// @param x1 domain upper bound
/// @param tol precision threshold
/// @param result value of the integral
/// @return error code
template <class Func>
inline int integrate_adaptive_bisection(const Func &func, double x0, double x1, double tol, double &result)
{
    DummyTracker tracker {};
    return integrate_adaptive_bisection(func, tracker, x0, x1, tol, result);
}

/// @brief Integrate the function and build quadrature-trustfully intervals
/// @tparam Func the definition of the integrand(function pointer or structure)
/// @param func the integrand
/// @param x0 domain lower bound
/// @param x1 domain upper bound
/// @param tol precision threshold
/// @param xs X coordinates of the integral [x0, ..., x1]
/// @param ys Y coordinates of the integral [0, ..., I]
/// @param n In: cap of coordinate array. Out: number of bisectional points
/// @return error code
template <class Func>
inline int integrate_adaptive_bisection(const Func &func, double x0, double x1, double tol, double *xs, double *ys, int &n)
{
    AdaptiveBisectionTracker tracker { xs, ys, n };
    double s; int err = integrate_adaptive_bisection(func, tracker, x0, x1, tol, s);
    n = tracker.n; xs[n] = x1; ys[n] = s; ++n;
    return err;
}

#endif