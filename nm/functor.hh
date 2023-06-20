#ifndef NM_FUNCTOR_HH
#define NM_FUNCTOR_HH

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

#endif