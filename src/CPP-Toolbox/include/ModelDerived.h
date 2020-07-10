#ifndef SRC_MODELCHILDREN_H_
#define SRC_MODELCHILDREN_H_

#include <cmath>

#include "Model.h"

#define CONCAT(arg) Model ## arg
#define CLASSNAME(arg) CONCAT(arg)

#ifdef DER_CPPAD
#include "cppad/cppad.hpp"
#endif

#if __cplusplus < 201103L
#define auto double
#define thread_local
#define constexpr const
#endif

// __restrict__ is called __restrict in VS
#if _MSC_VER
#define __restrict__ __restrict
#endif

using std::exp;
using std::log;
using std::pow;
#define min(x,y) ( x < y ? x : y )
#define max(x,y) ( x > y ? x : y )
constexpr double NaN=NAN;      // default MoBi output
constexpr double Inf=INFINITY; // default MoBi output

struct CLASSNAME(MODELNAME) : Model {

    CLASSNAME(MODELNAME)();

    virtual void ODEOptions(double &absTol, double &relTol, double &hInit, double &hMin, double &hMax, long int &maxSteps, bool &useJac, const double *__restrict__ P);
    virtual void ODEInitialValues(double* y, double* P, unsigned int *S);
    virtual void ODEInitialParameters(double* P);
    virtual void ODERHSFunction(const double Time, const double* y, const double* P, const unsigned int *S, double* dy);
    virtual void ODEObservers(const double Time, const double *y, const double *P, const unsigned int *S, double *obs);
    virtual bool ODEExplicitSwitch(const double Time, double *y, double *P, unsigned int *S);
    virtual void ODEImplicitSwitch(const double Time, double* y, double* P, unsigned int *S, double* gout);
    virtual void ODEJacDense    (const double Time, const double* y, const double* P, const unsigned int *S, double* J);
    virtual void ODEJacDense_fp (const double Time, const double* y, const double* P, const unsigned int *S, double* J);
    virtual void ODEJacDense_ox (const double Time, const double* y, const double* P, const unsigned int *S, double* J);
    virtual void ODEJacDense_op (const double Time, const double* y, const double* P, const unsigned int *S, double* J);
    virtual void ODEJacSparse   (const double Time, const double* y, const double *P, const unsigned int *S, double *J);
    virtual void ODEJacSparse_fp(const double Time, const double* y, const double *P, const unsigned int *S, double *J);
    virtual void ODEJacSparse_ox(const double Time, const double* y, const double *P, const unsigned int *S, double *J);
    virtual void ODEJacSparse_op(const double Time, const double* y, const double *P, const unsigned int *S, double *J);

#ifdef DER_CPPAD
    virtual void ODERHSFunction(CppAD::AD<double> &Time, CppAD::AD<double> *y, const double *P, CppAD::AD<double> *dy);
#endif
};

#endif // SRC_MODELCHILDREN_H_
