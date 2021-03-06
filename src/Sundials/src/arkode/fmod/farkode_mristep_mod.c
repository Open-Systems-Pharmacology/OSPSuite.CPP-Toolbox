/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.0
 *
 * This file is not intended to be easily readable and contains a number of
 * coding conventions designed to improve portability and efficiency. Do not make
 * changes to this file unless you know what you are doing--modify the SWIG
 * interface file instead.
 * ----------------------------------------------------------------------------- */

/* ---------------------------------------------------------------
 * Programmer(s): Auto-generated by swig.
 * ---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -------------------------------------------------------------*/

/* -----------------------------------------------------------------------------
 *  This section contains generic SWIG labels for method/variable
 *  declarations/attributes, and other compiler dependent labels.
 * ----------------------------------------------------------------------------- */

/* template workaround for compilers that cannot correctly implement the C++ standard */
#ifndef SWIGTEMPLATEDISAMBIGUATOR
# if defined(__SUNPRO_CC) && (__SUNPRO_CC <= 0x560)
#  define SWIGTEMPLATEDISAMBIGUATOR template
# elif defined(__HP_aCC)
/* Needed even with `aCC -AA' when `aCC -V' reports HP ANSI C++ B3910B A.03.55 */
/* If we find a maximum version that requires this, the test would be __HP_aCC <= 35500 for A.03.55 */
#  define SWIGTEMPLATEDISAMBIGUATOR template
# else
#  define SWIGTEMPLATEDISAMBIGUATOR
# endif
#endif

/* inline attribute */
#ifndef SWIGINLINE
# if defined(__cplusplus) || (defined(__GNUC__) && !defined(__STRICT_ANSI__))
#   define SWIGINLINE inline
# else
#   define SWIGINLINE
# endif
#endif

/* attribute recognised by some compilers to avoid 'unused' warnings */
#ifndef SWIGUNUSED
# if defined(__GNUC__)
#   if !(defined(__cplusplus)) || (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#     define SWIGUNUSED __attribute__ ((__unused__))
#   else
#     define SWIGUNUSED
#   endif
# elif defined(__ICC)
#   define SWIGUNUSED __attribute__ ((__unused__))
# else
#   define SWIGUNUSED
# endif
#endif

#ifndef SWIG_MSC_UNSUPPRESS_4505
# if defined(_MSC_VER)
#   pragma warning(disable : 4505) /* unreferenced local function has been removed */
# endif
#endif

#ifndef SWIGUNUSEDPARM
# ifdef __cplusplus
#   define SWIGUNUSEDPARM(p)
# else
#   define SWIGUNUSEDPARM(p) p SWIGUNUSED
# endif
#endif

/* internal SWIG method */
#ifndef SWIGINTERN
# define SWIGINTERN static SWIGUNUSED
#endif

/* internal inline SWIG method */
#ifndef SWIGINTERNINLINE
# define SWIGINTERNINLINE SWIGINTERN SWIGINLINE
#endif

/* qualifier for exported *const* global data variables*/
#ifndef SWIGEXTERN
# ifdef __cplusplus
#   define SWIGEXTERN extern
# else
#   define SWIGEXTERN
# endif
#endif

/* exporting methods */
#if defined(__GNUC__)
#  if (__GNUC__ >= 4) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#    ifndef GCC_HASCLASSVISIBILITY
#      define GCC_HASCLASSVISIBILITY
#    endif
#  endif
#endif

#ifndef SWIGEXPORT
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   if defined(STATIC_LINKED)
#     define SWIGEXPORT
#   else
#     define SWIGEXPORT __declspec(dllexport)
#   endif
# else
#   if defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#     define SWIGEXPORT __attribute__ ((visibility("default")))
#   else
#     define SWIGEXPORT
#   endif
# endif
#endif

/* calling conventions for Windows */
#ifndef SWIGSTDCALL
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   define SWIGSTDCALL __stdcall
# else
#   define SWIGSTDCALL
# endif
#endif

/* Deal with Microsoft's attempt at deprecating C standard runtime functions */
#if !defined(SWIG_NO_CRT_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_CRT_SECURE_NO_DEPRECATE)
# define _CRT_SECURE_NO_DEPRECATE
#endif

/* Deal with Microsoft's attempt at deprecating methods in the standard C++ library */
#if !defined(SWIG_NO_SCL_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_SCL_SECURE_NO_DEPRECATE)
# define _SCL_SECURE_NO_DEPRECATE
#endif

/* Deal with Apple's deprecated 'AssertMacros.h' from Carbon-framework */
#if defined(__APPLE__) && !defined(__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES)
# define __ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES 0
#endif

/* Intel's compiler complains if a variable which was never initialised is
 * cast to void, which is a common idiom which we use to indicate that we
 * are aware a variable isn't used.  So we just silence that warning.
 * See: https://github.com/swig/swig/issues/192 for more discussion.
 */
#ifdef __INTEL_COMPILER
# pragma warning disable 592
#endif

/*  Errors in SWIG */
#define  SWIG_UnknownError    	   -1
#define  SWIG_IOError        	   -2
#define  SWIG_RuntimeError   	   -3
#define  SWIG_IndexError     	   -4
#define  SWIG_TypeError      	   -5
#define  SWIG_DivisionByZero 	   -6
#define  SWIG_OverflowError  	   -7
#define  SWIG_SyntaxError    	   -8
#define  SWIG_ValueError     	   -9
#define  SWIG_SystemError    	   -10
#define  SWIG_AttributeError 	   -11
#define  SWIG_MemoryError    	   -12
#define  SWIG_NullReferenceError   -13




#include <assert.h>
#define SWIG_exception_impl(DECL, CODE, MSG, RETURNNULL) \
 { printf("In " DECL ": " MSG); assert(0); RETURNNULL; }


#include <stdio.h>
#if defined(_MSC_VER) || defined(__BORLANDC__) || defined(_WATCOM)
# ifndef snprintf
#  define snprintf _snprintf
# endif
#endif


/* Support for the `contract` feature.
 *
 * Note that RETURNNULL is first because it's inserted via a 'Replaceall' in
 * the fortran.cxx file.
 */
#define SWIG_contract_assert(RETURNNULL, EXPR, MSG) \
 if (!(EXPR)) { SWIG_exception_impl("$decl", SWIG_ValueError, MSG, RETURNNULL); } 


#define SWIGVERSION 0x040000 
#define SWIG_VERSION SWIGVERSION


#define SWIG_as_voidptr(a) (void *)((const void *)(a)) 
#define SWIG_as_voidptrptr(a) ((void)SWIG_as_voidptr(*a),(void**)(a)) 


#include "arkode/arkode_mristep.h"


#include <stdlib.h>
#ifdef _MSC_VER
# ifndef strtoull
#  define strtoull _strtoui64
# endif
# ifndef strtoll
#  define strtoll _strtoi64
# endif
#endif


typedef struct {
    void* data;
    size_t size;
} SwigArrayWrapper;


SWIGINTERN SwigArrayWrapper SwigArrayWrapper_uninitialized() {
  SwigArrayWrapper result;
  result.data = NULL;
  result.size = 0;
  return result;
}


#include <string.h>

SWIGEXPORT void * _wrap_FMRIStepCreate(ARKRhsFn farg1, double const *farg2, N_Vector farg3, int const *farg4, void *farg5) {
  void * fresult ;
  ARKRhsFn arg1 = (ARKRhsFn) 0 ;
  realtype arg2 ;
  N_Vector arg3 = (N_Vector) 0 ;
  MRISTEP_ID arg4 ;
  void *arg5 = (void *) 0 ;
  void *result = 0 ;
  
  arg1 = (ARKRhsFn)(farg1);
  arg2 = (realtype)(*farg2);
  arg3 = (N_Vector)(farg3);
  arg4 = (MRISTEP_ID)(*farg4);
  arg5 = (void *)(farg5);
  result = (void *)MRIStepCreate(arg1,arg2,arg3,arg4,arg5);
  fresult = result;
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepResize(void *farg1, N_Vector farg2, double const *farg3, ARKVecResizeFn farg4, void *farg5) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  N_Vector arg2 = (N_Vector) 0 ;
  realtype arg3 ;
  ARKVecResizeFn arg4 = (ARKVecResizeFn) 0 ;
  void *arg5 = (void *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (N_Vector)(farg2);
  arg3 = (realtype)(*farg3);
  arg4 = (ARKVecResizeFn)(farg4);
  arg5 = (void *)(farg5);
  result = (int)MRIStepResize(arg1,arg2,arg3,arg4,arg5);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepReInit(void *farg1, ARKRhsFn farg2, double const *farg3, N_Vector farg4) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  ARKRhsFn arg2 = (ARKRhsFn) 0 ;
  realtype arg3 ;
  N_Vector arg4 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (ARKRhsFn)(farg2);
  arg3 = (realtype)(*farg3);
  arg4 = (N_Vector)(farg4);
  result = (int)MRIStepReInit(arg1,arg2,arg3,arg4);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepRootInit(void *farg1, int const *farg2, ARKRootFn farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  ARKRootFn arg3 = (ARKRootFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  arg3 = (ARKRootFn)(farg3);
  result = (int)MRIStepRootInit(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetDefaults(void *farg1) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  result = (int)MRIStepSetDefaults(arg1);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetInterpolantType(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)MRIStepSetInterpolantType(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetInterpolantDegree(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)MRIStepSetInterpolantDegree(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetDenseOrder(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)MRIStepSetDenseOrder(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetTable(void *farg1, int const *farg2, void *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  ARKodeButcherTable arg3 = (ARKodeButcherTable) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  arg3 = (ARKodeButcherTable)(farg3);
  result = (int)MRIStepSetTable(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetTableNum(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)MRIStepSetTableNum(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetMaxNumSteps(void *farg1, long const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long)(*farg2);
  result = (int)MRIStepSetMaxNumSteps(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetMaxHnilWarns(void *farg1, int const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int)(*farg2);
  result = (int)MRIStepSetMaxHnilWarns(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetStopTime(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  result = (int)MRIStepSetStopTime(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetFixedStep(void *farg1, double const *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  result = (int)MRIStepSetFixedStep(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetRootDirection(void *farg1, int *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int *arg2 = (int *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int *)(farg2);
  result = (int)MRIStepSetRootDirection(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetNoInactiveRootWarn(void *farg1) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  result = (int)MRIStepSetNoInactiveRootWarn(arg1);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetErrHandlerFn(void *farg1, ARKErrHandlerFn farg2, void *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  ARKErrHandlerFn arg2 = (ARKErrHandlerFn) 0 ;
  void *arg3 = (void *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (ARKErrHandlerFn)(farg2);
  arg3 = (void *)(farg3);
  result = (int)MRIStepSetErrHandlerFn(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetErrFile(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  FILE *arg2 = (FILE *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (FILE *)(farg2);
  result = (int)MRIStepSetErrFile(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetUserData(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  void *arg2 = (void *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (void *)(farg2);
  result = (int)MRIStepSetUserData(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetDiagnostics(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  FILE *arg2 = (FILE *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (FILE *)(farg2);
  result = (int)MRIStepSetDiagnostics(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetPostprocessStepFn(void *farg1, ARKPostProcessFn farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  ARKPostProcessFn arg2 = (ARKPostProcessFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (ARKPostProcessFn)(farg2);
  result = (int)MRIStepSetPostprocessStepFn(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetPostprocessStageFn(void *farg1, ARKPostProcessFn farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  ARKPostProcessFn arg2 = (ARKPostProcessFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (ARKPostProcessFn)(farg2);
  result = (int)MRIStepSetPostprocessStageFn(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetPreInnerFn(void *farg1, MRIStepPreInnerFn farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  MRIStepPreInnerFn arg2 = (MRIStepPreInnerFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (MRIStepPreInnerFn)(farg2);
  result = (int)MRIStepSetPreInnerFn(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepSetPostInnerFn(void *farg1, MRIStepPostInnerFn farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  MRIStepPostInnerFn arg2 = (MRIStepPostInnerFn) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (MRIStepPostInnerFn)(farg2);
  result = (int)MRIStepSetPostInnerFn(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepEvolve(void *farg1, double const *farg2, N_Vector farg3, double *farg4, int const *farg5) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  N_Vector arg3 = (N_Vector) 0 ;
  realtype *arg4 = (realtype *) 0 ;
  int arg5 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  arg3 = (N_Vector)(farg3);
  arg4 = (realtype *)(farg4);
  arg5 = (int)(*farg5);
  result = (int)MRIStepEvolve(arg1,arg2,arg3,arg4,arg5);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepGetDky(void *farg1, double const *farg2, int const *farg3, N_Vector farg4) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype arg2 ;
  int arg3 ;
  N_Vector arg4 = (N_Vector) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype)(*farg2);
  arg3 = (int)(*farg3);
  arg4 = (N_Vector)(farg4);
  result = (int)MRIStepGetDky(arg1,arg2,arg3,arg4);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepGetNumRhsEvals(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)MRIStepGetNumRhsEvals(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepGetCurrentButcherTables(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  ARKodeButcherTable *arg2 = (ARKodeButcherTable *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (ARKodeButcherTable *)(farg2);
  result = (int)MRIStepGetCurrentButcherTables(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepGetWorkSpace(void *farg1, long *farg2, long *farg3) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  long *arg3 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  arg3 = (long *)(farg3);
  result = (int)MRIStepGetWorkSpace(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepGetNumSteps(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)MRIStepGetNumSteps(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepGetLastStep(void *farg1, double *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype *arg2 = (realtype *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype *)(farg2);
  result = (int)MRIStepGetLastStep(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepGetCurrentTime(void *farg1, double *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  realtype *arg2 = (realtype *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (realtype *)(farg2);
  result = (int)MRIStepGetCurrentTime(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepGetCurrentState(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  N_Vector *arg2 = (N_Vector *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (N_Vector *)(farg2);
  result = (int)MRIStepGetCurrentState(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepGetNumGEvals(void *farg1, long *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  long *arg2 = (long *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (long *)(farg2);
  result = (int)MRIStepGetNumGEvals(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepGetRootInfo(void *farg1, int *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int *arg2 = (int *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int *)(farg2);
  result = (int)MRIStepGetRootInfo(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepGetLastInnerStepFlag(void *farg1, int *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  int *arg2 = (int *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (int *)(farg2);
  result = (int)MRIStepGetLastInnerStepFlag(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT SwigArrayWrapper _wrap_FMRIStepGetReturnFlagName(long const *farg1) {
  SwigArrayWrapper fresult ;
  long arg1 ;
  char *result = 0 ;
  
  arg1 = (long)(*farg1);
  result = (char *)MRIStepGetReturnFlagName(arg1);
  fresult.size = strlen((const char*)(result));
  fresult.data = (char *)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepWriteParameters(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  FILE *arg2 = (FILE *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (FILE *)(farg2);
  result = (int)MRIStepWriteParameters(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FMRIStepWriteButcher(void *farg1, void *farg2) {
  int fresult ;
  void *arg1 = (void *) 0 ;
  FILE *arg2 = (FILE *) 0 ;
  int result;
  
  arg1 = (void *)(farg1);
  arg2 = (FILE *)(farg2);
  result = (int)MRIStepWriteButcher(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT void _wrap_FMRIStepFree(void *farg1) {
  void **arg1 = (void **) 0 ;
  
  arg1 = (void **)(farg1);
  MRIStepFree(arg1);
}


SWIGEXPORT void _wrap_FMRIStepPrintMem(void *farg1, void *farg2) {
  void *arg1 = (void *) 0 ;
  FILE *arg2 = (FILE *) 0 ;
  
  arg1 = (void *)(farg1);
  arg2 = (FILE *)(farg2);
  MRIStepPrintMem(arg1,arg2);
}



