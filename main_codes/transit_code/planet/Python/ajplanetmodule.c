/* File: ajplanetmodule.c
 * This file is auto-generated with f2py (version:2).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * See http://cens.ioc.ee/projects/f2py2e/
 * Generation date: Thu Jul 15 14:56:58 2010
 * $Revision:$
 * $Date:$
 * Do not edit this file directly unless you know what you are doing!!!
 */
#ifdef __cplusplus
extern "C" {
#endif

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "Python.h"
#include "fortranobject.h"
#include <math.h>

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *ajplanet_error;
static PyObject *ajplanet_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
/*need_typedefs*/

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif

#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
  PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
  fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif

#define rank(var) var ## _Rank
#define shape(var,dim) var ## _Dims[dim]
#define old_rank(var) (((PyArrayObject *)(capi_ ## var ## _tmp))->nd)
#define old_shape(var,dim) (((PyArrayObject *)(capi_ ## var ## _tmp))->dimensions[dim])
#define fshape(var,dim) shape(var,rank(var)-dim-1)
#define len(var) shape(var,0)
#define flen(var) fshape(var,0)
#define size(var) PyArray_SIZE((PyArrayObject *)(capi_ ## var ## _tmp))
/* #define index(i) capi_i ## i */
#define slen(var) capi_ ## var ## _len


/************************ See f2py2e/cfuncs.py: cfuncs ************************/
static int double_from_pyobj(double* v,PyObject *obj,const char *errmess) {
  PyObject* tmp = NULL;
  if (PyFloat_Check(obj)) {
#ifdef __sgi
    *v = PyFloat_AsDouble(obj);
#else
    *v = PyFloat_AS_DOUBLE(obj);
#endif
    return 1;
  }
  tmp = PyNumber_Float(obj);
  if (tmp) {
#ifdef __sgi
    *v = PyFloat_AsDouble(tmp);
#else
    *v = PyFloat_AS_DOUBLE(tmp);
#endif
    Py_DECREF(tmp);
    return 1;
  }
  if (PyComplex_Check(obj))
    tmp = PyObject_GetAttrString(obj,"real");
  else if (PyString_Check(obj))
    /*pass*/;
  else if (PySequence_Check(obj))
    tmp = PySequence_GetItem(obj,0);
  if (tmp) {
    PyErr_Clear();
    if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
    Py_DECREF(tmp);
  }
  {
    PyObject* err = PyErr_Occurred();
    if (err==NULL) err = ajplanet_error;
    PyErr_SetString(err,errmess);
  }
  return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
extern double pl_rv(double,double,double,double,double,double,double);
extern void pl_rv_array(double*,double,double,double,double,double,double,int,double*);
extern void pl_Orbit_array(double*,double,double,double,double,double,double,double,double*,double*,double*,int);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/*********************************** pl_rv ***********************************/
static char doc_f2py_rout_ajplanet_pl_rv[] = "\
Function signature:\n\
  pl_rv = pl_rv(t,v0,K,w,e,t0,P)\n\
Required arguments:\n"
"  t : input float\n"
"  v0 : input float\n"
"  K : input float\n"
"  w : input float\n"
"  e : input float\n"
"  t0 : input float\n"
"  P : input float\n"
"Return objects:\n"
"  pl_rv : float";
/* extern double pl_rv(double,double,double,double,double,double,double); */
static PyObject *f2py_rout_ajplanet_pl_rv(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           double (*f2py_func)(double,double,double,double,double,double,double)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  double pl_rv_return_value=0;
  double t = 0;
  PyObject *t_capi = Py_None;
  double v0 = 0;
  PyObject *v0_capi = Py_None;
  double K = 0;
  PyObject *K_capi = Py_None;
  double w = 0;
  PyObject *w_capi = Py_None;
  double e = 0;
  PyObject *e_capi = Py_None;
  double t0 = 0;
  PyObject *t0_capi = Py_None;
  double P = 0;
  PyObject *P_capi = Py_None;
  static char *capi_kwlist[] = {"t","v0","K","w","e","t0","P",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOOO:ajplanet.pl_rv",\
    capi_kwlist,&t_capi,&v0_capi,&K_capi,&w_capi,&e_capi,&t0_capi,&P_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable e */
    f2py_success = double_from_pyobj(&e,e_capi,"ajplanet.pl_rv() 5th argument (e) can't be converted to double");
  if (f2py_success) {
  /* Processing variable t0 */
    f2py_success = double_from_pyobj(&t0,t0_capi,"ajplanet.pl_rv() 6th argument (t0) can't be converted to double");
  if (f2py_success) {
  /* Processing variable v0 */
    f2py_success = double_from_pyobj(&v0,v0_capi,"ajplanet.pl_rv() 2nd argument (v0) can't be converted to double");
  if (f2py_success) {
  /* Processing variable P */
    f2py_success = double_from_pyobj(&P,P_capi,"ajplanet.pl_rv() 7th argument (P) can't be converted to double");
  if (f2py_success) {
  /* Processing variable t */
    f2py_success = double_from_pyobj(&t,t_capi,"ajplanet.pl_rv() 1st argument (t) can't be converted to double");
  if (f2py_success) {
  /* Processing variable w */
    f2py_success = double_from_pyobj(&w,w_capi,"ajplanet.pl_rv() 4th argument (w) can't be converted to double");
  if (f2py_success) {
  /* Processing variable K */
    f2py_success = double_from_pyobj(&K,K_capi,"ajplanet.pl_rv() 3rd argument (K) can't be converted to double");
  if (f2py_success) {
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
  pl_rv_return_value = (*f2py_func)(t,v0,K,w,e,t0,P);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("d",pl_rv_return_value);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  } /*if (f2py_success) of K*/
  /* End of cleaning variable K */
  } /*if (f2py_success) of w*/
  /* End of cleaning variable w */
  } /*if (f2py_success) of t*/
  /* End of cleaning variable t */
  } /*if (f2py_success) of P*/
  /* End of cleaning variable P */
  } /*if (f2py_success) of v0*/
  /* End of cleaning variable v0 */
  } /*if (f2py_success) of t0*/
  /* End of cleaning variable t0 */
  } /*if (f2py_success) of e*/
  /* End of cleaning variable e */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/******************************** end of pl_rv ********************************/

/******************************** pl_rv_array ********************************/
static char doc_f2py_rout_ajplanet_pl_rv_array[] = "\
Function signature:\n\
  res = pl_rv_array(t,v0,K,w,e,t0,P)\n\
Required arguments:\n"
"  t : input rank-1 array('d') with bounds (n)\n"
"  v0 : input float\n"
"  K : input float\n"
"  w : input float\n"
"  e : input float\n"
"  t0 : input float\n"
"  P : input float\n"
"Return objects:\n"
"  res : rank-1 array('d') with bounds (n)";
/* extern void pl_rv_array(double*,double,double,double,double,double,double,int,double*); */
static PyObject *f2py_rout_ajplanet_pl_rv_array(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(double*,double,double,double,double,double,double,int,double*)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  double *t = NULL;
  npy_intp t_Dims[1] = {-1};
  const int t_Rank = 1;
  PyArrayObject *capi_t_tmp = NULL;
  int capi_t_intent = 0;
  PyObject *t_capi = Py_None;
  double v0 = 0;
  PyObject *v0_capi = Py_None;
  double K = 0;
  PyObject *K_capi = Py_None;
  double w = 0;
  PyObject *w_capi = Py_None;
  double e = 0;
  PyObject *e_capi = Py_None;
  double t0 = 0;
  PyObject *t0_capi = Py_None;
  double P = 0;
  PyObject *P_capi = Py_None;
  int n = 0;
  double *res = NULL;
  npy_intp res_Dims[1] = {-1};
  const int res_Rank = 1;
  PyArrayObject *capi_res_tmp = NULL;
  int capi_res_intent = 0;
  static char *capi_kwlist[] = {"t","v0","K","w","e","t0","P",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOOO:ajplanet.pl_rv_array",\
    capi_kwlist,&t_capi,&v0_capi,&K_capi,&w_capi,&e_capi,&t0_capi,&P_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable e */
    f2py_success = double_from_pyobj(&e,e_capi,"ajplanet.pl_rv_array() 5th argument (e) can't be converted to double");
  if (f2py_success) {
  /* Processing variable K */
    f2py_success = double_from_pyobj(&K,K_capi,"ajplanet.pl_rv_array() 3rd argument (K) can't be converted to double");
  if (f2py_success) {
  /* Processing variable t0 */
    f2py_success = double_from_pyobj(&t0,t0_capi,"ajplanet.pl_rv_array() 6th argument (t0) can't be converted to double");
  if (f2py_success) {
  /* Processing variable v0 */
    f2py_success = double_from_pyobj(&v0,v0_capi,"ajplanet.pl_rv_array() 2nd argument (v0) can't be converted to double");
  if (f2py_success) {
  /* Processing variable P */
    f2py_success = double_from_pyobj(&P,P_capi,"ajplanet.pl_rv_array() 7th argument (P) can't be converted to double");
  if (f2py_success) {
  /* Processing variable t */
  ;
  capi_t_intent |= F2PY_INTENT_C|F2PY_INTENT_IN;
  capi_t_tmp = array_from_pyobj(PyArray_DOUBLE,t_Dims,t_Rank,capi_t_intent,t_capi);
  if (capi_t_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(ajplanet_error,"failed in converting 1st argument `t' of ajplanet.pl_rv_array to C/Fortran array" );
  } else {
    t = (double *)(capi_t_tmp->data);

  /* Processing variable w */
    f2py_success = double_from_pyobj(&w,w_capi,"ajplanet.pl_rv_array() 4th argument (w) can't be converted to double");
  if (f2py_success) {
  /* Processing variable n */
  n = len(t);
  /* Processing variable res */
  res_Dims[0]=n;
  capi_res_intent |= F2PY_INTENT_C|F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_res_tmp = array_from_pyobj(PyArray_DOUBLE,res_Dims,res_Rank,capi_res_intent,Py_None);
  if (capi_res_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(ajplanet_error,"failed in converting hidden `res' of ajplanet.pl_rv_array to C/Fortran array" );
  } else {
    res = (double *)(capi_res_tmp->data);

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
        (*f2py_func)(t,v0,K,w,e,t0,P,n,res);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("N",capi_res_tmp);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  }  /*if (capi_res_tmp == NULL) ... else of res*/
  /* End of cleaning variable res */
  /* End of cleaning variable n */
  } /*if (f2py_success) of w*/
  /* End of cleaning variable w */
  if((PyObject *)capi_t_tmp!=t_capi) {
    Py_XDECREF(capi_t_tmp); }
  }  /*if (capi_t_tmp == NULL) ... else of t*/
  /* End of cleaning variable t */
  } /*if (f2py_success) of P*/
  /* End of cleaning variable P */
  } /*if (f2py_success) of v0*/
  /* End of cleaning variable v0 */
  } /*if (f2py_success) of t0*/
  /* End of cleaning variable t0 */
  } /*if (f2py_success) of K*/
  /* End of cleaning variable K */
  } /*if (f2py_success) of e*/
  /* End of cleaning variable e */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/***************************** end of pl_rv_array *****************************/

/******************************* pl_Orbit_array *******************************/
static char doc_f2py_rout_ajplanet_pl_Orbit_array[] = "\
Function signature:\n\
  x,y,z = pl_Orbit_array(t,t0,P,a,e,w,inc,Omega)\n\
Required arguments:\n"
"  t : input rank-1 array('d') with bounds (n)\n"
"  t0 : input float\n"
"  P : input float\n"
"  a : input float\n"
"  e : input float\n"
"  w : input float\n"
"  inc : input float\n"
"  Omega : input float\n"
"Return objects:\n"
"  x : rank-1 array('d') with bounds (n)\n"
"  y : rank-1 array('d') with bounds (n)\n"
"  z : rank-1 array('d') with bounds (n)";
/* extern void pl_Orbit_array(double*,double,double,double,double,double,double,double,double*,double*,double*,int); */
static PyObject *f2py_rout_ajplanet_pl_Orbit_array(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(double*,double,double,double,double,double,double,double,double*,double*,double*,int)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  double *t = NULL;
  npy_intp t_Dims[1] = {-1};
  const int t_Rank = 1;
  PyArrayObject *capi_t_tmp = NULL;
  int capi_t_intent = 0;
  PyObject *t_capi = Py_None;
  double t0 = 0;
  PyObject *t0_capi = Py_None;
  double P = 0;
  PyObject *P_capi = Py_None;
  double a = 0;
  PyObject *a_capi = Py_None;
  double e = 0;
  PyObject *e_capi = Py_None;
  double w = 0;
  PyObject *w_capi = Py_None;
  double inc = 0;
  PyObject *inc_capi = Py_None;
  double Omega = 0;
  PyObject *Omega_capi = Py_None;
  double *x = NULL;
  npy_intp x_Dims[1] = {-1};
  const int x_Rank = 1;
  PyArrayObject *capi_x_tmp = NULL;
  int capi_x_intent = 0;
  double *y = NULL;
  npy_intp y_Dims[1] = {-1};
  const int y_Rank = 1;
  PyArrayObject *capi_y_tmp = NULL;
  int capi_y_intent = 0;
  double *z = NULL;
  npy_intp z_Dims[1] = {-1};
  const int z_Rank = 1;
  PyArrayObject *capi_z_tmp = NULL;
  int capi_z_intent = 0;
  int n = 0;
  static char *capi_kwlist[] = {"t","t0","P","a","e","w","inc","Omega",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOOOO:ajplanet.pl_Orbit_array",\
    capi_kwlist,&t_capi,&t0_capi,&P_capi,&a_capi,&e_capi,&w_capi,&inc_capi,&Omega_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable a */
    f2py_success = double_from_pyobj(&a,a_capi,"ajplanet.pl_Orbit_array() 4th argument (a) can't be converted to double");
  if (f2py_success) {
  /* Processing variable e */
    f2py_success = double_from_pyobj(&e,e_capi,"ajplanet.pl_Orbit_array() 5th argument (e) can't be converted to double");
  if (f2py_success) {
  /* Processing variable t0 */
    f2py_success = double_from_pyobj(&t0,t0_capi,"ajplanet.pl_Orbit_array() 2nd argument (t0) can't be converted to double");
  if (f2py_success) {
  /* Processing variable P */
    f2py_success = double_from_pyobj(&P,P_capi,"ajplanet.pl_Orbit_array() 3rd argument (P) can't be converted to double");
  if (f2py_success) {
  /* Processing variable t */
  ;
  capi_t_intent |= F2PY_INTENT_C|F2PY_INTENT_IN;
  capi_t_tmp = array_from_pyobj(PyArray_DOUBLE,t_Dims,t_Rank,capi_t_intent,t_capi);
  if (capi_t_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(ajplanet_error,"failed in converting 1st argument `t' of ajplanet.pl_Orbit_array to C/Fortran array" );
  } else {
    t = (double *)(capi_t_tmp->data);

  /* Processing variable w */
    f2py_success = double_from_pyobj(&w,w_capi,"ajplanet.pl_Orbit_array() 6th argument (w) can't be converted to double");
  if (f2py_success) {
  /* Processing variable Omega */
    f2py_success = double_from_pyobj(&Omega,Omega_capi,"ajplanet.pl_Orbit_array() 8th argument (Omega) can't be converted to double");
  if (f2py_success) {
  /* Processing variable inc */
    f2py_success = double_from_pyobj(&inc,inc_capi,"ajplanet.pl_Orbit_array() 7th argument (inc) can't be converted to double");
  if (f2py_success) {
  /* Processing variable n */
  n = len(t);
  /* Processing variable y */
  y_Dims[0]=n;
  capi_y_intent |= F2PY_INTENT_C|F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_y_tmp = array_from_pyobj(PyArray_DOUBLE,y_Dims,y_Rank,capi_y_intent,Py_None);
  if (capi_y_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(ajplanet_error,"failed in converting hidden `y' of ajplanet.pl_Orbit_array to C/Fortran array" );
  } else {
    y = (double *)(capi_y_tmp->data);

  /* Processing variable x */
  x_Dims[0]=n;
  capi_x_intent |= F2PY_INTENT_C|F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_x_tmp = array_from_pyobj(PyArray_DOUBLE,x_Dims,x_Rank,capi_x_intent,Py_None);
  if (capi_x_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(ajplanet_error,"failed in converting hidden `x' of ajplanet.pl_Orbit_array to C/Fortran array" );
  } else {
    x = (double *)(capi_x_tmp->data);

  /* Processing variable z */
  z_Dims[0]=n;
  capi_z_intent |= F2PY_INTENT_C|F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_z_tmp = array_from_pyobj(PyArray_DOUBLE,z_Dims,z_Rank,capi_z_intent,Py_None);
  if (capi_z_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(ajplanet_error,"failed in converting hidden `z' of ajplanet.pl_Orbit_array to C/Fortran array" );
  } else {
    z = (double *)(capi_z_tmp->data);

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
        (*f2py_func)(t,t0,P,a,e,w,inc,Omega,x,y,z,n);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("NNN",capi_x_tmp,capi_y_tmp,capi_z_tmp);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  }  /*if (capi_z_tmp == NULL) ... else of z*/
  /* End of cleaning variable z */
  }  /*if (capi_x_tmp == NULL) ... else of x*/
  /* End of cleaning variable x */
  }  /*if (capi_y_tmp == NULL) ... else of y*/
  /* End of cleaning variable y */
  /* End of cleaning variable n */
  } /*if (f2py_success) of inc*/
  /* End of cleaning variable inc */
  } /*if (f2py_success) of Omega*/
  /* End of cleaning variable Omega */
  } /*if (f2py_success) of w*/
  /* End of cleaning variable w */
  if((PyObject *)capi_t_tmp!=t_capi) {
    Py_XDECREF(capi_t_tmp); }
  }  /*if (capi_t_tmp == NULL) ... else of t*/
  /* End of cleaning variable t */
  } /*if (f2py_success) of P*/
  /* End of cleaning variable P */
  } /*if (f2py_success) of t0*/
  /* End of cleaning variable t0 */
  } /*if (f2py_success) of e*/
  /* End of cleaning variable e */
  } /*if (f2py_success) of a*/
  /* End of cleaning variable a */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/*************************** end of pl_Orbit_array ***************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
  {"pl_rv",-1,{{-1}},0,(char *)pl_rv,(f2py_init_func)f2py_rout_ajplanet_pl_rv,doc_f2py_rout_ajplanet_pl_rv},
  {"pl_rv_array",-1,{{-1}},0,(char *)pl_rv_array,(f2py_init_func)f2py_rout_ajplanet_pl_rv_array,doc_f2py_rout_ajplanet_pl_rv_array},
  {"pl_Orbit_array",-1,{{-1}},0,(char *)pl_Orbit_array,(f2py_init_func)f2py_rout_ajplanet_pl_Orbit_array,doc_f2py_rout_ajplanet_pl_Orbit_array},

/*eof routine_defs*/
  {NULL}
};

static PyMethodDef f2py_module_methods[] = {

  {NULL,NULL}
};

PyMODINIT_FUNC initajplanet(void) {
  int i;
  PyObject *m,*d, *s;
  m = ajplanet_module = Py_InitModule("ajplanet", f2py_module_methods);
  PyFortran_Type.ob_type = &PyType_Type;
  import_array();
  if (PyErr_Occurred())
    {PyErr_SetString(PyExc_ImportError, "can't initialize module ajplanet (failed to import numpy)"); return;}
  d = PyModule_GetDict(m);
  s = PyString_FromString("$Revision: $");
  PyDict_SetItemString(d, "__version__", s);
  s = PyString_FromString("This module 'ajplanet' is auto-generated with f2py (version:2).\nFunctions:\n"
"  pl_rv = pl_rv(t,v0,K,w,e,t0,P)\n"
"  res = pl_rv_array(t,v0,K,w,e,t0,P)\n"
"  x,y,z = pl_Orbit_array(t,t0,P,a,e,w,inc,Omega)\n"
".");
  PyDict_SetItemString(d, "__doc__", s);
  ajplanet_error = PyErr_NewException ("ajplanet.error", NULL, NULL);
  Py_DECREF(s);
  for(i=0;f2py_routine_defs[i].name!=NULL;i++)
    PyDict_SetItemString(d, f2py_routine_defs[i].name,PyFortranObject_NewAsAttr(&f2py_routine_defs[i]));



/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
  if (! PyErr_Occurred())
    on_exit(f2py_report_on_exit,(void*)"ajplanet");
#endif

}
#ifdef __cplusplus
}
#endif
