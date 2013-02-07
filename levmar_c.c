#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "levmar-2.6/levmar.h"
#include "levmar_h.h"

void callback_func(double *p, double *x, int m, int n, void *data) {
  Callback_func(p,x,data);
}

void callback_jacfunc(double *p, double *jac, int m, int n, void *data) {
  Callback_jacfunc(p,jac,data);
}

void levmar_der( double* ygiven, double* p, const int m, const int n, void* data ) {
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

  // optimization control parameters; passing to levmar NULL instead of opts reverts to defaults
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

  // invoke the optimization function
  dlevmar_der(callback_func, callback_jacfunc, p, ygiven, m, n, 1000, opts, info, NULL, NULL, data); // with analytic Jacobian
}

void levmar_dif( double* ygiven, double* p, const int m, const int n, void* data ) {
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

  // optimization control parameters; passing to levmar NULL instead of opts reverts to defaults
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

  // invoke the optimization function
  dlevmar_dif(callback_func, p, ygiven, m, n, 1000, opts, info, NULL, NULL, data); // without Jacobian
}


