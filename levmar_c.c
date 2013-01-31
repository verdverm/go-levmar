#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "levmar-2.6/levmar.h"
#include "levmar_h.h"

typedef struct {
	double* data;
	int pos;
} D_Stack;

typedef struct {
	int* data;
	int pos;
} I_Stack;

double stack_eval(double *c_in, double *x_in, I_Stack *i_stack, D_Stack *d_stack, StackExpr expr );

void func(double *p, double *x, int m, int n, void *data) {
  Callback_func(p,x,data);
}

void jacfunc(double *p, double *jac, int m, int n, void *data) {
  Callback_jacfunc(p,jac,data);
}

void levmar( double* ygiven, double* p, const int m, const int n, void* data ) {
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

  // optimization control parameters; passing to levmar NULL instead of opts reverts to defaults
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

  // invoke the optimization function
  dlevmar_der(func, jacfunc, p, ygiven, m, n, 1000, opts, info, NULL, NULL, data); // with analytic Jacobian
  //   dlevmar_dif(f1, p, x, m, n, 1000, opts, info, NULL, NULL, data); // without Jacobian
}


void stack_func( double *p, double *x, int m, int n, void *data) {
	StackData *sdata = (StackData*)data;

	int x_len = sdata->x_len;
	int x_dim = sdata->x_dim;
	double* x_data = sdata->x_data;
	D_Stack d_stack;
	d_stack.data = (double*) malloc(sizeof(double)*64);
	d_stack.pos = 0;
	I_Stack i_stack;
	i_stack.data = (int*) malloc(sizeof(int)*64);
	i_stack.pos = 0;

	int i;
	for( i=0; i <  sdata->x_len; i++ ) {
		x[i] = stack_eval(p,&(x_data[i*x_dim]),&i_stack,&d_stack,sdata->expr);
	}

}

void stack_jacfunc( double *p, double *jac, int m, int n, void *data) {
	StackData *sdata = (StackData*)data;

	int x_len = sdata->x_len;
	int x_dim = sdata->x_dim;
	double* x_data = sdata->x_data;
	D_Stack d_stack;
	d_stack.data = (double*) malloc(sizeof(double)*64);
	d_stack.pos = 0;
	I_Stack i_stack;
	i_stack.data = (int*) malloc(sizeof(int)*64);
	i_stack.pos = 0;

	int i,j;
	for( i=0; i < sdata->x_len; i++ ) {
		for( j=0; j < sdata->d_len; j++ ) {
			jac[i*x_len+j] = stack_eval(p,&(x_data[i*x_dim]),&i_stack,&d_stack,sdata->derivs[j]);
		}
	}

}

void stack_levmar( double* ygiven, double* p, const int m, const int n, void* data ) {
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

  // optimization control parameters; passing to levmar NULL instead of opts reverts to defaults
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

  // invoke the optimization function
  dlevmar_der(stack_func, stack_jacfunc, p, ygiven, m, n, 1000, opts, info, NULL, NULL, data); // with analytic Jacobian
  //   dlevmar_dif(f1, p, x, m, n, 1000, opts, info, NULL, NULL, data); // without Jacobian
}


/*
ExprTypes:
---------------
NULL:      0
STARTLEAF: 1
CONSTANT:  2
TIME:      4
SYSTEM:    5
VAR:       6
LASTLEAF:  7
STARTFUNC: 8
NEG:       9
ABS:       10
SQRT:      11
SIN:       12
COS:       13
TAN:       14
EXP:       15
LOG:       16
LASTFUNC:  17
POWI:      18
POWF:      19
POWE:      20
DIV:       21
ADD:       22
MUL:       23
EXPR_MAX:  24
STARTVAR:  25
*/


double stack_eval(double *c_in, double *x_in, I_Stack *i_stack, D_Stack *d_stack, StackExpr expr ) {

	int s_len = expr.s_len;
	int *serial = expr.serial;

	int s,i,d;
	// fill i_stack with cmds and d_stack with leaves
	for( s=0; s < s_len; s++ ) {
		int val = serial[s];
		switch (val) {
			
			// CONSTANT:  2
			case 2:
				s++;
				int c = serial[s];
				d_stack->data[d_stack->pos] = c_in[c];
				d_stack->pos++;
				break; 
			// TIME:      4
			case 4:
				d_stack->data[d_stack->pos] = x_in[0];;
				d_stack->pos++;
				break;
			// SYSTEM:    5
			// VAR:       6

			// NEG:       9
			// ABS:       10
			// SQRT:      11
			// SIN:       12
			// COS:       13
			// TAN:       14
			// EXP:       15
			// LOG:       16
			// POWF:      19
			// POWE:      20
			// DIV:       21
			case 9:
			case 10:
			case 11:
			case 12:
			case 13:
			case 14:
			case 15:
			case 16:
			case 17:
			case 19:
			case 20:
			case 21:
				i_stack->data[i_stack->pos] = val;
				i_stack->pos++;
				break;

			// POWI:      18
			case 18:
				s++;
				i_stack->data[i_stack->pos] = serial[s];
				i_stack->pos++;
				i_stack->data[i_stack->pos] = val;
				i_stack->pos++;
				break;

			// ADD:       22
			case 22:
				s++;
				i_stack->data[i_stack->pos] = serial[s];
				i_stack->pos++;
				i_stack->data[i_stack->pos] = val;
				i_stack->pos++;
				break;

			// MUL:       23
			case 23:
				s++;
				i_stack->data[i_stack->pos] = serial[s];
				i_stack->pos++;
				i_stack->data[i_stack->pos] = val;
				i_stack->pos++;
				break;

			// STARTVAR:  25

			case 0:
			default:
				if (serial[s] > 25) {
					d_stack->data[d_stack->pos] = x_in[val-25];
					d_stack->pos++;
				} else {
					abort();
				}
		}
	}

	// unwind i_stack and d_stack; d_stack->data[0] should have result after this loop
	while( i_stack->pos > 0 ) {
		int val = i_stack->data[i_stack->pos];
		switch (val) {
			
			// SYSTEM:    5
			// VAR:       6

			// NEG:       9
			case 9:
				d_stack->data[d_stack->pos] *= -1.0;
				break;
			// ABS:       10
			case 10:
				d_stack->data[d_stack->pos] = fabs(d_stack->data[d_stack->pos]);
				break;
			// SQRT:      11
			case 11:
				d_stack->data[d_stack->pos] = sqrt(d_stack->data[d_stack->pos]);
				break;
			// SIN:       12
			case 12:
				d_stack->data[d_stack->pos] = sin(d_stack->data[d_stack->pos]);
				break;
			// COS:       13
			case 13:
				d_stack->data[d_stack->pos] = cos(d_stack->data[d_stack->pos]);
				break;
			// TAN:       14
			case 14:
				d_stack->data[d_stack->pos] = tan(d_stack->data[d_stack->pos]);
				break;
			// EXP:       15
			case 15:
				d_stack->data[d_stack->pos] = exp(d_stack->data[d_stack->pos]);
				break;
			// LOG:       16
			case 16:
				d_stack->data[d_stack->pos] = log(d_stack->data[d_stack->pos]);
				break;

			// POWI:      18
			case 18:
				i_stack->pos--;
				d_stack->data[d_stack->pos] = pow(d_stack->data[d_stack->pos],i_stack->data[i_stack->pos]);
				break;
			// POWE:      20
			case 20: {
					double top = d_stack->data[d_stack->pos];
					d_stack->pos--;
					d_stack->data[d_stack->pos] = pow(d_stack->data[d_stack->pos],top);
				}
				break;
			// DIV:       21
			case 21: {
					double top = d_stack->data[d_stack->pos];
					d_stack->pos--;
					d_stack->data[d_stack->pos] = pow(d_stack->data[d_stack->pos],top);
				}
				break;
			// ADD:       22
			case 22: {
					i_stack->pos--;
					int nc = i_stack->data[i_stack->pos];
					double sum = 0;
					while( nc > 0 ) {
						sum += d_stack->data[d_stack->pos];
						d_stack->pos--;
					}
					d_stack->data[d_stack->pos] = sum
				}
				break;
			// MUL:       23
			case 23: {
					i_stack->pos--;
					int nc = i_stack->data[i_stack->pos];
					double sum = 1;
					while( nc > 0 ) {
						sum *= d_stack->data[d_stack->pos];
						d_stack->pos--;
					}
					d_stack->data[d_stack->pos] = sum
				}
				break;

			// STARTVAR:  25

			default:
				abort();
		}
	}

	return d_stack->data[0];

}
