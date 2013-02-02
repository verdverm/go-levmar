#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../levmar-2.6/levmar.h"
#include "levmar_h.h"
#include "stack.h"

double low = -2;
double hgh =  2;
double step = 0.1;
int len = -1;



int sF1[] = {2,0,2,1,25,23,22,2,2,25,25,23,23,22};
int sF1_len = 14;
int sF1_c0[] = {3,1};  // 1 ~ x/x
int sF1_c0_len = 2;
double F1 (double in) {
	return 1.0 + 2.0*in + 3.0*in*in;
}

// int sF1[] = {2,0,25,23,25,23,25,23,2,1,25,23,22,2,2,25,25,23,23,22};
// int sF1_len = 20;
// int sF1_c0[] = {25,25,23,25,23};  // x*x*x
// int sF1_c0_len = 5;
// double F1 (double in) {
// 	return 1.0*in*in*in + 2.0*in + 3.0*in*in;
// }

int sF1_c1[] = {25};		// x
int sF1_c1_len = 1;
int sF1_c2[] = {25,25,23};  // x*x
int sF1_c2_len = 3;


double pF1( double in, double *p) {
	return p[0]*in*in*in + p[1]*in + p[2]*in*in;
}

double dF1_c0(double in, double *p) {
	return in*in*in;
}

double dF1_c1(double in, double *p) {
	return in;
}

double dF1_c2(double in, double *p) {
	return in*in;
}

void levmar_func(double *p, double *x, int m, int n, void *data) {
	double* input = (double*)data;

	int i;
	for( i=0; i<n; i++ ) {
		x[i] = pF1(input[i],p);
	}
}

void levmar_jacfunc(double *p, double *jac, int m, int n, void *data) {
	double* input = (double*)data;
	int i;
	for( i=0; i<n; i++ ) {
		jac[i*m] = dF1_c0(input[i],p);
		jac[i*m+1] = dF1_c1(input[i],p);
		jac[i*m+2] = dF1_c2(input[i],p);
	}
}

double test_stack_eval(double t, double *c_in, double *x_in, D_Stack *d_stack, StackExpr expr );
void test_stack_func( double *p, double *x, int m, int n, void *data);
void test_stack_jacfunc( double *p, double *jac, int m, int n, void *data);


void test_levmar_der( double* ygiven, double* p, const int m, const int n, void* data ) {
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

  // optimization control parameters; passing to levmar NULL instead of opts reverts to defaults
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

  // invoke the optimization function
  // dlevmar_der(levmar_func, levmar_jacfunc, p, ygiven, m, n, 1000, opts, info, NULL, NULL, data); // with analytic Jacobian
  dlevmar_der(test_stack_func, test_stack_jacfunc, p, ygiven, m, n, 1000, opts, info, NULL, NULL, data); // with analytic Jacobian
}
void test_levmar_dif( double* ygiven, double* p, const int m, const int n, void* data ) {
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

  // optimization control parameters; passing to levmar NULL instead of opts reverts to defaults
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

  // invoke the optimization function
  // dlevmar_dif(levmar_func, p, ygiven, m, n, 1000, opts, info, NULL, NULL, data); // without Jacobian
  dlevmar_dif(test_stack_func, p, ygiven, m, n, 1000, opts, info, NULL, NULL, data); // without Jacobian
}


double* makeInput(int l, int h, double s) {
	double rng = h-l+s;
	len = (rng/s)+1;

	double* input = (double*) malloc(sizeof(double)*len);
	int i = 0;
	double curr = l;
	while( curr < h ) {
		input[i] = curr;
		curr += s;
		i++;
	}
	return input;
}
double* makeOutput( double* input, int len ) {
	double* out = (double*)malloc(sizeof(double)*len);

	int i;
	for( i=0; i<len; i++ )
		out[i] = F1(input[i]);
	return out;
}


int main(int argc, char const *argv[])
{
	double* input = makeInput(low,hgh,step);
	double* output = makeOutput(input,len);
	
	int i;
	for( i=0; i<len; i++) {
		printf( "%d %f -> %f\n", i,input[i],output[i]); 
	}


	StackData sdata;
	sdata.x_len = len;
	sdata.x_dim = 1;
	sdata.x_data = input;
	sdata.expr.serial = sF1;
	sdata.expr.s_len = sF1_len;

	sdata.d_len = 3;
	sdata.derivs = (StackExpr*)malloc(sizeof(StackExpr)*sdata.d_len);
	sdata.derivs[0].serial = sF1_c0;
	sdata.derivs[0].s_len = sF1_c0_len;
	sdata.derivs[1].serial = sF1_c1;
	sdata.derivs[1].s_len = sF1_c1_len;
	sdata.derivs[2].serial = sF1_c2;
	sdata.derivs[2].s_len = sF1_c2_len;


	D_Stack *d_stack = new_dstack();


	double *cs = (double*) malloc(sizeof(double)*3);
	cs[0] = 1;
	cs[1] = 2;
	cs[2] = 3;

	// for( i=0; i<1; i++) {
	// 	double out = pF1(input[i],cs);
	// 	// double out = dF1_c2(input[i],cs);
	// 	printf( "basic: %d %f -> %f  ~ %f\n", i,input[i],out,output[i]); 
	// 	double out2 = test_stack_eval(0,cs,&input[i],d_stack,sdata.expr);
	// 	// double out2 = test_stack_eval(0,cs,&input[i],d_stack,sdata.derivs[2]);
	// 	printf( "stack: %d %f -> %f  ~ %f\n\n", i,input[i],out2,output[i]); 
	// }
	


	// return 0;


	cs[0] = 1;
	cs[1] = 1;
	cs[2] = 1;
	test_levmar_der(output,cs,3,len,&sdata);
	// test_levmar_der(output,cs,3,len,input);
	printf( "deriv = %f %f %f\n", cs[0], cs[1], cs[2]);

	cs[0] = 1;
	cs[1] = 1;
	cs[2] = 1;
	test_levmar_dif(output,cs,3,len,&sdata);
	// test_levmar_dif(output,cs,3,len,input);
	printf( "diff  = %f %f %f\n", cs[0], cs[1], cs[2]);

	return 0;
}



void test_stack_func( double *p, double *x, int m, int n, void *data) {
	StackData *sdata = (StackData*)data;

	int x_len = sdata->x_len;
	int x_dim = sdata->x_dim;
	double* x_data = sdata->x_data;
	D_Stack d_stack;

	int i;
	for( i=0; i < n; i++ ) {
		x[i] = test_stack_eval(0,p,&(x_data[i*x_dim]),&d_stack,sdata->expr);
	}

}

void test_stack_jacfunc( double *p, double *jac, int m, int n, void *data) {
	StackData *sdata = (StackData*)data;

	int x_len = sdata->x_len;
	int x_dim = sdata->x_dim;
	double* x_data = sdata->x_data;
	D_Stack d_stack;

	int i,j;
	for( i=0; i < n; i++ ) {
		for( j=0; j < m; j++ ) {
			jac[i*m+j] = test_stack_eval(0,p,&(x_data[i*x_dim]),&d_stack,sdata->derivs[j]);
		}
	}

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


// #define print_stack_eval 1

double test_stack_eval(double t, double *c_in, double *x_in, D_Stack *d_stack, StackExpr expr ) {

	int i,s;

	I_Stack *serial = new_istack();
	#ifdef print_stack_eval
	printf("Serial: ");
	#endif
	for( s=0; s < expr.s_len; s++ ) {
	#ifdef print_stack_eval
		printf( "%d ", expr.serial[s]);
	#endif
		push_istack(serial,expr.serial[expr.s_len-s-1]);
	}
	#ifdef print_stack_eval
	printf("\n");
	#endif
	
	clear_dstack(d_stack);

	// fill i_stack with cmds and d_stack with leaves
	while ( !is_empty_istack(serial) ) {
		// printf( "processing serial\n");
		int val = top_istack(serial);
		pop_istack(serial);

	#ifdef print_stack_eval
		int dlen = len_dstack(d_stack);
		int slen = len_istack(serial);
		printf( "S: %d    val:  %d   \n", slen, val );
		printf( "serial(%d): [ ", slen); 
		for( i=0; i < slen; i++ )
			printf( "%d ", get_istack(serial,i) );
		printf(" ]\n");
		printf( "d_stack(%d): [ ", dlen); 
		for( i=0; i < dlen; i++ )
			printf( "%.2f ", get_dstack(d_stack,i) );
		printf(" ]\n");
	#endif

		switch (val) {
			
			// CONSTANT:  2
			case 2: {
				int p = top_istack(serial);
				pop_istack(serial);
				push_dstack(d_stack,c_in[p]);
			}
				break; 

			// HACK***
			// CONSTANT:  3
			case 3: {
				int p = top_istack(serial);
				pop_istack(serial);
				push_dstack(d_stack,p);
			}
				break; 
			// TIME:      4
			case 4:
				push_dstack(d_stack,t);
				break;
			// SYSTEM:    5
			// case 5:
				// s++;
				// push_dstack(d_stack,sys_in[serial[s]]);
			// VAR:       6   should already be transformed, but just in case
			case 6: {
				int p = top_istack(serial);
				pop_istack(serial);
				push_dstack(d_stack,x_in[p]);
				break;
			}
			// NEG:       9
			case 9: {
				double top = top_dstack(d_stack);
				pop_dstack(d_stack);
				push_dstack(d_stack, -top);
			}
				break;
			// ABS:       10
			case 10: {
				double top = top_dstack(d_stack);
				pop_dstack(d_stack);
				push_dstack(d_stack, fabs(top));
			}
				break;
			// SQRT:      11
			case 11: {
				double top = top_dstack(d_stack);
				pop_dstack(d_stack);
				push_dstack(d_stack, sqrt(top));
			}
				break;
			// SIN:       12
			case 12: {
				double top = top_dstack(d_stack);
				pop_dstack(d_stack);
				push_dstack(d_stack, sin(top));
			}
				break;
			// COS:       13
			case 13: {
				double top = top_dstack(d_stack);
				pop_dstack(d_stack);
				push_dstack(d_stack, cos(top));
			}
				break;
			// TAN:       14
			case 14: {
				double top = top_dstack(d_stack);
				pop_dstack(d_stack);
				push_dstack(d_stack, tan(top));
			}
				break;
			// EXP:       15
			case 15: {
				double top = top_dstack(d_stack);
				pop_dstack(d_stack);
				push_dstack(d_stack, exp(top));
			}
				break;
			// LOG:       16
			case 16: {
				double top = top_dstack(d_stack);
				pop_dstack(d_stack);
				push_dstack(d_stack, log(top));
			}
				break;

						// POWI:      18
			case 18: {
				double top = top_dstack(d_stack);
				pop_dstack(d_stack);
				int pwr = top_istack(serial);
				pop_istack(serial);
				push_dstack(d_stack, pow(top,pwr));
			}
				break;
			// POWE:      20
			case 20: {
					double top = top_dstack(d_stack);
					pop_dstack(d_stack);
					double pwr = top_dstack(d_stack);
					pop_dstack(d_stack);
					push_dstack(d_stack, pow(top,pwr));
				}
				break;
			// DIV:       21
			case 21: {
					double denom = top_dstack(d_stack);
					pop_dstack(d_stack);
					double numer = top_dstack(d_stack);
					pop_dstack(d_stack);
					push_dstack(d_stack, numer / denom);
				}
				break;
			// ADD:       22
			case 22: {
					double lhs = top_dstack(d_stack);
					pop_dstack(d_stack);
					double rhs = top_dstack(d_stack);
					pop_dstack(d_stack);
					push_dstack(d_stack, lhs + rhs);
					// printf( "%f + %f = %f\n", lhs, rhs, top_dstack(d_stack));
				}
				break;
			// MUL:       23
			case 23: {
					double lhs = top_dstack(d_stack);
					pop_dstack(d_stack);
					double rhs = top_dstack(d_stack);
					pop_dstack(d_stack);
					push_dstack(d_stack, lhs * rhs);
					// printf( "%f * %f = %f\n", lhs, rhs, top_dstack(d_stack));
				}
				break;




			case 0:
			default:
				// STARTVAR:  25
				if (val >= 25) {
					push_dstack(d_stack,x_in[val-25]); // x_dim_of_var = val - STARTVAR
				} else {
					printf("pushing unknown cmd %d\n", val);
					abort();
				}
		}

	#ifdef print_stack_eval
		dlen = len_dstack(d_stack);
		printf( "S: %d val: %d\n", s, val);
		printf( "d_stack(%d): [ ", dlen); 
		for( i=0; i < dlen; i++ )
			printf( "%.2f ", get_dstack(d_stack,i) );
		printf(" ]\n\n");
	#endif
	}

	return top_dstack(d_stack);

}


