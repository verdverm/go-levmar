#ifndef PGE_LEVMAR_H
#define PGE_LEVMAR_H


typedef struct {
	int *serial;
	int  s_len;
} StackExpr;

typedef struct {
	int  x_len;
	int  x_dim;
	double *x_data;  // dim 0 == time

	StackExpr  expr;
	StackExpr *derivs;
	int  d_len;
} StackData;



void func(double *p, double *x, int m, int n, void *data);
void jacfunc(double *p, double *x, int m, int n, void *data);
void levmar(double* ygiven, double* p, const int m, const int n, void* e );
void stack_func(double *p, double *x, int m, int n, void *data);
void stack_jacfunc(double *p, double *x, int m, int n, void *data);
void stack_levmar(double* ygiven, double* p, const int m, const int n, void* e );


#endif