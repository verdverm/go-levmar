#ifndef PGE_STACK_H
#define PGE_STACK_H


#define STACK_BUF_LEN 128

typedef struct {
	double data[STACK_BUF_LEN];
	int pos;
} D_Stack;

typedef struct {
	int data[STACK_BUF_LEN];
	int pos;
} I_Stack;

typedef struct {
	double *data;
	int pos;
} D_StackD;

typedef struct {
	int *data;
	int pos;
} I_StackD;

I_Stack* new_istack();
int push_istack(I_Stack* stack, int i);
void pop_istack(I_Stack* stack);
int  top_istack(I_Stack* stack);
int  len_istack(I_Stack* stack);
int  get_istack(I_Stack* stack, int pos);
int  is_empty_istack(I_Stack* stack);
void clear_istack(I_Stack* stack);
void free_istack(I_Stack* stack);

D_Stack* new_dstack();
int push_dstack(D_Stack* stack, double f);
void pop_dstack(D_Stack* stack);
double  top_dstack(D_Stack* stack);
int  len_dstack(D_Stack* stack);
double  get_dstack(D_Stack* stack, int pos);
int  is_empty_dstack(D_Stack* stack);
void clear_dstack(D_Stack* stack);
void free_dstack(D_Stack* stack);




#endif