#include <stdlib.h>

#include "stack.h"

I_Stack* new_istack() { 
	I_Stack* stack = (I_Stack*) malloc(sizeof(I_Stack));
	stack->pos = -1;
	return stack;
}
int push_istack(I_Stack* stack, int i) { 
	if (stack->pos +1 >= STACK_BUF_LEN)
		return 0;
	stack->pos++;
	stack->data[stack->pos] = i;
	return 1;
}
void pop_istack(I_Stack* stack) { 
	if (stack->pos >= 0) stack->pos--;
}
int  top_istack(I_Stack* stack) { 
	return stack->data[stack->pos];
}
int  len_istack(I_Stack* stack) { 
	return stack->pos + 1;
}
int  get_istack(I_Stack* stack, int pos) { 
	return stack->data[pos];
}
int  is_empty_istack(I_Stack* stack) { 
	return (stack->pos < 0) ? 1 : 0;
}
void clear_istack(I_Stack* stack) { 
	stack->pos = -1;
}
void free_istack(I_Stack* stack) { 
	free(stack);
}

D_Stack* new_dstack() { 
	D_Stack* stack = (D_Stack*) malloc(sizeof(D_Stack));
	stack->pos = -1;
	return stack;
}
int push_dstack(D_Stack* stack, double f) { 
	if (stack->pos +1 >= STACK_BUF_LEN)
		return 0;
	stack->pos++;
	stack->data[stack->pos] = f;
	return 1;
}
void pop_dstack(D_Stack* stack) { 
	if (stack->pos >= 0) stack->pos--;
}
double  top_dstack(D_Stack* stack) { 
	return stack->data[stack->pos];
}
int  len_dstack(D_Stack* stack) { 
	return stack->pos + 1;
}
double  get_dstack(D_Stack* stack, int pos) { 
	return stack->data[pos];
}
int  is_empty_dstack(D_Stack* stack) { 
	return (stack->pos < 0) ? 1 : 0;
}
void clear_dstack(D_Stack* stack) { 
	stack->pos = -1;
}
void free_dstack(D_Stack* stack) { 
	free(stack);
}

