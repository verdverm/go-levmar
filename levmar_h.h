#ifndef PGE_LEVMAR_H
#define PGE_LEVMAR_H


void levmar_der(double* ygiven, double* p, int m, int n, void* data );
void levmar_dif(double* ygiven, double* p, int m, int n, void* data );


#endif