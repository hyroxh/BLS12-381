/*
* Fq2 = Fq[u]/(u^2 + 1).
*
*/

#ifndef BLS12_381_fq2_H
#define BLS12_381_fq2_H

#include "fq/bls_params.h"
#include "fq/fq.h"
#include <gmp.h>

typedef struct
{

	Fq val0;
	Fq val1;

} Fq2;

// init
void Fq2_init(Fq2* a);
void Fq2_clear(Fq2* a);
void Fq2_set(Fq2* a, Fq* num0, Fq* num1);
void Fq2_assign(Fq2* rop, Fq2* a);
void Fq2_set_str(Fq2* a, const char* str0, const char* str1);
void Fq2_set_ui(Fq2* a, unsigned long ui0, unsigned long ui1);

// comparison
int Fq2_cmp(Fq2* a, Fq2* b); // NB! return 0, if a==b; returns 1, otherwise.
							 // Return values are in accordance with GNU GMP

// Arithmetic
void Fq2_add(Fq2* rop, Fq2* a, Fq2* b);
void Fq2_sub(Fq2* rop, Fq2* a, Fq2* b);
void Fq2_mul(Fq2* rop, Fq2* a, Fq2* b);
void Fq2_add_inv(Fq2* rop, Fq2* a);
void Fq2_mul_inv(Fq2* rop, Fq2* a);

void Fq2_transform(Fq2* rop, Fq2* a); // rop<-a*u

// print
void Fq2_print(Fq2* rop);

#endif
