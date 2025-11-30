/*
* Fq6 = Fq2[v]/(v^3 - (u + 1)).
*
*/

#ifndef BLS12_381_fq6_H
#define BLS12_381_fq6_H

#include "fq/bls_params.h"
#include "fq/fq.h"
#include "fq2.h"
#include <gmp.h>
#include <stdio.h>

typedef struct
{

	Fq2 val0;
	Fq2 val1;
	Fq2 val2;

} Fq6;

extern Fq6 v1inv, v2inv;

// init
void Fq6_init(Fq6* a);
void Fq6_clear(Fq6* a);
void Fq6_set(Fq6* a, Fq2* num0, Fq2* num1, Fq2* num2);
void Fq6_assign(Fq6* rop, Fq6* a);
void Fq6_zero(Fq6* rop); // set rop to zero.
void Fq6_id(Fq6* rop); // set rop to 1.

// conversion
void Fq2_to_Fq6(Fq6* rop, Fq2* a);

// comparison
int Fq6_cmp(Fq6* a, Fq6* b); // NB! return 0, if a=b; return 1, otherwise.

// Arithmetic
void Fq6_add(Fq6* rop, Fq6* a, Fq6* b);
void Fq6_sub(Fq6* rop, Fq6* a, Fq6* b);
void Fq6_mul(Fq6* rop, Fq6* a, Fq6* b);
void Fq6_add_inv(Fq6* rop, Fq6* a);
void Fq6_mul_inv(Fq6* rop, Fq6* a);

// Determinant
void Fq2_determinant(Fq2* rop, Fq6* a); // a helper function to evalute the multiplicative inverse.

// Transformation
void Fq6_transform(Fq6* rop, Fq6* a); // rop <- a*v

// print
void Fq6_print(Fq6* rop);

#endif
