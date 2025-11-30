/*
* Fq12 = Fq6[w]/(w^2 - v).
*
*/

#ifndef BLS12_381_fq12_H
#define BLS12_381_fq12_H

#include "fq/bls_params.h"
#include "fq/fq.h"
#include "fq2.h"
#include "fq6.h"
#include <gmp.h>
#include <stdio.h>

typedef struct
{

	Fq6 val0;
	Fq6 val1;

} Fq12;

// init
void Fq12_init(Fq12* a);
void Fq12_clear(Fq12* a);
void Fq12_set(Fq12* a, Fq6* num0, Fq6* num1);
void Fq12_assign(Fq12* rop, Fq12* a);
void Fq12_zero(Fq12* rop); // set rop to zero.
void Fq12_id(Fq12* rop); // set rop to 1.

// conversion
void Fq6_to_Fq12(Fq12* rop, Fq6* a);
void Fq2_to_Fq12(Fq12* rop, Fq2* a);

// comparison
int Fq12_cmp(Fq12* a, Fq12* b); // NB! return 0, if a=b; return 1, otherwise

// Arithmetic
void Fq12_add(Fq12* rop, Fq12* a, Fq12* b);
void Fq12_sub(Fq12* rop, Fq12* a, Fq12* b);
void Fq12_mul(Fq12* rop, Fq12* a, Fq12* b);
void Fq12_add_inv(Fq12* rop, Fq12* a);
void Fq12_mul_inv(Fq12* rop, Fq12* a);

// print
void Fq12_print(Fq12* rop);

#endif
