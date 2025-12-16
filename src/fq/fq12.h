/*
* Fq12 = Fq6[w]/(w^2 - v).
*
*/

#ifndef BLS12_381_fq12_H
#define BLS12_381_fq12_H

#include "fq.h"
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
void Fq12_set(Fq12* a, Fq6* num0, Fq6* num1);
void Fq12_set_zero(Fq12* rop); // set rop to zero.
void Fq12_set_id(Fq12* rop); // set rop to 1.
void Fq12_assign(Fq12* rop, Fq12* a);

// conversion
void Fq6_to_Fq12(Fq12* rop, Fq6* a);
void Fq2_to_Fq12(Fq12* rop, Fq2* a);

// Montgomery representation
void Fq12_mont_rep(Fq12* aR, Fq12* a);
void Fq12_mont_rep_inv(Fq12* a, Fq12* aR);

//Conjugation
void Fq12_conj(Fq12* rop, Fq12* a);

// comparison
uint32_t Fq12_cmp(Fq12* a, Fq12* b); // if a=b, return 1; otherwise, return 0

// Arithmetic
static inline void Fq12_add(Fq12* rop, Fq12* a, Fq12* b)
{
	Fq6_add(&(rop->val0), &(a->val0), &(b->val0));
	Fq6_add(&(rop->val1), &(a->val1), &(b->val1));
}

static inline void Fq12_sub(Fq12* rop, Fq12* a, Fq12* b)
{
	Fq6_sub(&(rop->val0), &(a->val0), &(b->val0));
	Fq6_sub(&(rop->val1), &(a->val1), &(b->val1));
}

static inline void Fq12_mont_mul(Fq12* rop, Fq12* a, Fq12* b)
{
	/*
	
	a*b = a0b0 + (a0b1 + a1b0)w + a1b1 w^2 = [w^2 = v] =

		= (a0b0 + a1b1 * v) + (a0b1 + a1b0)w.
	
	*/
	
	// auxiliary init
	Fq6 temp00, temp01, temp10, temp11, temp11v, c0, c1;

	// auxiliary computation
	Fq6_mont_mul(&temp00, &(a->val0), &(b->val0));
	Fq6_mont_mul(&temp01, &(a->val0), &(b->val1));
	Fq6_mont_mul(&temp10, &(a->val1), &(b->val0));
	Fq6_mont_mul(&temp11, &(a->val1), &(b->val1));
	Fq6_transform(&temp11v, &temp11);

	// computation
	Fq6_add(&c0, &temp00, &temp11v); // a0b0 + a1b1 * v
	Fq6_add(&c1, &temp01, &temp10);

	// storing
	Fq6_assign(&(rop->val0), &c0);
	Fq6_assign(&(rop->val1), &c1);
}

static inline void Fq12_add_inv(Fq12* rop, Fq12* a)
{
	Fq6_add_inv(&(rop->val0), &(a->val0));
	Fq6_add_inv(&(rop->val1), &(a->val1));
}
static inline void Fq12_mont_mul_inv(Fq12* rop, Fq12* a)
{
	/*
	
	a^(-1) = (a0^2 - a1^2 * v)^(-1) * a0 + (a0^2 - a1^2 * v)^(-1) * (-a1) * w.	
	
	*/
	
	// auxiliary init
	Fq6 temp00, temp11, temp11v, norm, norminv, neg1, c0, c1;

	// auxiliary computation
	Fq6_mont_mul(&temp00, &(a->val0), &(a->val0));
	Fq6_mont_mul(&temp11, &(a->val1), &(a->val1));
	Fq6_transform(&temp11v, &temp11);

	Fq6_sub(&norm, &temp00, &temp11v);
	Fq6_mont_mul_inv(&norminv, &norm);

	Fq6_add_inv(&neg1, &(a->val1));

	// computation
	Fq6_mont_mul(&c0, &norminv, &(a->val0));
	Fq6_mont_mul(&c1, &norminv, &neg1);

	// storing
	Fq6_assign(&(rop->val0), &c0);
	Fq6_assign(&(rop->val1), &c1);
}

// print
void Fq12_print(Fq12* rop);
void Fq12_print_dec(Fq12* rop);

#endif
