/*
* Fq2 = Fq[u]/(u^2 + 1).
*
*/

#ifndef BLS12_381_fq2_H
#define BLS12_381_fq2_H

#include "fq.h"
#include <gmp.h>

typedef struct
{

	Fq val0;
	Fq val1;

} Fq2;

// init
void Fq2_set(Fq2* rop, Fq* num0, Fq* num1);
void Fq2_set_ui(Fq2* rop, uint32_t num0, uint32_t num1);
void Fq2_set_zero(Fq2* a);
void Fq2_set_id(Fq2* a);
void Fq2_assign(Fq2* rop, Fq2* a);
void Fq2_set_hex_str(Fq2* rop, char* hex_str0, char* hex_str1);

// comparison
uint32_t Fq2_cmp(Fq2* a, Fq2* b); // NB! return 0, if a==b; returns 1, otherwise.
							 // Return values are in accordance with GNU GMP

// Montgomery representation
void Fq2_mont_rep(Fq2* aR, Fq2* a);
void Fq2_mont_rep_inv(Fq2* a, Fq2* aR);

// Arithmetic
static inline void Fq2_add(Fq2* rop, Fq2* a, Fq2* b)
{
	Fq_add(&(rop->val0), &(a->val0), &(b->val0));
	Fq_add(&(rop->val1), &(a->val1), &(b->val1));
}
static inline void Fq2_sub(Fq2* rop, Fq2* a, Fq2* b)
{
	Fq_sub(&(rop->val0), &(a->val0), &(b->val0));
	Fq_sub(&(rop->val1), &(a->val1), &(b->val1));
}
static inline void Fq2_mont_mul(Fq2* rop, Fq2* a, Fq2* b)
{
	Fq temp00, temp11, temp01, temp10;

	Fq_mont_mul(&temp00, &(a->val0), &(b->val0));
	Fq_mont_mul(&temp11, &(a->val1), &(b->val1));
	Fq_mont_mul(&temp01, &(a->val0), &(b->val1));
	Fq_mont_mul(&temp10, &(a->val1), &(b->val0));

	Fq_sub(&(rop->val0), &temp00, &temp11);
	Fq_add(&(rop->val1), &temp01, &temp10);
}
static inline void Fq2_add_inv(Fq2* rop, Fq2* a)
{
	Fq_add_inv(&(rop->val0), &(a->val0));
	Fq_add_inv(&(rop->val1), &(a->val1));
}
static inline void Fq2_mont_mul_inv(Fq2* rop, Fq2* a)
{
	/*
	
	a^-1 = (a0 * a0 + a1 * a1)^(-1) * a0 + (a0 * a0 + a1 * a1)^(-1) * (-a1) * u.
	
	*/
	
	Fq temp00, temp11, abs, inv, neg1;

	Fq_mont_mul(&temp00, &(a->val0), &(a->val0));
	Fq_mont_mul(&temp11, &(a->val1), &(a->val1));
	Fq_add(&abs, &temp00, &temp11);
	Fq_mont_mul_inv(&inv, &abs);
	Fq_add_inv(&neg1, &(a->val1));

	Fq_mont_mul(&(rop->val0), &inv, &(a->val0));
	Fq_mont_mul(&(rop->val1), &inv, &neg1);
}

static inline void Fq2_transform(Fq2* rop, Fq2* a) // rop<-a*u
{
	Fq neg1, temp;

	Fq_assign(&temp, &(a->val0));
	Fq_add_inv(&neg1, &(a->val1));

	Fq_assign(&(rop->val0), &neg1);
	Fq_assign(&(rop->val1), &temp);
}

// print
void Fq2_print(Fq2* rop);
void Fq2_print_dec(Fq2* a);

#endif
