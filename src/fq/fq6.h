/*
* Fq6 = Fq2[v]/(v^3 - (u + 1)).
*
*/

#ifndef BLS12_381_fq6_H
#define BLS12_381_fq6_H

#include "fq.h"
#include "fq2.h"
#include <gmp.h>
#include <stdio.h>

typedef struct
{

	Fq2 val0;
	Fq2 val1;
	Fq2 val2;

} Fq6;

// init
void Fq6_set(Fq6* a, Fq2* num0, Fq2* num1, Fq2* num2);
void Fq6_set_zero(Fq6* a); // set a to zero.
void Fq6_set_id(Fq6* a); // set a to 1.
void Fq6_assign(Fq6* rop, Fq6* a);

// conversion
void Fq2_to_Fq6(Fq6* rop, Fq2* a);

// Montgomery representation
void Fq6_mont_rep(Fq6* aR, Fq6* a);
void Fq6_mont_rep_inv(Fq6* a, Fq6* aR);

// comparison
uint32_t Fq6_cmp(Fq6* a, Fq6* b); // if a=b, return 1; otherwise, return 0.

// Determinant
static inline void Fq2_determinant(Fq2* rop, Fq6* a) // a helper function to evalute the multiplicative inverse.
{
	/*
	
	det M = a0^3 - 3a0a1a2(u+1) + a1^3(u+1) + 2a2^3 u = Norm(a) = a * a^(q^2) *
		* a^(q^4).
	
	*/
	
	// auxiliary init
	Fq2 temp000, temp111, temp222, temp012, u1, c3temp012u1, temp111u1,
		c2temp222u, det;

	// auxiliary computation
	Fq2_mont_mul(&temp000, &(a->val0), &(a->val0));
	Fq2_mont_mul(&temp000, &temp000, &(a->val0));

	Fq2_mont_mul(&temp111, &(a->val1), &(a->val1));
	Fq2_mont_mul(&temp111, &temp111, &(a->val1));

	Fq2_mont_mul(&temp222, &(a->val2), &(a->val2));
	Fq2_mont_mul(&temp222, &temp222, &(a->val2));

	Fq2_mont_mul(&temp012, &(a->val0), &(a->val1));
	Fq2_mont_mul(&temp012, &temp012, &(a->val2));

	Fq2_set_ui(&u1, 1, 1);
	Fq2_mont_rep(&u1, &u1);

	Fq2_add(&c3temp012u1, &temp012, &temp012);
	Fq2_add(&c3temp012u1, &c3temp012u1, &temp012);
	Fq2_mont_mul(&c3temp012u1, &c3temp012u1, &u1);

	Fq2_mont_mul(&temp111u1, &temp111, &u1);

	Fq2_add(&c2temp222u, &temp222, &temp222);
	Fq2_transform(&c2temp222u, &c2temp222u);

	// computation
	Fq2_sub(&det, &temp000, &c3temp012u1);
	Fq2_add(&det, &det, &temp111u1);
	Fq2_add(&det, &det, &c2temp222u);

	// storing
	Fq2_assign(rop, &det);
}

// Arithmetic
static inline void Fq6_add(Fq6* rop, Fq6* a, Fq6* b)
{
	Fq2_add(&(rop->val0), &(a->val0), &(b->val0));
	Fq2_add(&(rop->val1), &(a->val1), &(b->val1));
	Fq2_add(&(rop->val2), &(a->val2), &(b->val2));
}
static inline void Fq6_sub(Fq6* rop, Fq6* a, Fq6* b)
{
	Fq2_sub(&(rop->val0), &(a->val0), &(b->val0));
	Fq2_sub(&(rop->val1), &(a->val1), &(b->val1));
	Fq2_sub(&(rop->val2), &(a->val2), &(b->val2));
}
static inline void Fq6_mont_mul(Fq6* rop, Fq6* a, Fq6* b)
{
	/*
	
	a*b = a0b0 + (a0b1 + a1b0)v + (a0b2 + a1b1 + a2b0)v^2 + (a1b2 + a2b1)v^3 +
		a2b2v^4 = [v^3 = u + 1] =

	= [a0b0 + (a1b2 + a2b1) + (a1b2 + a2b1)u] + [a0b1 + a1b0 + (a2b2) + (a2b2)u]v +
		[a0b2 + a1b1 + a2b0]v^2.
	
	*/
	
	// auxiliary init
	Fq2 temp00, temp01, temp02, temp10, temp11, temp12, temp20, temp21, temp22,
		u1, c0, c1, c2;

	// auxiliary computation
	Fq2_mont_mul(&temp00, &(a->val0), &(b->val0));
	Fq2_mont_mul(&temp01, &(a->val0), &(b->val1));
	Fq2_mont_mul(&temp02, &(a->val0), &(b->val2));
	Fq2_mont_mul(&temp10, &(a->val1), &(b->val0));
	Fq2_mont_mul(&temp11, &(a->val1), &(b->val1));
	Fq2_mont_mul(&temp12, &(a->val1), &(b->val2));
	Fq2_mont_mul(&temp20, &(a->val2), &(b->val0));
	Fq2_mont_mul(&temp21, &(a->val2), &(b->val1));
	Fq2_mont_mul(&temp22, &(a->val2), &(b->val2));
	
	Fq2_set_ui(&u1, 1, 1);
	Fq2_mont_rep(&u1, &u1);

	// computation
	Fq2_add(&c0, &temp12, &temp21);
	Fq2_mont_mul(&c0, &c0, &u1);
	Fq2_add(&c0, &c0, &temp00);

	Fq2_mont_mul(&c1, &temp22, &u1);
	Fq2_add(&c1, &c1, &temp10);
	Fq2_add(&c1, &c1, &temp01);

	Fq2_add(&c2, &temp02, &temp11);
	Fq2_add(&c2, &c2, &temp20);

	// storing
	Fq2_assign(&(rop->val0), &c0);
	Fq2_assign(&(rop->val1), &c1);
	Fq2_assign(&(rop->val2), &c2);
}
static inline void Fq6_add_inv(Fq6* rop, Fq6* a)
{

	Fq2_add_inv(&(rop->val0), &(a->val0));
	Fq2_add_inv(&(rop->val1), &(a->val1));
	Fq2_add_inv(&(rop->val2), &(a->val2));
}
static inline void Fq6_mont_mul_inv(Fq6* rop, Fq6* a)
{
	/*
	
	a^(-1) = detM ^(-1) * ((a0^2 - a1a2(u+1)) + (a2^2(u+1) - a0a1)v + (a1^2 -
		- a0a2)v^2).
	
	*/
	
	// auxiliary init
	Fq2 temp00, temp01, temp02, temp11, temp12, temp22, u1, temp12u1, temp22u1,
		c0, c1, c2;

	// determinant
	Fq2 detM, detMinv;
	Fq2_determinant(&detM, a);
	Fq2_mont_mul_inv(&detMinv, &detM);

	// auxiliary computation
	Fq2_mont_mul(&temp00, &(a->val0), &(a->val0));
	Fq2_mont_mul(&temp01, &(a->val0), &(a->val1));
	Fq2_mont_mul(&temp02, &(a->val0), &(a->val2));
	Fq2_mont_mul(&temp11, &(a->val1), &(a->val1));
	Fq2_mont_mul(&temp12, &(a->val1), &(a->val2));
	Fq2_mont_mul(&temp22, &(a->val2), &(a->val2));
	
	Fq2_set_ui(&u1, 1, 1);
	Fq2_mont_rep(&u1, &u1);
	
	Fq2_mont_mul(&temp12u1, &temp12, &u1);
	Fq2_mont_mul(&temp22u1, &temp22, &u1);

	// computatioon
	Fq2_sub(&c0, &temp00, &temp12u1);
	Fq2_sub(&c1, &temp22u1, &temp01);
	Fq2_sub(&c2, &temp11, &temp02);

	// normalization
	Fq2_mont_mul(&c0, &c0, &detMinv);
	Fq2_mont_mul(&c1, &c1, &detMinv);
	Fq2_mont_mul(&c2, &c2, &detMinv);

	// storing
	Fq2_assign(&(rop->val0), &c0);
	Fq2_assign(&(rop->val1), &c1);
	Fq2_assign(&(rop->val2), &c2);
}


// Transformation
static inline void Fq6_transform(Fq6* rop, Fq6* a) // rop <- a*v
{
	Fq2 temp0, temp1, u1, temp2u1;

	Fq2_assign(&temp0, &(a->val0));
	Fq2_assign(&temp1, &(a->val1));
	
	Fq2_set_ui(&u1, 1, 1);
	Fq2_mont_rep(&u1, &u1);
	
	Fq2_mont_mul(&temp2u1, &(a->val2), &u1);

	Fq2_assign(&(rop->val0), &temp2u1);
	Fq2_assign(&(rop->val1), &temp0);
	Fq2_assign(&(rop->val2), &temp1);
}

// print
void Fq6_print(Fq6* rop);
void Fq6_print_dec(Fq6* rop);

#endif
