/*
* G2 in E(Fq) - a group on the elliptic curve E(Fq): y^2 = x^3 + 4(1+u).
* It is specifically chosen to make pairings defined and computationally efficient.
*      Source: https://eth2book.info/latest/part2/building_blocks/bls12-381/
*
* The functions and notations used follow those of G2. 
*
*/

#ifndef BLS12_381_G2_H
#define BLS12_381_G2_H

#include "../fq/fq.h"
#include "../fq/fq2.h"
#include "../fr/fr.h"
#include <gmp.h>
#include <stdio.h>

typedef struct
{
	Fq2 valx;
	Fq2 valy;
	uint32_t inf;

} G2;

typedef struct
{

	Fq2 valX;
	Fq2 valY;
	Fq2 valZ;

} G2_Z;

extern Fr orderG2;
extern G2 genG2;
extern Fr cofactorG2;

// Group Init
void G2_group_init();

// Init
void G2_set(G2* a, Fq2* x, Fq2* y);
void G2_set_inf(G2* a);
void G2_assign(G2* rop, G2* a);

// Init projective
void G2_Z_set(G2_Z* a, Fq2* X, Fq2* Y, Fq2* Z);
void G2_Z_set_inf(G2_Z* a);
void G2_Z_assign(G2_Z* rop, G2_Z* a);

// Maps
void G2_Z_to_group(G2_Z* rop, G2_Z* a);
void G2_to_affine(G2* rop, G2_Z* a);
void G2_to_project(G2_Z* rop, G2* a);

// Curve check projective
int G2_Z_curve_check(G2_Z* a);
int G2_Z_group_check(G2_Z* a);

// Montgomery representation
void G2_mont_rep(G2* aR, G2* a);
void G2_mont_rep_inv(G2* a, G2* aR);
void G2_Z_mont_rep(G2_Z* aR, G2_Z* a);
void G2_Z_mont_rep_inv(G2_Z* a, G2_Z* aR);

//correction
void G2_Z_correction(G2_Z* a);

// Group Arithmetic Projective
void G2_Z_add_inv(G2_Z* rop, G2_Z* a);
uint32_t G2_Z_cmp(G2_Z* p, G2_Z* q);
void G2_Z_adding(G2_Z* rop, G2_Z* p, G2_Z* q);
void G2_Z_doubling(G2_Z* rop, G2_Z* p); 
void G2_Z_add(G2_Z* rop, G2_Z* p, G2_Z* q);
void G2_Z_mul(G2_Z* rop, Fr* n, G2_Z* p);

// Print
void G2_print(G2* a);
void G2_Z_print(G2_Z* a);
void G2_print_dec(G2* a);
void G2_Z_print_dec(G2_Z* a);

#endif
