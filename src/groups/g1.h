/*
* G1 in E(Fq) - a group on the elliptic curve E(Fq): y^2 = x^3 + 4.
* It is specifically chosen to make pairings defined and computationally efficient.
*      Source: https://eth2book.info/latest/part2/building_blocks/bls12-381/
*
* ===CONTENTS==================================================================
*
* -- Variable and Struct Descriptions -----------------------------------------
*
* struct G1 - a struct to contain coordinates.
*      Fq valx, valy - coordinates.
*      char inf - '0', if not a point at infinity. '1', otherwise.
*
* struct G1_Z - same, projective space.
*      Fq valX, valY, valZ - coordinates. Infinity corresponds to valZ = 0.
*
* Fr orderG1 - the order of the group G1.
*
* G1 genG1 - the generator of the group G1.
*
* Fr cofactorG1 - the cofactor of the group G1.
*
* char* bORD_STR - the order of G1, written as a binary string.
*
* int len_bORD_STR - the length of that string.
*
* -- Function Descriptions ----------------------------------------------------
*
* Initialisation functions
* ---------------
* void G1_group_init() - set the global parameters for
*      the group G1. See the source for reference.
*
* void G1_set(G1* a, Fq* x, Fq* y) - set a to (x, y).
* void G1_set_inf(G1* a) - set a to the infinity point.
* void G1_assign(G1* rop, G1* a) - set rop to a.
*
* void G1_Z_set(G1_Z* a, Fq* X, Fq* Y, Fq* Z),
* void G1_Z_set_inf(G1_Z* a),
* void G1_Z_assign(G1_Z* rop, G1_Z* a),
*      - same, projective space.
*
* Mapping functions
* ---------------
* void G1_Z_to_group(G1_Z* rop, G1_Z* a) - set rop to [cofactorG1]a, projective space.
*      This operation maps a curve element to the group G1.
*
* void G1_to_affine(G1* rop, G1_Z* a) - map a to the affine coordinates.
*      Set rop to the result.
*
* void G1_to_project(G1_Z* rop, G1* a) - map a to the projective coordinates.
*      Set rop to the result.
*
* Checking functions
* ---------------
* int G1_Z_curve_check(G1_Z* a) - return 1, if a is in E(Fq); return 0, otherwise.
*
* int G1_Z_group_check(G1_Z* a) - return 1, if a is in G1 as specified; return 0, otherwise.
*      To check if a is in G1, it suffices to verify that [orderG1]a is infinity.
*
* Montgomery representation
* ---------------
* void G1_mont_rep(G1* aR, G1* a) - set aR to the Montgomery representation of a.
*
* void G1_mont_rep_inv(G1* a, G1* aR) - revert the Montgomery representation
*      and set the result to a.
*
* void G1_Z_mont_rep(G1_Z* aR, G1_Z* a), void G1_Z_mont_rep_inv(G1_Z* a, G1_Z* aR)
*      - same, projective space.
*
* Correction
* ---------------
* void G1_Z_correction(G1_Z* a) - if a is inf, set it to (0:1:0).
*
* Arithmetic functions
* ---------------
* void G1_Z_add_inv(G1_Z* rop, G1_Z* a) - set rop to (-a)
*
* uint32_t G1_Z_cmp(G1_Z* p, G1_Z* q) - return 1, if p == q; return 0, otherwise
*
* void G1_Z_adding(G1_Z* rop, G1_Z* p, G1_Z* q) - set rop to p + q, if p!=q is known.
*      This is a helper function.
*
* void G1_Z_doubling(G1_Z* rop, G1_Z* p - set rop to [2]p.
*      This is a helper function.
*
* void G1_Z_add(G1_Z* rop, G1_Z* p, G1_Z* q) - set rop to p + q.
*
* void G1_Z_mul(G1_Z* rop, Fr* n, G1_Z* p) - set rop to [n]p = p + p + ...
*
* Print 
* ---------------
* void G1_print(G1* a) - print the element a in the limb form.
* void G1_Z_print(G1_Z* a) - same, projective space.
* void G1_print_dec(G1* a) - print the element a in its decimal form.
* void G1_Z_print_dec(G1_Z* a) - same, projective space.
*
* -- Note --------------------------------------------------------------------
* All arithmetic operations are constant-time.
*
*/


#ifndef BLS12_381_G1_H
#define BLS12_381_G1_H

#include "../fq/fq.h"
#include "../fr/fr.h"
#include <gmp.h>
#include <stdio.h>

typedef struct
{ 
	Fq valx;
	Fq valy;
	uint32_t inf;

} G1;

typedef struct
{

	Fq valX;
	Fq valY;
	Fq valZ;

} G1_Z;

extern Fr orderG1;
extern G1 genG1;
extern Fr cofactorG1;
extern char* bORD_STR;
extern int len_bORD_STR;

// Group Init
void G1_group_init();

// Init
void G1_set(G1* a, Fq* x, Fq* y);
void G1_set_inf(G1* a);
void G1_assign(G1* rop, G1* a);

// Init projective
void G1_Z_set(G1_Z* a, Fq* X, Fq* Y, Fq* Z);
void G1_Z_set_inf(G1_Z* a);
void G1_Z_assign(G1_Z* rop, G1_Z* a);

// Maps
void G1_Z_to_group(G1_Z* rop, G1_Z* a);
void G1_to_affine(G1* rop, G1_Z* a);
void G1_to_project(G1_Z* rop, G1* a);

// Curve check projective
int G1_Z_curve_check(G1_Z* a);
int G1_Z_group_check(G1_Z* a);

// Montgomery representation
void G1_mont_rep(G1* aR, G1* a);
void G1_mont_rep_inv(G1* a, G1* aR);
void G1_Z_mont_rep(G1_Z* aR, G1_Z* a);
void G1_Z_mont_rep_inv(G1_Z* a, G1_Z* aR);

//correction
void G1_Z_correction(G1_Z* a);

// Group Arithmetic Projective
void G1_Z_add_inv(G1_Z* rop, G1_Z* a);
uint32_t G1_Z_cmp(G1_Z* p, G1_Z* q);
void G1_Z_adding(G1_Z* rop, G1_Z* p, G1_Z* q);
void G1_Z_doubling(G1_Z* rop, G1_Z* p); 
void G1_Z_add(G1_Z* rop, G1_Z* p, G1_Z* q);
void G1_Z_mul(G1_Z* rop, Fr* n, G1_Z* p);

// Print
void G1_print(G1* a);
void G1_Z_print(G1_Z* a);
void G1_print_dec(G1* a);
void G1_Z_print_dec(G1_Z* a);

#endif
