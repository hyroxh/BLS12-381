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
* mpz_t orderG1 - the order of the group G1.
*
* G1 genG1 - the generator of the group G1.
*
* mpz_t cofactorG1 - the cofactor of the group G1.
*
* char* bORD_STR - the order of G1, written as a binary string.
*
* int len_bORD_STR - the length of that string.
*
* -- Function Descriptions ----------------------------------------------------
*
* Initialisation functions
* ---------------
* void G1_group_init(), G1_group_clear() - set/clear the global parameters for
*      the group G1. See the source for reference.
*
* void G1_init(G1* a) - initialize a, if not infinity.
* void G1_set(G1* a, Fq* x, Fq* y) - set a to (x, y).
* void G1_init_inf(G1* a) - initialize a, if infinity.
* void G1_set_inf(G1* a) - set a to the infinity point.
* void G1_assign(G1* rop, G1* a) - set rop to a.
* void G1_clear(G1* a) - free the space occupied by a.
*
* void G1_Z_init(G1_Z* a),
* void G1_Z_set(G1_Z* a, Fq* X, Fq* Y, Fq* Z),
* void G1_Z_set_inf(G1_Z* a),
* void G1_Z_assign(G1_Z* rop, G1_Z* a),
* void G1_Z_clear(G1_Z* a),
*      - same, projective space.
*
* Mapping functions
* ---------------
* void G1_to_group(G1* rop, G1* a) - set rop to [cofactorG1]a.
*      This operation maps a curve element to the group G1.
*
* void G1_Z_to_group(G1_Z* rop, G1_Z* a) - same, projective space.
*
* void G1_to_affine(G1* rop, G1_Z* a) - map a to the affine coordinates.
*      Set rop to the result.
*
* void G1_to_project(G1_Z* rop, G1* a) - map a to the projective coordinates.
*      Set rop to the result.
*
* Checking functions
* ---------------
* int G1_curve_check(G1* a) - return 1, if a is in E(Fq); return 0, otherwise.
*
* int G1_group_check(G1* a) - return 1, if a is in G1 as specified; return 0, otherwise.
*      To check if a is in G1, it suffices to verify that [orderG1]a is infinity.
*
* int G1_Z_curve_check(G1_Z* a), G1_Z_group_check(G1_Z* a) - same, projective space.
*
* Arithmetic functions
* ---------------
*
* void G1_add_inv(G1* rop, G1* a) - set rop to (-a)
*
* int G1_Z_cmp(G1_Z* p, G1_Z* q) - return 0, if p == q; return 1, otherwise
*
* void G1_add(G1* rop, G1* p, G1* q) - set rop to p + q.
*
* void G1_mul(G1* rop, mpz_t n, G1* p) - set rop to [n]p = p + p + ...
*
* void G1_Z_add_inv(G1_Z* rop, G1_Z* a),
* void G1_Z_add(G1_Z* rop, G1_Z* p, G1_Z* q),
* void G1_Z_mul(G1_Z* rop, mpz_t n, G1_Z* p), - same, projective space.
*
* Print 
* ---------------
* void G1_print(G1* a) - print the element a.
* void G1_Z_print(G1_Z* a) - same, projective space.
*
* -- Note --------------------------------------------------------------------
* For efficiency, it's recommended to use the projective space functions. 
*
*/


#ifndef BLS12_381_G1_H
#define BLS12_381_G1_H

#include "../fq/bls_params.h"
#include "../fq/fq.h"
#include <gmp.h>
#include <stdio.h>

typedef struct
{ 
	Fq valx;
	Fq valy;
	char inf;

} G1;

typedef struct
{

	Fq valX;
	Fq valY;
	Fq valZ;

} G1_Z;

extern mpz_t orderG1;
extern G1 genG1;
extern mpz_t cofactorG1;
extern char* bORD_STR;
extern int len_bORD_STR;

// Group Init
void G1_group_init();
void G1_group_clear();

// Init
void G1_init(G1* a);
void G1_set(G1* a, Fq* x, Fq* y);
void G1_init_inf(G1* a);
void G1_set_inf(G1* a);
void G1_assign(G1* rop, G1* a);
void G1_clear(G1* a);

// Init projective
void G1_Z_init(G1_Z* a);
void G1_Z_set(G1_Z* a, Fq* X, Fq* Y, Fq* Z);
void G1_Z_set_inf(G1_Z* a);
void G1_Z_assign(G1_Z* rop, G1_Z* a);
void G1_Z_clear(G1_Z* a);

// Maps
void G1_to_group(G1* rop, G1* a);
void G1_Z_to_group(G1_Z* rop, G1_Z* a);
void G1_to_affine(G1* rop, G1_Z* a);
void G1_to_project(G1_Z* rop, G1* a);

// Curve check
int G1_curve_check(G1* a);
int G1_group_check(G1* a);

// Curve check projective
int G1_Z_curve_check(G1_Z* a);
int G1_Z_group_check(G1_Z* a);

// Group Arithmetic
void G1_add_inv(G1* rop, G1* a);
int G1_Z_cmp(G1_Z* p, G1_Z* q); 
void G1_add(G1* rop, G1* p, G1* q);
void G1_mul(G1* rop, mpz_t n, G1* p);

// Group Arithmetic Projective
void G1_Z_add_inv(G1_Z* rop, G1_Z* a);
void G1_Z_add(G1_Z* rop, G1_Z* p, G1_Z* q);
void G1_Z_mul(G1_Z* rop, mpz_t n, G1_Z* p);

// Print
void G1_print(G1* a);
void G1_Z_print(G1_Z* a);

#endif
