/* 
 * The computation of the Tate and ate pairings on the curve BLS12-381.
 * 
 * For two points P and Q known to be in their respective groups G1 and G2,
 *      please use the following functions to compute the desired pairing:
 *
 * void tate_pairing(Fq12* rop, G1* P, G2* Q) - compute the Tate pairing of points
 *      P and Q using Miller's algorithm and store the result in rop.
 *
 * void ate_pairing(Fq12* rop, G1* P, G2* Q) - compute the ate pairing of points
 *      P and Q using Miller's algorithm and store the result in rop.
 *
 * ===CONTENTS==================================================================
 *
 * -- Variable and Struct Descriptions -----------------------------------------
 * 
 * struct G12 - an auxiliary struct to store coordinates of a point from G2
 *      after performing the sextic untwisting.
 *
 * struct G12_Z - same, projective space
 *
 * char* xPOW_STR - a string containting the hexadecimal representation of 
 *      (q^6 + 1)/r.
 *
 * char* bPOW_STR - same, binary.
 *
 * int len_bPOW_STR - the length of bPOW_STR;
 *
 * Fq2 sextic_const - a multiplier (1/2, -1/2), for the sextic untwisting.
 *
 * -- Function Descriptions ----------------------------------------------------
 *
 * Initialisation functions
 * ---------------
 * void G12_init(G12* a), void G12_clear(G12* a) - init./clear the G12 element a.
 * void G12_Z_init(G12_Z* a), void G12_Z_clear(G12_Z* a) - same, projective space.
 * 
 * void sextic_precompute() - compute of Fq2 sextic_const.
 * void sextic_clear() - free the space occupied by Fq2 sextic_const.
 * void tate_precompute() - compute xPOW_STR, bPOW_STR, and len_bPOW_STR.
 * void tate_clear() - free the space ocupied by xPOW_STR and bPOW_STR. 
 * 
 * Sextic untwist
 * -----------------
 * void sextic_untwist(G12* rop, G2* a) - perform the untwisting and store the result 
 *      in rop.
 *
 * void sextic_untwist_project(G12_Z* rop, G2_Z* a) - same, projective space
 *
 * Tate pairing
 * -----------------
 * void tate_line_doubling(Fq12* rop, G1_Z* RZ, G12* Q) - find a tangent line 
 *      to the point RZ and evaluate it at Q. Store the result in rop. This is 
 *      a helper function for Miller's algorithm.
 *
 * void tate_line_adding(Fq12* rop, G1_Z* RZ, G1_Z* PZ, G12* Q) - find a line
 *      passing through points RZ and PZ and evaluate it at Q. Store the result
 *      in rop. This is a helper function for Miller's algorithm.
 *
 * void tate(Fq12* rop, G1* P, G2* Q) - execute Miller's algorithm for the points
 *      P and Q prior to the final exponentiation and store the result in rop. 
 *      This is a helper function.
 *
 * void tate_exp(Fq12* rop, Fq12* t) - compute t to the power (q^12 - 1)/r and
 *      store the result in rop. We optimize exponentiation, since (q^12 - 1)/r =
 *      = (q^6 - 1) * (q^6 + 1)/r, and raising to the power (q^6 - 1) is easy.
 *      This is a helper function.
 *
 * void tate_pairing(Fq12* rop, G1* P, G2* Q) - compute the Tate pairing of points
 *      P and Q using Miller's algorithm and store the result in rop.
 *
 * Ate pairing
 * -----------------
 * void ate_line_doubling(Fq12* rop, G2_Z* RZ, G1* P) - find a tangent line 
 *      to the point RZ and evaluate it at P. Store the result in rop. This is 
 *      a helper function for Miller's algorithm.
 *
 * void ate_line_adding(Fq12* rop, G2_Z* RZ, G2_Z* QZ, G1* P) - find a line
 *      passing through points RZ and QZ and evaluate it at P. Store the result
 *      in rop. This is a helper function for Miller's algorithm.
 *
 * void ate(Fq12* rop, G1* P, G2* Q) - execute Miller's algorithm for the points
 *      P and Q prior to the final exponentiation and store the result in rop. 
 *      This is a helper function.
 *
 * void ate_exp(Fq12* rop, Fq12* a) - compute t to the power (q^12 - 1)/r and
 *      store the result in rop. We optimize exponentiation, since (q^12 - 1)/r =
 *      = (q^6 - 1) * (q^6 + 1)/r, and raising to the power (q^6 - 1) is easy.
 *      This is a helper function.
 *
 * void ate_pairing(Fq12* rop, G1* P, G2* Q) - compute the ate pairing of points
 *      P and Q using Miller's algorithm and store the result in rop.
 *
 * -- Note --------------------------------------------------------------------
 * At first glance, one would need to test whether two given points belong to their 
 * respective groups G1 and G2. However, in practice, it is often guarunteed that 
 * points belong to proper domains.
 *
 * For any additional info about Miller's Algorithm, please address
 *      "Computing the Tate pairing" by M. Scott and
 *      "Pairing for beginners" by C. Costello, Alg. 7.1 and 7.2
 *
*/ 

#ifndef BLS12_381_MILLER
#define BLS12_381_MILLER

#include "../extensions/fq12.h"
#include "../extensions/fq2.h"
#include "../extensions/fq6.h"
#include "../fq/bls_params.h"
#include "../fq/fq.h"
#include "../groups/g1.h"
#include "../groups/g2.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{

	Fq12 valx;
	Fq12 valy;

} G12;

typedef struct
{

	Fq12 valX;
	Fq12 valY;
	Fq12 valZ;

} G12_Z;

extern char* xPOW_STR;
extern char* bPOW_STR;
extern int len_bPOW_STR;

extern Fq2 sextic_const;

// init
void G12_init(G12* a);
void G12_clear(G12* a);
void G12_Z_init(G12_Z* a);
void G12_Z_clear(G12_Z* a);
void sextic_precompute();
void sextic_clear();
void tate_precompute();
void tate_clear();

// sextic untwist
void sextic_untwist(G12* rop, G2* a);
void sextic_untwist_project(G12_Z* rop, G2_Z* a);

// Tate pairing
void tate_line_doubling(Fq12* rop, G1_Z* RZ, G12* Q);
void tate_line_adding(Fq12* rop, G1_Z* RZ, G1_Z* PZ, G12* Q);
void tate(Fq12* rop, G1* P, G2* Q);
void tate_exp(Fq12* rop, Fq12* t);
void tate_pairing(Fq12* rop, G1* P, G2* Q);

// Ate pairing
void ate_line_doubling(Fq12* rop, G2_Z* RZ, G1* P);
void ate_line_adding(Fq12* rop, G2_Z* RZ, G2_Z* QZ, G1* P);
void ate(Fq12* rop, G1* P, G2* Q);
void ate_exp(Fq12* rop, Fq12* a);
void ate_pairing(Fq12* rop, G1* P, G2* Q);

#endif
