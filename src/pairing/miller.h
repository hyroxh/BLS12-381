/* 
 * The computation of the ate pairing on the curve BLS12-381.
 * 
 * For two points P and Q known to be in their respective groups G1 and G2,
 *      please use the following function to compute the ate pairing:
 *
 * void ate_pairing(Fq12* rop, G1* P, G2_Z* Q) - compute the ate pairing of points
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
 * char* bPOW_STR - a string containting the binary representation of 
 *      (q^6 + 1)/r.
 *
 * int len_bPOW_STR - the length of bPOW_STR;
 *
 * char* bATE_STR - a string containing the binary representation of 
 *      the curve parameter x.
 *
 * int len_bATE_STr - its length.
 *
 * -- Function Descriptions ----------------------------------------------------
 *
 * Initialisation functions
 * ---------------
 * void exp_precompute() - set bPOW_STR and len_bPOW_STR.
 * void ate_precompute() - set bATE_STR and len_bATE_STR.
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
 * void ate(Fq12* rop, G1* P, G2_Z* QZ) - execute Miller's algorithm for the points
 *      P and Q prior to the final exponentiation and store the result in rop. 
 *      This is a helper function.
 *
 * void ate_exp(Fq12* rop, Fq12* a) - compute t to the power (q^12 - 1)/r and
 *      store the result in rop. We optimize exponentiation, since (q^12 - 1)/r =
 *      = (q^6 - 1) * (q^6 + 1)/r, and raising to the power (q^6 - 1) is easy.
 *      This is a helper function.
 *
 * void ate_pairing(Fq12* rop, G1* P, G2_Z* QZ) - compute the ate pairing of points
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

#include "../fq/fq12.h"
#include "../fq/fq6.h"
#include "../fq/fq2.h"
#include "../fq/fq.h"
#include "../fr/fr.h"
#include "../groups/g1.h"
#include "../groups/g2.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern char* bPOW_STR;
extern int len_bPOW_STR;

extern char* bATE_STR;
extern int len_bATE_STR;

// init
void exp_precompute();
void ate_precompute();

// Ate pairing
void ate_line_doubling(Fq12* rop, G2_Z* RZ, G1* P);
void ate_line_adding(Fq12* rop, G2_Z* RZ, G2_Z* QZ, G1* P);
void ate(Fq12* rop, G1* P, G2_Z* QZ);
void ate_exp(Fq12* rop, Fq12* a);
void ate_pairing(Fq12* rop, G1* P, G2_Z* QZ);

#endif
