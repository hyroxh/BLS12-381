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

char* bPOW_STR;

int len_bPOW_STR;

char* bATE_STR;

int len_bATE_STR;

void exp_precompute()
{
	bPOW_STR = "10100010110011000101001000011101110101000000110111101101101"
	"1110010001101011100010101011100101001110101011001101101101111111000010"
	"1101011100110011001001100111101011011101100110110010101111001101011101"
	"0101000001111000100100011000001110110101110000011101100100100000011000"
	"1000101111001101111011110110011001010110101110011011101011010001101110"
	"1100011101111011111011111001100111101010110111110110001010101110011101"
	"1111010101010100011000101010010001100111000001000000010011011110001011"
	"1110110000110101111101101000110111000011001011110111101001011111010010"
	"0100010011001111100001100010100000001001001101100100000000010111011101"
	"1001000010110100010111001110000011111110000100001000001100001010101010"
	"0011101000100010010010111111110001011001011110010100100100010100101100"
	"1111000011110001111111010111100101110010101110100011111000100110000010"
	"0001001110101111011110100100110011101111111111011000100101101011010101"
	"0000111010011010010000110110111001111011010001010111000001000101111010"
	"1000101001010100111000100111101101111111011100110101010010101110110110"
	"0010001111001100011111101010110010111101000100001101100100101001100010"
	"0111110000010001110000100110000111011010111100010111101010101011111000"
	"0001011000101011111001001111101011110111101100100001001001101010000001"
	"0000111000011111100000000011111000000000111100111111010111110001110101"
	"1111100100000010110000100000001110111011101000001110110000100010111110"
	"1000111010110000101011101010010011110001001110110010011010111101111000"
	"0111011010111101111010010110001100010001010001000001011000000111000111"
	"0111000011100110101000111011110001110000011011110011100110011100000100"
	"1011001110000100010101000010010101111010100000000100010010001110100000"
	"1101100110110111111101100011101110110000000101101011100100111000101010"
	"1100011100010010000111100010011001100111010000010011100010001001001011"
	"1100100000011111101110110111010011100111100001111011100001010011000011"
	"1000111100100011110001000001001101001010010010101101101111000000011100"
	"000011010000101101000011100111001111000011100110111000000011100000101110101101010"; //(q^6 + 1)/r
	len_bPOW_STR = 2030;
}

void ate_precompute()
{
	bATE_STR = "1101001000000001000000000000000000000000000000010000000000000000"; //T = |x|, with x = t - 1, where t is the Frobenius trace
	
	len_bATE_STR = 64;
}


//-------Ate pairing--------------------------------

void ate_line_doubling(Fq12* rop, G2_Z* RZ, G1* P) // Lr,r (P)
{
	/*

	Computing L_R_R (P)

	L_R_R (P) = (2*YrZr^2) * Yp - 3*Xr^2Zr * Xp + 3Xr^3 - 2Yr^2Zr ->

	-> (2*YrZr^2) * Yp * const * vw - 3*Xr^2Zr * Xp * const * v + 3Xr^3 * const
	- 2Yr^2Zr * const,

	where const = sextic_const = (1-u)/2 <--- no need to multiply due to Frobenius endomorphism

	*/

	// init
	Fq zero;
	Fq_set_ui(&zero, 0);

	Fq2 P2x, P2y;
	Fq2_set(&P2x, &(P->valx), &zero);
	Fq2_set(&P2y, &(P->valy), &zero);

	// aux init
	Fq2 aux1, aux2, aux3, aux4;

	// aux comp
	Fq2_add(&aux1, &(RZ->valY), &(RZ->valY));
	Fq2_mont_mul(&aux1, &aux1, &(RZ->valZ));
	Fq2_mont_mul(&aux1, &aux1, &(RZ->valZ));

	Fq2_add(&aux2, &(RZ->valX), &(RZ->valX));
	Fq2_add(&aux2, &aux2, &(RZ->valX));
	Fq2_mont_mul(&aux2, &aux2, &(RZ->valX));
	Fq2_mont_mul(&aux2, &aux2, &(RZ->valZ));

	Fq2_add(&aux3, &(RZ->valX), &(RZ->valX));
	Fq2_add(&aux3, &aux3, &(RZ->valX));
	Fq2_mont_mul(&aux3, &aux3, &(RZ->valX));
	Fq2_mont_mul(&aux3, &aux3, &(RZ->valX));

	Fq2_add(&aux4, &(RZ->valY), &(RZ->valY));
	Fq2_mont_mul(&aux4, &aux4, &(RZ->valY));
	Fq2_mont_mul(&aux4, &aux4, &(RZ->valZ));

	Fq2_mont_mul(&aux1, &aux1, &P2y);
	Fq2_mont_mul(&aux2, &aux2, &P2x);
	Fq2_sub(&aux3, &aux3, &aux4);

	// untwisting

	Fq2 zero2;
	Fq2_set_ui(&zero2, 0, 0);

	Fq6 temp1, temp2, temp3, zero6;

	Fq6_set_zero(&zero6);

	Fq6_set(&temp1, &zero2, &aux1, &zero2);
	Fq6_set(&temp2, &zero2, &aux2, &zero2);
	Fq6_set(&temp3, &aux3, &zero2, &zero2);

	Fq12 ut1, ut2, ut3, L;

	Fq12_set(&ut1, &zero6, &temp1);
	Fq12_set(&ut2, &temp2, &zero6);
	Fq12_set(&ut3, &temp3, &zero6);

	Fq12_sub(&L, &ut1, &ut2);
	Fq12_add(&L, &L, &ut3);

	Fq12_assign(rop, &L);
}

void ate_line_adding(Fq12* rop, G2_Z* RZ, G2_Z* QZ, G1* P) // L_R,Q (P)
{
	/*

	Computing L_R,Q (P)

	L: (Xq*Zr - Xr*Zq) * Yp + (YrZq - YqZr) * Xp + (XrYq - XqYr) ->

	(Xq*Zr - Xr*Zq) * const * Yp * v^2 + (YrZq - YqZr) * const * Xp * v * w +
	(XrYq - XqYr) * const * w,

	where const = (1-u)/2 = sextic_const <--- no need to multiply due to Frobenius endomorphism

	*/

	// init
	Fq zero;
	Fq_set_ui(&zero, 0);

	Fq2 P2x, P2y;
	Fq2_set(&P2x, &(P->valx), &zero);
	Fq2_set(&P2y, &(P->valy), &zero);

	// aux init
	Fq2 aux1, aux2, aux3, aux4, aux5, aux6, c1, c2, c3;

	// aux computation
	Fq2_mont_mul(&aux1, &(QZ->valX), &(RZ->valZ));
	Fq2_mont_mul(&aux2, &(RZ->valX), &(QZ->valZ));
	Fq2_mont_mul(&aux3, &(RZ->valY), &(QZ->valZ));
	Fq2_mont_mul(&aux4, &(QZ->valY), &(RZ->valZ));
	Fq2_mont_mul(&aux5, &(RZ->valX), &(QZ->valY));
	Fq2_mont_mul(&aux6, &(QZ->valX), &(RZ->valY));

	Fq2_sub(&c1, &aux1, &aux2);
	Fq2_sub(&c2, &aux3, &aux4);
	Fq2_sub(&c3, &aux5, &aux6);

	Fq2_mont_mul(&c1, &c1, &P2y);
	Fq2_mont_mul(&c2, &c2, &P2x);

	// untwisting
	Fq2 zero2;
	Fq2_set_ui(&zero2, 0, 0);

	Fq6 temp1, temp2, temp3, zero6;
	Fq6_set_zero(&zero6);

	Fq6_set(&temp1, &zero2, &zero2, &c1);
	Fq6_set(&temp2, &zero2, &c2, &zero2);
	Fq6_set(&temp3, &c3, &zero2, &zero2);

	Fq12 ut1, ut2, ut3, L;
	
	Fq12_set(&ut1, &temp1, &zero6);
	Fq12_set(&ut2, &zero6, &temp2);
	Fq12_set(&ut3, &zero6, &temp3);

	Fq12_add(&L, &ut1, &ut2);
	Fq12_add(&L, &L, &ut3);

	Fq12_assign(rop, &L);
}

void ate(Fq12* rop, G1* P, G2_Z* QZ) // pre-exponent execution
{
	
	Fq zero;
	Fq_set_zero(&zero);
	
	uint32_t mask = (P->inf == 1) | ( (Fq_cmp(&QZ->valZ.val0, &zero) == 1) & (Fq_cmp(&QZ->valZ.val1, &zero) == 1) );
	mask = - mask; //0xFFFFFFFF, if either P or Q are inf; 0, otherwise
	
	G2_Z RZ;
	G2_Z_assign(&RZ, QZ);

	Fq12 f;
	Fq12_set_id(&f);
	Fq12_mont_rep(&f, &f);

	Fq12 L;

	for (int i = 1; i < len_bATE_STR; i++)
	{
		ate_line_doubling(&L, &RZ, P);
		G2_Z_add(&RZ, &RZ, &RZ);
		Fq12_mont_mul(&f, &f, &f);
		Fq12_mont_mul(&f, &f, &L);

		if (bATE_STR[i] == '1')
		{
			ate_line_adding(&L, &RZ, QZ, P);
			G2_Z_add(&RZ, &RZ, QZ);
			Fq12_mont_mul(&f, &f, &L);
		}
	}

	Fq12_assign(rop, &f);
	
	rop->val0.val0.val0.limbs[0] = (rop->val0.val0.val0.limbs[0] & ~mask) | (1 & mask); //if either P or Q are inf, return 1. 
										//(This case, the rest is guaranteed to be 0).
}

void ate_exp(Fq12* rop, Fq12* a)
{
	Fq6 val0, val1;

	Fq12 temp, inv, aux;

	Fq6_assign(&val0, &(a->val0));
	Fq6_assign(&val1, &(a->val1));
	Fq6_add_inv(&val1, &val1);
	Fq12_mont_mul_inv(&inv, a);

	Fq12_set(&temp, &val0, &val1);
	Fq12_mont_mul(&temp, &temp, &inv);

	Fq12_assign(&aux, &temp);

	for (int i = 1; i < len_bPOW_STR; i++)
	{
		Fq12_mont_mul(&aux, &aux, &aux);
		if (bPOW_STR[i] == '1')
		{
			Fq12_mont_mul(&aux, &aux, &temp);
		}
	}
	
	//Fq12_mont_mul_inv(&aux, &aux);
	Fq12_assign(rop, &aux);
}

void ate_pairing(Fq12* rop, G1* P, G2_Z* QZ)
{
	Fq12 temp;

	ate(&temp, P, QZ);
	ate_exp(&temp, &temp);

	Fq12_assign(rop, &temp);
}

