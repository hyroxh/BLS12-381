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

char* bATE_STR;

int len_bATE_STR;

Fq2 g11, g12, g13, g14, g15;
Fq2 g21, g22, g23, g24, g25;

char* bX_STR;
char* bX2_STR;
char* bXP1_STR;
char* bXP1_DIV3_STR;

int len_bX_STR;
int len_bX2_STR;
int len_bXP1_STR;
int len_bXP1_DIV3_STR;

void exp_precompute()
{	
	bX_STR = "1101001000000001000000000000000000000000000000010000000000000000"; // |x|
	bXP1_STR = "1101001000000001000000000000000000000000000000010000000000000001";// |x| + 1
	bXP1_DIV3_STR = "100011000000000010101010101010101010101010101011010101010101011"; // (1 + |x|)/3
	bX2_STR = "10101100010001011010010000000001000000000000000110100100000000100000000000000000000000000000000100000000000000000000000000000000"; //x^2
	
	len_bX_STR = 64;
	len_bXP1_STR = 64;
	len_bXP1_DIV3_STR = 63;
	len_bX2_STR = 128;
	
	/*
	
	If f = h + gw, then f = g0 + h0 * w + g1 * w^2 + h1 * w^3 + g2 * w^4 + h2 * w^5.
	
	f^q = ~g0 + ~h0 * g11 * w + ~g1 * g12 *  w^2 + ~h1 * g13 * w^3 + ~g2 * g14 * w^4 + ~h2 * g15 * w^5, where ~a is conjugation of a.
	
	f^(q^2) = f = g0 + h0 * g21 * w + g1 * g22 * w^2 + h1 * g23 * w^3 + g2 * g24 * w^4 + h2 * g25 * w^5.
	
	Coefficients g11, ..., g15, g21, ..., g25 are given below. 
	
	Notice that g2j = g1j * ~g1j.	
	
	*/
	
	Fq2_set_hex_str(&g11, "1904d3bf02bb0667c231beb4202c0d1f0fd603fd3cbd5f4f7b2443d784bab9c4f67ea53d63e7813d8d0775ed92235fb8",
	"fc3e2b36c4e03288e9e902231f9fb854a14787b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec22cf78a126ddc4af3");
	
	Fq2_set_hex_str(&g12, "0", "1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac");
	
	Fq2_set_hex_str(&g13, "6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09",
	"6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09");
	
	Fq2_set_hex_str(&g14, "1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad", "0");
	
	Fq2_set_hex_str(&g15,  "5b2cfd9013a5fd8df47fa6b48b1e045f39816240c0b8fee8beadf4d8e9c0566c63a3e6e257f87329b18fae980078116",
	"144e4211384586c16bd3ad4afa99cc9170df3560e77982d0db45f3536814f0bd5871c1908bd478cd1ee605167ff82995");
	
	Fq2_set_hex_str(&g21, "5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffeffff", "0");
	Fq2_set_hex_str(&g22, "5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe", "0");
	Fq2_set_hex_str(&g23, "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaaa", "0");
	Fq2_set_hex_str(&g24, "1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac", "0");
	Fq2_set_hex_str(&g25, "1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad", "0");
	
	//Mont representation
	Fq2_mont_rep(&g11, &g11);
	Fq2_mont_rep(&g12, &g12);
	Fq2_mont_rep(&g13, &g13);
	Fq2_mont_rep(&g14, &g14);
	Fq2_mont_rep(&g15, &g15);
	
	Fq2_mont_rep(&g21, &g21);
	Fq2_mont_rep(&g22, &g22);
	Fq2_mont_rep(&g23, &g23);
	Fq2_mont_rep(&g24, &g24);
	Fq2_mont_rep(&g25, &g25);
	
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

void ate_exp_q(Fq12* rop, Fq12* a) //rop = a^q
{
	/*
	
	f^q = ~g0 + ~h0 * g11 * w + ~g1 * g12 *  w^2 + ~h1 * g13 * w^3 + ~g2 * g14 * w^4 + ~h2 * g15 * w^5
	
	*/
	
	Fq2 g0, g1, g2, h0, h1, h2;
	Fq2_assign(&g0, &(a->val0.val0));
	Fq2_assign(&g1, &(a->val0.val1));
	Fq2_assign(&g2, &(a->val0.val2));
	Fq2_assign(&h0, &(a->val1.val0));
	Fq2_assign(&h1, &(a->val1.val1));
	Fq2_assign(&h2, &(a->val1.val2));
	
	Fq2_conj(&g0, &g0);
	Fq2_conj(&g1, &g1);
	Fq2_conj(&g2, &g2);
	Fq2_conj(&h0, &h0);
	Fq2_conj(&h1, &h1);
	Fq2_conj(&h2, &h2);
	
	Fq2_mont_mul(&g1, &g1, &g12);
	Fq2_mont_mul(&g2, &g2, &g14);
	Fq2_mont_mul(&h0, &h0, &g11);
	Fq2_mont_mul(&h1, &h1, &g13);
	Fq2_mont_mul(&h2, &h2, &g15);
	
	Fq2_assign(&(rop->val0.val0), &g0);
	Fq2_assign(&(rop->val0.val1), &g1);
	Fq2_assign(&(rop->val0.val2), &g2);
	Fq2_assign(&(rop->val1.val0), &h0);
	Fq2_assign(&(rop->val1.val1), &h1);
	Fq2_assign(&(rop->val1.val2), &h2);
	
}

void ate_exp_q2(Fq12* rop, Fq12* a) //rop = a^(q^2).
{
	/*
	
	f^(q^2) = f = g0 + h0 * g21 * w + g1 * g22 * w^2 + h1 * g23 * w^3 + g2 * g24 * w^4 + h2 * g25 * w^5.
	
	*/
	Fq2 g0, g1, g2, h0, h1, h2;
	Fq2_assign(&g0, &(a->val0.val0));
	Fq2_assign(&g1, &(a->val0.val1));
	Fq2_assign(&g2, &(a->val0.val2));
	Fq2_assign(&h0, &(a->val1.val0));
	Fq2_assign(&h1, &(a->val1.val1));
	Fq2_assign(&h2, &(a->val1.val2));
	
	Fq2_mont_mul(&g1, &g1, &g22);
	Fq2_mont_mul(&g2, &g2, &g24);
	Fq2_mont_mul(&h0, &h0, &g21);
	Fq2_mont_mul(&h1, &h1, &g23);
	Fq2_mont_mul(&h2, &h2, &g25);
	
	Fq2_assign(&(rop->val0.val0), &g0);
	Fq2_assign(&(rop->val0.val1), &g1);
	Fq2_assign(&(rop->val0.val2), &g2);
	Fq2_assign(&(rop->val1.val0), &h0);
	Fq2_assign(&(rop->val1.val1), &h1);
	Fq2_assign(&(rop->val1.val2), &h2);
}

void ate_exp_q6(Fq12* rop, Fq12* a) //rop = a^(q^6)
{
	/*
	
	f = g + hw
	
	f^(q^6) = g - hw
	
	*/
	
	Fq6 g, h;
	Fq6_assign(&g, &(a->val0));
	Fq6_assign(&h, &(a->val1));
	
	Fq6_add_inv(&h, &h);
	
	Fq6_assign(&(rop->val0), &g);
	Fq6_assign(&(rop->val1), &h);
}

void ate_exp_x(Fq12* rop, Fq12* a) // rop = a ^ x
{
	Fq12 aux;
	Fq12_assign(&aux, a);
	
	for(int i = 1; i < len_bX_STR; i++)
	{
		Fq12_mont_mul(&aux, &aux, &aux);
		if(bX_STR[i] == '1')
		{
			Fq12_mont_mul(&aux, &aux, a);
		}
	}
	
	Fq12_assign(rop, &aux);
}

void ate_exp_xp1(Fq12* rop, Fq12* a) // rop = a ^ (x+1)
{
	Fq12 aux;
	Fq12_assign(&aux, a);
	
	for(int i = 1; i < len_bXP1_STR; i++)
	{
		Fq12_mont_mul(&aux, &aux, &aux);
		if(bXP1_STR[i] == '1')
		{
			Fq12_mont_mul(&aux, &aux, a);
		}
	}
	
	Fq12_assign(rop, &aux);
}

void ate_exp_xp1div3(Fq12* rop, Fq12* a) // rop = a ^ (x+1)/3
{
	Fq12 aux;
	Fq12_assign(&aux, a);
	
	for(int i = 1; i < len_bXP1_DIV3_STR; i++)
	{
		Fq12_mont_mul(&aux, &aux, &aux);
		if(bXP1_DIV3_STR[i] == '1')
		{
			Fq12_mont_mul(&aux, &aux, a);
		}
	}
	
	Fq12_assign(rop, &aux);
}

void ate_exp_hard(Fq12* rop, Fq12* f) //rop = f^((q^4 - q^2 + 1)/r)
{
	/*
	
	Note the following identity.
	
	(q^4 - q^2 + 1)/r = (x+1)^2 /3 * (q-x)(x^2 + q^2 - 1) + 1, where x = 0xd201000000010000 (x is taken positive here).
	
	Also, f^(-1) = ~f, because of the consequent raising to the power (q^6-1).
	
	*/
	
	Fq12 a, b1, b2, b, c1, c2, c3, c; //auxiliary variables
	
	// I. powering to (x+1)^2/3
	
	ate_exp_xp1div3(&a, f);
	ate_exp_xp1(&a, &a);
	
	// II. powering to q - x
	
	ate_exp_q(&b1, &a);
	
	Fq12_conj(&b2, &a);
	ate_exp_x(&b2, &b2); //b2 = a^(-x).
	Fq12_mont_mul(&b, &b1, &b2);
	
	// III. powering to x^2 + q^2 - 1
	
	ate_exp_x(&c1, &b);
	ate_exp_x(&c1, &c1);
	
	ate_exp_q2(&c2, &b);
	
	Fq12_conj(&c3, &b);
	
	Fq12_mont_mul(&c, &c1, &c2);
	Fq12_mont_mul(&c, &c, &c3);
	
	// IV. multiplying by f and returning the result
	
	Fq12_mont_mul(rop, &c, f);
}

void ate_exp(Fq12* rop, Fq12* a)
{	
	Fq12 aux1, aux2, inv2;
	
	// Computing a^(q^2 + 1)
	
	Fq12_assign(&aux1, a);
	
	ate_exp_q2(&aux1, &aux1);
	Fq12_mont_mul(&aux1, &aux1, a);
	
	// Computing (a^(q^2 + 1)) ^ (q^6 - 1)
	
	Fq12_assign(&aux2, &aux1);
	Fq12_mont_mul_inv(&inv2, &aux2);
	
	ate_exp_q6(&aux2, &aux2);
	Fq12_mont_mul(&aux2, &aux2, &inv2); 
	
	// Powering to the (q^4 - q^2 + 1)/r
	
	Fq12 result;
	Fq12_assign(&result, &aux2);
	ate_exp_hard(rop, &result);
}

void ate_pairing(Fq12* rop, G1* P, G2_Z* QZ)
{
	Fq12 temp;

	ate(&temp, P, QZ);
	ate_exp(&temp, &temp);

	Fq12_assign(rop, &temp);
}

