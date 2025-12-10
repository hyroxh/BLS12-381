#include "g1.h"
#include "../fq/fq.h"
#include "../fr/fr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Fr orderG1;
G1 genG1;
Fr cofactorG1;
char* bORD_STR;
int len_bORD_STR;

// Group init
void G1_group_init()
{

	// auxiliary
	Fq xgen, ygen;

	Fq_set_hex_str(&xgen, "17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a"
					  "3f171bac586c55e83ff97a1aeffb3af00adb22c6bb");
	Fq_set_hex_str(&ygen, "8b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18"
					  "cb2c04b3edd03cc744a2888ae40caa232946c5e7e1");

	// Defining
	Fr_set_hex_str(
		&orderG1,
		"73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001");
	G1_set(&genG1, &xgen, &ygen);
	Fr_set_hex_str(&cofactorG1, "396c8c005555e1568c00aaab0000aaab");

	bORD_STR = "1110011111011011010011101010011001010011001110101111101010010000011001100111001110110000000"
		"1000000010011010000111011000000001010101001110111101101001000000001011111111111111100101101111"
		"1111101111111111111111111111111111111100000000000000000000000000000001";
	len_bORD_STR = 255;
}


// Init

void G1_set(G1* a, Fq* x, Fq* y)
{
	Fq_assign(&(a->valx), x);
	Fq_assign(&(a->valy), y);

	a->inf = 0;
}

void G1_set_inf(G1* a)
{
	Fq_set_ui(&(a->valx), 0);
	Fq_set_ui(&(a->valy), 0);

	a->inf = 1;
}

void G1_assign(G1* rop, G1* a)
{
	Fq_assign(&(rop->valx), &(a->valx));
	Fq_assign(&(rop->valy), &(a->valy));

	rop->inf = a->inf;
}

// Init projective

void G1_Z_set(G1_Z* a, Fq* X, Fq* Y, Fq* Z)
{
	Fq_assign(&(a->valX), X);
	Fq_assign(&(a->valY), Y);
	Fq_assign(&(a->valZ), Z);
}
void G1_Z_set_inf(G1_Z* a)
{
	Fq_set_ui(&(a->valX), 0);
	Fq_set_ui(&(a->valY), 1);
	Fq_set_ui(&(a->valZ), 0);
}
void G1_Z_assign(G1_Z* rop, G1_Z* a)
{
	Fq_assign(&(rop->valX), &(a->valX));
	Fq_assign(&(rop->valY), &(a->valY));
	Fq_assign(&(rop->valZ), &(a->valZ));
}

// Maps

void G1_Z_to_group(G1_Z* rop, G1_Z* a)
{
	if(G1_Z_curve_check(a) != 1)
	{
		errno = EDOM;
		perror("Domain error occured");
		exit(EXIT_FAILURE);
	}
	
	G1_Z_mul(rop, &cofactorG1, a);

}

void G1_to_affine(G1* rop, G1_Z* a)
{
	Fq zero;
	Fq_set_zero(&zero);
	
	uint32_t mask = -Fq_cmp(&(a->valZ), &zero); //0xFFFFFFFF, if aZ == 0; 0, otherwise
	
	G1 res1, res2;
	
	G1_set_inf(&res1);
	
	// auxiliary init
	Fq tempZ, tempZinv;

	// auxiliary comp
	Fq_assign(&tempZ, &(a->valZ));
	Fq_mont_mul_inv(&tempZinv, &tempZ);

	// computation
	Fq_mont_mul(&(res2.valx), &(a->valX), &tempZinv);
	Fq_mont_mul(&(res2.valy), &(a->valY), &tempZinv);
	
	res2.inf = 0;
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->valx.limbs[i] = (res1.valx.limbs[i] & mask) | (res2.valx.limbs[i] & ~mask);
		rop->valy.limbs[i] = (res1.valy.limbs[i] & mask) | (res2.valy.limbs[i] & ~mask);
		rop->inf = (res1.inf & mask) | (res2.inf & ~mask);
	}

}
void G1_to_project(G1_Z* rop, G1* a)
{
	
	uint32_t mask = -(a->inf == 1); //0xFFFFFFFF if a=O; 0, otherwise
	
	G1_Z inf;
	G1_Z_set_inf(&inf);
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->valX.limbs[i] = ( a->valx.limbs[i] & ~mask );
		rop->valY.limbs[i] = ( a->valy.limbs[i] & ~mask ) | ( inf.valY.limbs[i] & mask );
		rop->valZ.limbs[i] = 0;
	}
	
	rop->valZ.limbs[0] = (1 & ~mask) | (inf.valZ.limbs[0] & mask);	
}

// Curve check projective
int G1_Z_curve_check(G1_Z* a)
{
	// E(Fq): Y^2 * Z = X^3 + 4*Z^3
	Fq zero;
	Fq_set_ui(&zero, 0);
	
	if(Fq_cmp(&(a->valZ), &zero) == 1) //inf
	{
		return 1;
	}
	
	int check = 0; 
	
	//init
	Fq tempX, tempY, tempZ, aux;
	
	Fq_assign(&tempX, &(a->valX));
	Fq_assign(&tempY, &(a->valY));
	Fq_assign(&tempZ, &(a->valZ));
	
	//computations
	Fq_mont_mul(&tempY, &tempY, &tempY);
	Fq_mont_mul(&tempY, &tempY, &(a->valZ));
	
	Fq_mont_mul(&tempX, &tempX, &tempX);
	Fq_mont_mul(&tempX, &tempX, &(a->valX));
	
	Fq_mont_mul(&tempZ, &tempZ, &tempZ);
	Fq_mont_mul(&tempZ, &tempZ, &(a->valZ));
	Fq_add(&tempZ, &tempZ, &tempZ);
	Fq_add(&tempZ, &tempZ, &tempZ);
	
	Fq_add(&aux, &tempX, &tempZ);
	
	if(Fq_cmp(&tempY, &aux) == 0)
	{
		check = 1;
	}
	
	return check;
	
}

int G1_Z_group_check(G1_Z* a)
{
	// If a is in G1, then [r]a = 0.
	// In the projective coordinates, Za = 0

	int inGroup = 1;

	Fq zero;
	Fq_set_ui(&zero, 0);

	G1_Z temp;

	G1_Z_mul(&temp, &R_MODULUS, a);
	
	uint32_t mask = -Fq_cmp(&(temp.valZ), &zero);
	
	inGroup = (1 & mask);
	
	return inGroup;
}


// Montgomery representation

void G1_mont_rep(G1* aR, G1* a)
{
	Fq_mont_rep(&(aR->valx), &(a->valx));
	Fq_mont_rep(&(aR->valy), &(a->valy));
}

void G1_mont_rep_inv(G1* a, G1* aR)
{
	Fq_mont_rep_inv(&(a->valx), &(aR->valx));
	Fq_mont_rep_inv(&(a->valy), &(aR->valy));	
}

void G1_Z_mont_rep(G1_Z* aR, G1_Z* a)
{
	Fq_mont_rep(&(aR->valX), &(a->valX));
	Fq_mont_rep(&(aR->valY), &(a->valY));
	Fq_mont_rep(&(aR->valZ), &(a->valZ));
}

void G1_Z_mont_rep_inv(G1_Z* a, G1_Z* aR)
{
	Fq_mont_rep_inv(&(a->valX), &(aR->valX));
	Fq_mont_rep_inv(&(a->valY), &(aR->valY));
	Fq_mont_rep_inv(&(a->valZ), &(aR->valZ));	
}

//Correction
void G1_Z_correction(G1_Z* a) //a -> (0:1:0), if a is inf; nothing changes, otherwise
{
	Fq zero;
	Fq_set_zero(&zero);
	
	uint32_t mask = -Fq_cmp(&(a->valZ), &zero); //0xFFFFFFFF, if a is inf; 0, otherwise
	
	for(int i = 0; i < N_LIMBS; i++)
	{	
		a->valX.limbs[i] = (a->valX.limbs[i] & ~mask);
		a->valY.limbs[i] = (a->valY.limbs[i] & ~mask);
		a->valZ.limbs[i] = (a->valZ.limbs[i] & ~mask);
	}
	
	a->valY.limbs[0] = (a->valY.limbs[0] & ~mask) | (1 & mask);
}

// Group Arithmetic Projective
void G1_Z_add_inv(G1_Z* rop, G1_Z* a)
{
	// aux init
	Fq tempY;

	Fq_add_inv(&tempY, &(a->valY));

	Fq_assign(&(rop->valX), &(a->valX));
	Fq_assign(&(rop->valY), &tempY);
	Fq_assign(&(rop->valZ), &(a->valZ));
}

uint32_t G1_Z_cmp(G1_Z* p, G1_Z* q) //returns 1, if P = Q. 0, otherwise
{
	//(Xp : Yp : Zp) ~ (Xq : Yq : Zq) <=>
	// <=> [XpYq - XqYp, YpZq - YqZp, ZpXq - ZqXp] = [0, 0, 0]
	
	uint32_t mask = 1;

	// aux init
	Fq tempXpYq, tempXqYp, tempYpZq, tempYqZp, tempZpXq, tempZqXp;

	// computation
	Fq_mont_mul(&tempXpYq, &(p->valX), &(q->valY));
	Fq_mont_mul(&tempXqYp, &(p->valY), &(q->valX));
	Fq_mont_mul(&tempYpZq, &(p->valY), &(q->valZ));
	Fq_mont_mul(&tempYqZp, &(p->valZ), &(q->valY));
	Fq_mont_mul(&tempZpXq, &(p->valZ), &(q->valX));
	Fq_mont_mul(&tempZqXp, &(p->valX), &(q->valZ));
	
	mask = mask & Fq_cmp(&tempXpYq, &tempXqYp);
	mask = mask & Fq_cmp(&tempYpZq, &tempYqZp);
	mask = mask & Fq_cmp(&tempZpXq, &tempZqXp);

	return mask;
}

void G1_Z_adding(G1_Z* rop, G1_Z* p, G1_Z* q) //if p!=q
{
	/*

	Zr = ZpZq * (XpZq - XqZp)^3.

	Xr = (XpZq - XqZp) * [ZpZq * (YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq +
	XqZp)]

	Yr = ZpZq * (XqYp - XpYq) * (XpZq - XqZp)^2 - (YpZq - YqZp) * [ZpZq *
	(YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq + XqZp)]

	*/

	// aux init
	Fq tempXpYq, tempXqYp, tempYpZq, tempYqZp, tempXqZp, tempXpZq, tempZpZq,
		diff1, diff2, diff3, sum, delta1, delta2, aux1, aux2, bracket, c1,
		c2, c3;

	// aux computation
	Fq_mont_mul(&tempXpYq, &(p->valX), &(q->valY));
	Fq_mont_mul(&tempXqYp, &(p->valY), &(q->valX));
	Fq_mont_mul(&tempYpZq, &(p->valY), &(q->valZ));
	Fq_mont_mul(&tempYqZp, &(p->valZ), &(q->valY));
	Fq_mont_mul(&tempXqZp, &(p->valZ), &(q->valX));
	Fq_mont_mul(&tempXpZq, &(p->valX), &(q->valZ));
	Fq_mont_mul(&tempZpZq, &(p->valZ), &(q->valZ));

	Fq_sub(&diff1, &tempXpZq, &tempXqZp);
	Fq_mont_mul(&diff2, &diff1, &diff1);
	Fq_mont_mul(&diff3, &diff2, &diff1);

	Fq_add(&sum, &tempXpZq, &tempXqZp);
	Fq_sub(&delta1, &tempYpZq, &tempYqZp);
	Fq_mont_mul(&delta2, &delta1, &delta1);

	Fq_mont_mul(&aux1, &tempZpZq, &delta2);
	Fq_mont_mul(&aux2, &diff2, &sum);
	Fq_sub(
		&bracket, &aux1,
		&aux2); //[ZpZq * (YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq + XqZp)]

	// Zr
	Fq_mont_mul(&(rop->valZ), &tempZpZq, &diff3);

	// Xr
	Fq_mont_mul(&(rop->valX), &diff1, &bracket);

	// Yr
	Fq_sub(&c3, &tempXqYp, &tempXpYq);
	Fq_mont_mul(&c1, &tempZpZq, &c3);
	Fq_mont_mul(&c1, &c1, &diff2);
	Fq_mont_mul(&c2, &delta1, &bracket);

	Fq_sub(&(rop->valY), &c1, &c2);
	
	G1_Z_correction(rop);
}

void G1_Z_doubling(G1_Z* rop, G1_Z* p)
{
	/*

	Zr = 8Yp^3Zp^3

	Xr = 18Xp^4YpZp - 16XpYp^3Zp^2

	Yr = -27Xp^6 + 36Xp^3Yp^2Zp - 8Yp^4Zp^2

	*/

	// aux init
	Fq c8, c16, c18, c27, c36, aux1, aux2, aux3, storeX, storeY, storeZ;

	// constants setting
	Fq_set_ui(&c8, 8);
	Fq_set_ui(&c16, 16);
	Fq_set_ui(&c18, 18);
	Fq_set_ui(&c27, 27);
	Fq_set_ui(&c36, 36);
	
	//constants to Montgomery rep
	Fq_mont_rep(&c8, &c8);
	Fq_mont_rep(&c16, &c16);
	Fq_mont_rep(&c18, &c18);
	Fq_mont_rep(&c27, &c27);
	Fq_mont_rep(&c36, &c36);

	//---NB! aux1,2,3 are fixed within the computation of the current
	// coordinate only
	// Zr
	Fq_mont_mul(&aux1, &(p->valY), &(p->valY));
	Fq_mont_mul(&aux1, &aux1, &(p->valY));
	Fq_mont_mul(&aux1, &aux1, &(p->valZ));
	Fq_mont_mul(&aux1, &aux1, &(p->valZ));
	Fq_mont_mul(&aux1, &aux1, &(p->valZ));
	Fq_mont_mul(&storeZ, &c8, &aux1);

	// Xr
	Fq_mont_mul(&aux1, &(p->valX), &(p->valX));
	Fq_mont_mul(&aux1, &aux1, &(p->valX));
	Fq_mont_mul(&aux1, &aux1, &(p->valX));
	Fq_mont_mul(&aux1, &aux1, &(p->valY));
	Fq_mont_mul(&aux1, &aux1, &(p->valZ));
	Fq_mont_mul(&aux1, &c18, &aux1);

	Fq_mont_mul(&aux2, &(p->valX), &(p->valY));
	Fq_mont_mul(&aux2, &aux2, &(p->valY));
	Fq_mont_mul(&aux2, &aux2, &(p->valY));
	Fq_mont_mul(&aux2, &aux2, &(p->valZ));
	Fq_mont_mul(&aux2, &aux2, &(p->valZ));
	Fq_mont_mul(&aux2, &c16, &aux2);

	Fq_sub(&storeX, &aux1, &aux2);

	// Yr
	Fq_mont_mul(&aux1, &(p->valX), &(p->valX));
	Fq_mont_mul(&aux1, &aux1, &(p->valX));
	Fq_mont_mul(&aux1, &aux1, &(p->valX));
	Fq_mont_mul(&aux1, &aux1, &(p->valX));
	Fq_mont_mul(&aux1, &aux1, &(p->valX));
	Fq_mont_mul(&aux1, &c27, &aux1);

	Fq_mont_mul(&aux2, &(p->valX), &(p->valX));
	Fq_mont_mul(&aux2, &aux2, &(p->valX));
	Fq_mont_mul(&aux2, &aux2, &(p->valY));
	Fq_mont_mul(&aux2, &aux2, &(p->valY));
	Fq_mont_mul(&aux2, &aux2, &(p->valZ));
	Fq_mont_mul(&aux2, &c36, &aux2);

	Fq_mont_mul(&aux3, &(p->valY), &(p->valY));
	Fq_mont_mul(&aux3, &aux3, &(p->valY));
	Fq_mont_mul(&aux3, &aux3, &(p->valY));
	Fq_mont_mul(&aux3, &aux3, &(p->valZ));
	Fq_mont_mul(&aux3, &aux3, &(p->valZ));
	Fq_mont_mul(&aux3, &c8, &aux3);

	Fq_sub(&aux2, &aux2, &aux1);
	Fq_sub(&storeY, &aux2, &aux3);

	// storing the result
	Fq_assign(&(rop->valX), &storeX);
	Fq_assign(&(rop->valY), &storeY);
	Fq_assign(&(rop->valZ), &storeZ);
	
	G1_Z_correction(rop);
}


void G1_Z_add(G1_Z* rop, G1_Z* p, G1_Z* q)
{
	
	Fq zero;
	Fq_set_zero(&zero);
	
	uint32_t mask_inf1 = -Fq_cmp(&(p->valZ), &zero); //0xFFFFFFFF if p=O; 0, otherwise
	uint32_t mask_inf2 = -Fq_cmp(&(q->valZ), &zero); //0xFFFFFFFF if q=O; 0, otherwise
	uint32_t mask_inf = ~mask_inf1 & ~mask_inf2; //0xFFFFFFFF, if p!=O & q!=O; 0, otherwise
	
	G1_Z res1, res2;
	G1_Z_adding(&res1, p, q);
	G1_Z_doubling(&res2, p);
	
	uint32_t mask = -G1_Z_cmp(p, q); //0xFFFFFFFF if p==q; 0, otherwise
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->valX.limbs[i] = ( ((res1.valX.limbs[i] & ~mask) | (res2.valX.limbs[i] & mask)) & mask_inf ) | ( p->valX.limbs[i] & mask_inf2 ) | ( q->valX.limbs[i] & mask_inf1 );
		rop->valY.limbs[i] = ( ((res1.valY.limbs[i] & ~mask) | (res2.valY.limbs[i] & mask)) & mask_inf ) | ( p->valY.limbs[i] & mask_inf2 ) | ( q->valY.limbs[i] & mask_inf1 );
		rop->valZ.limbs[i] = ( ((res1.valZ.limbs[i] & ~mask) | (res2.valZ.limbs[i] & mask)) & mask_inf ) | ( p->valZ.limbs[i] & mask_inf2 ) | ( q->valZ.limbs[i] & mask_inf1 );
	}
	
}

void G1_Z_mul(G1_Z* rop, Fr* n, G1_Z* p) // rop = [n]P, Montgomery ladder
{
	G1_Z R0, R1, T_swap;
	G1_Z_set_inf(&R0);
	G1_Z_assign(&R1, p);
	
	unsigned int binstr[4*LEN_MAX_R];
	
	for(int i = 0; i < FR_LIMBS; i++)
	{
		for(int j = 0; j < 32; j++)
		{
			binstr[32 * i + j ] = (n->limbs[i] >> j) & 1;
		}
	}
	
	for(int i = 4*LEN_MAX_R - 1; i >=0; i--)
	{		
		uint32_t mask = -(binstr[i]); //0xFFFFFFFF, if 1; 0, otherwise
		
		//swapping
		for(int j = 0; j < N_LIMBS; j++)
		{
			T_swap.valX.limbs[j] = (R0.valX.limbs[j] ^ R1.valX.limbs[j]) & mask;
			T_swap.valY.limbs[j] = (R0.valY.limbs[j] ^ R1.valY.limbs[j]) & mask;
			T_swap.valZ.limbs[j] = (R0.valZ.limbs[j] ^ R1.valZ.limbs[j]) & mask;
			
			R0.valX.limbs[j] ^= T_swap.valX.limbs[j];
			R0.valY.limbs[j] ^= T_swap.valY.limbs[j];
			R0.valZ.limbs[j] ^= T_swap.valZ.limbs[j];
			
			R1.valX.limbs[j] ^= T_swap.valX.limbs[j];
			R1.valY.limbs[j] ^= T_swap.valY.limbs[j];
			R1.valZ.limbs[j] ^= T_swap.valZ.limbs[j];
		}
		
		G1_Z_add(&R1, &R0, &R1);
		G1_Z_add(&R0, &R0, &R0);
		
		//swapping back
		for(int j = 0; j < N_LIMBS; j++)
		{
			T_swap.valX.limbs[j] = (R0.valX.limbs[j] ^ R1.valX.limbs[j]) & mask;
			T_swap.valY.limbs[j] = (R0.valY.limbs[j] ^ R1.valY.limbs[j]) & mask;
			T_swap.valZ.limbs[j] = (R0.valZ.limbs[j] ^ R1.valZ.limbs[j]) & mask;
			
			R0.valX.limbs[j] ^= T_swap.valX.limbs[j];
			R0.valY.limbs[j] ^= T_swap.valY.limbs[j];
			R0.valZ.limbs[j] ^= T_swap.valZ.limbs[j];
			
			R1.valX.limbs[j] ^= T_swap.valX.limbs[j];
			R1.valY.limbs[j] ^= T_swap.valY.limbs[j];
			R1.valZ.limbs[j] ^= T_swap.valZ.limbs[j];
		} 
		
	}
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->valX.limbs[i] = R0.valX.limbs[i];
		rop->valY.limbs[i] = R0.valY.limbs[i];
		rop->valZ.limbs[i] = R0.valZ.limbs[i];
	}
}

// Print
void G1_print(G1* a)
{
	if (a->inf == 1)
	{
		printf("inf");
	}

	if (a->inf == 0)
	{
		printf("(");
		Fq_print(&(a->valx));
		printf(", ");
		Fq_print(&(a->valy));
		printf(")");
	}
}

void G1_Z_print(G1_Z* a)
{
	printf("(");
	Fq_print(&(a->valX));
	printf(" : ");
	Fq_print(&(a->valY));
	printf(" : ");
	Fq_print(&(a->valZ));
	printf(")");
}

void G1_print_dec(G1* a)
{
	if (a->inf == 1)
	{
		printf("inf");
	}

	if (a->inf == 0)
	{
		printf("(");
		Fq_print_dec(&(a->valx));
		printf(", ");
		Fq_print_dec(&(a->valy));
		printf(")");
	}
}

void G1_Z_print_dec(G1_Z* a)
{
	printf("(");
	Fq_print_dec(&(a->valX));
	printf(" : ");
	Fq_print_dec(&(a->valY));
	printf(" : ");
	Fq_print_dec(&(a->valZ));
	printf(")");
}
