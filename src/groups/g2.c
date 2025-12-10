#include "g2.h"
#include "../fr/fr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Fr orderG2;
G2 genG2;
Fr cofactorG2;

// Group init
void G2_group_init()
{
	// auxiliary
	Fq2 xgen, ygen;

	Fq2_set_hex_str(&xgen,
				"24aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3"
				"d1770bac0326a805bbefd48056c8c121bdb8",
				"13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f"
				"5049334cf11213945d57e5ac7d055d042b7e");
	Fq2_set_hex_str(&ygen,
				"ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160"
				"d12c923ac9cc3baca289e193548608b82801",
				"606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e"
				"99ab3f370d275cec1da1aaa9075ff05f79be");

	// Defining
	Fr_set_hex_str(
		&orderG2,
		"73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001");
	G2_set(&genG2, &xgen, &ygen);
	Fr_set_hex_str(
		&cofactorG2,
		"5d543a95414e7f1091d50792876a202cd91de4547085abaa68a205b2e5a7ddfa628f"
		"1cb4d9e82ef21537e293a6691ae1616ec6e786f0c70cf1c38e31c7238e5");
}

// Init

void G2_set(G2* a, Fq2* x, Fq2* y)
{
	Fq2_assign(&(a->valx), x);
	Fq2_assign(&(a->valy), y);

	a->inf = 0;
}

void G2_set_inf(G2* a)
{
	Fq2_set_ui(&(a->valx), 0, 0);
	Fq2_set_ui(&(a->valy), 0, 0);

	a->inf = 1;
}

void G2_assign(G2* rop, G2* a)
{
	Fq2_assign(&(rop->valx), &(a->valx));
	Fq2_assign(&(rop->valy), &(a->valy));

	rop->inf = a->inf;
}

// Init projective

void G2_Z_set(G2_Z* a, Fq2* X, Fq2* Y, Fq2* Z)
{
	Fq2_assign(&(a->valX), X);
	Fq2_assign(&(a->valY), Y);
	Fq2_assign(&(a->valZ), Z);
}
void G2_Z_set_inf(G2_Z* a)
{
	Fq2_set_ui(&(a->valX), 0, 0);
	Fq2_set_ui(&(a->valY), 1, 0);
	Fq2_set_ui(&(a->valZ), 0, 0);
}
void G2_Z_assign(G2_Z* rop, G2_Z* a)
{
	Fq2_assign(&(rop->valX), &(a->valX));
	Fq2_assign(&(rop->valY), &(a->valY));
	Fq2_assign(&(rop->valZ), &(a->valZ));
}

// Maps

void G2_Z_to_group(G2_Z* rop, G2_Z* a)
{
	if(G2_Z_curve_check(a) != 1)
	{
		errno = EDOM;
		perror("Domain error occured");
		exit(EXIT_FAILURE);
	}
	
	G2_Z_mul(rop, &cofactorG2, a);

}

void G2_to_affine(G2* rop, G2_Z* a)
{
	Fq2 zero;
	Fq2_set_zero(&zero);
	
	uint32_t mask = -Fq2_cmp(&(a->valZ), &zero); //0xFFFFFFFF, if aZ == 0; 0, otherwise
	
	G2 res1, res2;
	
	G2_set_inf(&res1);
	
	// auxiliary init
	Fq2 tempZ, tempZinv;

	// auxiliary comp
	Fq2_assign(&tempZ, &(a->valZ));
	Fq2_mont_mul_inv(&tempZinv, &tempZ);

	// computation
	Fq2_mont_mul(&(res2.valx), &(a->valX), &tempZinv);
	Fq2_mont_mul(&(res2.valy), &(a->valY), &tempZinv);
	
	res2.inf = 0;
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->valx.val0.limbs[i] = (res1.valx.val0.limbs[i] & mask) | (res2.valx.val0.limbs[i] & ~mask);
		rop->valx.val1.limbs[i] = (res1.valx.val1.limbs[i] & mask) | (res2.valx.val1.limbs[i] & ~mask);
		
		rop->valy.val0.limbs[i] = (res1.valy.val0.limbs[i] & mask) | (res2.valy.val0.limbs[i] & ~mask);
		rop->valy.val1.limbs[i] = (res1.valy.val1.limbs[i] & mask) | (res2.valy.val1.limbs[i] & ~mask);
		
		rop->inf = (res1.inf & mask) | (res2.inf & ~mask);
	}

}
void G2_to_project(G2_Z* rop, G2* a)
{
	
	uint32_t mask = -(a->inf == 1); //0xFFFFFFFF if a=O; 0, otherwise
	
	G2_Z inf;
	G2_Z_set_inf(&inf);
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->valX.val0.limbs[i] = ( a->valx.val0.limbs[i] & ~mask );
		rop->valX.val1.limbs[i] = ( a->valx.val1.limbs[i] & ~mask );
		
		rop->valY.val0.limbs[i] = ( a->valy.val0.limbs[i] & ~mask ) | ( inf.valY.val0.limbs[i] & mask );
		rop->valY.val1.limbs[i] = ( a->valy.val1.limbs[i] & ~mask ) | ( inf.valY.val1.limbs[i] & mask );
		
		rop->valZ.val0.limbs[i] = 0;
		rop->valZ.val1.limbs[i] = 0;
	}
	
	rop->valZ.val0.limbs[0] = (1 & ~mask) | (inf.valZ.val0.limbs[0] & mask);	
}

// Curve check projective
int G2_Z_curve_check(G2_Z* a)
{
	// E(Fq2): Y^2 * Z = X^3 + (4 + 4u)*Z^3
	Fq2 zero;
	Fq2_set_zero(&zero);
	
	if(Fq2_cmp(&(a->valZ), &zero) == 1) //inf
	{
		return 1;
	}
	
	int check = 0; 
	
	//init
	Fq2 tempX, tempY, tempZ, aux;
	
	Fq2_assign(&tempX, &(a->valX));
	Fq2_assign(&tempY, &(a->valY));
	Fq2_assign(&tempZ, &(a->valZ));
	
	//computations
	Fq2_mont_mul(&tempY, &tempY, &tempY);
	Fq2_mont_mul(&tempY, &tempY, &(a->valZ));
	
	Fq2_mont_mul(&tempX, &tempX, &tempX);
	Fq2_mont_mul(&tempX, &tempX, &(a->valX));
	
	Fq2_mont_mul(&tempZ, &tempZ, &tempZ);
	Fq2_mont_mul(&tempZ, &tempZ, &(a->valZ));
	Fq2 c44;
	Fq2_set_ui(&c44, 4, 4);
	Fq2_mont_rep(&c44, &c44);
	Fq2_mont_mul(&tempZ, &tempZ, &c44);
	
	Fq2_add(&aux, &tempX, &tempZ);
	
	if(Fq2_cmp(&tempY, &aux) == 0)
	{
		check = 1;
	}
	
	return check;
	
}

int G2_Z_group_check(G2_Z* a)
{
	// If a is in G2, then [r]a = 0.
	// In the projective coordinates, Za = 0

	int inGroup = 1;

	Fq2 zero;
	Fq2_set_zero(&zero);

	G2_Z temp;

	G2_Z_mul(&temp, &R_MODULUS, a);
	
	uint32_t mask = -Fq2_cmp(&(temp.valZ), &zero);
	
	inGroup = (1 & mask);
	
	return inGroup;
}


// Montgomery representation

void G2_mont_rep(G2* aR, G2* a)
{
	Fq2_mont_rep(&(aR->valx), &(a->valx));
	Fq2_mont_rep(&(aR->valy), &(a->valy));
}

void G2_mont_rep_inv(G2* a, G2* aR)
{
	Fq2_mont_rep_inv(&(a->valx), &(aR->valx));
	Fq2_mont_rep_inv(&(a->valy), &(aR->valy));	
}

void G2_Z_mont_rep(G2_Z* aR, G2_Z* a)
{
	Fq2_mont_rep(&(aR->valX), &(a->valX));
	Fq2_mont_rep(&(aR->valY), &(a->valY));
	Fq2_mont_rep(&(aR->valZ), &(a->valZ));
}

void G2_Z_mont_rep_inv(G2_Z* a, G2_Z* aR)
{
	Fq2_mont_rep_inv(&(a->valX), &(aR->valX));
	Fq2_mont_rep_inv(&(a->valY), &(aR->valY));
	Fq2_mont_rep_inv(&(a->valZ), &(aR->valZ));	
}

//Correction
void G2_Z_correction(G2_Z* a) //a -> (0:1:0), if a is inf; nothing changes, otherwise
{
	Fq2 zero;
	Fq2_set_zero(&zero);
	
	uint32_t mask = -Fq2_cmp(&(a->valZ), &zero); //0xFFFFFFFF, if a is inf; 0, otherwise
	
	for(int i = 0; i < N_LIMBS; i++)
	{	
		a->valX.val0.limbs[i] = (a->valX.val0.limbs[i] & ~mask);
		a->valX.val1.limbs[i] = (a->valX.val1.limbs[i] & ~mask);
		
		a->valY.val0.limbs[i] = (a->valY.val0.limbs[i] & ~mask);
		a->valY.val1.limbs[i] = (a->valY.val1.limbs[i] & ~mask);
		
		a->valZ.val0.limbs[i] = (a->valZ.val0.limbs[i] & ~mask);
		a->valZ.val1.limbs[i] = (a->valZ.val1.limbs[i] & ~mask);
	}
	
	a->valY.val0.limbs[0] = (a->valY.val0.limbs[0] & ~mask) | (1 & mask);
}

// Group Arithmetic Projective
void G2_Z_add_inv(G2_Z* rop, G2_Z* a)
{
	// aux init
	Fq2 tempY;

	Fq2_add_inv(&tempY, &(a->valY));

	Fq2_assign(&(rop->valX), &(a->valX));
	Fq2_assign(&(rop->valY), &tempY);
	Fq2_assign(&(rop->valZ), &(a->valZ));
}

uint32_t G2_Z_cmp(G2_Z* p, G2_Z* q) //returns 1, if P = Q. 0, otherwise
{
	//(Xp : Yp : Zp) ~ (Xq : Yq : Zq) <=>
	// <=> [XpYq - XqYp, YpZq - YqZp, ZpXq - ZqXp] = [0, 0, 0]
	
	uint32_t mask = 1;

	// aux init
	Fq2 tempXpYq, tempXqYp, tempYpZq, tempYqZp, tempZpXq, tempZqXp;

	// computation
	Fq2_mont_mul(&tempXpYq, &(p->valX), &(q->valY));
	Fq2_mont_mul(&tempXqYp, &(p->valY), &(q->valX));
	Fq2_mont_mul(&tempYpZq, &(p->valY), &(q->valZ));
	Fq2_mont_mul(&tempYqZp, &(p->valZ), &(q->valY));
	Fq2_mont_mul(&tempZpXq, &(p->valZ), &(q->valX));
	Fq2_mont_mul(&tempZqXp, &(p->valX), &(q->valZ));
	
	mask = mask & Fq2_cmp(&tempXpYq, &tempXqYp);
	mask = mask & Fq2_cmp(&tempYpZq, &tempYqZp);
	mask = mask & Fq2_cmp(&tempZpXq, &tempZqXp);

	return mask;
}

void G2_Z_adding(G2_Z* rop, G2_Z* p, G2_Z* q) //if p!=q
{
	/*

	Zr = ZpZq * (XpZq - XqZp)^3.

	Xr = (XpZq - XqZp) * [ZpZq * (YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq +
	XqZp)]

	Yr = ZpZq * (XqYp - XpYq) * (XpZq - XqZp)^2 - (YpZq - YqZp) * [ZpZq *
	(YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq + XqZp)]

	*/

	// aux init
	Fq2 tempXpYq, tempXqYp, tempYpZq, tempYqZp, tempXqZp, tempXpZq, tempZpZq,
		diff1, diff2, diff3, sum, delta1, delta2, aux1, aux2, bracket, c1,
		c2, c3;

	// aux computation
	Fq2_mont_mul(&tempXpYq, &(p->valX), &(q->valY));
	Fq2_mont_mul(&tempXqYp, &(p->valY), &(q->valX));
	Fq2_mont_mul(&tempYpZq, &(p->valY), &(q->valZ));
	Fq2_mont_mul(&tempYqZp, &(p->valZ), &(q->valY));
	Fq2_mont_mul(&tempXqZp, &(p->valZ), &(q->valX));
	Fq2_mont_mul(&tempXpZq, &(p->valX), &(q->valZ));
	Fq2_mont_mul(&tempZpZq, &(p->valZ), &(q->valZ));

	Fq2_sub(&diff1, &tempXpZq, &tempXqZp);
	Fq2_mont_mul(&diff2, &diff1, &diff1);
	Fq2_mont_mul(&diff3, &diff2, &diff1);

	Fq2_add(&sum, &tempXpZq, &tempXqZp);
	Fq2_sub(&delta1, &tempYpZq, &tempYqZp);
	Fq2_mont_mul(&delta2, &delta1, &delta1);

	Fq2_mont_mul(&aux1, &tempZpZq, &delta2);
	Fq2_mont_mul(&aux2, &diff2, &sum);
	Fq2_sub(
		&bracket, &aux1,
		&aux2); //[ZpZq * (YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq + XqZp)]

	// Zr
	Fq2_mont_mul(&(rop->valZ), &tempZpZq, &diff3);

	// Xr
	Fq2_mont_mul(&(rop->valX), &diff1, &bracket);

	// Yr
	Fq2_sub(&c3, &tempXqYp, &tempXpYq);
	Fq2_mont_mul(&c1, &tempZpZq, &c3);
	Fq2_mont_mul(&c1, &c1, &diff2);
	Fq2_mont_mul(&c2, &delta1, &bracket);

	Fq2_sub(&(rop->valY), &c1, &c2);
	
	G2_Z_correction(rop);
}

void G2_Z_doubling(G2_Z* rop, G2_Z* p)
{
	/*

	Zr = 8Yp^3Zp^3

	Xr = 18Xp^4YpZp - 16XpYp^3Zp^2

	Yr = -27Xp^6 + 36Xp^3Yp^2Zp - 8Yp^4Zp^2

	*/

	// aux init
	Fq2 c8, c16, c18, c27, c36, aux1, aux2, aux3, storeX, storeY, storeZ;

	// constants setting
	Fq2_set_ui(&c8, 8, 0);
	Fq2_set_ui(&c16, 16, 0);
	Fq2_set_ui(&c18, 18, 0);
	Fq2_set_ui(&c27, 27, 0);
	Fq2_set_ui(&c36, 36, 0);
	
	//constants to Montgomery rep
	Fq2_mont_rep(&c8, &c8);
	Fq2_mont_rep(&c16, &c16);
	Fq2_mont_rep(&c18, &c18);
	Fq2_mont_rep(&c27, &c27);
	Fq2_mont_rep(&c36, &c36);

	//---NB! aux1,2,3 are fixed within the computation of the current
	// coordinate only
	// Zr
	Fq2_mont_mul(&aux1, &(p->valY), &(p->valY));
	Fq2_mont_mul(&aux1, &aux1, &(p->valY));
	Fq2_mont_mul(&aux1, &aux1, &(p->valZ));
	Fq2_mont_mul(&aux1, &aux1, &(p->valZ));
	Fq2_mont_mul(&aux1, &aux1, &(p->valZ));
	Fq2_mont_mul(&storeZ, &c8, &aux1);

	// Xr
	Fq2_mont_mul(&aux1, &(p->valX), &(p->valX));
	Fq2_mont_mul(&aux1, &aux1, &(p->valX));
	Fq2_mont_mul(&aux1, &aux1, &(p->valX));
	Fq2_mont_mul(&aux1, &aux1, &(p->valY));
	Fq2_mont_mul(&aux1, &aux1, &(p->valZ));
	Fq2_mont_mul(&aux1, &c18, &aux1);

	Fq2_mont_mul(&aux2, &(p->valX), &(p->valY));
	Fq2_mont_mul(&aux2, &aux2, &(p->valY));
	Fq2_mont_mul(&aux2, &aux2, &(p->valY));
	Fq2_mont_mul(&aux2, &aux2, &(p->valZ));
	Fq2_mont_mul(&aux2, &aux2, &(p->valZ));
	Fq2_mont_mul(&aux2, &c16, &aux2);

	Fq2_sub(&storeX, &aux1, &aux2);

	// Yr
	Fq2_mont_mul(&aux1, &(p->valX), &(p->valX));
	Fq2_mont_mul(&aux1, &aux1, &(p->valX));
	Fq2_mont_mul(&aux1, &aux1, &(p->valX));
	Fq2_mont_mul(&aux1, &aux1, &(p->valX));
	Fq2_mont_mul(&aux1, &aux1, &(p->valX));
	Fq2_mont_mul(&aux1, &c27, &aux1);

	Fq2_mont_mul(&aux2, &(p->valX), &(p->valX));
	Fq2_mont_mul(&aux2, &aux2, &(p->valX));
	Fq2_mont_mul(&aux2, &aux2, &(p->valY));
	Fq2_mont_mul(&aux2, &aux2, &(p->valY));
	Fq2_mont_mul(&aux2, &aux2, &(p->valZ));
	Fq2_mont_mul(&aux2, &c36, &aux2);

	Fq2_mont_mul(&aux3, &(p->valY), &(p->valY));
	Fq2_mont_mul(&aux3, &aux3, &(p->valY));
	Fq2_mont_mul(&aux3, &aux3, &(p->valY));
	Fq2_mont_mul(&aux3, &aux3, &(p->valZ));
	Fq2_mont_mul(&aux3, &aux3, &(p->valZ));
	Fq2_mont_mul(&aux3, &c8, &aux3);

	Fq2_sub(&aux2, &aux2, &aux1);
	Fq2_sub(&storeY, &aux2, &aux3);

	// storing the result
	Fq2_assign(&(rop->valX), &storeX);
	Fq2_assign(&(rop->valY), &storeY);
	Fq2_assign(&(rop->valZ), &storeZ);
	
	G2_Z_correction(rop);
}


void G2_Z_add(G2_Z* rop, G2_Z* p, G2_Z* q)
{
	
	Fq2 zero;
	Fq2_set_zero(&zero);
	
	uint32_t mask_inf1 = -Fq2_cmp(&(p->valZ), &zero); //0xFFFFFFFF if p=O; 0, otherwise
	uint32_t mask_inf2 = -Fq2_cmp(&(q->valZ), &zero); //0xFFFFFFFF if q=O; 0, otherwise
	uint32_t mask_inf = ~mask_inf1 & ~mask_inf2; //0xFFFFFFFF, if p!=O & q!=O; 0, otherwise
	
	G2_Z res1, res2;
	G2_Z_adding(&res1, p, q);
	G2_Z_doubling(&res2, p);
	
	uint32_t mask = -G2_Z_cmp(p, q); //0xFFFFFFFF if p==q; 0, otherwise
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->valX.val0.limbs[i] = ( ((res1.valX.val0.limbs[i] & ~mask) | (res2.valX.val0.limbs[i] & mask)) & mask_inf ) | ( p->valX.val0.limbs[i] & mask_inf2 ) | ( q->valX.val0.limbs[i] & mask_inf1 );
		rop->valX.val1.limbs[i] = ( ((res1.valX.val1.limbs[i] & ~mask) | (res2.valX.val1.limbs[i] & mask)) & mask_inf ) | ( p->valX.val1.limbs[i] & mask_inf2 ) | ( q->valX.val1.limbs[i] & mask_inf1 );
		
		rop->valY.val0.limbs[i] = ( ((res1.valY.val0.limbs[i] & ~mask) | (res2.valY.val0.limbs[i] & mask)) & mask_inf ) | ( p->valY.val0.limbs[i] & mask_inf2 ) | ( q->valY.val0.limbs[i] & mask_inf1 );
		rop->valY.val1.limbs[i] = ( ((res1.valY.val1.limbs[i] & ~mask) | (res2.valY.val1.limbs[i] & mask)) & mask_inf ) | ( p->valY.val1.limbs[i] & mask_inf2 ) | ( q->valY.val1.limbs[i] & mask_inf1 );
		
		rop->valZ.val0.limbs[i] = ( ((res1.valZ.val0.limbs[i] & ~mask) | (res2.valZ.val0.limbs[i] & mask)) & mask_inf ) | ( p->valZ.val0.limbs[i] & mask_inf2 ) | ( q->valZ.val0.limbs[i] & mask_inf1 );
		rop->valZ.val1.limbs[i] = ( ((res1.valZ.val1.limbs[i] & ~mask) | (res2.valZ.val1.limbs[i] & mask)) & mask_inf ) | ( p->valZ.val1.limbs[i] & mask_inf2 ) | ( q->valZ.val1.limbs[i] & mask_inf1 );
	}
	
}

void G2_Z_mul(G2_Z* rop, Fr* n, G2_Z* p) // rop = [n]P, Montgomery ladder
{
	G2_Z R0, R1, T_swap;
	G2_Z_set_inf(&R0);
	G2_Z_assign(&R1, p);
	
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
			T_swap.valX.val0.limbs[j] = (R0.valX.val0.limbs[j] ^ R1.valX.val0.limbs[j]) & mask;
			T_swap.valX.val1.limbs[j] = (R0.valX.val1.limbs[j] ^ R1.valX.val1.limbs[j]) & mask;
			T_swap.valY.val0.limbs[j] = (R0.valY.val0.limbs[j] ^ R1.valY.val0.limbs[j]) & mask;
			T_swap.valY.val1.limbs[j] = (R0.valY.val1.limbs[j] ^ R1.valY.val1.limbs[j]) & mask;
			T_swap.valZ.val0.limbs[j] = (R0.valZ.val0.limbs[j] ^ R1.valZ.val0.limbs[j]) & mask;
			T_swap.valZ.val1.limbs[j] = (R0.valZ.val1.limbs[j] ^ R1.valZ.val1.limbs[j]) & mask;
			
			R0.valX.val0.limbs[j] ^= T_swap.valX.val0.limbs[j];
			R0.valX.val1.limbs[j] ^= T_swap.valX.val1.limbs[j];
			R0.valY.val0.limbs[j] ^= T_swap.valY.val0.limbs[j];
			R0.valY.val1.limbs[j] ^= T_swap.valY.val1.limbs[j];
			R0.valZ.val0.limbs[j] ^= T_swap.valZ.val0.limbs[j];
			R0.valZ.val1.limbs[j] ^= T_swap.valZ.val1.limbs[j];
			
			R1.valX.val0.limbs[j] ^= T_swap.valX.val0.limbs[j];
			R1.valX.val1.limbs[j] ^= T_swap.valX.val1.limbs[j];
			R1.valY.val0.limbs[j] ^= T_swap.valY.val0.limbs[j];
			R1.valY.val1.limbs[j] ^= T_swap.valY.val1.limbs[j];
			R1.valZ.val0.limbs[j] ^= T_swap.valZ.val0.limbs[j];
			R1.valZ.val1.limbs[j] ^= T_swap.valZ.val1.limbs[j];
		}
		
		G2_Z_add(&R1, &R0, &R1);
		G2_Z_add(&R0, &R0, &R0);
		
		//swapping back
		for(int j = 0; j < N_LIMBS; j++)
		{
			T_swap.valX.val0.limbs[j] = (R0.valX.val0.limbs[j] ^ R1.valX.val0.limbs[j]) & mask;
			T_swap.valX.val1.limbs[j] = (R0.valX.val1.limbs[j] ^ R1.valX.val1.limbs[j]) & mask;
			T_swap.valY.val0.limbs[j] = (R0.valY.val0.limbs[j] ^ R1.valY.val0.limbs[j]) & mask;
			T_swap.valY.val1.limbs[j] = (R0.valY.val1.limbs[j] ^ R1.valY.val1.limbs[j]) & mask;
			T_swap.valZ.val0.limbs[j] = (R0.valZ.val0.limbs[j] ^ R1.valZ.val0.limbs[j]) & mask;
			T_swap.valZ.val1.limbs[j] = (R0.valZ.val1.limbs[j] ^ R1.valZ.val1.limbs[j]) & mask;
			
			R0.valX.val0.limbs[j] ^= T_swap.valX.val0.limbs[j];
			R0.valX.val1.limbs[j] ^= T_swap.valX.val1.limbs[j];
			R0.valY.val0.limbs[j] ^= T_swap.valY.val0.limbs[j];
			R0.valY.val1.limbs[j] ^= T_swap.valY.val1.limbs[j];
			R0.valZ.val0.limbs[j] ^= T_swap.valZ.val0.limbs[j];
			R0.valZ.val1.limbs[j] ^= T_swap.valZ.val1.limbs[j];
			
			R1.valX.val0.limbs[j] ^= T_swap.valX.val0.limbs[j];
			R1.valX.val1.limbs[j] ^= T_swap.valX.val1.limbs[j];
			R1.valY.val0.limbs[j] ^= T_swap.valY.val0.limbs[j];
			R1.valY.val1.limbs[j] ^= T_swap.valY.val1.limbs[j];
			R1.valZ.val0.limbs[j] ^= T_swap.valZ.val0.limbs[j];
			R1.valZ.val1.limbs[j] ^= T_swap.valZ.val1.limbs[j];
		}
		
	}
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->valX.val0.limbs[i] = R0.valX.val0.limbs[i];
		rop->valX.val1.limbs[i] = R0.valX.val1.limbs[i];
		rop->valY.val0.limbs[i] = R0.valY.val0.limbs[i];
		rop->valY.val1.limbs[i] = R0.valY.val1.limbs[i];
		rop->valZ.val0.limbs[i] = R0.valZ.val0.limbs[i];
		rop->valZ.val1.limbs[i] = R0.valZ.val1.limbs[i];
	}
}

// Print
void G2_print(G2* a)
{
	if (a->inf == 1)
	{
		printf("inf");
	}

	if (a->inf == 0)
	{
		printf("(");
		Fq2_print(&(a->valx));
		printf(", ");
		Fq2_print(&(a->valy));
		printf(")");
	}
}

void G2_Z_print(G2_Z* a)
{
	printf("(");
	Fq2_print(&(a->valX));
	printf(" : ");
	Fq2_print(&(a->valY));
	printf(" : ");
	Fq2_print(&(a->valZ));
	printf(")");
}

void G2_print_dec(G2* a)
{
	if (a->inf == 1)
	{
		printf("inf");
	}

	if (a->inf == 0)
	{
		printf("(");
		Fq2_print_dec(&(a->valx));
		printf(", ");
		Fq2_print_dec(&(a->valy));
		printf(")");
	}
}

void G2_Z_print_dec(G2_Z* a)
{
	printf("(");
	Fq2_print_dec(&(a->valX));
	printf(" : ");
	Fq2_print_dec(&(a->valY));
	printf(" : ");
	Fq2_print_dec(&(a->valZ));
	printf(")");
}
/*
int main()
{
	Fq_field_init();
	Fr_field_init();
	G2_group_init();
	G2_print_dec(&genG2);
	printf("\n");
	
	G2_Z p;
	G2_to_project(&p, &genG2);
	G2_Z_mont_rep(&p, &p);
	
	Fr num;
	
	Fr_set_hex_str(&num, "aa1345fffffff");
	
	Fr_print_dec(&num);
	printf("\n");
	
	Fr_print(&num);
	printf("\n");
	
	G2_Z rop;
	G2_Z_mul(&rop, &num, &p);
	
	G2 final;
	G2_to_affine(&final ,&rop);
	G2_mont_rep_inv(&final, &final);
	
	G2_print_dec(&final);
	printf("\n");
	
	
	return 0;
}
*/
