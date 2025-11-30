#include "g1.h"
#include "../fq/bls_params.h"
#include "../fq/fq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

mpz_t orderG1;
G1 genG1;
mpz_t cofactorG1;
char* bORD_STR;
int len_bORD_STR;

// Group init
void G1_group_init()
{
	mpz_init(orderG1);
	G1_init(&genG1);
	mpz_init(cofactorG1);

	// auxiliary
	Fq xgen, ygen;
	Fq_init(&xgen);
	Fq_init(&ygen);

	Fq_set_str(&xgen, "0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a"
					  "3f171bac586c55e83ff97a1aeffb3af00adb22c6bb");
	Fq_set_str(&ygen, "0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18"
					  "cb2c04b3edd03cc744a2888ae40caa232946c5e7e1");

	// Defining
	mpz_set_str(
		orderG1,
		"0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
		0);
	G1_set(&genG1, &xgen, &ygen);
	mpz_set_str(cofactorG1, "0x396c8c005555e1568c00aaab0000aaab", 0);

	bORD_STR = mpz_get_str(NULL, 2, orderG1);
	len_bORD_STR = strlen(bORD_STR);

	// clearing
	Fq_clear(&xgen);
	Fq_clear(&ygen);
}

void G1_group_clear()
{
	mpz_clears(orderG1, cofactorG1, NULL);
	G1_clear(&genG1);
	free(bORD_STR);
}

// Init
void G1_init(G1* a)
{
	Fq_init(&(a->valx));
	Fq_init(&(a->valy));

	a->inf = '0';
}

void G1_set(G1* a, Fq* x, Fq* y)
{
	Fq_set(&(a->valx), x->val);
	Fq_set(&(a->valy), y->val);

	a->inf = '0';
}

void G1_init_inf(G1* a)
{
	Fq_init(&(a->valx));
	Fq_init(&(a->valy));

	a->inf = '1';
}

void G1_set_inf(G1* a)
{
	Fq_set_ui(&(a->valx), 0);
	Fq_set_ui(&(a->valy), 0);

	a->inf = '1';
}

void G1_assign(G1* rop, G1* a)
{
	Fq_assign(&(rop->valx), &(a->valx));
	Fq_assign(&(rop->valy), &(a->valy));

	rop->inf = a->inf;
}

void G1_clear(G1* a)
{
	Fq_clear(&(a->valx));
	Fq_clear(&(a->valy));

	a->inf = '0';
}

// Init projective
void G1_Z_init(G1_Z* a)
{
	Fq_init(&(a->valX));
	Fq_init(&(a->valY));
	Fq_init(&(a->valZ));
}
void G1_Z_set(G1_Z* a, Fq* X, Fq* Y, Fq* Z)
{
	Fq_set(&(a->valX), X->val);
	Fq_set(&(a->valY), Y->val);
	Fq_set(&(a->valZ), Z->val);
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
void G1_Z_clear(G1_Z* a)
{
	Fq_clear(&(a->valX));
	Fq_clear(&(a->valY));
	Fq_clear(&(a->valZ));
}

// Maps
void G1_to_group(G1* rop, G1* a)
{
	if(G1_curve_check(a) != 1)
	{
		errno = EDOM;
		perror("Domain error occured");
		exit(EXIT_FAILURE);
	}
	
	G1_mul(rop, cofactorG1, a);
}
void G1_Z_to_group(G1_Z* rop, G1_Z* a)
{
	if(G1_Z_curve_check(a) != 1)
	{
		errno = EDOM;
		perror("Domain error occured");
		exit(EXIT_FAILURE);
	}
	
	G1_Z_mul(rop, cofactorG1, a);

}
void G1_to_affine(G1* rop, G1_Z* a)
{
	Fq zero;
	Fq_init(&zero);
	Fq_set_ui(&zero, 0);

	if (Fq_cmp(&(a->valZ), &zero) == 0) // if Z == 0
	{
		G1_set_inf(rop);
	}

	else // if Z != 0
	{
		// auxiliary init
		Fq tempZ, tempZinv;
		Fq_init(&tempZ);
		Fq_init(&tempZinv);

		// auxiliary comp
		Fq_assign(&tempZ, &(a->valZ));
		Fq_mul_inv(&tempZinv, &tempZ);

		// computation
		Fq_mul(&(rop->valx), &(a->valX), &tempZinv);
		Fq_mul(&(rop->valy), &(a->valY), &tempZinv);

		// aux clear
		Fq_clear(&tempZ);
		Fq_clear(&tempZinv);
	}

	Fq_clear(&zero);
}
void G1_to_project(G1_Z* rop, G1* a)
{
	if (a->inf == '1')
	{
		G1_Z_set_inf(rop);
	}
	else
	{
		Fq_assign(&(rop->valX), &(a->valx));
		Fq_assign(&(rop->valY), &(a->valy));
		Fq_set_ui(&(rop->valZ), 1);
	}
}

// Curve check
int G1_curve_check(G1* a)
{

	if ((a->inf) == '1')
	{
		return 1;
	}
	
	int check = 0;	
	
	// E(Fq): y^2 = x^3 + 4
	// auxiliary init
	Fq tempx3, tempy2, c4temp, temp;
	Fq_init(&tempx3);
	Fq_init(&tempy2);
	Fq_init(&c4temp);
	Fq_init(&temp);

	// auxiliary comp
	Fq_mul(&tempx3, &(a->valx), &(a->valx));
	Fq_mul(&tempx3, &tempx3, &(a->valx));

	Fq_mul(&tempy2, &(a->valy), &(a->valy));

	Fq_set_ui(&c4temp, 4);

	// Computation
	Fq_add(&temp, &tempx3, &c4temp);

	if (Fq_cmp(&tempy2, &temp) == 0)
	{
		check = 1;
	}

	// clear
	Fq_clear(&tempx3);
	Fq_clear(&tempy2);
	Fq_clear(&c4temp);
	Fq_clear(&temp);

	return check;
}

// Group check
int G1_group_check(G1* a)
{
	// If a is in G1, then [r]a = O
	G1 temp;
	G1_init(&temp);
	G1_assign(&temp, a);

	G1_mul(&temp, orderG1, a);

	if (temp.inf == '1')
	{
		G1_clear(&temp);
		return 1;
	}

	G1_clear(&temp);
	return 0;
}

// Curve check projective
int G1_Z_curve_check(G1_Z* a)
{
	// E(Fq): Y^2 * Z = X^3 + 4*Z^3
	Fq zero;
	Fq_init(&zero);
	Fq_set_ui(&zero, 0);
	
	if(Fq_cmp(&(a->valZ), &zero) == 0) //inf
	{
		Fq_clear(&zero);
		return 1;
	}
	
	int check = 0; 
	
	//init
	Fq tempX, tempY, tempZ, aux;
	Fq_init(&tempX);
	Fq_init(&tempY);
	Fq_init(&tempZ);
	Fq_init(&aux);
	
	Fq_assign(&tempX, &(a->valX));
	Fq_assign(&tempY, &(a->valY));
	Fq_assign(&tempZ, &(a->valZ));
	
	//computations
	Fq_mul(&tempY, &tempY, &tempY);
	Fq_mul(&tempY, &tempY, &(a->valZ));
	
	Fq_mul(&tempX, &tempX, &tempX);
	Fq_mul(&tempX, &tempX, &(a->valX));
	
	Fq_mul(&tempZ, &tempZ, &tempZ);
	Fq_mul(&tempZ, &tempZ, &(a->valZ));
	Fq_add(&tempZ, &tempZ, &tempZ);
	Fq_add(&tempZ, &tempZ, &tempZ);
	
	Fq_add(&aux, &tempX, &tempZ);
	
	if(Fq_cmp(&tempY, &aux) == 0)
	{
		check = 1;
	}
	
	Fq_clear(&tempX);
	Fq_clear(&tempY);
	Fq_clear(&tempZ);
	Fq_clear(&aux);
	Fq_clear(&zero);
	
	return check;
	
}

int G1_Z_group_check(G1_Z* a)
{
	// If a is in G1, then [r]a = 0.
	// In the projective coordinates, Za = 0

	int inGroup = 0;

	Fq zero;
	Fq_init(&zero);
	Fq_set_ui(&zero, 0);

	G1_Z temp;
	G1_Z_init(&temp);

	G1_Z_mul(&temp, orderG1, a);

	if (Fq_cmp(&(temp.valZ), &zero) == 0)
	{
		inGroup = 1;
	}

	Fq_clear(&zero);
	G1_Z_clear(&temp);

	return inGroup;
}

// Arithmetic
void G1_add_inv(G1* rop, G1* a)
{
	if (a->inf == '1')
	{
		G1_assign(rop, a);
	}

	if (a->inf == '0')
	{
		Fq negy;
		Fq_init(&negy);

		Fq_add_inv(&negy, &(a->valy));

		Fq_assign(&(rop->valx), &(a->valx));
		Fq_assign(&(rop->valy), &negy);
		rop->inf = '0';

		Fq_clear(&negy);
	}
}

void G1_add(G1* rop, G1* p, G1* q)
{
	if (p->inf == '1') // if P == O, rop = P + Q = Q
	{
		G1_assign(rop, q);
		return;
	}

	if (q->inf == '1') // if Q == O, rop = P + Q = P
	{
		G1_assign(rop, p);
		return;
	}

	// if P!= O && Q!= O

	// Case 1: xp != xq

	if (Fq_cmp(&(p->valx), &(q->valx)) != 0)
	{
		// auxiliary init
		Fq tempL, tempL2, tempV, tempXpXqInv, tempXr, tempYr;
		Fq_init(&tempL);
		Fq_init(&tempL2);
		Fq_init(&tempV);
		Fq_init(&tempXpXqInv);
		Fq_init(&tempXr);
		Fq_init(&tempYr);

		// auxiliary computation
		Fq_sub(&tempXpXqInv, &(p->valx), &(q->valx));
		Fq_mul_inv(&tempXpXqInv, &tempXpXqInv);

		Fq_sub(&tempL, &(p->valy), &(q->valy));
		Fq_mul(&tempL, &tempL, &tempXpXqInv);

		Fq_mul(&tempL2, &tempL, &tempL);

		Fq_mul(&tempV, &tempL, &(p->valx));
		Fq_sub(&tempV, &(p->valy), &tempV);

		// computation
		// xr
		Fq_sub(&tempXr, &tempL2, &(p->valx));
		Fq_sub(&tempXr, &tempXr, &(q->valx));

		// yr
		Fq_mul(&tempYr, &tempL, &tempXr);
		Fq_add(&tempYr, &tempYr, &tempV);
		Fq_add_inv(&tempYr, &tempYr);

		// storing
		Fq_assign(&(rop->valx), &tempXr);
		Fq_assign(&(rop->valy), &tempYr);
		rop->inf = '0';

		// clear
		Fq_clear(&tempL);
		Fq_clear(&tempL2);
		Fq_clear(&tempV);
		Fq_clear(&tempXpXqInv);
		Fq_clear(&tempXr);
		Fq_clear(&tempYr);

		return;
	}

	Fq tempNull, tempYpAddInv;
	Fq_init(&tempNull);
	Fq_init(&tempYpAddInv);

	Fq_set_ui(&tempNull, 0);
	Fq_add_inv(&tempYpAddInv, &(p->valy));

	// Case 2: xp = xq, yp = yq, and yp!=0 (Doubling)
	if (Fq_cmp(&(p->valx), &(q->valx)) == 0 &&
		Fq_cmp(&(p->valy), &(q->valy)) == 0 &&
		Fq_cmp(&(p->valy), &tempNull) != 0)
	{
		// auxiliary init
		Fq tempL, tempL2, tempV, c2tempYpMulInv, tempXp2, tempXr, tempYr;
		Fq_init(&tempL);
		Fq_init(&tempL2);
		Fq_init(&tempV);
		Fq_init(&c2tempYpMulInv);
		Fq_init(&tempXp2);
		Fq_init(&tempXr);
		Fq_init(&tempYr);

		// auxiliary computation
		Fq_add(&c2tempYpMulInv, &(p->valy), &(p->valy));
		Fq_mul_inv(&c2tempYpMulInv, &c2tempYpMulInv);

		Fq_mul(&tempXp2, &(p->valx), &(p->valx));
		Fq_add(&tempL, &tempXp2, &tempXp2);
		Fq_add(&tempL, &tempL, &tempXp2);

		Fq_mul(&tempL, &tempL, &c2tempYpMulInv);

		Fq_mul(&tempL2, &tempL, &tempL);

		Fq_mul(&tempV, &tempL, &(p->valx));
		Fq_sub(&tempV, &(p->valy), &tempV);

		// computation

		// xr
		Fq_sub(&tempXr, &tempL2, &(p->valx));
		Fq_sub(&tempXr, &tempXr, &(p->valx));

		// yr
		Fq_mul(&tempYr, &tempL, &tempXr);
		Fq_add(&tempYr, &tempYr, &tempV);
		Fq_add_inv(&tempYr, &tempYr);

		// storing
		Fq_assign(&(rop->valx), &tempXr);
		Fq_assign(&(rop->valy), &tempYr);
		rop->inf = '0';

		// clear
		Fq_clear(&tempL);
		Fq_clear(&tempL2);
		Fq_clear(&tempV);
		Fq_clear(&c2tempYpMulInv);
		Fq_clear(&tempXp2);
		Fq_clear(&tempNull);
		Fq_clear(&tempYpAddInv);
		Fq_clear(&tempXr);
		Fq_clear(&tempYr);

		return;
	}

	// Case 3: Q = -P. xp = xq, yq = -yp.
	if (Fq_cmp(&(p->valx), &(q->valx)) == 0 &&
		Fq_cmp(&(q->valy), &tempYpAddInv) == 0)
	{
		G1_set_inf(rop);

		Fq_clear(&tempNull);
		Fq_clear(&tempYpAddInv);

		return;
	}

	// Safeguard
	Fq_clear(&tempNull);
	Fq_clear(&tempYpAddInv);
}

void G1_mul(G1* rop, mpz_t n, G1* p) 
{
	if (mpz_cmp_ui(n, 0) == 0)
	{
		G1_set_inf(rop);
		return;
	}

	// auxiliary init
	G1 temp;
	G1_init(&temp);
	G1_assign(&temp, p);

	char* bin = mpz_get_str(NULL, 2, n);
	int len = strlen(bin);

	for (int i = 1; i < len; i++)
	{
		G1_add(&temp, &temp, &temp);

		if (bin[i] == '1')
		{
			G1_add(&temp, &temp, p);
		}
	}

	G1_assign(rop, &temp);

	// clear
	G1_clear(&temp);
	free(bin);
}

// Group Arithmetic Projective
void G1_Z_add_inv(G1_Z* rop, G1_Z* a)
{
	// aux init
	Fq tempY;
	Fq_init(&tempY);

	Fq_add_inv(&tempY, &(a->valY));

	Fq_assign(&(rop->valX), &(a->valX));
	Fq_assign(&(rop->valY), &tempY);
	Fq_assign(&(rop->valZ), &(a->valZ));

	// aux clear
	Fq_clear(&tempY);
}

int G1_Z_cmp(G1_Z* p, G1_Z* q)
{
	//(Xp : Yp : Zp) ~ (Xq : Yq : Zq) <=>
	// <=> [XpYq - XqYp, YpZq - YqZp, ZpXq - ZqXp] = [0, 0, 0]

	int cmp = 0;

	// aux init
	Fq tempXpYq, tempXqYp, tempYpZq, tempYqZp, tempZpXq, tempZqXp;
	Fq_init(&tempXpYq);
	Fq_init(&tempXqYp);
	Fq_init(&tempYpZq);
	Fq_init(&tempYqZp);
	Fq_init(&tempZpXq);
	Fq_init(&tempZqXp);

	// computation
	Fq_mul(&tempXpYq, &(p->valX), &(q->valY));
	Fq_mul(&tempXqYp, &(p->valY), &(q->valX));
	Fq_mul(&tempYpZq, &(p->valY), &(q->valZ));
	Fq_mul(&tempYqZp, &(p->valZ), &(q->valY));
	Fq_mul(&tempZpXq, &(p->valZ), &(q->valX));
	Fq_mul(&tempZqXp, &(p->valX), &(q->valZ));

	if (Fq_cmp(&tempXpYq, &tempXqYp) != 0)
	{
		cmp = 1;
	}

	if (Fq_cmp(&tempYpZq, &tempYqZp) != 0)
	{
		cmp = 1;
	}

	if (Fq_cmp(&tempZpXq, &tempZqXp) != 0)
	{
		cmp = 1;
	}

	// aux clear
	Fq_clear(&tempXpYq);
	Fq_clear(&tempXqYp);
	Fq_clear(&tempYpZq);
	Fq_clear(&tempYqZp);
	Fq_clear(&tempZpXq);
	Fq_clear(&tempZqXp);

	return cmp;
}

void G1_Z_add(G1_Z* rop, G1_Z* p, G1_Z* q)
{
	Fq zero;
	Fq_init(&zero);
	Fq_set_ui(&zero, 0);

	// if P == O, P + Q = Q
	if (Fq_cmp(&(p->valZ), &zero) == 0)
	{
		G1_Z_assign(rop, q);
		Fq_clear(&zero);
		return;
	}

	// if Q == O, P + Q = P
	if (Fq_cmp(&(q->valZ), &zero) == 0)
	{
		G1_Z_assign(rop, p);
		Fq_clear(&zero);
		return;
	}

	// General case. P!=Q
	if (G1_Z_cmp(p, q) != 0)
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
		Fq_init(&tempXpYq);
		Fq_init(&tempXqYp);
		Fq_init(&tempYpZq);
		Fq_init(&tempYqZp);
		Fq_init(&tempXqZp);
		Fq_init(&tempXpZq);
		Fq_init(&tempZpZq);
		Fq_init(&diff1);
		Fq_init(&diff2);
		Fq_init(&diff3);
		Fq_init(&sum);
		Fq_init(&delta1);
		Fq_init(&delta2);
		Fq_init(&aux1);
		Fq_init(&aux2);
		Fq_init(&bracket);
		Fq_init(&c1);
		Fq_init(&c2);
		Fq_init(&c3);

		// aux computation
		Fq_mul(&tempXpYq, &(p->valX), &(q->valY));
		Fq_mul(&tempXqYp, &(p->valY), &(q->valX));
		Fq_mul(&tempYpZq, &(p->valY), &(q->valZ));
		Fq_mul(&tempYqZp, &(p->valZ), &(q->valY));
		Fq_mul(&tempXqZp, &(p->valZ), &(q->valX));
		Fq_mul(&tempXpZq, &(p->valX), &(q->valZ));
		Fq_mul(&tempZpZq, &(p->valZ), &(q->valZ));

		Fq_sub(&diff1, &tempXpZq, &tempXqZp);
		Fq_mul(&diff2, &diff1, &diff1);
		Fq_mul(&diff3, &diff2, &diff1);

		Fq_add(&sum, &tempXpZq, &tempXqZp);
		Fq_sub(&delta1, &tempYpZq, &tempYqZp);
		Fq_mul(&delta2, &delta1, &delta1);

		Fq_mul(&aux1, &tempZpZq, &delta2);
		Fq_mul(&aux2, &diff2, &sum);
		Fq_sub(
			&bracket, &aux1,
			&aux2); //[ZpZq * (YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq + XqZp)]

		// Zr
		Fq_mul(&(rop->valZ), &tempZpZq, &diff3);

		// Xr
		Fq_mul(&(rop->valX), &diff1, &bracket);

		// Yr
		Fq_sub(&c3, &tempXqYp, &tempXpYq);
		Fq_mul(&c1, &tempZpZq, &c3);
		Fq_mul(&c1, &c1, &diff2);
		Fq_mul(&c2, &delta1, &bracket);

		Fq_sub(&(rop->valY), &c1, &c2);

		// aux clear
		Fq_clear(&tempXpYq);
		Fq_clear(&tempXqYp);
		Fq_clear(&tempYpZq);
		Fq_clear(&tempYqZp);
		Fq_clear(&tempXqZp);
		Fq_clear(&tempXpZq);
		Fq_clear(&tempZpZq);
		Fq_clear(&diff1);
		Fq_clear(&diff2);
		Fq_clear(&diff3);
		Fq_clear(&sum);
		Fq_clear(&delta1);
		Fq_clear(&delta2);
		Fq_clear(&aux1);
		Fq_clear(&aux2);
		Fq_clear(&bracket);
		Fq_clear(&c1);
		Fq_clear(&c2);
		Fq_clear(&c3);
	}

	else // P == Q, doubling
	{
		/*

		Zr = 8Yp^3Zp^3

		Xr = 18Xp^4YpZp - 16XpYp^3Zp^2

		Yr = -27Xp^6 + 36Xp^3Yp^2Zp - 8Yp^4Zp^2

		*/

		// aux init
		Fq c8, c16, c18, c27, c36, aux1, aux2, aux3, storeX, storeY, storeZ;
		Fq_init(&c8);
		Fq_init(&c16);
		Fq_init(&c18);
		Fq_init(&c27);
		Fq_init(&c36);
		Fq_init(&aux1);
		Fq_init(&aux2);
		Fq_init(&aux3);
		Fq_init(&storeX);
		Fq_init(&storeY);
		Fq_init(&storeZ);

		// constants setting
		Fq_set_ui(&c8, 8);
		Fq_set_ui(&c16, 16);
		Fq_set_ui(&c18, 18);
		Fq_set_ui(&c27, 27);
		Fq_set_ui(&c36, 36);

		//---NB! aux1,2,3 are fixed within the computation of the current
		// coordinate only
		// Zr
		Fq_mul(&aux1, &(p->valY), &(p->valY));
		Fq_mul(&aux1, &aux1, &(p->valY));
		Fq_mul(&aux1, &aux1, &(p->valZ));
		Fq_mul(&aux1, &aux1, &(p->valZ));
		Fq_mul(&aux1, &aux1, &(p->valZ));
		Fq_mul(&storeZ, &c8, &aux1);

		// Xr
		Fq_mul(&aux1, &(p->valX), &(p->valX));
		Fq_mul(&aux1, &aux1, &(p->valX));
		Fq_mul(&aux1, &aux1, &(p->valX));
		Fq_mul(&aux1, &aux1, &(p->valY));
		Fq_mul(&aux1, &aux1, &(p->valZ));
		Fq_mul(&aux1, &c18, &aux1);

		Fq_mul(&aux2, &(p->valX), &(p->valY));
		Fq_mul(&aux2, &aux2, &(p->valY));
		Fq_mul(&aux2, &aux2, &(p->valY));
		Fq_mul(&aux2, &aux2, &(p->valZ));
		Fq_mul(&aux2, &aux2, &(p->valZ));
		Fq_mul(&aux2, &c16, &aux2);

		Fq_sub(&storeX, &aux1, &aux2);

		// Yr
		Fq_mul(&aux1, &(p->valX), &(p->valX));
		Fq_mul(&aux1, &aux1, &(p->valX));
		Fq_mul(&aux1, &aux1, &(p->valX));
		Fq_mul(&aux1, &aux1, &(p->valX));
		Fq_mul(&aux1, &aux1, &(p->valX));
		Fq_mul(&aux1, &c27, &aux1);

		Fq_mul(&aux2, &(p->valX), &(p->valX));
		Fq_mul(&aux2, &aux2, &(p->valX));
		Fq_mul(&aux2, &aux2, &(p->valY));
		Fq_mul(&aux2, &aux2, &(p->valY));
		Fq_mul(&aux2, &aux2, &(p->valZ));
		Fq_mul(&aux2, &c36, &aux2);

		Fq_mul(&aux3, &(p->valY), &(p->valY));
		Fq_mul(&aux3, &aux3, &(p->valY));
		Fq_mul(&aux3, &aux3, &(p->valY));
		Fq_mul(&aux3, &aux3, &(p->valZ));
		Fq_mul(&aux3, &aux3, &(p->valZ));
		Fq_mul(&aux3, &c8, &aux3);

		Fq_sub(&aux2, &aux2, &aux1);
		Fq_sub(&storeY, &aux2, &aux3);

		// storing the result
		Fq_assign(&(rop->valX), &storeX);
		Fq_assign(&(rop->valY), &storeY);
		Fq_assign(&(rop->valZ), &storeZ);

		// aux clear
		Fq_clear(&c8);
		Fq_clear(&c16);
		Fq_clear(&c18);
		Fq_clear(&c27);
		Fq_clear(&c36);
		Fq_clear(&aux1);
		Fq_clear(&aux2);
		Fq_clear(&aux3);
		Fq_clear(&storeX);
		Fq_clear(&storeY);
		Fq_clear(&storeZ);
	}

	// Safeguard
	Fq_clear(&zero);
}
void G1_Z_mul(G1_Z* rop, mpz_t n, G1_Z* p) // rop = [n]P
{
	G1_Z temp;
	G1_Z_init(&temp);
	G1_Z_assign(&temp, p);

	char* bin = mpz_get_str(NULL, 2, n);
	int len = strlen(bin);

	for (int i = 1; i < len; i++)
	{
		G1_Z_add(&temp, &temp, &temp);

		if (bin[i] == '1')
		{
			G1_Z_add(&temp, &temp, p);
		}
	}

	G1_Z_assign(rop, &temp);

	G1_Z_clear(&temp);
	free(bin);
}

// Print
void G1_print(G1* a)
{
	if (a->inf == '1')
	{
		printf("inf");
	}

	if (a->inf == '0')
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
