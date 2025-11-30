#include "g2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

mpz_t orderG2;
G2 genG2;
mpz_t cofactorG2;

// Group init
void G2_group_init()
{
	mpz_init(orderG2);
	G2_init(&genG2);
	mpz_init(cofactorG2);

	// auxiliary
	Fq2 xgen, ygen;
	Fq2_init(&xgen);
	Fq2_init(&ygen);

	Fq2_set_str(&xgen,
				"0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3"
				"d1770bac0326a805bbefd48056c8c121bdb8",
				"0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f"
				"5049334cf11213945d57e5ac7d055d042b7e");
	Fq2_set_str(&ygen,
				"0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160"
				"d12c923ac9cc3baca289e193548608b82801",
				"0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e"
				"99ab3f370d275cec1da1aaa9075ff05f79be");

	// Defining
	mpz_set_str(
		orderG2,
		"0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
		0);
	G2_set(&genG2, &xgen, &ygen);
	mpz_set_str(
		cofactorG2,
		"0x5d543a95414e7f1091d50792876a202cd91de4547085abaa68a205b2e5a7ddfa628f"
		"1cb4d9e82ef21537e293a6691ae1616ec6e786f0c70cf1c38e31c7238e5",
		0);

	// clearing
	Fq2_clear(&xgen);
	Fq2_clear(&ygen);
}
void G2_group_clear()
{
	mpz_clears(orderG2, cofactorG2, NULL);
	G2_clear(&genG2);
}

// Init
void G2_init(G2* a)
{
	Fq2_init(&(a->valx));
	Fq2_init(&(a->valy));

	a->inf = '0';
}

void G2_set(G2* a, Fq2* x, Fq2* y)
{
	Fq2_assign(&(a->valx), x);
	Fq2_assign(&(a->valy), y);

	a->inf = '0';
}

void G2_init_inf(G2* a)
{
	Fq2_init(&(a->valx));
	Fq2_init(&(a->valy));

	a->inf = '1';
}

void G2_set_inf(G2* a)
{
	Fq2_set_ui(&(a->valx), 0, 0);
	Fq2_set_ui(&(a->valy), 0, 0);

	a->inf = '1';
}

void G2_assign(G2* rop, G2* a)
{
	Fq2_assign(&(rop->valx), &(a->valx));
	Fq2_assign(&(rop->valy), &(a->valy));

	rop->inf = a->inf;
}

void G2_clear(G2* a)
{
	Fq2_clear(&(a->valx));
	Fq2_clear(&(a->valy));

	a->inf = '0';
}

void G2_Z_init(G2_Z* a)
{
	Fq2_init(&(a->valX));
	Fq2_init(&(a->valY));
	Fq2_init(&(a->valZ));
}
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
void G2_Z_clear(G2_Z* a)
{
	Fq2_clear(&(a->valX));
	Fq2_clear(&(a->valY));
	Fq2_clear(&(a->valZ));
}

// Maps
void G2_to_group(G2* rop, G2* a)
{
	if(G2_curve_check(a) != 1)
	{
		errno = EDOM;
		perror("Domain error occured");
		exit(EXIT_FAILURE);
	}
	
	G2_mul(rop, cofactorG2, a);
}

void G2_Z_to_group(G2_Z* rop, G2_Z* a)
{
	if(G2_Z_curve_check(a) != 1)
	{
		errno = EDOM;
		perror("Domain error occured");
		exit(EXIT_FAILURE);
	}
	
	G2_Z_mul(rop, cofactorG2, a);
}

void G2_to_affine(G2* rop, G2_Z* a)
{
	Fq2 zero;
	Fq2_init(&zero);
	Fq2_set_ui(&zero, 0, 0);

	if (Fq2_cmp(&(a->valZ), &zero) == 0) // if Z == 0
	{
		G2_set_inf(rop);
	}

	else // if Z != 0
	{
		// auxiliary init
		Fq2 tempZ, tempZinv;
		Fq2_init(&tempZ);
		Fq2_init(&tempZinv);

		// auxiliary comp
		Fq2_assign(&tempZ, &(a->valZ));
		Fq2_mul_inv(&tempZinv, &tempZ);

		// computation
		Fq2_mul(&(rop->valx), &(a->valX), &tempZinv);
		Fq2_mul(&(rop->valy), &(a->valY), &tempZinv);

		// aux clear
		Fq2_clear(&tempZ);
		Fq2_clear(&tempZinv);
	}

	Fq2_clear(&zero);
}
void G2_to_project(G2_Z* rop, G2* a)
{
	if (a->inf == '1')
	{
		G2_Z_set_inf(rop);
	}
	else
	{
		Fq2_assign(&(rop->valX), &(a->valx));
		Fq2_assign(&(rop->valY), &(a->valy));
		Fq2_set_ui(&(rop->valZ), 1, 0);
	}
}

// Curve check

int G2_curve_check(G2* a)
{
	int check = 0;

	if ((a->inf) == '1')
	{
		return 1;
	}

	// E(Fq2): y^2 = x^3 + 4(1 + u)
	// auxiliary init
	Fq2 tempx3, tempy2, c4temp, temp;
	Fq2_init(&tempx3);
	Fq2_init(&tempy2);
	Fq2_init(&c4temp);
	Fq2_init(&temp);

	// auxiliary comp
	Fq2_mul(&tempx3, &(a->valx), &(a->valx));
	Fq2_mul(&tempx3, &tempx3, &(a->valx));

	Fq2_mul(&tempy2, &(a->valy), &(a->valy));

	Fq2_set_ui(&c4temp, 4, 4);

	// Computation
	Fq2_add(&temp, &tempx3, &c4temp);

	if (Fq2_cmp(&tempy2, &temp) == 0)
	{
		check = 1;
	}

	// clear
	Fq2_clear(&tempx3);
	Fq2_clear(&tempy2);
	Fq2_clear(&c4temp);
	Fq2_clear(&temp);

	return check;
}

// Group check
int G2_group_check(G2* a)
{
	// If a is in G2, then [r]a = O
	G2 temp;
	G2_init(&temp);
	G2_assign(&temp, a);

	G2_mul(&temp, orderG2, a);

	if (temp.inf == '1')
	{
		G2_clear(&temp);
		return 1;
	}

	G2_clear(&temp);
	return 0;
}
// Curve check projective
int G2_Z_curve_check(G2_Z* a)
{
	// E(Fq): Y^2 * Z = X^3 + 4*(1+u)*Z^3
	Fq2 zero2;
	Fq2_init(&zero2);
	Fq2_set_ui(&zero2, 0, 0);
	
	if(Fq2_cmp(&(a->valZ), &zero2) == 0) //inf
	{
		Fq2_clear(&zero2);
		return 1;
	}
	
	int check = 0;
	
	Fq2 tempX, tempY, tempZ, c44, aux;
	Fq2_init(&tempX);
	Fq2_init(&tempY);
	Fq2_init(&tempZ);
	Fq2_init(&c44);
	Fq2_init(&aux);
	
	Fq2_assign(&tempX, &(a->valX));
	Fq2_assign(&tempY, &(a->valY));
	Fq2_assign(&tempZ, &(a->valZ));
	Fq2_set_ui(&c44, 4, 4);
	
	Fq2_mul(&tempY, &tempY, &tempY);
	Fq2_mul(&tempY, &tempY, &(a->valZ));
	
	Fq2_mul(&tempX, &tempX, &tempX);
	Fq2_mul(&tempX, &tempX, &(a->valX));
	
	Fq2_mul(&tempZ, &tempZ, &tempZ);
	Fq2_mul(&tempZ, &tempZ, &(a->valZ));
	Fq2_mul(&tempZ, &tempZ, &c44);
	
	Fq2_add(&aux, &tempX, &tempZ);
	
	if(Fq2_cmp(&aux, &tempY) == 0)
	{
		check = 1;
	}
	
	Fq2_clear(&zero2);
	Fq2_clear(&tempX);
	Fq2_clear(&tempY);
	Fq2_clear(&tempZ);
	Fq2_clear(&c44);
	Fq2_clear(&aux);
	
	return check;
}

int G2_Z_group_check(G2_Z* a)
{
	// If a is in G1, then [r]a = 0.
	// In the projective coordinates, Za = 0

	int inGroup = 0;

	Fq2 zero;
	Fq2_init(&zero);
	Fq2_set_ui(&zero, 0, 0);

	G2_Z temp;
	G2_Z_init(&temp);

	G2_Z_mul(&temp, orderG2, a);

	if (Fq2_cmp(&(temp.valZ), &zero) == 0)
	{
		inGroup = 1;
	}

	Fq2_clear(&zero);
	G2_Z_clear(&temp);

	return inGroup;
}
// Arithmetic
void G2_add_inv(G2* rop, G2* a)
{
	if (a->inf == '1')
	{
		G2_assign(rop, a);
	}

	if (a->inf == '0')
	{
		Fq2 negy;
		Fq2_init(&negy);

		Fq2_add_inv(&negy, &(a->valy));

		Fq2_assign(&(rop->valx), &(a->valx));
		Fq2_assign(&(rop->valy), &negy);
		rop->inf = '0';

		Fq2_clear(&negy);
	}
}

void G2_add(G2* rop, G2* p, G2* q)
{
	if (p->inf == '1') // if P == O, rop = P + Q = Q
	{
		G2_assign(rop, q);
		return;
	}

	if (q->inf == '1') // if Q == O, rop = P + Q = P
	{
		G2_assign(rop, p);
		return;
	}

	// if P!= O && Q!= O

	// Case 1: xp != xq

	if (Fq2_cmp(&(p->valx), &(q->valx)) != 0)
	{
		// auxiliary init
		Fq2 tempL, tempL2, tempV, tempXpXqInv, tempXr, tempYr;
		Fq2_init(&tempL);
		Fq2_init(&tempL2);
		Fq2_init(&tempV);
		Fq2_init(&tempXpXqInv);
		Fq2_init(&tempXr);
		Fq2_init(&tempYr);

		// auxiliary computation
		Fq2_sub(&tempXpXqInv, &(p->valx), &(q->valx));
		Fq2_mul_inv(&tempXpXqInv, &tempXpXqInv); // (xp - xq)^(-1);

		Fq2_sub(&tempL, &(p->valy), &(q->valy));
		Fq2_mul(&tempL, &tempL, &tempXpXqInv); // L = (yp - yq)/(xp - xq)

		Fq2_mul(&tempL2, &tempL, &tempL); // L^2 = L * L

		Fq2_mul(&tempV, &tempL, &(p->valx));
		Fq2_sub(&tempV, &(p->valy), &tempV); // V = yp - L * xp

		// computation
		// xr
		Fq2_sub(&tempXr, &tempL2, &(p->valx));
		Fq2_sub(&tempXr, &tempXr, &(q->valx)); // xr = L^2 - xp - xq

		// yr
		Fq2_mul(&tempYr, &tempL, &tempXr);
		Fq2_add(&tempYr, &tempYr, &tempV);
		Fq2_add_inv(&tempYr, &tempYr); // yr = -(L * xp + V)

		// storing
		Fq2_assign(&(rop->valx), &tempXr);
		Fq2_assign(&(rop->valy), &tempYr);
		rop->inf = '0'; // not infinity point

		// clear
		Fq2_clear(&tempL);
		Fq2_clear(&tempL2);
		Fq2_clear(&tempV);
		Fq2_clear(&tempXpXqInv);
		Fq2_clear(&tempXr);
		Fq2_clear(&tempYr);

		return;
	}

	Fq2 tempNull, tempYpAddInv;
	Fq2_init(&tempNull);
	Fq2_init(&tempYpAddInv);

	Fq2_set_ui(&tempNull, 0, 0);
	Fq2_add_inv(&tempYpAddInv, &(p->valy)); //-yp

	// Case 2: xp = xq, yp = yq, and yp!=0 (Doubling)
	if (Fq2_cmp(&(p->valx), &(q->valx)) == 0 &&
		Fq2_cmp(&(p->valy), &(q->valy)) == 0 &&
		Fq2_cmp(&(p->valy), &tempNull) != 0)
	{
		// auxiliary init
		Fq2 tempL, tempL2, tempV, c2tempYpMulInv, tempXp2, tempXr, tempYr;
		Fq2_init(&tempL);
		Fq2_init(&tempL2);
		Fq2_init(&tempV);
		Fq2_init(&c2tempYpMulInv);
		Fq2_init(&tempXp2);
		Fq2_init(&tempXr);
		Fq2_init(&tempYr);

		// auxiliary computation
		Fq2_add(&c2tempYpMulInv, &(p->valy), &(p->valy)); // 2yp
		Fq2_mul_inv(&c2tempYpMulInv, &c2tempYpMulInv);	  // (2yp)^(-1)

		Fq2_mul(&tempXp2, &(p->valx), &(p->valx));
		Fq2_add(&tempL, &tempXp2, &tempXp2);
		Fq2_add(&tempL, &tempL, &tempXp2); // 3xp^2

		Fq2_mul(&tempL, &tempL, &c2tempYpMulInv); // L = (3xp^2)/(2yp)

		Fq2_mul(&tempL2, &tempL, &tempL); // L^2

		Fq2_mul(&tempV, &tempL, &(p->valx));
		Fq2_sub(&tempV, &(p->valy), &tempV); // V = yp - L * xp

		// computation

		// xr
		Fq2_sub(&tempXr, &tempL2, &(p->valx));
		Fq2_sub(&tempXr, &tempXr, &(p->valx)); // xr = L^2 - 2xp

		// yr
		Fq2_mul(&tempYr, &tempL, &tempXr);
		Fq2_add(&tempYr, &tempYr, &tempV);
		Fq2_add_inv(&tempYr, &tempYr); // yr = -(L*xr + V)

		// storing
		Fq2_assign(&(rop->valx), &tempXr);
		Fq2_assign(&(rop->valy), &tempYr);
		rop->inf = '0'; // not infinity point

		// clear
		Fq2_clear(&tempL);
		Fq2_clear(&tempL2);
		Fq2_clear(&tempV);
		Fq2_clear(&c2tempYpMulInv);
		Fq2_clear(&tempXp2);
		Fq2_clear(&tempNull);
		Fq2_clear(&tempYpAddInv);
		Fq2_clear(&tempXr);
		Fq2_clear(&tempYr);

		return;
	}

	// Case 3: Q = -P. xp = xq, yq = -yp.
	if (Fq2_cmp(&(p->valx), &(q->valx)) == 0 &&
		Fq2_cmp(&(q->valy), &tempYpAddInv) == 0)
	{
		G2_set_inf(rop);

		Fq2_clear(&tempNull);
		Fq2_clear(&tempYpAddInv);

		return;
	}

	// Safeguard
	Fq2_clear(&tempNull);
	Fq2_clear(&tempYpAddInv);
}

void G2_mul(G2* rop, mpz_t n, G2* p) // rop = [n]P
{

	if (mpz_cmp_ui(n, 0) == 0)
	{
		G2_set_inf(rop);
		return;
	}

	// auxiliary init
	G2 temp;
	G2_init(&temp);
	G2_assign(&temp, p);

	char* bin = mpz_get_str(NULL, 2, n);
	int len = strlen(bin);

	for (int i = 1; i < len; i++)
	{
		G2_add(&temp, &temp, &temp);

		if (bin[i] == '1')
		{
			G2_add(&temp, &temp, p);
		}
	}

	G2_assign(rop, &temp);

	// clear
	G2_clear(&temp);
	free(bin);
}

// Group Arithmetic Projective
void G2_Z_add_inv(G2_Z* rop, G2_Z* a)
{
	// aux init
	Fq2 tempY;
	Fq2_init(&tempY);

	Fq2_add_inv(&tempY, &(a->valY));

	Fq2_assign(&(rop->valX), &(a->valX));
	Fq2_assign(&(rop->valY), &tempY);
	Fq2_assign(&(rop->valZ), &(a->valZ));

	// aux clear
	Fq2_clear(&tempY);
}

int G2_Z_cmp(G2_Z* p, G2_Z* q)
{
	//(Xp : Yp : Zp) = L * (Xq : Yq : Zq) <=>
	// <=> [XpYq - XqYp, YpZq - YqZp, ZpXq - ZqXp] = [0, 0, 0]

	int cmp = 0;

	// aux init
	Fq2 tempXpYq, tempXqYp, tempYpZq, tempYqZp, tempZpXq, tempZqXp;
	Fq2_init(&tempXpYq);
	Fq2_init(&tempXqYp);
	Fq2_init(&tempYpZq);
	Fq2_init(&tempYqZp);
	Fq2_init(&tempZpXq);
	Fq2_init(&tempZqXp);

	// computation
	Fq2_mul(&tempXpYq, &(p->valX), &(q->valY));
	Fq2_mul(&tempXqYp, &(p->valY), &(q->valX));
	Fq2_mul(&tempYpZq, &(p->valY), &(q->valZ));
	Fq2_mul(&tempYqZp, &(p->valZ), &(q->valY));
	Fq2_mul(&tempZpXq, &(p->valZ), &(q->valX));
	Fq2_mul(&tempZqXp, &(p->valX), &(q->valZ));

	if (Fq2_cmp(&tempXpYq, &tempXqYp) != 0)
	{
		cmp = 1;
	}

	if (Fq2_cmp(&tempYpZq, &tempYqZp) != 0)
	{
		cmp = 1;
	}

	if (Fq2_cmp(&tempZpXq, &tempZqXp) != 0)
	{
		cmp = 1;
	}

	// aux clear
	Fq2_clear(&tempXpYq);
	Fq2_clear(&tempXqYp);
	Fq2_clear(&tempYpZq);
	Fq2_clear(&tempYqZp);
	Fq2_clear(&tempZpXq);
	Fq2_clear(&tempZqXp);

	return cmp;
}

void G2_Z_add(G2_Z* rop, G2_Z* p, G2_Z* q)
{
	Fq2 zero;
	Fq2_init(&zero);
	Fq2_set_ui(&zero, 0, 0);

	// if P == O, P + Q = Q
	if (Fq2_cmp(&(p->valZ), &zero) == 0)
	{
		G2_Z_assign(rop, q);
		Fq2_clear(&zero);
		return;
	}

	// if Q == O, P + Q = P
	if (Fq2_cmp(&(q->valZ), &zero) == 0)
	{
		G2_Z_assign(rop, p);
		Fq2_clear(&zero);
		return;
	}

	// General case. P!=Q
	if (G2_Z_cmp(p, q) != 0)
	{
		/*

		Zr = ZpZq * (XpZq - XqZp)^3.

		Xr = (XpZq - XqZp) * [ZpZq * (YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq +
		XqZp)]

		Yr = ZpZq * (XqYp - XpYq) * (XpZq - XqZp)^2 - (YpZq - YqZp) * [ZpZq *
		(YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq + XqZp)]

		*/

		// aux init
		Fq2 tempXpYq, tempXqYp, tempYpZq, tempYqZp, tempXqZp, tempXpZq,
			tempZpZq, diff1, diff2, diff3, sum, delta1, delta2, aux1, aux2,
			bracket, c1, c2, c3;
		Fq2_init(&tempXpYq);
		Fq2_init(&tempXqYp);
		Fq2_init(&tempYpZq);
		Fq2_init(&tempYqZp);
		Fq2_init(&tempXqZp);
		Fq2_init(&tempXpZq);
		Fq2_init(&tempZpZq);
		Fq2_init(&diff1);
		Fq2_init(&diff2);
		Fq2_init(&diff3);
		Fq2_init(&sum);
		Fq2_init(&delta1);
		Fq2_init(&delta2);
		Fq2_init(&aux1);
		Fq2_init(&aux2);
		Fq2_init(&bracket);
		Fq2_init(&c1);
		Fq2_init(&c2);
		Fq2_init(&c3);

		// aux computation
		Fq2_mul(&tempXpYq, &(p->valX), &(q->valY));
		Fq2_mul(&tempXqYp, &(p->valY), &(q->valX));
		Fq2_mul(&tempYpZq, &(p->valY), &(q->valZ));
		Fq2_mul(&tempYqZp, &(p->valZ), &(q->valY));
		Fq2_mul(&tempXqZp, &(p->valZ), &(q->valX));
		Fq2_mul(&tempXpZq, &(p->valX), &(q->valZ));
		Fq2_mul(&tempZpZq, &(p->valZ), &(q->valZ));

		Fq2_sub(&diff1, &tempXpZq, &tempXqZp);
		Fq2_mul(&diff2, &diff1, &diff1);
		Fq2_mul(&diff3, &diff2, &diff1);

		Fq2_add(&sum, &tempXpZq, &tempXqZp);
		Fq2_sub(&delta1, &tempYpZq, &tempYqZp);
		Fq2_mul(&delta2, &delta1, &delta1);

		Fq2_mul(&aux1, &tempZpZq, &delta2);
		Fq2_mul(&aux2, &diff2, &sum);
		Fq2_sub(
			&bracket, &aux1,
			&aux2); //[ZpZq * (YpZq - YqZp)^2 - (XpZq - XqZp)^2 * (XpZq + XqZp)]

		// Zr
		Fq2_mul(&(rop->valZ), &tempZpZq, &diff3);

		// Xr
		Fq2_mul(&(rop->valX), &diff1, &bracket);

		// Yr
		Fq2_sub(&c3, &tempXqYp, &tempXpYq);
		Fq2_mul(&c1, &tempZpZq, &c3);
		Fq2_mul(&c1, &c1, &diff2);
		Fq2_mul(&c2, &delta1, &bracket);

		Fq2_sub(&(rop->valY), &c1, &c2);

		// aux clear
		Fq2_clear(&tempXpYq);
		Fq2_clear(&tempXqYp);
		Fq2_clear(&tempYpZq);
		Fq2_clear(&tempYqZp);
		Fq2_clear(&tempXqZp);
		Fq2_clear(&tempXpZq);
		Fq2_clear(&tempZpZq);
		Fq2_clear(&diff1);
		Fq2_clear(&diff2);
		Fq2_clear(&diff3);
		Fq2_clear(&sum);
		Fq2_clear(&delta1);
		Fq2_clear(&delta2);
		Fq2_clear(&aux1);
		Fq2_clear(&aux2);
		Fq2_clear(&bracket);
		Fq2_clear(&c1);
		Fq2_clear(&c2);
		Fq2_clear(&c3);
	}

	else // P == Q, doubling
	{
		/*

		Zr = 8Yp^3 Zp^3

		Xr = 18Xp^4YpZp - 16XpYp^3Zp^2

		Yr = -27Xp^6 + 36Xp^3Yp^2Zp - 8Yp^4Zp^2

		*/

		// aux init
		Fq2 c8, c16, c18, c27, c36, aux1, aux2, aux3, storeX, storeY, storeZ;
		Fq2_init(&c8);
		Fq2_init(&c16);
		Fq2_init(&c18);
		Fq2_init(&c27);
		Fq2_init(&c36);
		Fq2_init(&aux1);
		Fq2_init(&aux2);
		Fq2_init(&aux3);
		Fq2_init(&storeX);
		Fq2_init(&storeY);
		Fq2_init(&storeZ);

		// constants setting
		Fq2_set_ui(&c8, 8, 0);
		Fq2_set_ui(&c16, 16, 0);
		Fq2_set_ui(&c18, 18, 0);
		Fq2_set_ui(&c27, 27, 0);
		Fq2_set_ui(&c36, 36, 0);

		//---NB! aux1,2,3 are not fixed and may take different values
		// Zr
		Fq2_mul(&aux1, &(p->valY), &(p->valY));
		Fq2_mul(&aux1, &aux1, &(p->valY));
		Fq2_mul(&aux1, &aux1, &(p->valZ));
		Fq2_mul(&aux1, &aux1, &(p->valZ));
		Fq2_mul(&aux1, &aux1, &(p->valZ));
		Fq2_mul(&storeZ, &c8, &aux1);

		// Xr
		Fq2_mul(&aux1, &(p->valX), &(p->valX));
		Fq2_mul(&aux1, &aux1, &(p->valX));
		Fq2_mul(&aux1, &aux1, &(p->valX));
		Fq2_mul(&aux1, &aux1, &(p->valY));
		Fq2_mul(&aux1, &aux1, &(p->valZ));
		Fq2_mul(&aux1, &c18, &aux1);

		Fq2_mul(&aux2, &(p->valX), &(p->valY));
		Fq2_mul(&aux2, &aux2, &(p->valY));
		Fq2_mul(&aux2, &aux2, &(p->valY));
		Fq2_mul(&aux2, &aux2, &(p->valZ));
		Fq2_mul(&aux2, &aux2, &(p->valZ));
		Fq2_mul(&aux2, &c16, &aux2);

		Fq2_sub(&storeX, &aux1, &aux2);

		// Yr
		Fq2_mul(&aux1, &(p->valX), &(p->valX));
		Fq2_mul(&aux1, &aux1, &(p->valX));
		Fq2_mul(&aux1, &aux1, &(p->valX));
		Fq2_mul(&aux1, &aux1, &(p->valX));
		Fq2_mul(&aux1, &aux1, &(p->valX));
		Fq2_mul(&aux1, &c27, &aux1);

		Fq2_mul(&aux2, &(p->valX), &(p->valX));
		Fq2_mul(&aux2, &aux2, &(p->valX));
		Fq2_mul(&aux2, &aux2, &(p->valY));
		Fq2_mul(&aux2, &aux2, &(p->valY));
		Fq2_mul(&aux2, &aux2, &(p->valZ));
		Fq2_mul(&aux2, &c36, &aux2);

		Fq2_mul(&aux3, &(p->valY), &(p->valY));
		Fq2_mul(&aux3, &aux3, &(p->valY));
		Fq2_mul(&aux3, &aux3, &(p->valY));
		Fq2_mul(&aux3, &aux3, &(p->valZ));
		Fq2_mul(&aux3, &aux3, &(p->valZ));
		Fq2_mul(&aux3, &c8, &aux3);

		Fq2_sub(&aux2, &aux2, &aux1);
		Fq2_sub(&storeY, &aux2, &aux3);

		// storing the result
		Fq2_assign(&(rop->valX), &storeX);
		Fq2_assign(&(rop->valY), &storeY);
		Fq2_assign(&(rop->valZ), &storeZ);

		// aux clear
		Fq2_clear(&c8);
		Fq2_clear(&c16);
		Fq2_clear(&c18);
		Fq2_clear(&c27);
		Fq2_clear(&c36);
		Fq2_clear(&aux1);
		Fq2_clear(&aux2);
		Fq2_clear(&aux3);
		Fq2_clear(&storeX);
		Fq2_clear(&storeY);
		Fq2_clear(&storeZ);
	}

	// Safeguard
	Fq2_clear(&zero);
}

void G2_Z_mul(G2_Z* rop, mpz_t n, G2_Z* p)
{
	// auxiliary init
	G2_Z temp;
	G2_Z_init(&temp);
	G2_Z_assign(&temp, p);

	char* bin = mpz_get_str(NULL, 2, n);
	int len = strlen(bin);

	for (int i = 1; i < len; i++)
	{
		G2_Z_add(&temp, &temp, &temp);

		if (bin[i] == '1')
		{
			G2_Z_add(&temp, &temp, p);
		}
	}

	G2_Z_assign(rop, &temp);

	// clear
	G2_Z_clear(&temp);
	free(bin);
}

// Print
void G2_print(G2* a)
{
	if (a->inf == '1')
	{
		printf("inf");
	}

	if (a->inf == '0')
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
