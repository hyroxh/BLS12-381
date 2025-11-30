#include "miller.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* xPOW_STR;
char* bPOW_STR;

int len_bPOW_STR;

Fq2 sextic_const;

// init
void G12_init(G12* a)
{
	Fq12_init(&(a->valx));
	Fq12_init(&(a->valy));
}

void G12_clear(G12* a)
{
	Fq12_clear(&(a->valx));
	Fq12_clear(&(a->valy));
}

void G12_Z_init(G12_Z* a)
{
	Fq12_init(&(a->valX));
	Fq12_init(&(a->valY));
	Fq12_init(&(a->valZ));
}

void G12_Z_clear(G12_Z* a)
{
	Fq12_clear(&(a->valX));
	Fq12_clear(&(a->valY));
	Fq12_clear(&(a->valZ));
}

void tate_precompute()
{
	mpz_t pow;
	mpz_init(pow);

	mpz_pow_ui(pow, Q_MODULUS, 6);
	mpz_add_ui(pow, pow, 1);

	mpz_cdiv_q(pow, pow, orderG1);

	xPOW_STR = mpz_get_str(NULL, 16, pow);
	bPOW_STR = mpz_get_str(NULL, 2, pow);
	len_bPOW_STR = strlen(bPOW_STR);

	mpz_clear(pow);
}

void tate_clear()
{
	free(xPOW_STR);
	free(bPOW_STR);
}

void sextic_precompute()
{
	Fq2_init(&sextic_const);

	Fq c1, c2; // 1/2 and -1/2
	Fq_init(&c1);
	Fq_init(&c2);
	Fq_set_ui(&c1, 2);
	Fq_set_ui(&c2, 2);

	Fq_mul_inv(&c1, &c1);
	Fq_add_inv(&c2, &c1);

	Fq2_set(&sextic_const, &c1, &c2);

	Fq_clear(&c1);
	Fq_clear(&c2);
}

void sextic_clear()
{
	Fq2_clear(&sextic_const);
}

// sextic untwist
void sextic_untwist(G12* rop, G2* a) // map E(Fq2) -> E(Fq12)
{
	/*

	Map G2 -> G12
	(x', y') -> (x'/w^2, y'/w^3)

	x = x'/w^2 = (x' * (1/2 - 1/2 * u)) * v^2 + 0 * w
	y = y'/w^3 = 0 + (y' * (1/2 - 1/2 * u)) * v * w

	*/

	Fq2 tempx, tempy, zero2;
	Fq2_init(&tempx);
	Fq2_init(&tempy);
	Fq2_init(&zero2);

	Fq2_assign(&tempx, &(a->valx));
	Fq2_assign(&tempy, &(a->valy));
	Fq2_set_ui(&zero2, 0, 0);

	Fq2_mul(&tempx, &tempx, &sextic_const);
	Fq2_mul(&tempy, &tempy, &sextic_const);

	Fq6 ax, ay, zero6;
	Fq6_init(&ax);
	Fq6_init(&ay);
	Fq6_init(&zero6);

	Fq6_zero(&zero6);

	Fq6_set(&ax, &zero2, &zero2, &tempx);
	Fq6_set(&ay, &zero2, &tempy, &zero2);

	Fq12_set(&(rop->valx), &ax, &zero6);
	Fq12_set(&(rop->valy), &zero6, &ay);

	// clear
	Fq2_clear(&tempx);
	Fq2_clear(&tempy);
	Fq2_clear(&zero2);
	Fq6_clear(&ax);
	Fq6_clear(&ay);
	Fq6_clear(&zero6);
}

// sextic untwist
void sextic_untwist_project(G12_Z* rop, G2_Z* a) // map E(Fq2) -> E(Fq12)
{
	/*

	Map G2 -> G12
	(X' : Y' : Z') -> (X'/w^2 : Y'/w^3 : Z')

	x = (x' * (1/2 - 1/2 * u)) * v^2 + 0 * w, x' in Fq2
	y = 0 + (y' * (1/2 - 1/2 * u)) * v * w, y' in Fq2
	z = z'

	*/

	Fq2 tempx, tempy, zero2;
	Fq2_init(&tempx);
	Fq2_init(&tempy);
	Fq2_init(&zero2);

	Fq2_assign(&tempx, &(a->valX));
	Fq2_assign(&tempy, &(a->valY));
	Fq2_set_ui(&zero2, 0, 0);

	Fq2_mul(&tempx, &tempx, &sextic_const);
	Fq2_mul(&tempy, &tempy, &sextic_const);

	Fq6 ax, ay, zero6;
	Fq6_init(&ax);
	Fq6_init(&ay);
	Fq6_init(&zero6);

	Fq6_zero(&zero6);

	Fq6_set(&ax, &zero2, &zero2, &tempx);
	Fq6_set(&ay, &zero2, &tempy, &zero2);

	Fq12_set(&(rop->valX), &ax, &zero6);
	Fq12_set(&(rop->valY), &zero6, &ay);
	Fq2_to_Fq12(&(rop->valZ), &(a->valZ));

	// clear
	Fq2_clear(&tempx);
	Fq2_clear(&tempy);
	Fq2_clear(&zero2);
	Fq6_clear(&ax);
	Fq6_clear(&ay);
	Fq6_clear(&zero6);
}

//------------Tate pairing-------------------------

void tate_line_doubling(Fq12* rop, G1_Z* RZ, G12* Q) // Lr,r (Q)
{
	/*Computing L_R_R (Q)

	L_R_R (Q) = (2*YrZr^2) * Yq - 3*Xr^2Zr * Xq + 3Xr^3 - 2Yr^2Zr

	*/

	// init
	Fq zero;
	Fq_init(&zero);
	Fq_set_ui(&zero, 0);

	// aux init
	Fq aux1, aux2, aux3, aux4;
	Fq_init(&aux1);
	Fq_init(&aux2);
	Fq_init(&aux3);
	Fq_init(&aux4);

	// aux comp
	Fq_add(&aux1, &(RZ->valY), &(RZ->valY));
	Fq_mul(&aux1, &aux1, &(RZ->valZ));
	Fq_mul(&aux1, &aux1, &(RZ->valZ));

	Fq_add(&aux2, &(RZ->valX), &(RZ->valX));
	Fq_add(&aux2, &aux2, &(RZ->valX));
	Fq_mul(&aux2, &aux2, &(RZ->valX));
	Fq_mul(&aux2, &aux2, &(RZ->valZ));

	Fq_add(&aux3, &(RZ->valX), &(RZ->valX));
	Fq_add(&aux3, &aux3, &(RZ->valX));
	Fq_mul(&aux3, &aux3, &(RZ->valX));
	Fq_mul(&aux3, &aux3, &(RZ->valX));

	Fq_add(&aux4, &(RZ->valY), &(RZ->valY));
	Fq_mul(&aux4, &aux4, &(RZ->valY));
	Fq_mul(&aux4, &aux4, &(RZ->valZ));

	// temp init
	Fq2 temp1, temp2, temp3, temp4;
	Fq2_init(&temp1);
	Fq2_init(&temp2);
	Fq2_init(&temp3);
	Fq2_init(&temp4);

	Fq2_set(&temp1, &aux1, &zero);
	Fq2_set(&temp2, &aux2, &zero);
	Fq2_set(&temp3, &aux3, &zero);
	Fq2_set(&temp4, &aux4, &zero);

	Fq12 t1y, t2x, t3, t4, L;
	Fq12_init(&t1y);
	Fq12_init(&t2x);
	Fq12_init(&t3);
	Fq12_init(&t4);
	Fq12_init(&L);

	Fq2_to_Fq12(&t1y, &temp1);
	Fq2_to_Fq12(&t2x, &temp2);
	Fq2_to_Fq12(&t3, &temp3);
	Fq2_to_Fq12(&t4, &temp4);

	Fq12_mul(&t1y, &t1y, &(Q->valy));
	Fq12_mul(&t2x, &t2x, &(Q->valx));

	Fq12_sub(&L, &t1y, &t2x);
	Fq12_add(&L, &L, &t3);
	Fq12_sub(&L, &L, &t4);

	// assigning
	Fq12_assign(rop, &L);

	// clear
	Fq_clear(&zero);
	Fq_clear(&aux1);
	Fq_clear(&aux2);
	Fq_clear(&aux3);
	Fq_clear(&aux4);

	Fq2_clear(&temp1);
	Fq2_clear(&temp2);
	Fq2_clear(&temp3);
	Fq2_clear(&temp4);

	Fq12_clear(&t1y);
	Fq12_clear(&t2x);
	Fq12_clear(&t3);
	Fq12_clear(&t4);
	Fq12_clear(&L);
}

void tate_line_adding(Fq12* rop, G1_Z* RZ, G1_Z* PZ, G12* Q) // L_R,P (Q)
{
	/*

	Computing L_R,P (Q)

	L: (Xp*Zr - Xr*Zp) * Yq + (YrZp - YpZr) * Xq + (XrYp - XpYr)

	*/

	// aux init
	Fq aux1, aux2, aux3, aux4, aux5, aux6, c1, c2, c3, zero;
	Fq_init(&aux1);
	Fq_init(&aux2);
	Fq_init(&aux3);
	Fq_init(&aux4);
	Fq_init(&aux5);
	Fq_init(&aux6);
	Fq_init(&c1);
	Fq_init(&c2);
	Fq_init(&c3);
	Fq_init(&zero);

	Fq_set_ui(&zero, 0);

	// aux computation
	Fq_mul(&aux1, &(PZ->valX), &(RZ->valZ));
	Fq_mul(&aux2, &(RZ->valX), &(PZ->valZ));
	Fq_mul(&aux3, &(RZ->valY), &(PZ->valZ));
	Fq_mul(&aux4, &(PZ->valY), &(RZ->valZ));
	Fq_mul(&aux5, &(RZ->valX), &(PZ->valY));
	Fq_mul(&aux6, &(PZ->valX), &(RZ->valY));

	Fq_sub(&c1, &aux1, &aux2);
	Fq_sub(&c2, &aux3, &aux4);
	Fq_sub(&c3, &aux5, &aux6);

	// temp init and conversion
	Fq2 temp1, temp2, temp3;
	Fq2_init(&temp1);
	Fq2_init(&temp2);
	Fq2_init(&temp3);

	Fq2_set(&temp1, &c1, &zero);
	Fq2_set(&temp2, &c2, &zero);
	Fq2_set(&temp3, &c3, &zero);

	Fq12 t1, t2, t3, L;
	Fq12_init(&t1);
	Fq12_init(&t2);
	Fq12_init(&t3);
	Fq12_init(&L);

	Fq2_to_Fq12(&t1, &temp1);
	Fq2_to_Fq12(&t2, &temp2);
	Fq2_to_Fq12(&t3, &temp3);

	Fq12_mul(&t1, &t1, &(Q->valy));
	Fq12_mul(&t2, &t2, &(Q->valx));

	Fq12_add(&L, &t1, &t2);
	Fq12_add(&L, &L, &t3);

	// assigning
	Fq12_assign(rop, &L);

	// clearing
	Fq_clear(&aux1);
	Fq_clear(&aux2);
	Fq_clear(&aux3);
	Fq_clear(&aux4);
	Fq_clear(&aux5);
	Fq_clear(&aux6);
	Fq_clear(&c1);
	Fq_clear(&c2);
	Fq_clear(&c3);
	Fq_clear(&zero);

	Fq2_clear(&temp1);
	Fq2_clear(&temp2);
	Fq2_clear(&temp3);

	Fq12_clear(&t1);
	Fq12_clear(&t2);
	Fq12_clear(&t3);
	Fq12_clear(&L);
}

void tate(Fq12* rop, G1* P, G2* Q) // pre-exponent execution
{
	if (P->inf == '1' || Q->inf == '1')
	{
		Fq12_id(rop);
		return;
	}

	// init
	G1_Z PZ;
	G1_Z_init(&PZ);
	G1_to_project(&PZ, P);

	G1_Z RZ;
	G1_Z_init(&RZ);
	G1_Z_assign(&RZ, &PZ);

	G12 Quntwist;
	G12_init(&Quntwist);
	sextic_untwist(&Quntwist, Q);

	Fq12 f;
	Fq12_init(&f);
	Fq12_id(&f);

	Fq12 L;
	Fq12_init(&L);

	for (int i = 1; i < len_bORD_STR; i++)
	{

		tate_line_doubling(&L, &RZ, &Quntwist);
		G1_Z_add(&RZ, &RZ, &RZ);
		Fq12_mul(&f, &f, &f);
		Fq12_mul(&f, &f, &L);

		if (bORD_STR[i] == '1')
		{
			tate_line_adding(&L, &RZ, &PZ, &Quntwist);
			G1_Z_add(&RZ, &RZ, &PZ);
			Fq12_mul(&f, &f, &L);
		}
	}

	Fq12_assign(rop, &f);

	G1_Z_clear(&PZ);
	G1_Z_clear(&RZ);
	G12_clear(&Quntwist);
	Fq12_clear(&f);
	Fq12_clear(&L);
}

void tate_exp(Fq12* rop, Fq12* t)
{
	Fq6 val0, val1;
	Fq6_init(&val0);
	Fq6_init(&val1);

	Fq12 temp, inv, aux;
	Fq12_init(&temp);
	Fq12_init(&inv);
	Fq12_init(&aux);

	Fq6_assign(&val0, &(t->val0));
	Fq6_assign(&val1, &(t->val1));
	Fq6_add_inv(&val1, &val1);
	Fq12_mul_inv(&inv, t);

	Fq12_set(&temp, &val0, &val1);
	Fq12_mul(&temp, &temp, &inv);

	Fq12_assign(&aux, &temp);

	for (int i = 1; i < len_bPOW_STR; i++)
	{
		Fq12_mul(&aux, &aux, &aux);
		if (bPOW_STR[i] == '1')
		{
			Fq12_mul(&aux, &aux, &temp);
		}
	}

	Fq12_assign(rop, &aux);

	Fq6_clear(&val0);
	Fq6_clear(&val1);
	Fq12_clear(&temp);
	Fq12_clear(&inv);
	Fq12_clear(&aux);
}

void tate_pairing(Fq12* rop, G1* P, G2* Q)
{
	Fq12 temp;
	Fq12_init(&temp);

	tate(&temp, P, Q);
	tate_exp(&temp, &temp);

	Fq12_assign(rop, &temp);

	Fq12_clear(&temp);
}

//-------Ate pairing--------------------------------

void ate_line_doubling(Fq12* rop, G2_Z* RZ, G1* P) // Lr,r (P)
{
	/*

	Computing L_R_R (P)

	L_R_R (P) = (2*YrZr^2) * Yp - 3*Xr^2Zr * Xp + 3Xr^3 - 2Yr^2Zr ->

	-> (2*YrZr^2) * Yp * const * vw - 3*Xr^2Zr * Xp * const * v + 3Xr^3 * const
	- 2Yr^2Zr * const,

	where const = sextic_const = (1-u)/2.

	*/

	// init
	Fq zero;
	Fq_init(&zero);
	Fq_set_ui(&zero, 0);

	Fq2 P2x, P2y;
	Fq2_init(&P2x);
	Fq2_init(&P2y);
	Fq2_set(&P2x, &(P->valx), &zero);
	Fq2_set(&P2y, &(P->valy), &zero);

	// aux init
	Fq2 aux1, aux2, aux3, aux4;
	Fq2_init(&aux1);
	Fq2_init(&aux2);
	Fq2_init(&aux3);
	Fq2_init(&aux4);

	// aux comp
	Fq2_add(&aux1, &(RZ->valY), &(RZ->valY));
	Fq2_mul(&aux1, &aux1, &(RZ->valZ));
	Fq2_mul(&aux1, &aux1, &(RZ->valZ));

	Fq2_add(&aux2, &(RZ->valX), &(RZ->valX));
	Fq2_add(&aux2, &aux2, &(RZ->valX));
	Fq2_mul(&aux2, &aux2, &(RZ->valX));
	Fq2_mul(&aux2, &aux2, &(RZ->valZ));

	Fq2_add(&aux3, &(RZ->valX), &(RZ->valX));
	Fq2_add(&aux3, &aux3, &(RZ->valX));
	Fq2_mul(&aux3, &aux3, &(RZ->valX));
	Fq2_mul(&aux3, &aux3, &(RZ->valX));

	Fq2_add(&aux4, &(RZ->valY), &(RZ->valY));
	Fq2_mul(&aux4, &aux4, &(RZ->valY));
	Fq2_mul(&aux4, &aux4, &(RZ->valZ));

	Fq2_mul(&aux1, &aux1, &P2y);
	Fq2_mul(&aux2, &aux2, &P2x);
	Fq2_sub(&aux3, &aux3, &aux4);

	// Fq2_mul(&aux1, &aux1, &sextic_const); <--- not needed due to Frobenius E.
	// Fq2_mul(&aux2, &aux2, &sextic_const); 
	// Fq2_mul(&aux3, &aux3, &sextic_const);

	// untwisting

	Fq2 zero2;
	Fq2_init(&zero2);
	Fq2_set_ui(&zero2, 0, 0);

	Fq6 temp1, temp2, temp3, zero6;
	Fq6_init(&temp1);
	Fq6_init(&temp2);
	Fq6_init(&temp3);
	Fq6_init(&zero6);

	Fq6_zero(&zero6);

	Fq6_set(&temp1, &zero2, &aux1, &zero2);
	Fq6_set(&temp2, &zero2, &aux2, &zero2);
	Fq6_set(&temp3, &aux3, &zero2, &zero2);

	Fq12 ut1, ut2, ut3, L;
	Fq12_init(&ut1);
	Fq12_init(&ut2);
	Fq12_init(&ut3);
	Fq12_init(&L);

	Fq12_set(&ut1, &zero6, &temp1);
	Fq12_set(&ut2, &temp2, &zero6);
	Fq12_set(&ut3, &temp3, &zero6);

	Fq12_sub(&L, &ut1, &ut2);
	Fq12_add(&L, &L, &ut3);

	Fq12_assign(rop, &L);

	// clear
	Fq_clear(&zero);
	Fq2_clear(&P2x);
	Fq2_clear(&P2y);
	Fq2_clear(&aux1);
	Fq2_clear(&aux2);
	Fq2_clear(&aux3);
	Fq2_clear(&aux4);
	Fq2_clear(&zero2);
	Fq6_clear(&temp1);
	Fq6_clear(&temp2);
	Fq6_clear(&temp3);
	Fq6_clear(&zero6);
	Fq12_clear(&ut1);
	Fq12_clear(&ut2);
	Fq12_clear(&ut3);
	Fq12_clear(&L);
}

void ate_line_adding(Fq12* rop, G2_Z* RZ, G2_Z* QZ, G1* P) // L_R,Q (P)
{
	/*

	Computing L_R,Q (P)

	L: (Xq*Zr - Xr*Zq) * Yp + (YrZq - YqZr) * Xp + (XrYq - XqYr) ->

	(Xq*Zr - Xr*Zq) * const * Yp * v^2 + (YrZq - YqZr) * const * Xp * v * w +
	(XrYq - XqYr) * const * w,

	where const = (1-u)/2 = sextic_const

	*/

	// init
	Fq zero;
	Fq_init(&zero);
	Fq_set_ui(&zero, 0);

	Fq2 P2x, P2y;
	Fq2_init(&P2x);
	Fq2_init(&P2y);
	Fq2_set(&P2x, &(P->valx), &zero);
	Fq2_set(&P2y, &(P->valy), &zero);

	// aux init
	Fq2 aux1, aux2, aux3, aux4, aux5, aux6, c1, c2, c3;
	Fq2_init(&aux1);
	Fq2_init(&aux2);
	Fq2_init(&aux3);
	Fq2_init(&aux4);
	Fq2_init(&aux5);
	Fq2_init(&aux6);
	Fq2_init(&c1);
	Fq2_init(&c2);
	Fq2_init(&c3);

	// aux computation
	Fq2_mul(&aux1, &(QZ->valX), &(RZ->valZ));
	Fq2_mul(&aux2, &(RZ->valX), &(QZ->valZ));
	Fq2_mul(&aux3, &(RZ->valY), &(QZ->valZ));
	Fq2_mul(&aux4, &(QZ->valY), &(RZ->valZ));
	Fq2_mul(&aux5, &(RZ->valX), &(QZ->valY));
	Fq2_mul(&aux6, &(QZ->valX), &(RZ->valY));

	Fq2_sub(&c1, &aux1, &aux2);
	Fq2_sub(&c2, &aux3, &aux4);
	Fq2_sub(&c3, &aux5, &aux6);

	Fq2_mul(&c1, &c1, &P2y);
	Fq2_mul(&c2, &c2, &P2x);

	// Fq2_mul(&c1, &c1, &sextic_const); <--- not needed due to Frobenius E.
	// Fq2_mul(&c2, &c2, &sextic_const); 
	// Fq2_mul(&c3, &c3, &sextic_const);

	// untwisting
	Fq2 zero2;
	Fq2_init(&zero2);
	Fq2_set_ui(&zero2, 0, 0);

	Fq6 temp1, temp2, temp3, zero6;
	Fq6_init(&temp1);
	Fq6_init(&temp2);
	Fq6_init(&temp3);
	Fq6_init(&zero6);
	Fq6_zero(&zero6);

	Fq6_set(&temp1, &zero2, &zero2, &c1);
	Fq6_set(&temp2, &zero2, &c2, &zero2);
	Fq6_set(&temp3, &c3, &zero2, &zero2);

	Fq12 ut1, ut2, ut3, L;
	Fq12_init(&ut1);
	Fq12_init(&ut2);
	Fq12_init(&ut3);
	Fq12_init(&L);

	Fq12_set(&ut1, &temp1, &zero6);
	Fq12_set(&ut2, &zero6, &temp2);
	Fq12_set(&ut3, &zero6, &temp3);

	Fq12_add(&L, &ut1, &ut2);
	Fq12_add(&L, &L, &ut3);

	Fq12_assign(rop, &L);

	// clearing
	Fq_clear(&zero);
	Fq2_clear(&P2x);
	Fq2_clear(&P2y);
	Fq2_clear(&aux1);
	Fq2_clear(&aux2);
	Fq2_clear(&aux3);
	Fq2_clear(&aux4);
	Fq2_clear(&aux5);
	Fq2_clear(&aux6);
	Fq2_clear(&c1);
	Fq2_clear(&c2);
	Fq2_clear(&c3);
	Fq2_clear(&zero2);
	Fq6_clear(&temp1);
	Fq6_clear(&temp2);
	Fq6_clear(&temp3);
	Fq6_clear(&zero6);
	Fq12_clear(&ut1);
	Fq12_clear(&ut2);
	Fq12_clear(&ut3);
	Fq12_clear(&L);
}

void ate(Fq12* rop, G1* P, G2* Q) // pre-exponent execution
{

	if (P->inf == '1' || Q->inf == '1')
	{
		Fq12_id(rop);
		return;
	}

	G2_Z QZ;
	G2_Z_init(&QZ);
	G2_to_project(&QZ, Q);

	G2_Z RZ;
	G2_Z_init(&RZ);
	G2_Z_assign(&RZ, &QZ);

	Fq12 f;
	Fq12_init(&f);
	Fq12_id(&f);

	Fq12 L;
	Fq12_init(&L);

	for (int i = 1; i < len_bORD_STR; i++)
	{
		ate_line_doubling(&L, &RZ, P);
		G2_Z_add(&RZ, &RZ, &RZ);
		Fq12_mul(&f, &f, &f);
		Fq12_mul(&f, &f, &L);

		if (bORD_STR[i] == '1')
		{
			ate_line_adding(&L, &RZ, &QZ, P);
			G2_Z_add(&RZ, &RZ, &QZ);
			Fq12_mul(&f, &f, &L);
		}
	}

	Fq12_assign(rop, &f);

	G2_Z_clear(&QZ);
	G2_Z_clear(&RZ);
	Fq12_clear(&f);
	Fq12_clear(&L);
}

void ate_exp(Fq12* rop, Fq12* a)
{
	Fq6 val0, val1;
	Fq6_init(&val0);
	Fq6_init(&val1);

	Fq12 temp, inv, aux;
	Fq12_init(&temp);
	Fq12_init(&inv);
	Fq12_init(&aux);

	Fq6_assign(&val0, &(a->val0));
	Fq6_assign(&val1, &(a->val1));
	Fq6_add_inv(&val1, &val1);
	Fq12_mul_inv(&inv, a);

	Fq12_set(&temp, &val0, &val1);
	Fq12_mul(&temp, &temp, &inv);

	Fq12_assign(&aux, &temp);

	for (int i = 1; i < len_bPOW_STR; i++)
	{
		Fq12_mul(&aux, &aux, &aux);
		if (bPOW_STR[i] == '1')
		{
			Fq12_mul(&aux, &aux, &temp);
		}
	}

	Fq12_assign(rop, &aux);

	Fq6_clear(&val0);
	Fq6_clear(&val1);
	Fq12_clear(&temp);
	Fq12_clear(&inv);
	Fq12_clear(&aux);
}

void ate_pairing(Fq12* rop, G1* P, G2* Q)
{
	Fq12 temp;
	Fq12_init(&temp);

	ate(&temp, P, Q);
	ate_exp(&temp, &temp);

	Fq12_assign(rop, &temp);

	Fq12_clear(&temp);
}
