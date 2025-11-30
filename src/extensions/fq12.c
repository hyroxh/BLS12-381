#include "fq12.h"

// Initialization
void Fq12_init(Fq12* a)
{
	Fq6_init(&(a->val0));
	Fq6_init(&(a->val1));
}

void Fq12_clear(Fq12* a)
{
	Fq6_clear(&(a->val0));
	Fq6_clear(&(a->val1));
}

void Fq12_set(Fq12* a, Fq6* num0, Fq6* num1)
{
	Fq6_assign(&(a->val0), num0);
	Fq6_assign(&(a->val1), num1);
}

void Fq12_assign(Fq12* rop, Fq12* a)
{
	Fq6_assign(&(rop->val0), &(a->val0));
	Fq6_assign(&(rop->val1), &(a->val1));
}

void Fq12_zero(Fq12* rop)
{
	Fq6 num0, num1;
	Fq6_init(&num0);
	Fq6_init(&num1);

	Fq6_zero(&num0);
	Fq6_zero(&num1);

	Fq12_set(rop, &num0, &num1);

	Fq6_clear(&num0);
	Fq6_clear(&num1);
}

void Fq12_id(Fq12* rop)
{
	Fq6 num0, num1;
	Fq6_init(&num0);
	Fq6_init(&num1);

	Fq6_id(&num0);
	Fq6_zero(&num1);

	Fq12_set(rop, &num0, &num1);

	Fq6_clear(&num0);
	Fq6_clear(&num1);
}

// conversion
void Fq6_to_Fq12(Fq12* rop, Fq6* a)
{
	Fq6 zero6;
	Fq6_init(&zero6);
	Fq6_zero(&zero6);

	Fq12_set(rop, a, &zero6);

	Fq6_clear(&zero6);
}
void Fq2_to_Fq12(Fq12* rop, Fq2* a)
{
	Fq6 aux;
	Fq6_init(&aux);

	Fq2_to_Fq6(&aux, a);
	Fq6_to_Fq12(rop, &aux);

	Fq6_clear(&aux);
}

// comparison
int Fq12_cmp(Fq12* a, Fq12* b)
{
	int cmp0, cmp1;
	cmp0 = Fq6_cmp(&(a->val0), &(b->val0));
	cmp1 = Fq6_cmp(&(a->val1), &(b->val1));

	if (cmp0 == 0 && cmp1 == 0)
	{
		return 0;
	}

	return 1;
}

// Arithmetic
void Fq12_add(Fq12* rop, Fq12* a, Fq12* b)
{
	Fq6_add(&(rop->val0), &(a->val0), &(b->val0));
	Fq6_add(&(rop->val1), &(a->val1), &(b->val1));
}

void Fq12_sub(Fq12* rop, Fq12* a, Fq12* b)
{
	Fq6_sub(&(rop->val0), &(a->val0), &(b->val0));
	Fq6_sub(&(rop->val1), &(a->val1), &(b->val1));
}

void Fq12_mul(Fq12* rop, Fq12* a, Fq12* b)
{
	/*
	
	a*b = a0b0 + (a0b1 + a1b0)w + a1b1 w^2 = [w^2 = v] =

		= (a0b0 + a1b1 * v) + (a0b1 + a1b0)w.
	
	*/
	
	// auxiliary init
	Fq6 temp00, temp01, temp10, temp11, temp11v, c0, c1;
	Fq6_init(&temp00);
	Fq6_init(&temp01);
	Fq6_init(&temp10);
	Fq6_init(&temp11);
	Fq6_init(&temp11v);
	Fq6_init(&c0);
	Fq6_init(&c1);

	// auxiliary computation
	Fq6_mul(&temp00, &(a->val0), &(b->val0));
	Fq6_mul(&temp01, &(a->val0), &(b->val1));
	Fq6_mul(&temp10, &(a->val1), &(b->val0));
	Fq6_mul(&temp11, &(a->val1), &(b->val1));
	Fq6_transform(&temp11v, &temp11);

	// computation
	Fq6_add(&c0, &temp00, &temp11v); // a0b0 + a1b1 * v
	Fq6_add(&c1, &temp01, &temp10);

	// storing
	Fq6_assign(&(rop->val0), &c0);
	Fq6_assign(&(rop->val1), &c1);

	// clear
	Fq6_clear(&temp00);
	Fq6_clear(&temp01);
	Fq6_clear(&temp10);
	Fq6_clear(&temp11);
	Fq6_clear(&temp11v);
	Fq6_clear(&c0);
	Fq6_clear(&c1);
}

void Fq12_add_inv(Fq12* rop, Fq12* a)
{
	Fq6_add_inv(&(rop->val0), &(a->val0));
	Fq6_add_inv(&(rop->val1), &(a->val1));
}

void Fq12_mul_inv(Fq12* rop, Fq12* a)
{
	/*
	
	a^(-1) = (a0^2 - a1^2 * v)^(-1) * a0 + (a0^2 - a1^2 * v)^(-1) * (-a1) * w.	
	
	*/
	
	// auxiliary init
	Fq6 temp00, temp11, temp11v, norm, norminv, neg1, c0, c1;
	Fq6_init(&temp00);
	Fq6_init(&temp11);
	Fq6_init(&temp11v);
	Fq6_init(&norm);
	Fq6_init(&norminv);
	Fq6_init(&neg1);
	Fq6_init(&c0);
	Fq6_init(&c1);

	// auxiliary computation
	Fq6_mul(&temp00, &(a->val0), &(a->val0));
	Fq6_mul(&temp11, &(a->val1), &(a->val1));
	Fq6_transform(&temp11v, &temp11);

	Fq6_sub(&norm, &temp00, &temp11v);
	Fq6_mul_inv(&norminv, &norm);

	Fq6_add_inv(&neg1, &(a->val1));

	// computation
	Fq6_mul(&c0, &norminv, &(a->val0));
	Fq6_mul(&c1, &norminv, &neg1);

	// storing
	Fq6_assign(&(rop->val0), &c0);
	Fq6_assign(&(rop->val1), &c1);

	// clear
	Fq6_clear(&temp00);
	Fq6_clear(&temp11);
	Fq6_clear(&temp11v);
	Fq6_clear(&norm);
	Fq6_clear(&norminv);
	Fq6_clear(&neg1);
	Fq6_clear(&c0);
	Fq6_clear(&c1);
}

void Fq12_print(Fq12* rop)
{
	printf("{");
	Fq6_print(&(rop->val0));
	printf(", ");
	Fq6_print(&(rop->val1));
	printf("}");
}
