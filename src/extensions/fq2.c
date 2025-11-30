#include "fq2.h"

// Init
void Fq2_init(Fq2* a)
{
	Fq_init(&(a->val0));
	Fq_init(&(a->val1));
}

void Fq2_clear(Fq2* a)
{
	Fq_clear(&(a->val0));
	Fq_clear(&(a->val1));
}

void Fq2_set(Fq2* a, Fq* num0, Fq* num1)
{
	Fq_set(&(a->val0), num0->val);
	Fq_set(&(a->val1), num1->val);
}

void Fq2_assign(Fq2* rop, Fq2* a)
{
	Fq_assign(&(rop->val0), &(a->val0));
	Fq_assign(&(rop->val1), &(a->val1));
}

void Fq2_set_str(Fq2* a, const char* str0, const char* str1)
{
	Fq_set_str(&(a->val0), str0);
	Fq_set_str(&(a->val1), str1);
}

void Fq2_set_ui(Fq2* a, unsigned long ui0, unsigned long ui1)
{
	Fq_set_ui(&(a->val0), ui0);
	Fq_set_ui(&(a->val1), ui1);
}

// Comparison
int Fq2_cmp(Fq2* a, Fq2* b)
{
	int cmp0, cmp1;
	cmp0 = Fq_cmp(&(a->val0), &(b->val0));
	cmp1 = Fq_cmp(&(a->val1), &(b->val1));

	if (cmp0 == 0 && cmp1 == 0)
	{
		return 0;
	}

	return 1;
}

// Arithmetic
void Fq2_add(Fq2* rop, Fq2* a, Fq2* b)
{
	Fq_add(&(rop->val0), &(a->val0), &(b->val0));
	Fq_add(&(rop->val1), &(a->val1), &(b->val1));
}

void Fq2_sub(Fq2* rop, Fq2* a, Fq2* b)
{
	Fq_sub(&(rop->val0), &(a->val0), &(b->val0));
	Fq_sub(&(rop->val1), &(a->val1), &(b->val1));
}

void Fq2_mul(Fq2* rop, Fq2* a, Fq2* b)
{
	/*
	
	a * b = (a0 * b0 - a1* b1) + (a0 * b1 + a1 * b0) * u
	
	*/
	
	Fq temp00, temp11, temp01, temp10;
	Fq_init(&temp00);
	Fq_init(&temp11);
	Fq_init(&temp01);
	Fq_init(&temp10);

	Fq_mul(&temp00, &(a->val0), &(b->val0));
	Fq_mul(&temp11, &(a->val1), &(b->val1));
	Fq_mul(&temp01, &(a->val0), &(b->val1));
	Fq_mul(&temp10, &(a->val1), &(b->val0));

	Fq_sub(&(rop->val0), &temp00, &temp11);
	Fq_add(&(rop->val1), &temp01, &temp10);

	Fq_clear(&temp00);
	Fq_clear(&temp11);
	Fq_clear(&temp01);
	Fq_clear(&temp10);
}

void Fq2_add_inv(Fq2* rop, Fq2* a)
{
	Fq_add_inv(&(rop->val0), &(a->val0));
	Fq_add_inv(&(rop->val1), &(a->val1));
}

void Fq2_mul_inv(Fq2* rop, Fq2* a)
{
	/*
	
	a^-1 = (a0 * a0 + a1 * a1)^(-1) * a0 + (a0 * a0 + a1 * a1)^(-1) * (-a1) * u.
	
	*/
	
	Fq temp00, temp11, abs, inv, neg1;
	Fq_init(&temp00);
	Fq_init(&temp11);
	Fq_init(&abs);
	Fq_init(&inv);
	Fq_init(&neg1);

	Fq_mul(&temp00, &(a->val0), &(a->val0));
	Fq_mul(&temp11, &(a->val1), &(a->val1));
	Fq_add(&abs, &temp00, &temp11);
	Fq_mul_inv(&inv, &abs);
	Fq_add_inv(&neg1, &(a->val1));

	Fq_mul(&(rop->val0), &inv, &(a->val0));
	Fq_mul(&(rop->val1), &inv, &neg1);

	Fq_clear(&temp00);
	Fq_clear(&temp11);
	Fq_clear(&abs);
	Fq_clear(&inv);
	Fq_clear(&neg1);
}

void Fq2_transform(Fq2* rop, Fq2* a)
{
	Fq neg1, temp;
	Fq_init(&neg1);
	Fq_init(&temp);

	Fq_assign(&temp, &(a->val0));
	Fq_add_inv(&neg1, &(a->val1));

	Fq_assign(&(rop->val0), &neg1);
	Fq_assign(&(rop->val1), &temp);

	Fq_clear(&neg1);
	Fq_clear(&temp);
}

// Print
void Fq2_print(Fq2* a)
{
	gmp_printf("(%Zd, %Zd)", a->val0.val, a->val1.val);
}
