#include "fq6.h"

// init
void Fq6_init(Fq6* a)
{
	Fq2_init(&(a->val0));
	Fq2_init(&(a->val1));
	Fq2_init(&(a->val2));
}

void Fq6_clear(Fq6* a)
{
	Fq2_clear(&(a->val0));
	Fq2_clear(&(a->val1));
	Fq2_clear(&(a->val2));
}

void Fq6_set(Fq6* a, Fq2* num0, Fq2* num1, Fq2* num2)
{
	Fq2_set(&(a->val0), &(num0->val0), &(num0->val1));
	Fq2_set(&(a->val1), &(num1->val0), &(num1->val1));
	Fq2_set(&(a->val2), &(num2->val0), &(num2->val1));
}

void Fq6_assign(Fq6* rop, Fq6* a)
{
	Fq2_assign(&(rop->val0), &(a->val0));
	Fq2_assign(&(rop->val1), &(a->val1));
	Fq2_assign(&(rop->val2), &(a->val2));
}

void Fq6_zero(Fq6* rop)
{
	Fq2 num0, num1, num2;
	Fq2_init(&num0);
	Fq2_init(&num1);
	Fq2_init(&num2);

	Fq2_set_ui(&num0, 0, 0);
	Fq2_set_ui(&num1, 0, 0);
	Fq2_set_ui(&num2, 0, 0);

	Fq6_set(rop, &num0, &num1, &num2);

	Fq2_clear(&num0);
	Fq2_clear(&num1);
	Fq2_clear(&num2);
}

void Fq6_id(Fq6* rop)
{
	Fq2 num0, num1, num2;
	Fq2_init(&num0);
	Fq2_init(&num1);
	Fq2_init(&num2);

	Fq2_set_ui(&num0, 1, 0);
	Fq2_set_ui(&num1, 0, 0);
	Fq2_set_ui(&num2, 0, 0);

	Fq6_set(rop, &num0, &num1, &num2);

	Fq2_clear(&num0);
	Fq2_clear(&num1);
	Fq2_clear(&num2);
}

// conversion
void Fq2_to_Fq6(Fq6* rop, Fq2* a)
{
	Fq2 zero2;
	Fq2_init(&zero2);
	Fq2_set_ui(&zero2, 0, 0);

	Fq6_set(rop, a, &zero2, &zero2);

	Fq2_clear(&zero2);
}

// comparison
int Fq6_cmp(Fq6* a, Fq6* b) // NB! return 0, if a=b; return 1, otherwise
{
	int cmp0, cmp1, cmp2;
	cmp0 = Fq2_cmp(&(a->val0), &(b->val0));
	cmp1 = Fq2_cmp(&(a->val1), &(b->val1));
	cmp2 = Fq2_cmp(&(a->val2), &(b->val2));

	if (cmp0 == 0 && cmp1 == 0 && cmp2 == 0)
	{
		return 0;
	}

	return 1;
}

// Arithmetic
void Fq6_add(Fq6* rop, Fq6* a, Fq6* b)
{
	Fq2_add(&(rop->val0), &(a->val0), &(b->val0));
	Fq2_add(&(rop->val1), &(a->val1), &(b->val1));
	Fq2_add(&(rop->val2), &(a->val2), &(b->val2));
}

void Fq6_sub(Fq6* rop, Fq6* a, Fq6* b)
{
	Fq2_sub(&(rop->val0), &(a->val0), &(b->val0));
	Fq2_sub(&(rop->val1), &(a->val1), &(b->val1));
	Fq2_sub(&(rop->val2), &(a->val2), &(b->val2));
}

void Fq6_mul(Fq6* rop, Fq6* a, Fq6* b)
{
	/*
	
	a*b = a0b0 + (a0b1 + a1b0)v + (a0b2 + a1b1 + a2b0)v^2 + (a1b2 + a2b1)v^3 +
		a2b2v^4 = [v^3 = u + 1] =

	= [a0b0 + (a1b2 + a2b1) + (a1b2 + a2b1)u] + [a0b1 + a1b0 + (a2b2) + (a2b2)u]v +
		[a0b2 + a1b1 + a2b0]v^2.
	
	*/
	
	// auxiliary init
	Fq2 temp00, temp01, temp02, temp10, temp11, temp12, temp20, temp21, temp22,
		u1, c0, c1, c2;
	Fq2_init(&temp00);
	Fq2_init(&temp01);
	Fq2_init(&temp02);
	Fq2_init(&temp10);
	Fq2_init(&temp11);
	Fq2_init(&temp12);
	Fq2_init(&temp20);
	Fq2_init(&temp21);
	Fq2_init(&temp22);
	Fq2_init(&u1);
	Fq2_init(&c0);
	Fq2_init(&c1);
	Fq2_init(&c2);

	// auxiliary computation
	Fq2_mul(&temp00, &(a->val0), &(b->val0));
	Fq2_mul(&temp01, &(a->val0), &(b->val1));
	Fq2_mul(&temp02, &(a->val0), &(b->val2));
	Fq2_mul(&temp10, &(a->val1), &(b->val0));
	Fq2_mul(&temp11, &(a->val1), &(b->val1));
	Fq2_mul(&temp12, &(a->val1), &(b->val2));
	Fq2_mul(&temp20, &(a->val2), &(b->val0));
	Fq2_mul(&temp21, &(a->val2), &(b->val1));
	Fq2_mul(&temp22, &(a->val2), &(b->val2));
	Fq2_set_str(&u1, "1", "1");

	// computation
	Fq2_add(&c0, &temp12, &temp21);
	Fq2_mul(&c0, &c0, &u1);
	Fq2_add(&c0, &c0, &temp00);

	Fq2_mul(&c1, &temp22, &u1);
	Fq2_add(&c1, &c1, &temp10);
	Fq2_add(&c1, &c1, &temp01);

	Fq2_add(&c2, &temp02, &temp11);
	Fq2_add(&c2, &c2, &temp20);

	// storing
	Fq2_assign(&(rop->val0), &c0);
	Fq2_assign(&(rop->val1), &c1);
	Fq2_assign(&(rop->val2), &c2);

	// clearing
	Fq2_clear(&temp00);
	Fq2_clear(&temp01);
	Fq2_clear(&temp02);
	Fq2_clear(&temp10);
	Fq2_clear(&temp11);
	Fq2_clear(&temp12);
	Fq2_clear(&temp20);
	Fq2_clear(&temp21);
	Fq2_clear(&temp22);
	Fq2_clear(&u1);
	Fq2_clear(&c0);
	Fq2_clear(&c1);
	Fq2_clear(&c2);
}

void Fq6_add_inv(Fq6* rop, Fq6* a)
{

	Fq2_add_inv(&(rop->val0), &(a->val0));
	Fq2_add_inv(&(rop->val1), &(a->val1));
	Fq2_add_inv(&(rop->val2), &(a->val2));
}

void Fq2_determinant(Fq2* rop, Fq6* a)
{
	/*
	
	det M = a0^3 - 3a0a1a2(u+1) + a1^3(u+1) + 2a2^3 u = Norm(a) = a * a^(q^2) *
		* a^(q^4).
	
	*/
	
	// auxiliary init
	Fq2 temp000, temp111, temp222, temp012, u1, c3temp012u1, temp111u1,
		c2temp222u, det;
	Fq2_init(&temp000);
	Fq2_init(&temp111);
	Fq2_init(&temp222);
	Fq2_init(&temp012);
	Fq2_init(&u1);
	Fq2_init(&c3temp012u1);
	Fq2_init(&temp111u1);
	Fq2_init(&c2temp222u);
	Fq2_init(&det);

	// auxiliary computation
	Fq2_mul(&temp000, &(a->val0), &(a->val0));
	Fq2_mul(&temp000, &temp000, &(a->val0));

	Fq2_mul(&temp111, &(a->val1), &(a->val1));
	Fq2_mul(&temp111, &temp111, &(a->val1));

	Fq2_mul(&temp222, &(a->val2), &(a->val2));
	Fq2_mul(&temp222, &temp222, &(a->val2));

	Fq2_mul(&temp012, &(a->val0), &(a->val1));
	Fq2_mul(&temp012, &temp012, &(a->val2));

	Fq2_set_str(&u1, "1", "1");

	Fq2_add(&c3temp012u1, &temp012, &temp012);
	Fq2_add(&c3temp012u1, &c3temp012u1, &temp012);
	Fq2_mul(&c3temp012u1, &c3temp012u1, &u1);

	Fq2_mul(&temp111u1, &temp111, &u1);

	Fq2_add(&c2temp222u, &temp222, &temp222);
	Fq2_transform(&c2temp222u, &c2temp222u);

	// computation
	Fq2_sub(&det, &temp000, &c3temp012u1);
	Fq2_add(&det, &det, &temp111u1);
	Fq2_add(&det, &det, &c2temp222u);

	// storing
	Fq2_assign(rop, &det);

	// clearing
	Fq2_clear(&temp000);
	Fq2_clear(&temp111);
	Fq2_clear(&temp222);
	Fq2_clear(&temp012);
	Fq2_clear(&u1);
	Fq2_clear(&c3temp012u1);
	Fq2_clear(&temp111u1);
	Fq2_clear(&c2temp222u);
	Fq2_clear(&det);
}

void Fq6_mul_inv(Fq6* rop, Fq6* a)
{
	/*
	
	a^(-1) = detM ^(-1) * ((a0^2 - a1a2(u+1)) + (a2^2(u+1) - a0a1)v + (a1^2 -
		- a0a2)v^2).
	
	*/
	
	// auxiliary init
	Fq2 temp00, temp01, temp02, temp11, temp12, temp22, u1, temp12u1, temp22u1,
		c0, c1, c2;
	Fq2_init(&temp00);
	Fq2_init(&temp01);
	Fq2_init(&temp02);
	Fq2_init(&temp11);
	Fq2_init(&temp12);
	Fq2_init(&temp22);
	Fq2_init(&u1);
	Fq2_init(&temp12u1);
	Fq2_init(&temp22u1);
	Fq2_init(&c0);
	Fq2_init(&c1);
	Fq2_init(&c2);

	// determinant
	Fq2 detM, detMinv;
	Fq2_init(&detM);
	Fq2_init(&detMinv);
	Fq2_determinant(&detM, a);
	Fq2_mul_inv(&detMinv, &detM);

	// auxiliary computation
	Fq2_mul(&temp00, &(a->val0), &(a->val0));
	Fq2_mul(&temp01, &(a->val0), &(a->val1));
	Fq2_mul(&temp02, &(a->val0), &(a->val2));
	Fq2_mul(&temp11, &(a->val1), &(a->val1));
	Fq2_mul(&temp12, &(a->val1), &(a->val2));
	Fq2_mul(&temp22, &(a->val2), &(a->val2));
	Fq2_set_str(&u1, "1", "1");
	Fq2_mul(&temp12u1, &temp12, &u1);
	Fq2_mul(&temp22u1, &temp22, &u1);

	// computatioon
	Fq2_sub(&c0, &temp00, &temp12u1);
	Fq2_sub(&c1, &temp22u1, &temp01);
	Fq2_sub(&c2, &temp11, &temp02);

	// normalization
	Fq2_mul(&c0, &c0, &detMinv);
	Fq2_mul(&c1, &c1, &detMinv);
	Fq2_mul(&c2, &c2, &detMinv);

	// storing
	Fq2_assign(&(rop->val0), &c0);
	Fq2_assign(&(rop->val1), &c1);
	Fq2_assign(&(rop->val2), &c2);

	// clearing
	Fq2_clear(&temp00);
	Fq2_clear(&temp01);
	Fq2_clear(&temp02);
	Fq2_clear(&temp11);
	Fq2_clear(&temp12);
	Fq2_clear(&temp22);
	Fq2_clear(&u1);
	Fq2_clear(&temp12u1);
	Fq2_clear(&temp22u1);
	Fq2_clear(&c0);
	Fq2_clear(&c1);
	Fq2_clear(&c2);
	Fq2_clear(&detM);
	Fq2_clear(&detMinv);
}

// Transformation
void Fq6_transform(Fq6* rop, Fq6* a) // rop <- a*v
{
	Fq2 temp0, temp1, u1, temp2u1;
	Fq2_init(&temp0);
	Fq2_init(&temp1);
	Fq2_init(&u1);
	Fq2_init(&temp2u1);

	Fq2_assign(&temp0, &(a->val0));
	Fq2_assign(&temp1, &(a->val1));
	Fq2_set_str(&u1, "1", "1");
	Fq2_mul(&temp2u1, &(a->val2), &u1);

	Fq2_assign(&(rop->val0), &temp2u1);
	Fq2_assign(&(rop->val1), &temp0);
	Fq2_assign(&(rop->val2), &temp1);

	Fq2_clear(&temp0);
	Fq2_clear(&temp1);
	Fq2_clear(&u1);
	Fq2_clear(&temp2u1);
}

// print
void Fq6_print(Fq6* rop)
{
	printf("[");
	Fq2_print(&(rop->val0));
	printf(", ");
	Fq2_print(&(rop->val1));
	printf(", ");
	Fq2_print(&(rop->val2));
	printf("]");
}
