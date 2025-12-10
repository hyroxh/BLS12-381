#include "fq12.h"

// Initialization

void Fq12_set(Fq12* a, Fq6* num0, Fq6* num1)
{
	Fq6_assign(&(a->val0), num0);
	Fq6_assign(&(a->val1), num1);
}

void Fq12_set_zero(Fq12* rop)
{
	Fq6 num0, num1;

	Fq6_set_zero(&num0);
	Fq6_set_zero(&num1);

	Fq12_set(rop, &num0, &num1);
}

void Fq12_set_id(Fq12* rop)
{
	Fq6 num0, num1;

	Fq6_set_id(&num0);
	Fq6_set_zero(&num1);

	Fq12_set(rop, &num0, &num1);
}

void Fq12_assign(Fq12* rop, Fq12* a)
{
	Fq6_assign(&(rop->val0), &(a->val0));
	Fq6_assign(&(rop->val1), &(a->val1));
}

// conversion
void Fq6_to_Fq12(Fq12* rop, Fq6* a)
{
	Fq6 zero6;
	Fq6_set_zero(&zero6);

	Fq12_set(rop, a, &zero6);
}
void Fq2_to_Fq12(Fq12* rop, Fq2* a)
{
	Fq6 aux;

	Fq2_to_Fq6(&aux, a);
	Fq6_to_Fq12(rop, &aux);
}

// Montgomery representation
void Fq12_mont_rep(Fq12* aR, Fq12* a)
{
	Fq6_mont_rep(&(aR->val0), &(a->val0));
	Fq6_mont_rep(&(aR->val1), &(a->val1));
}
void Fq12_mont_rep_inv(Fq12* a, Fq12* aR)
{
	Fq6_mont_rep_inv(&(a->val0), &(aR->val0));
	Fq6_mont_rep_inv(&(a->val1), &(aR->val1));
}

uint32_t Fq12_cmp(Fq12* a, Fq12* b)
{
	return Fq6_cmp(&(a->val0), &(b->val0)) & Fq6_cmp(&(a->val1), &(b->val1));
}

void Fq12_print(Fq12* rop)
{
	printf("{");
	Fq6_print(&(rop->val0));
	printf(", ");
	Fq6_print(&(rop->val1));
	printf("}");
}

void Fq12_print_dec(Fq12* rop)
{
	printf("{");
	Fq6_print_dec(&(rop->val0));
	printf(", ");
	Fq6_print_dec(&(rop->val1));
	printf("}");
}
