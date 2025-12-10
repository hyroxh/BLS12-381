#include "fq6.h"

// init
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

void Fq6_set_zero(Fq6* rop)
{
	Fq2 zero2;
	Fq2_set_zero(&zero2);
	Fq6_set(rop, &zero2, &zero2, &zero2);
}

void Fq6_set_id(Fq6* rop)
{
	Fq2 id2, zero2;

	Fq2_set_id(&id2);
	Fq2_set_zero(&zero2);

	Fq6_set(rop, &id2, &zero2, &zero2);
}

// conversion
void Fq2_to_Fq6(Fq6* rop, Fq2* a)
{
	Fq2 zero2;
	Fq2_set_zero(&zero2);

	Fq6_set(rop, a, &zero2, &zero2);
}

// Montgomery representation
void Fq6_mont_rep(Fq6* aR, Fq6* a)
{
	Fq2_mont_rep(&(aR->val0), &(a->val0));
	Fq2_mont_rep(&(aR->val1), &(a->val1));
	Fq2_mont_rep(&(aR->val2), &(a->val2));
}

void Fq6_mont_rep_inv(Fq6* a, Fq6* aR)
{
	Fq2_mont_rep_inv(&(a->val0), &(aR->val0));
	Fq2_mont_rep_inv(&(a->val1), &(aR->val1));
	Fq2_mont_rep_inv(&(a->val2), &(aR->val2));
}

/*
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
*/

uint32_t Fq6_cmp(Fq6* a, Fq6* b)
{
	return Fq2_cmp(&(a->val0), &(b->val0)) & Fq2_cmp(&(a->val1), &(b->val1)) & Fq2_cmp(&(a->val2), &(b->val2));
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

void Fq6_print_dec(Fq6* rop)
{
	printf("[");
	Fq2_print_dec(&(rop->val0));
	printf(", ");
	Fq2_print_dec(&(rop->val1));
	printf(", ");
	Fq2_print_dec(&(rop->val2));
	printf("]");
}

