#include "fq2.h"

// Init

void Fq2_set(Fq2* a, Fq* num0, Fq* num1)
{
	Fq_assign(&(a->val0), num0);
	Fq_assign(&(a->val1), num1);
}

void Fq2_set_ui(Fq2* a, uint32_t num0, uint32_t num1)
{
	Fq_set_ui(&(a->val0), num0);
	Fq_set_ui(&(a->val1), num1);
}

void Fq2_set_zero(Fq2* a)
{
	Fq_set_ui(&(a->val0), 0);
	Fq_set_ui(&(a->val1), 0);
}

void Fq2_set_id(Fq2* a)
{
	Fq_set_ui(&(a->val0), 1);
	Fq_set_ui(&(a->val1), 0);
}

void Fq2_assign(Fq2* rop, Fq2* a)
{
	Fq_assign(&(rop->val0), &(a->val0));
	Fq_assign(&(rop->val1), &(a->val1));
}

void Fq2_set_hex_str(Fq2* rop, char* hex_str0, char* hex_str1)
{
	Fq_set_hex_str(&(rop->val0), hex_str0);
	Fq_set_hex_str(&(rop->val1), hex_str1);
}


// Comparison

uint32_t Fq2_cmp(Fq2* a, Fq2* b)
{
	return (Fq_cmp(&(a->val0), &(b->val0)) & Fq_cmp(&(a->val1), &(b->val1)));
}


//Montgomery representation
void Fq2_mont_rep(Fq2* aR, Fq2* a)
{
	Fq_mont_rep(&(aR->val0), &(a->val0));
	Fq_mont_rep(&(aR->val1), &(a->val1));
}

void Fq2_mont_rep_inv(Fq2* a, Fq2* aR)
{
	Fq_mont_rep_inv(&(a->val0), &(aR->val0));
	Fq_mont_rep_inv(&(a->val1), &(aR->val1));
}



// Print
void Fq2_print(Fq2* a)
{
	printf("(");
	Fq_print(&(a->val0));
	printf(",");
	Fq_print(&(a->val1));
	printf(")");
}

void Fq2_print_dec(Fq2* a)
{
	printf("(");
	Fq_print_dec(&(a->val0));
	printf(", ");
	Fq_print_dec(&(a->val1));
	printf(")");
}

