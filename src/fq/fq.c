#include "fq.h"

// Initialization functions
void Fq_init(Fq* a)
{
	mpz_init(a->val);
}

void Fq_clear(Fq* a)
{
	mpz_clear(a->val);
}

void Fq_set(Fq* a, mpz_t num) // a <- num
{
	mpz_set(a->val, num);
}

void Fq_assign(Fq* rop, Fq* a)
{
	mpz_set(rop->val, a->val);
}

void Fq_set_str(Fq* a, const char* str)
{
	mpz_set_str(a->val, str, 0);
}

void Fq_set_ui(Fq* a, unsigned long ui)
{
	mpz_set_ui(a->val, ui);
}

// comparison
int Fq_cmp(Fq* a, Fq* b)
{
	return mpz_cmp(a->val, b->val);
}

// Operations
void Fq_add(Fq* rop, Fq* a, Fq* b)
{
	mpz_add(rop->val, a->val, b->val);
	mpz_mod(rop->val, rop->val, Q_MODULUS);
}

void Fq_sub(Fq* rop, Fq* a, Fq* b)
{
	mpz_sub(rop->val, a->val, b->val);
	mpz_mod(rop->val, rop->val, Q_MODULUS);
}

void Fq_mul(Fq* rop, Fq* a, Fq* b)
{
	mpz_mul(rop->val, a->val, b->val);
	mpz_mod(rop->val, rop->val, Q_MODULUS);
}

void Fq_add_inv(Fq* rop, Fq* a)
{
	mpz_sub(rop->val, Q_MODULUS, a->val);
	mpz_mod(rop->val, rop->val, Q_MODULUS);
}

void Fq_mul_inv(Fq* rop, Fq* a)
{
	if(mpz_cmp_ui(a->val, 0) == 0)
	{
		errno = EDOM;
		perror("Domain error occured");
		exit(EXIT_FAILURE);
	}
	mpz_t exp;
	mpz_init(exp);
	mpz_sub_ui(exp, Q_MODULUS, 2);

	mpz_powm_sec(rop->val, a->val, exp,
				 Q_MODULUS); // sec for cryptographic purposes

	mpz_clear(exp);
}

// Printing
void Fq_print(Fq* a)
{
	gmp_printf("%Zd", a->val);
}
