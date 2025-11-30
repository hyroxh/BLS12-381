/*
 * For testing purposes, several helper methods are introduced. They are primarily
 * used to verify properties of the pairings.
 *
 * -- Function Descriptions ----------------------------------------------------
 *
 * Helper functions
 * ---------------
 * void tate_pairing_mul(Fq12* rop, mpz_t a, mpz_t b)
 *      Compute the Tate pairing e([a]genG1, [b]genG2) and store the result in *rop.
 *
 * void tate_power(Fq12* rop, mpz_t a, mpz_t b)
 *      Compute e(genG1, genG2) and raise it to the power a*b.
 *      This function is not used in practice.
 *
 * void ate_pairing_mul(Fq12* rop, mpz_t a, mpz_t b)
 *      Compute the ate pairing e([a]genG1, [b]genG2) and store the result in *rop.
 *
 * void ate_power(Fq12* rop, mpz_t a, mpz_t b)
 *      Compute e(genG1, genG2) and raise it to the power a*b.
 *      This function is not used in practice.
 *
 *
 * Testing functions
 * -----------------
 * int tate_test(mpz_t a, mpz_t b, int print)
 *      Compute the Tate pairings e([a]genG1, [b]genG2) and e([b]genG1, [a]genG2),
 *      and compare them.
 *      Return 1 if they are equal, and 0 otherwise.
 *      If 'print' is set to 1, the function prints e([a]genG1, [b]genG2).
 * int ate_test(mpz_t a, mpz_t b, int print)
 *      Compute the ate pairings e([a]genG1, [b]genG2) and e([b]genG1, [a]genG2),
 *      and compare them.
 *      Return 1 if they are equal, and 0 otherwise.
 *      If 'print' is set to 1, the function prints e([a]genG1, [b]genG2).
 *
 *
 * -- Note --------------------------------------------------------------------
 * There is no need to test underlying units, such as field arithmetic, individually. 
 * Because of the nature of the task, the probability that an error in one of the 
 * submodules would somehow cancel out and still produce a correct, non-trivial 
 * bilinear pairing is negligible.
 *
 *
*/ 



#include "../src/pairing/miller.h"
#include <stdio.h>

//-------------Tate pairing

void tate_pairing_mul(Fq12* rop, mpz_t a, mpz_t b)
{
	G1 P;
	G1_init(&P);

	G2 Q;
	G2_init(&Q);

	G1_mul(&P, a, &genG1);
	G2_mul(&Q, b, &genG2);

	tate_pairing(rop, &P, &Q);

	G1_clear(&P);
	G2_clear(&Q);
}

void tate_power(Fq12* rop, mpz_t a, mpz_t b)
{
	if (mpz_cmp_ui(a, 0) == 0 || mpz_cmp_ui(b, 0) == 0)
	{
		Fq12_id(rop);
		return;
	}

	Fq12 pairing;
	Fq12_init(&pairing);

	tate_pairing(&pairing, &genG1, &genG2);

	mpz_t c;
	mpz_init(c);
	mpz_mul(c, a, b);

	char* r = mpz_get_str(NULL, 2, c);
	int len = strlen(r);

	Fq12 aux;
	Fq12_init(&aux);
	Fq12_assign(&aux, &pairing);

	for (int i = 1; i < len; i++)
	{
		Fq12_mul(&aux, &aux, &aux);
		if (r[i] == '1')
		{
			Fq12_mul(&aux, &aux, &pairing);
		}
	}

	Fq12_assign(rop, &aux);

	Fq12_clear(&pairing);
	mpz_clear(c);
	Fq12_clear(&aux);
}

int tate_test(mpz_t a, mpz_t b, int print)
{
	int check = 1;
	Fq12 rop1, rop2;
	Fq12_init(&rop1);
	Fq12_init(&rop2);

	tate_pairing_mul(&rop1, a, b);
	tate_pairing_mul(&rop2, b, a);

	if (Fq12_cmp(&rop1, &rop2) != 0)
	{
		errno = ERANGE;
		perror("Range error occured");
		//exit(EXIT_FAILURE);
		check = 0;
	}

	if (print == 1)
	{
		Fq12_print(&rop1);
		printf("\n");
	}

	Fq12_clear(&rop1);
	Fq12_clear(&rop2);
	
	return check;
}


//-------------Ate pairing

void ate_pairing_mul(Fq12* rop, mpz_t a, mpz_t b)
{
	G1 P;
	G1_init(&P);

	G2 Q;
	G2_init(&Q);

	G1_mul(&P, a, &genG1);
	G2_mul(&Q, b, &genG2);

	ate_pairing(rop, &P, &Q);

	G1_clear(&P);
	G2_clear(&Q);
}

void ate_power(Fq12* rop, mpz_t a, mpz_t b)
{
	if (mpz_cmp_ui(a, 0) == 0 || mpz_cmp_ui(b, 0) == 0)
	{
		Fq12_id(rop);
		return;
	}

	Fq12 pairing;
	Fq12_init(&pairing);

	ate_pairing(&pairing, &genG1, &genG2);

	mpz_t c;
	mpz_init(c);
	mpz_mul(c, a, b);

	char* r = mpz_get_str(NULL, 2, c);
	int len = strlen(r);

	Fq12 aux;
	Fq12_init(&aux);
	Fq12_assign(&aux, &pairing);

	for (int i = 1; i < len; i++)
	{
		Fq12_mul(&aux, &aux, &aux);
		if (r[i] == '1')
		{
			Fq12_mul(&aux, &aux, &pairing);
		}
	}

	Fq12_assign(rop, &aux);

	Fq12_clear(&pairing);
	mpz_clear(c);
	Fq12_clear(&aux);
}

int ate_test(mpz_t a, mpz_t b, int print)
{
	int check = 1;
	Fq12 rop1, rop2;
	Fq12_init(&rop1);
	Fq12_init(&rop2);

	ate_pairing_mul(&rop1, a, b);
	ate_pairing_mul(&rop2, b, a);

	if (Fq12_cmp(&rop1, &rop2) != 0)
	{
		errno = ERANGE;
		perror("Range error occured");
		//exit(EXIT_FAILURE);
		check = 0;
	}

	if (print == 1)
	{
		Fq12_print(&rop1);
		printf("\n");
	}

	Fq12_clear(&rop1);
	Fq12_clear(&rop2);
	
	return check;
}

int main(void)
{
	//---init
	bls_params_init();
	G1_group_init();
	G2_group_init();
	sextic_precompute();
	tate_precompute();

	//----code
	
	mpz_t a, b;
	mpz_inits(a, b, NULL);
	
	//----Test1

	mpz_set_str(a, "52435875175126190479447740508185965837690938581184513", 10);
	mpz_set_str(b, "25256", 10);
	
	if(tate_test(a, b, 0) == 1 && ate_test(a, b, 0) == 1)
	{
		printf("Success!\n");
	}
	
	//----Test2
	
	mpz_set_str(a, "0", 10);
	mpz_set_str(b, "0", 10);
	
	if(tate_test(a, b, 0) == 1 && ate_test(a, b, 0) == 1)
	{
		printf("Success!\n");
	}
	
	//----Test3
	
	mpz_set_str(a, "3198563", 10);
	
	if(tate_test(a, orderG1, 0) == 1 && ate_test(a, orderG1, 0) == 1)
	{
		printf("Success!\n");
	}
	
	//----Test4
	
	mpz_set_str(b, "984765894268363", 10);
	
	if(tate_test(orderG2, b, 0) == 1 && ate_test(orderG2, b, 0) == 1)
	{
		printf("Success!\n");
	}
	
	//-------
	
	mpz_clears(a, b, NULL);
	
	//---clear
	bls_params_clear();
	G1_group_clear();
	G2_group_clear();
	sextic_clear();
	tate_clear();

	return 0;
}
