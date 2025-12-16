/*
 * For testing purposes, several helper methods are introduced. They are primarily
 * used to verify properties of the pairings.
 *
 * -- Function Descriptions ----------------------------------------------------
 *
 * Helper functions
 * ---------------
 * void ate_pairing_mul(Fq12* rop, Fr* a, Fr* b)
 *      Compute the ate pairing e([a]genG1, [b]genG2) and store the result in *rop.
 *
 *
 * Testing functions
 * -----------------
 * int ate_test(Fr* a, Fr* b, int print)
 *      @Compute the ate pairings e([a]genG1, [b]genG2) and e([b]genG1, [a]genG2),
 *      and compare them.
 *      @Return 1 if they are equal, and 0 otherwise.
 *      @If 'print' is set to 1, and the return value is 1, the function 
 *      prints e([a]genG1, [b]genG2).
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

//-------------Ate pairing

void ate_pairing_mul(Fq12* rop, Fr* a, Fr* b)
{
	G1_Z PZ;
	G1_to_project(&PZ, &genG1);
	G1_Z_mont_rep(&PZ, &PZ);
	
	G1_Z_mul(&PZ, a, &PZ);
	
	G1 P;
	G1_to_affine(&P, &PZ);

	G2_Z QZ;
	G2_to_project(&QZ, &genG2);
	G2_Z_mont_rep(&QZ, &QZ);
	
	G2_Z_mul(&QZ, b, &QZ);

	ate_pairing(rop, &P, &QZ);
	
	Fq12_mont_rep_inv(rop, rop);
}

int ate_test(Fr* a, Fr* b, int print)
{
	int check = 1;
	Fq12 rop1, rop2;

	ate_pairing_mul(&rop1, a, b);
	ate_pairing_mul(&rop2, b, a);
	
	
	if (Fq12_cmp(&rop1, &rop2) != 1)
	{
		errno = ERANGE;
		perror("Range error occured");
		//exit(EXIT_FAILURE);
		check = 0;
	}
	
	
	if (print == 1)
	{
		Fq12_print_dec(&rop1);
		printf("\n");
	}
	
	return check;
}

int main(void)
{
	//---init
	Fq_field_init();
	Fr_field_init();
	G1_group_init();
	G2_group_init();
	exp_precompute();
	ate_precompute();
	
	
	//----code
	
	Fr a, b;
	
	
	//----Test1

	Fr_set_hex_str(&a, "1aaaaaaaaaaffffffffbf");
	Fr_set_hex_str(&b, "12345ccccddf");
	
	if(ate_test(&a, &b, 1) == 1)
	{
		printf("Test 1:\tSuccess!\n");
		printf("----------\n");
	}
	
	//----Test2
	
	Fr_set_hex_str(&a, "0");
	Fr_set_hex_str(&b, "0");
	
	if(ate_test(&a, &b, 0) == 1)
	{
		printf("Test 2:\tSuccess!\n");
		printf("----------\n");
	}
	
	//----Test3
	
	Fr_set_hex_str(&a, "aaaaaaaaaa33333333333333333333");
	
	if(ate_test(&a, &orderG1, 0) == 1)
	{
		printf("Test 3:\tSuccess!\n");
		printf("----------\n");
	}
	
	//----Test4
	
	Fr_set_hex_str(&b, "ffffffffffffffffffff7777");
	
	if(ate_test(&orderG2, &b, 0) == 1)
	{
		printf("Test 4:\tSuccess!\n");
		printf("----------\n");
	}
	
	//-------
	
	return 0;
}
