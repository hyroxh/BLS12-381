/*
* Field arithmetic over a given prime q. 
*
* ===CONTENTS=================================================================
*
* struct Fq - a struct to store a field element.
*
* -- Function Descriptions ----------------------------------------------------
*
* Initialisation functions
* ---------------
*
* void Fq_init(Fq* a), void Fq_clear(Fq* a) - init./clear the Fq element a.
* void Fq_set(Fq* a, mpz_t num) - set the value of a to num.
* void Fq_set_str(Fq* a, const char* str) - set the value a to str.
*      The base depends on the leading characters - 0x, 0b, or 0 for
*      hexadecimal, binary, or octal, respectively; decimal, otherwise.
*
* void Fq_set_ui(Fq* a, unsigned long ui) - set the value a to ui.
*      By construction, ui is guaranteed to be smaller than q.
*
* void Fq_assign(Fq* rop, Fq* a) - set the value of rop to a.
*
* Comparison functions
* ---------------
*
* int Fq_cmp(Fq* a, Fq* b) - return 0, if a == b, and +-1, otherwise.
*      Return values are chosen in accordance with GNU GMP.
*
* Arithmetic functions
* ---------------
*
* void Fq_add(Fq* rop, Fq* a, Fq* b) - set rop to a + b.
* void Fq_sub(Fq* rop, Fq* a, Fq* b) - set rop to a - b.
* void Fq_mul(Fq* rop, Fq* a, Fq* b) - set rop to a * b.
* void Fq_add_inv(Fq* rop, Fq* a) - set rop to (-a).
* void Fq_mul_inv(Fq* rop, Fq* a) - set rop to 1/a.
*      If a == 0, output an error message and terminate
*      the program.
*
* Print 
* ---------------
* void Fq_print(Fq* a) - print the field element a in the decimal base.
*
* -- Note --------------------------------------------------------------------
* The similar functions and notations are used for the field extensions. 
*
*/


#ifndef BLS12_381_Fq_H
#define BLS12_381_Fq_H

#include "bls_params.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <gmp.h>

typedef struct
{

	mpz_t val;

} Fq;

// init
void Fq_init(Fq* a);
void Fq_clear(Fq* a);
void Fq_set(Fq* a, mpz_t num);
void Fq_set_str(Fq* a, const char* str);
void Fq_set_ui(Fq* a, unsigned long ui);
void Fq_assign(Fq* rop, Fq* a);

// comparison
int Fq_cmp(Fq* a, Fq* b);

// operations
void Fq_add(Fq* rop, Fq* a, Fq* b);
void Fq_sub(Fq* rop, Fq* a, Fq* b);
void Fq_mul(Fq* rop, Fq* a, Fq* b);
void Fq_add_inv(Fq* rop, Fq* a);
void Fq_mul_inv(Fq* rop, Fq* a);

// printing
void Fq_print(Fq* a);

#endif
