/*
* Field arithmetic over a given prime q. 
*
* ===CONTENTS=================================================================
*
* struct Fq - a struct to store a field element.
*	@ We use N_LIMBS number of limbs of the type uint32_t.
*       @ N_LIMBS = 12
*
* -- Function Descriptions ----------------------------------------------------
*
* Initialisation functions
* ---------------
*
* void Fq_set(Fq* a, mpz_t num) - set the value of a to num.
* void Fq_set_hex_str(Fq* a, const char* str) - set the value a to str.
*      @The base is assumed to be hexadecimal. Don't use '0x' as leading characters
*      @The length is limited to LEN_MAX = 96.
*      @Constant-time. Its processing doesn't depend on the actual string length. 
*
* void Fq_set_ui(Fq* a, uint32_t ui) - set LS limb of a to ui.
* void Fq_set_zero(Fq* a) - set a to 0.
* void Fq_set_id(Fq* a) - set a to 1.
* uint8_t Fq_hex_to_nibble (char c) - returns the decimal form of the hex char c.
*      This is a helper function.
*
* void Fq_assign(Fq* rop, Fq* a) - set the value of rop to a.
*
* Montgomery representation
* ---------------
* void Fq_mont_rep(Fq* aR, Fq* a) - transform a and set the result to aR.
* void Fq_mont_rep_inv(Fq* a, Fq* aR) - transform aR back and set the result to aR.
*
* Comparison functions
* ---------------
* uint32_t Fq_cmp(Fq* a, Fq* b) - return 1, if a == b, and 0, otherwise.
*
* Arithmetic functions
* ---------------
*
* void Fq_add(Fq* rop, Fq* a, Fq* b) - set rop to a + b.
* void Fq_sub(Fq* rop, Fq* a, Fq* b) - set rop to a - b.
* void Fq_mont_mul(Fq* rop, Fq* a, Fq* b) - set rop to aR * bR mod q = abR mod q.
* void Fq_add_inv(Fq* rop, Fq* a) - set rop to (-a).
* void Fq_mont_mul_inv(Fq* rop, Fq* a) - set rop to (1/a)R mod q.
*      If a == 0, the output is 0.
*
* Print 
* ---------------
* void Fq_print(Fq* a) - print the limbs of a in the 2^32 base.
* void Fq_print_dec(Fq* a) - print a in the decimal form
*      @Uses GNU GMP for convience 
*
* -- Note --------------------------------------------------------------------
* The similar functions and notations are used for the field extensions. 
*
*/


#ifndef BLS12_381_Fq_H
#define BLS12_381_Fq_H

#define N_LIMBS 12
#define LEN_MAX 96

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <inttypes.h>
#include <string.h>
#include <gmp.h> //to print in the decimal format only

typedef struct
{
	uint32_t limbs[N_LIMBS]; //Base = 2^32
} Fq;
extern uint32_t Q_LSLIMB_INV;
extern Fq Q_MODULUS;
//extern uint64_t BASE;
extern char* bQ_2_STR;
extern int bQ_2_STR_len;

void Fq_field_init();

// init
void Fq_set_ui(Fq* rop, uint32_t a);
void Fq_set_zero(Fq* a);
void Fq_set_id(Fq* a);
uint8_t Fq_hex_to_nibble (char c);
void Fq_set_hex_str(Fq* rop, char* hex_str);
void Fq_assign(Fq* rop, Fq* a);

static inline void Fq_correction(Fq* a) //if a==q, make a = 0; otherwise, leave untouched
{
	uint32_t mask = 1;
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		mask = mask & (a->limbs[i] == Q_MODULUS.limbs[i]);
	}
	
	mask = -mask;
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		a->limbs[i] = (a->limbs[i] & ~mask);
	}
}

//Montgomery representation
void Fq_mont_rep(Fq* aR, Fq* a);
void Fq_mont_rep_inv(Fq* a, Fq* aR);

//Comparison
uint32_t Fq_cmp(Fq* a, Fq* b); //returns 1, if a == b; 0, otherwise.

//Arithmetic operations
static inline void Fq_add(Fq* rop, Fq* a, Fq* b)
{
	Fq d, e;
	
	uint64_t c1 = 0;
	uint64_t c2 = 0;
	
	//d = a + b
	for(int i = 0; i < N_LIMBS; i++)
	{
		uint64_t sum = (uint64_t)a->limbs[i] + b->limbs[i] + c1;
		d.limbs[i] = (uint32_t)sum;
		c1 = sum >> 32;
	}
	
	//e = d - Q_MODULUS
	for(int i = 0; i < N_LIMBS; i++)
	{
		uint64_t diff = (uint64_t)d.limbs[i] - Q_MODULUS.limbs[i] - c2;
		e.limbs[i] = (uint32_t)diff;
		c2 = diff >> 63;		
	}
	
	uint32_t mask = -(uint32_t) (c2 == 0); 
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->limbs[i] = (e.limbs[i] & mask) | (d.limbs[i] & ~mask);
	}
	
	Fq_correction(rop);	
}

static inline void Fq_sub(Fq* rop, Fq* a, Fq* b)
{
	Fq d, e;
	
	uint64_t c1 = 0;
	uint64_t c2 = 0;
	
	//d = a - b
	for(int i = 0; i < N_LIMBS; i++)
	{
		uint64_t diff = (uint64_t)a->limbs[i] - b->limbs[i] - c1;
		d.limbs[i] = (uint32_t)diff;
		c1 = diff >> 63;
	}
	
	//e = d + Q_MODULUS
	for(int i = 0; i < N_LIMBS; i++)
	{
		uint64_t sum = (uint64_t)d.limbs[i] + Q_MODULUS.limbs[i] + c2;
		e.limbs[i] = (uint32_t)sum;
		c2 = sum >> 32;
	}
	
	uint32_t mask = -(uint32_t) (c1 == 1); //(c2 >= c1)
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->limbs[i] = (e.limbs[i] & mask) | (d.limbs[i] & ~mask);
	}
	
	Fq_correction(rop);	
}
static inline void Fq_add_inv(Fq* rop, Fq* a)
{
	Fq_sub(rop, &Q_MODULUS, a);
	Fq_correction(rop);
}
static inline void Fq_mont_mul(Fq* rop, Fq* a, Fq* b)
{
	uint32_t g = -Q_LSLIMB_INV;
	Fq d;
	Fq_set_zero(&d);
	uint32_t dh = 0;
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		uint32_t f = (d.limbs[0] + a->limbs[i] * b->limbs[0]) * g;
		
		uint64_t c = 0;
		
		for(int j = 0; j < N_LIMBS; j++)
		{
			__uint128_t z = (__uint128_t)d.limbs[j] + (__uint128_t)a->limbs[i] * b->limbs[j] + (__uint128_t)f * Q_MODULUS.limbs[j] + (__uint128_t)c;
			if(j > 0)
			{
				d.limbs[j-1] = (uint32_t)z;
			}
			c = (uint64_t)(z >> 32);
		} 
		
		uint64_t z = (uint64_t)dh + c;
		d.limbs[N_LIMBS - 1] = (uint32_t)z;
		dh = (uint32_t)(z >> 32); 
	}
	
	Fq d_sub_Q;
	uint64_t borrow = 0;
	for(int i = 0; i < N_LIMBS; i++)
	{
		uint64_t diff = (uint64_t)d.limbs[i] - Q_MODULUS.limbs[i] - borrow;
		d_sub_Q.limbs[i] = (uint32_t)diff;
		borrow = diff >> 63;
	}
	
	//if dh >= borrow, it's d_sub_Q; it's d, otherwise
	uint32_t mask = -(uint32_t) (dh >= borrow);
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->limbs[i] = (d_sub_Q.limbs[i] & mask) | (d.limbs[i] & ~mask);
	}
	
	Fq_correction(rop);
}

static inline void Fq_mont_mul_inv(Fq* rop, Fq* a)
{
	Fq aux;
	Fq_assign(&aux, a);
	
	for(int i = 1; i < bQ_2_STR_len; i++)
	{
		Fq_mont_mul(&aux, &aux, &aux);
		
		if(bQ_2_STR[i] == '1')
		{
			Fq_mont_mul(&aux, &aux, a);
		}
	}
	
	Fq_assign(rop, &aux);	
}

// printing
void Fq_print(Fq* a);
void Fq_print_dec(Fq* a);

#endif
