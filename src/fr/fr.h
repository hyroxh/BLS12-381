/*
* Field of scalars. According to the setting, the field is defined over r,
*	where r = orderG1 = orderG2
*
* The functions are similar to those of Fq, adjusted to the size of r.
*
*/


#ifndef FR_H
#define FR_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <inttypes.h>
#include <string.h>
#include <gmp.h>

#define FR_LIMBS 8
#define LEN_MAX_R 64

typedef struct
{
	uint32_t limbs[FR_LIMBS];
} Fr;

extern Fr R_MODULUS;
extern char* bR_STR;
extern int bR_STR_len;

void Fr_field_init();
uint8_t Fr_hex_to_nibble (char c);
void Fr_set_hex_str(Fr* rop, char* hex_str);
void Fr_print(Fr* a);
void Fr_print_dec(Fr* a);

#endif
