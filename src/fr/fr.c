#include "fr.h"

Fr R_MODULUS;
char* bR_STR;
int bR_STR_len;

void Fr_field_init()
{
	Fr_set_hex_str(&R_MODULUS, "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001");
	bR_STR_len = 256;
}

uint8_t Fr_hex_to_nibble (char c)
{
	uint8_t result  = 0;
	uint8_t x = (uint8_t)c;
	
	//masking
	//digit
	uint8_t ge_0 = ((uint8_t)('0' - 1 - x)) >> 7;
	uint8_t le_9 = ((uint8_t)(x - '9' - 1)) >> 7;
	uint8_t mask_digit = ge_0 & le_9;
	
	//letter
	uint8_t ge_a = ((uint8_t)('a' - 1 - x)) >> 7;
	uint8_t le_f = ((uint8_t)(x - 'f' - 1)) >> 7;
	uint8_t mask_letter = ge_a & le_f;
	
	//Letter
	uint8_t ge_A = ((uint8_t)('A' - 1 - x)) >> 7;
	uint8_t le_F = ((uint8_t)(x - 'F' - 1)) >> 7;
	uint8_t mask_Letter = ge_A & le_F;
	
	result |= mask_digit * (uint8_t)(x - '0');
	result |= mask_letter * (uint8_t)(x - 'a' + 10);
	result |= mask_Letter * (uint8_t)(x - 'A' + 10);
	
	return result;
}

void Fr_set_hex_str(Fr* rop, char* hex_str) //NB! len_str <= LEN_MAX_R
{
	char buff[LEN_MAX_R];
	for(int i = 0; i < LEN_MAX_R; i++)
	{
		buff[i] = '0';
	}
	
	//evaluating len_str
	int len_str = 0;
	int check = 1;
	
	for(int i = 0; i < LEN_MAX_R; i++)
	{
		int mask = (hex_str[i] != '\0');
		
		check = check & mask; //1, until it's the end of the string; 0, afterwards 
		len_str = len_str + check; 
	}
	
	//filling the buffer
	
	for(int i = LEN_MAX_R - 1; i >= 0; i--)
	{
		int j = i - (LEN_MAX_R - len_str);
		uint32_t mask = (uint32_t)((int32_t)j >> 31); // 0xFFFFFFFF, if j < 0; 0, otherwise
		mask = ~mask; // 0xFFFFFFFF, if j >= 0; 0, otherwise
		
		buff[i] = (char)((hex_str[j] & mask) | (buff[i] & ~mask));
	}
	
	for(int ind_limb = FR_LIMBS - 1; ind_limb >= 0; ind_limb--)
	{
		int j = 8 * (FR_LIMBS - 1 - ind_limb);
		uint32_t val = 0;
		
		for(int i = j; i < j + 8; i++)
		{
			val = val << 4;
			uint32_t temp = (uint32_t)Fr_hex_to_nibble(buff[i]);
			val = val | temp;
		}
		
		rop->limbs[ind_limb] = val;
	}
	
}

void Fr_print(Fr* a)
{	
	for(int i = FR_LIMBS - 1; i>0; i--)
	{
		printf("%"PRIu32"|", a->limbs[i]);
	}
	
	printf("%"PRIu32, a->limbs[0]);
}

void Fr_print_dec(Fr* a)
{
	mpz_t rop, aux, base, deg;
	mpz_inits(rop, aux, base, deg, NULL);
	
	mpz_set_ui(base, 2);
	mpz_pow_ui(base, base, 32);
	
	for(int i = 0; i < FR_LIMBS; i++)
	{
		mpz_set_ui(aux, a->limbs[i]);
		mpz_pow_ui(deg, base, i);
		
		mpz_mul(aux, aux, deg);
		
		mpz_add(rop, rop, aux);
	}
	
	gmp_printf("%Zd", rop);
	
	mpz_clears(rop, aux, base, deg, NULL);
}
