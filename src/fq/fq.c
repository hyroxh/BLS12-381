#include "fq.h"

#define N_LIMBS 12
#define LEN_MAX 96
//LEN_MAX = 8 * N_LIMBS

Fq Q_MODULUS;
uint32_t Q_LSLIMB_INV;
char* bQ_2_STR;
int bQ_2_STR_len;

//uint64_t BASE = 4294967296;

void Fq_field_init()
{
	Fq_set_hex_str(&Q_MODULUS, "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab");
	Q_LSLIMB_INV = 196611; // Q_LSLIMB_INV = (Q_MODULUS.limbs[0]) ^ (-1) mod 2^32
	bQ_2_STR = "110100000000100010001111010100011100101111111111"
	"0011010011010010010110001101110100111101101100100001101001011"
	"10101100110101110110010001110111010010111000010011110011100001"
	"01000100101011111101100111001100001101001010100000111101101011"
	"00001111011000100100000111101010101111111111111111101011000101"
	"01001111111111111111111011100111111110111111111111111111111111"
	"111111111010101010101001"; //Q - 2
	bQ_2_STR_len = 381;
}

//Init 

void Fq_set_ui(Fq* rop, uint32_t a)
{
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->limbs[i] = 0;
	}
	
	rop->limbs[0] = a;
}

void Fq_set_zero(Fq* a)
{
	for(int i = 0; i < N_LIMBS; i++)
	{
		a->limbs[i] = 0;
	}
}

void Fq_set_id(Fq* a)
{
	for(int i = 0; i < N_LIMBS; i++)
	{
		a->limbs[i] = 0;
	}
	
	a->limbs[0] = 1;
}

void Fq_assign(Fq* rop, Fq* a)
{
	for(int i = 0; i < N_LIMBS; i++)
	{
		rop->limbs[i] = a->limbs[i];
	}
}
/*
void Fq_correction(Fq* a) //if a==q, make a = 0; otherwise, leave untouched
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
*/
uint8_t Fq_hex_to_nibble (char c)
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

void Fq_set_hex_str(Fq* rop, char* hex_str) //NB! len_str <= LEN_MAX
{
	char buff[LEN_MAX];
	for(int i = 0; i < LEN_MAX; i++)
	{
		buff[i] = '0';
	}
	
	//evaluating len_str
	int len_str = 0;
	int check = 1;
	
	for(int i = 0; i < LEN_MAX; i++)
	{
		int mask = (hex_str[i] != '\0');
		
		check = check & mask; //1, until it's the end of the string; 0, afterwards 
		len_str = len_str + check; 
	}
	
	//filling the buffer
	
	for(int i = LEN_MAX - 1; i >= 0; i--)
	{
		int j = i - (LEN_MAX - len_str);
		uint32_t mask = (uint32_t)((int32_t)j >> 31); // 0xFFFFFFFF, if j < 0; 0, otherwise
		mask = ~mask; // 0xFFFFFFFF, if j >= 0; 0, otherwise
		
		buff[i] = (char)((hex_str[j] & mask) | (buff[i] & ~mask));
	}
	
	for(int ind_limb = N_LIMBS - 1; ind_limb >= 0; ind_limb--)
	{
		int j = 8 * (N_LIMBS - 1 - ind_limb);
		uint32_t val = 0;
		
		for(int i = j; i < j + 8; i++)
		{
			val = val << 4;
			uint32_t temp = (uint32_t)Fq_hex_to_nibble(buff[i]);
			val = val | temp;
		}
		
		rop->limbs[ind_limb] = val;
	}
	
}

//Montgomery representation

void Fq_mont_rep(Fq* aR, Fq* a) //aR mod q, where R = W^N, with W = 2^32.
{
	Fq aux;
	Fq_assign(&aux, a);
	
	for(int i = 0; i < 32 * N_LIMBS; i++)
	{
		Fq_add(&aux, &aux, &aux);
	}
	
	Fq_assign(aR, &aux);
	Fq_correction(aR);
}



void Fq_mont_rep_inv(Fq* a, Fq* aR)
{
	Fq id;
	Fq_set_id(&id);
	
	Fq_mont_mul(a, aR, &id);
	Fq_correction(a);
}

uint32_t Fq_cmp(Fq* a, Fq* b)
{
	uint32_t mask = 1;
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		mask = mask & (a->limbs[i] == b->limbs[i]);
	}
	
	return mask;
}

//Arithmetic

void Fq_mul(Fq* rop, Fq* a, Fq* b) //NB! Not adivsed to use. Use Mont_mul instead
{
	Fq aR, bR, cR;
	Fq_mont_rep(&aR, a);
	Fq_mont_rep(&bR, b);
	
	Fq_mont_mul(&cR, &aR, &bR);
	
	Fq_mont_rep_inv(rop, &cR);
}

//print
void Fq_print(Fq* a) //not const time
{
	int leading_limb = N_LIMBS - 1;
	while(leading_limb > 0)
	{
		if(a->limbs[leading_limb] == 0)
		{
			leading_limb--;
		}
		else
		{
			break;
		}
	}
	
	for(int i = leading_limb; i>0; i--)
	{
		printf("%"PRIu32"|", a->limbs[i]);
	}
	
	printf("%"PRIu32, a->limbs[0]);
}

void Fq_print_dec(Fq* a)
{
	mpz_t rop, aux, base, deg;
	mpz_inits(rop, aux, base, deg, NULL);
	
	mpz_set_ui(base, 2);
	mpz_pow_ui(base, base, 32);
	
	for(int i = 0; i < N_LIMBS; i++)
	{
		mpz_set_ui(aux, a->limbs[i]);
		mpz_pow_ui(deg, base, i);
		
		mpz_mul(aux, aux, deg);
		
		mpz_add(rop, rop, aux);
	}
	
	gmp_printf("%Zd", rop);
	
	mpz_clears(rop, aux, base, deg, NULL);
}

