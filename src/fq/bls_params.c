#include "bls_params.h"

/*

Source: https://electriccoin.co/blog/new-snark-curve/

*/

const char* Q_STR = "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0"
					"f6b0f6241eabfffeb153ffffb9feffffffffaaab";

mpz_t Q_MODULUS;

void bls_params_init()
{
	mpz_init_set_str(Q_MODULUS, Q_STR, 0);
}

void bls_params_clear()
{
	mpz_clear(Q_MODULUS);
}
