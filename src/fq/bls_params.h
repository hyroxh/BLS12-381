#ifndef BLS12_381_PARAMS_H
#define BLS12_381_PARAMS_H

#include <gmp.h>

// Declare the modulus as an external constant.
// Its actual value will be defined in bls_params.c
extern mpz_t Q_MODULUS;

// The actual BLS12-381 q is very large (381 bits)

void bls_params_init();

void bls_params_clear();

#endif // BLS12_381_PARAMS_H
