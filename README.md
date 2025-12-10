# TatePairing

- [Project overview](#project-overview)
- [Setup](#setup)
- [Build Instructions](#build)
- [Testing](#testing)
- [License](#license)

## Project overview

This repository provides a proof-of-concept implementation of the ate pairing on the BLS12-381 curve, utilizing constant-time operations for its evaluation. This implementation is for demonstration purposes only and not intended for use as a cryptographic library. (NB! GMP library is used for printing functions only as a matter of convenience. It is not utilized in computations).

Function specifications and API are fully documented in the source headers (`.h` files).  
Each function includes a comment block describing:

- Purpose
- Parameters
- Return values
- Special notes or behavior

Please refer to the headers for detailed specifications.

## Setup

You will need to install [GNU GMP](https://gmplib.org/) first. Then follow the standard 'cmake' procedure. For example:

## Build instructions

mkdir build
cd build
cmake ..
make
../tests/test

## Testing

See the '/test' directory for examples of tests. 

## License

MIT License
