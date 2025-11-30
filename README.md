# TatePairing

- [Project overview](#project-overview)
- [Setup](#setup)
- [Build Instructions](#build)
- [Testing](#testing)
- [License](#license)

## Project overview

This repository provides a proof-of-concept implementation of the Tate and ate pairings for the BLS12-381 curve. These implementations are for demonstration purposes only and not intended for use as a cryptographic library.

This function specifications and API are fully documented in the source headers (`.h` files).  
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
./tests/test

## Testing

See the '/test' directory for examples of tests. 

## License

MIT License
