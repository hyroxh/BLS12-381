cmake_minimum_required(VERSION 3.16)
project(BLS_Project C)

# Set C standard

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Set optimization flags (similar to -O2 -march=native)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -march=native")

# Include directories

include_directories(${PROJECT_SOURCE_DIR}/src)

# Source files

set(SOURCES
    src/fq/fq.c
    src/fq/fq2.c
    src/fq/fq6.c
    src/fq/fq12.c
    src/fr/fr.c
    src/groups/g1.c
    src/groups/g2.c
    src/pairing/miller.c
    tests/test.c
)

include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/Fq
    ${PROJECT_SOURCE_DIR}/src/Fr
    ${PROJECT_SOURCE_DIR}/src/groups
    ${PROJECT_SOURCE_DIR}/src/pairing
)

# Executable

add_executable(test_exec ${SOURCES})

# Link GMP manually

target_link_libraries(test_exec PRIVATE gmp)

# Set output directory for executable

set_target_properties(test_exec PROPERTIES
RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests
OUTPUT_NAME "test"
)

