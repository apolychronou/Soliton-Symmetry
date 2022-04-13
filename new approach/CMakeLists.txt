cmake_minimum_required(VERSION 3.7.2)
project(SolitonSymmetry
        DESCRIPTION "Solution of LL equation"
        LANGUAGES C)
add_executable(axialSym
                src/functions.c
                src/rk4.c
                src/axial_sym.c)

target_compile_features(axialSym PRIVATE)
target_compile_options(myApp PRIVATE -Wall -Wextra -Wpedantic)
set(CMAKE_C_COMPILER gcc)
