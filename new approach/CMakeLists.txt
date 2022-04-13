cmake_minimum_required(VERSION 3.7.2)
project(SolitonSymmetry
        DESCRIPTION "Solution of LL equation"
        LANGUAGES C)
add_executable(axialSym
                src/functions.c
                src/rk4.c
                src/axial_sym.c)
target_link_libraries(axialSym m)
target_compile_features(axialSym PRIVATE)
target_compile_options(axialSym PRIVATE  -O3 -x c -g -Wall)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_C_COMPILER gcc)
