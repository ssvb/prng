cmake_minimum_required(VERSION 2.8)

project(prng)

add_executable(chi chi.c)
target_link_libraries(chi m)

add_executable(countx countx.c prng.c prng.h)
target_link_libraries(countx m)

add_executable(freq freq.c prng.c prng.h)

add_executable(rngav rngav.c prng.c prng.h)
target_link_libraries(rngav m)
