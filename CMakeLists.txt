cmake_minimum_required(VERSION 2.8)

project(prng)
set(CMAKE_C_FLAGS "-O3 -Wno-unused-result")

add_executable(chi chi.c prng.c prng.h)
target_link_libraries(chi m)

add_executable(countx countx.c prng.c prng.h)
target_link_libraries(countx m)

add_executable(freq freq.c prng.c prng.h)

add_executable(rngav rngav.c prng.c prng.h)
target_link_libraries(rngav m)
