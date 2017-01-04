cmake_minimum_required (VERSION 2.8.11)

set(SRC "../lodepng/lodepng.h" "../lodepng/lodepng.cpp")

add_subdirectory(1-matrixless)
add_subdirectory(2-better-advection)
add_subdirectory(3-conjugate-gradients)
add_subdirectory(4-solid-boundaries)
add_subdirectory(5-curved-boundaries)
add_subdirectory(6-heat)
add_subdirectory(7-variable-density)
add_subdirectory(8-flip)