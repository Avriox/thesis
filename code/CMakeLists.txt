cmake_minimum_required(VERSION 3.10)
project(code)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

# Special settings to show full command line
set(CMAKE_VERBOSE_MAKEFILE ON)

# Explicitly find BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(Armadillo REQUIRED)

# Print Armadillo config details
message(STATUS "Armadillo includes: ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "Armadillo libs: ${ARMADILLO_LIBRARIES}")
message(STATUS "BLAS libs: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK libs: ${LAPACK_LIBRARIES}")


add_executable(code main.cpp
        timer.h
        onion_fixed.h
        isserlis_theorem.h
        isserlis_theorem.h)

target_include_directories(code PRIVATE ${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(code PRIVATE
        ${ARMADILLO_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
)