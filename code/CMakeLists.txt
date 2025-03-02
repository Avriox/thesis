cmake_minimum_required(VERSION 3.10)
project(code)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

# Force compiler and print info
set(CMAKE_CXX_COMPILER "g++")
message(STATUS "Using compiler: ${CMAKE_CXX_COMPILER}")

# Special settings to show full command line
set(CMAKE_CXX_COMPILER_LAUNCHER "echo")  # Prints commands instead of executing
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(Armadillo REQUIRED)

# Print Armadillo config details
message(STATUS "Armadillo includes: ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "Armadillo libs: ${ARMADILLO_LIBRARIES}")

add_executable(code main.cpp)

target_include_directories(code PRIVATE ${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(code PRIVATE
        ${ARMADILLO_LIBRARIES}
        openblas
        lapack
        cblas
        gfortran
        m
)

# For even more verbosity during linking
set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_COMPILER} <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")