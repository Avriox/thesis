cmake_minimum_required(VERSION 3.10)
project(code)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

# Special settings to show full command line
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(Armadillo REQUIRED)

# Print Armadillo config details
message(STATUS "Armadillo includes: ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "Armadillo libs: ${ARMADILLO_LIBRARIES}")

add_executable(code main.cpp
        onion_monte_carlo.h
        monte_carlo.h
        onion.h
        timer.h
        onion_fixed.h)

target_include_directories(code PRIVATE ${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(code PRIVATE
        ${ARMADILLO_LIBRARIES}
)
