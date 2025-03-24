cmake_minimum_required(VERSION 3.10)
project(bisam_v2)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

# Special settings to show full command line
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(Armadillo REQUIRED)

# Print Armadillo config details
message(STATUS "Armadillo includes: ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "Armadillo libs: ${ARMADILLO_LIBRARIES}")

add_executable(bisam_v2 main.cpp
        controlled-simulation.h
        mombf-bridge.h)

target_include_directories(bisam_v2 PRIVATE ${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(bisam_v2 PRIVATE
        ${ARMADILLO_LIBRARIES}
)
