cmake_minimum_required(VERSION 3.10)
project(bisam_v2)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

# Special settings to show full command line
set(CMAKE_VERBOSE_MAKEFILE ON)

# Find Rscript (more suitable for scripting than R)
find_program(RSCRIPT_EXECUTABLE Rscript)
if (NOT RSCRIPT_EXECUTABLE)
    message(FATAL_ERROR "Rscript executable not found")
endif ()

# Get R include directory
execute_process(
        COMMAND ${RSCRIPT_EXECUTABLE} -e "cat(R.home('include'))"
        OUTPUT_VARIABLE R_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get R library directory
execute_process(
        COMMAND ${RSCRIPT_EXECUTABLE} -e "cat(R.home('lib'))"
        OUTPUT_VARIABLE R_LIB_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get RcppArmadillo include directory
execute_process(
        COMMAND ${RSCRIPT_EXECUTABLE} -e "cat(system.file('include', package='RcppArmadillo'))"
        OUTPUT_VARIABLE RCPPARMADILLO_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get Rcpp include directory
execute_process(
        COMMAND ${RSCRIPT_EXECUTABLE} -e "cat(system.file('include', package='Rcpp'))"
        OUTPUT_VARIABLE RCPP_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Find regular Armadillo (might still be needed)
find_package(Armadillo REQUIRED)

# Print config details
message(STATUS "R include: ${R_INCLUDE_DIR}")
message(STATUS "R lib: ${R_LIB_DIR}")
message(STATUS "RcppArmadillo include: ${RCPPARMADILLO_INCLUDE_DIR}")
message(STATUS "Rcpp include: ${RCPP_INCLUDE_DIR}")
message(STATUS "Armadillo includes: ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "Armadillo libs: ${ARMADILLO_LIBRARIES}")

add_executable(bisam_v2 main.cpp
        controlled-simulation.h
        mombf-bridge.h)

# Include all necessary directories
target_include_directories(bisam_v2 PRIVATE
        ${R_INCLUDE_DIR}
        ${RCPPARMADILLO_INCLUDE_DIR}
        ${RCPP_INCLUDE_DIR}
        ${ARMADILLO_INCLUDE_DIRS}
)

# Link against libraries
target_link_libraries(bisam_v2 PRIVATE
        ${ARMADILLO_LIBRARIES}
        -L${R_LIB_DIR} -lR
)