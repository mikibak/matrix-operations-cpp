cmake_minimum_required(VERSION 3.25)
project(matrix_operations_cpp)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(matrix_operations_cpp
        Matrix.cpp
        Matrix_Operations.cpp)
