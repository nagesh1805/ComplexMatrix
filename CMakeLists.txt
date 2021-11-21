cmake_minimum_required(VERSION 3.20)
project(ComplexMatrix)

set(CMAKE_CXX_STANDARD 14)

add_executable(ComplexMatrix main.cpp ComplexNumber.h ComplexVector.h RationalNumber.h Matrix.h SquareMatrix.h boolean.h ComplexNumber.cpp ComplexVector.cpp Matrix.cpp SquareMatrix.cpp RationalNumber.cpp)
