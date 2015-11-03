cmake_minimum_required(VERSION 3.3)
project(MinErrMethod)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp Matrix.cpp Matrix.h Generator.cpp Generator.h MinimalErrorsMethod.cpp MinimalErrorsMethod.h)
add_executable(MinErrMethod ${SOURCE_FILES})