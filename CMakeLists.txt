cmake_minimum_required(VERSION 3.3)
project(RichardsonMethod)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp Matrix.cpp Matrix.h Generator.cpp Generator.h RichardsonMethod.cpp RichardsonMethod.h Vector.cpp Vector.h)
add_executable(RichardsonMethod ${SOURCE_FILES})