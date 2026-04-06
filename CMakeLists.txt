cmake_minimum_required(VERSION 3.10)
project(Four_Thermique)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

set(EIGEN_INCLUDE_DIR "/usr/include/eigen3")
include_directories(${EIGEN_INCLUDE_DIR})

add_executable(Four_Thermique 
    src/main.cpp 
    src/utils_mef.cpp
)