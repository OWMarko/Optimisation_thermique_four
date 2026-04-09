cmake_minimum_required(VERSION 3.10)
project(Four)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /O2 /arch:AVX2")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -march=native")
endif()

set(EIGEN_INCLUDE_DIR "/usr/include/eigen3")
include_directories(${EIGEN_INCLUDE_DIR})
include_directories(include)

file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(${PROJECT_NAME} ${SOURCES})

message(STATUS "Fichiers sources détectés : ${SOURCES}")