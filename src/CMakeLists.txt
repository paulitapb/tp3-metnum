project(tp3)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8...3.13)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -std=c++17 -Wall -Wpedantic -Wextra -Wshadow -Wpointer-arith -Wcast-qual")
set(SOURCE_FILES src/main.cpp src/utils.cpp src/eigen.cpp)

add_executable(tp3 ${SOURCE_FILES} src/main/cpp/main.cpp)