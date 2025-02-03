cmake_minimum_required(VERSION 3.0.0)
project(float-analysis VERSION 1.0.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
include_directories(include)

# Add compiler optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

# Add executable
add_executable(float-analysis src/main.cpp)
add_executable(zstd-test src/zstd-test.cpp)
add_executable(op-test src/op-test.cpp)

# Link Zstandard library to your executable
target_link_libraries(float-analysis PRIVATE zstd)
target_link_libraries(zstd-test PRIVATE zstd)