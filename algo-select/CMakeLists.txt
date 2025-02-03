cmake_minimum_required(VERSION 3.5.0)
project(algo-select VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(LibPressio REQUIRED)

# Add compiler optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

get_filename_component(DATADIR "${CMAKE_CURRENT_SOURCE_DIR}/../huffman/real-datasets/" ABSOLUTE)

add_executable(algo-select src/main.cpp src/precise_time.cpp)
target_include_directories(algo-select PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../float-analysis/include")
target_link_libraries(algo-select PRIVATE LibPressio::libpressio zstd)
target_compile_definitions(algo-select PRIVATE DATADIR="${DATADIR}/")