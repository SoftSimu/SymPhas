cmake_minimum_required(VERSION 3.14)
project(MyFirstSim)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../" "out/build")
add_executable(symphas_impl "main.cpp")
target_link_libraries(symphas_impl SymPhas)