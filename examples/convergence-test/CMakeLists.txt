cmake_minimum_required(VERSION 3.14)

project(SymPhas_Implementation 
	VERSION 1.0
	DESCRIPTION "Implementation of Symbolic Phase-Field Solver"
	LANGUAGES CXX)


add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../" "out/build")

add_executable(symphas_impl "main.cpp")
target_link_libraries(symphas_impl SymPhas)
target_compile_definitions(symphas_impl PUBLIC SOLVER_TYPE=SolverFT)

