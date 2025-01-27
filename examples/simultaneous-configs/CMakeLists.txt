cmake_minimum_required(VERSION 3.14)

project(SymPhas_Implementation
	VERSION 1.0
	DESCRIPTION "Implementation of Symbolic Phase-Field Solver"
	LANGUAGES CXX CUDA)

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED True)
set(CMAKE_CUDA_ARCHITECTURES "75")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../" "out/build")

add_executable(symphas_impl "main.cu")
target_link_libraries(symphas_impl SymPhas)
set_target_properties(symphas_impl PROPERTIES CUDA_ARCHITECTURES "${SYMPHAS_CUDA_ARCHITECTURES}")

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	# target_compile_options(symphas_impl PRIVATE -v -fmacro-backtrace-limit=0)
	else()
		target_compile_options(symphas_impl PRIVATE -v -ftemplate-backtrace-limit=0)
	endif()
endif()
