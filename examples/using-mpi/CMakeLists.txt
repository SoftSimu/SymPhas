cmake_minimum_required(VERSION 3.14)

project(SymPhas_Implementation 
	VERSION 1.0
	DESCRIPTION "Testing MPI in SymPhas"
	LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../" "out/build")

add_executable(symphas_impl "main.cpp")
target_link_libraries(symphas_impl SymPhas)

if (WIN32)
	target_compile_options(symphas_impl PRIVATE -v -fmacro-backtrace-limit=0)
else()
	target_compile_options(symphas_impl PRIVATE -v -ftemplate-backtrace-limit=0)
endif()
