
# ***************************************************************************
# This file is part of the SymPhas library, a framework for implementing
# solvers for phase-field problems with compile-time symbolic algebra.
# 
# Copyright (c) 2018-2021 by Steven A. Silber and Mikko Karttunen
# 
# SymPhas is free software, which can be redistributed or modified under
# the terms of the GNU Lesser General Public License (LGPL) as published
# by the Free Software Foundation; LGPL version 3, or later versions at
# your choice.
# 
# SymPhas is distributed with the faith that it will be helpful and
# practical but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
# General Public License for more details.
# 
# ***************************************************************************


find_package(Threads REQUIRED)
find_package(FFTW3)
find_package(TBB)
find_package(MPI)

include_directories(SYSTEM ${FFTW_INCLUDE_PATH})

file(GLOB_RECURSE SOURCES_LIB "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(THREAD_COUNT 2)
set(VECTOR_TYPE_NAME VectorValue)
set(ENABLE_TIMINGS ON CACHE BOOL "enable collecting the elapsed time information of the simulation")
set(ENABLE_INLINE_TIMINGS OFF CACHE BOOL "enable collecting the elapsed time of inline timings")

#add_library(fftw3_imported SHARED IMPORTED)
#set_target_properties(fftw3_imported PROPERTIES
#	IMPORTED_LOCATION "${FFTW_DIR}/bin/fftw3.dll"
#	IMPORTED_IMPLIB "${FFTW3_LIBRARY_DIRS}/fftw3.lib")

configure_file("inc/names.h.in" "${CMAKE_SOURCE_DIR}/inc/names.h")
	
if(NOT WIN32)
	list(APPEND ADD_LINUX_LIBS m ${CMAKE_THREAD_LIBS_INIT})
endif()

add_library(symphas_lib_o OBJECT ${SOURCES_LIB})
target_include_directories(symphas_lib_o PRIVATE ${FFTW3_INCLUDE_DIRS} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
	$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/inc>
	$<INSTALL_INTERFACE:include>)
target_link_directories(symphas_lib_o PUBLIC ${FFTW3_LIBRARY_DIRS})
target_link_libraries(symphas_lib_o PUBLIC ${ADD_LINUX_LIBS} PRIVATE ${FFTW3_LIBRARIES})
target_compile_definitions(symphas_lib_o PRIVATE LIB_EXPORTS)

if(ENABLE_TIMINGS)
	target_compile_definitions(symphas_lib_o PUBLIC PRINT_TIMINGS)
endif()

if(FFTW3_FOUND)
	target_compile_definitions(symphas_lib_o PUBLIC USING_FFTW)
endif()

if(TBB_FOUND)
	target_link_libraries(symphas_lib_o PUBLIC TBB::tbb)
	message(STATUS "using TBB v${TBB_VERSION}")
endif()

if(MPI_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_FLAGS} ${MPI_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_EXE_LINKER_FLAGS}")
	target_compile_definitions(symphas_lib_o PUBLIC USING_MPI)
	target_link_libraries(symphas_lib_o PRIVATE ${MPI_CXX_LIBRARIES})
	target_include_directories(symphas_lib_o PUBLIC ${MPI_CXX_INCLUDE_DIRS})
	message(STATUS "using MPI v${MPI_CXX_VERSION}")
endif()

if(ENABLE_INLINE_TIMINGS)
	target_compile_definitions(symphas_lib_o PUBLIC PRINT_TIMINGS_INLINE)
endif()

message(STATUS "'CMake Build Type' = ${CMAKE_BUILD_TYPE}")
if("${CMAKE_BUILD_TYPE}" MATCHES "Debug" OR "${CMAKE_CONFIGURATION_TYPES}" MATCHES "Debug")
	target_compile_definitions(symphas_lib_o PUBLIC DEBUG)
	message(STATUS "using Debug mode")
endif()

if(MSVC)
	find_package(OpenMP)
	if(OpenMP_CXX_FOUND)
		target_link_libraries(symphas_lib_o PUBLIC OpenMP::OpenMP_CXX)
	endif()
	target_compile_options(symphas_lib_o PUBLIC /openmp:experimental)
else()
	find_package(OpenMP)
	if(OpenMP_CXX_FOUND)
		target_link_libraries(symphas_lib_o PUBLIC OpenMP::OpenMP_CXX)
		if (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
			target_compile_options(symphas_lib_o PUBLIC "-fopenmp")
		endif()
	endif()
endif()

if(OpenMP_CXX_FOUND)
	message(STATUS "Using OpenMP " ${OpenMP_CXX_VERSION})
else()
	message(STATUS "No OpenMP")
endif()

if(MSVC)
	target_compile_options(symphas_lib_o PUBLIC /bigobj -F85894592)
else()
	target_compile_options(symphas_lib_o PUBLIC -Wno-format-security)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	target_compile_options(symphas_lib_o PUBLIC 
		-Wno-unused-variable -Wno-unused-function -ftime-trace -mavx)
endif()

if(WIN32)
	target_compile_definitions(symphas_lib_o PUBLIC _CRT_SECURE_NO_WARNINGS)
endif()


add_library(symphas_lib SHARED)
target_link_libraries(symphas_lib PUBLIC symphas_lib_o)



