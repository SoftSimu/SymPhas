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

include_directories(SYSTEM ${FFTW_INCLUDE_PATH})

file(GLOB_RECURSE SOURCES_LIB "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
	set(DISABLE_EXCEPTIONS ON)
else()
	unset(DISABLE_EXCEPTIONS)
endif()

set(THREAD_COUNT 2)
set(VECTOR_TYPE_NAME VectorValue)
set(ENABLE_TIMINGS ON CACHE BOOL "enable collecting the elapsed time information of the simulation")
set(ENABLE_INLINE_TIMINGS OFF CACHE BOOL "enable collecting the elapsed time of inline timings")

# add_library(fftw3_imported SHARED IMPORTED)
# set_target_properties(fftw3_imported PROPERTIES
# IMPORTED_LOCATION "${FFTW_DIR}/bin/fftw3.dll"
# IMPORTED_IMPLIB "${FFTW3_LIBRARY_DIRS}/fftw3.lib")
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
	set(PRINT_TIMINGS ON)
endif()

if(FFTW3_FOUND)
	set(USING_FFTW ON)
else()
	unset(USING_FFTW)
endif()

if(TBB_FOUND)
	target_link_libraries(symphas_lib_o PUBLIC TBB::tbb)
	message(STATUS "using TBB v${TBB_VERSION}")
endif()

if(USING_MPI)
	find_package(MPI)

	if(MPI_FOUND)
		set(CMAKE_C_FLAGS "${CMAKE_FLAGS} ${MPI_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_EXE_LINKER_FLAGS}")

		set(USING_MPI ON)
		target_link_libraries(symphas_lib_o PRIVATE ${MPI_CXX_LIBRARIES})
		target_include_directories(symphas_lib_o PUBLIC ${MPI_CXX_INCLUDE_DIRS})
		message(STATUS "using MPI v${MPI_CXX_VERSION}")
	else()
		message(STATUS "MPI found but not enabled")
	endif()
else()
	unset(USING_MPI)
	message(STATUS "MPI not found")
endif()

if(CMAKE_CUDA_COMPILER)
	if(USING_CUDA)
		find_package(CUDA REQUIRED)
		target_include_directories(symphas_lib_o PUBLIC ${CUDA_INCLUDE_DIRS})

		message(STATUS "using CUDA architectures: '${CMAKE_CUDA_ARCHITECTURES}'")
		message(STATUS "using CUDA v${CUDA_VERSION}")
	else()
		message(STATUS "CUDA found but not enabled")
	endif()
else()
	unset(USING_CUDA)
	message(STATUS "CUDA not found")
endif()

if(ENABLE_INLINE_TIMINGS)
	set(PRINT_TIMINGS_INLINE ON)
endif()

message(STATUS "'CMake Build Type' = ${CMAKE_BUILD_TYPE} and compiler = '${CMAKE_CXX_COMPILER_ID}'")

if("${CMAKE_BUILD_TYPE}" MATCHES "Debug" OR
	("${CMAKE_CONFIGURATION_TYPES}" MATCHES "Debug" AND NOT "${CMAKE_CONFIGURATION_TYPES}" MATCHES "Release"))
	set(DEBUG ON)
	message(STATUS "using Debug mode")
else()
	set(DEBUG OFF)
endif()

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
	message(STATUS "Using OpenMP " ${OpenMP_CXX_VERSION})
	target_link_libraries(symphas_lib_o PUBLIC OpenMP::OpenMP_CXX)

	if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
		if(WIN32)
			target_compile_options(symphas_lib_o PUBLIC "-openmp")
		else()
			target_compile_options(symphas_lib_o PUBLIC "-fopenmp")
		endif()
	endif()

	if(MSVC)
		if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
			target_compile_options(symphas_lib_o PUBLIC /Qiopenmp)
		else()
			# target_compile_options(symphas_lib_o PUBLIC /openmp:experimental)
		endif()
	else()
	endif()
else()
	message(STATUS "No OpenMP")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	message(STATUS "compiling with MSVC")
	target_compile_options(symphas_lib_o PUBLIC
		$<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:/bigobj>
		$<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:/F85894592>
		$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler /bigobj>
	)

	if(DISABLE_EXCEPTIONS)
		# target_compile_options(symphas_lib_o PUBLIC /EHs-c-)
	endif()
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	message(STATUS "compiling with Clang")
	target_compile_options(symphas_lib_o PUBLIC
		$<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:-Wno-unused-variable>
		$<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:-Wno-unused-function>
		$<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:-Wno-format-security>
		$<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:-ftime-trace>
		$<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:-mavx>
		$<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:-fno-limit-debug-info>
	)

	if(MSVC)
		target_compile_options(symphas_lib_o PUBLIC
			$<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:-fmacro-backtrace-limit=0>
			$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler /bigobj>
		)
	else()
	endif()

	if(DISABLE_EXCEPTIONS)
		if(MSVC)
			target_compile_options(symphas_lib_o PUBLIC "/EHs-c-")
		else()
			target_compile_options(symphas_lib_o PUBLIC "-fno-exceptions")
		endif()
	endif()
endif()

if(WIN32)
	target_compile_definitions(symphas_lib_o PUBLIC _CRT_SECURE_NO_WARNINGS)
endif()

configure_file("lib-defs.h.in" "${CMAKE_SOURCE_DIR}/inc/lib-defs.h")

add_library(symphas_lib SHARED)
target_link_libraries(symphas_lib PUBLIC symphas_lib_o)
