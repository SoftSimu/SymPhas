
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
find_package(FFTW3 REQUIRED)

include_directories(SYSTEM ${FFTW_INCLUDE_PATH})


file(GLOB_RECURSE SOURCES_LIB "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")




set(THREAD_COUNT 2)
set(VECTOR_TYPE_NAME VectorValue)

#add_library(fftw3_imported SHARED IMPORTED)
#set_target_properties(fftw3_imported PROPERTIES
#	IMPORTED_LOCATION "${FFTW_DIR}/bin/fftw3.dll"
#	IMPORTED_IMPLIB "${FFTW3_LIBRARY_DIRS}/fftw3.lib")

configure_file("inc/names.h.in" "${CMAKE_CURRENT_SOURCE_DIR}/inc/names.h")
	
if(NOT MSVC)
	list(APPEND ADD_LINUX_LIBS tbb m ${CMAKE_THREAD_LIBS_INIT})
endif()

add_library(symphas_lib_o OBJECT ${SOURCES_LIB})
target_include_directories(symphas_lib_o PRIVATE ${FFTW3_INCLUDE_DIRS} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
	$<INSTALL_INTERFACE:include>)
target_link_directories(symphas_lib_o PUBLIC ${FFTW3_LIBRARY_DIRS})
target_link_libraries(symphas_lib_o PUBLIC ${ADD_LINUX_LIBS} PRIVATE ${FFTW3_LIBRARIES})
target_compile_definitions(symphas_lib_o PRIVATE LIB_EXPORTS)

if(MSVC)
	target_compile_options(symphas_lib_o PUBLIC /openmp:experimental)
else()
	find_package(OpenMP)
	if(OpenMP_CXX_FOUND)
		target_link_libraries(symphas_lib_o PUBLIC OpenMP::OpenMP_CXX)
	endif()
endif()

add_library(symphas_lib SHARED)
target_link_libraries(symphas_lib PUBLIC symphas_lib_o)



