
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

set(PROCESSES_INCLUDE_H "procinclude.h")

set(DATAPROC_INC "${CMAKE_SOURCE_DIR}/dataprocesses" CACHE STRING "the folder containing '${PROCESSES_INCLUDE_H}'")
set(SOURCES_PROC "${CMAKE_CURRENT_SOURCE_DIR}/src/data.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/processtypes.cpp")


if(DATAPROC_INC)
	if(NOT IS_ABSOLUTE ${DATAPROC_INC})
		get_filename_component(DATAPROC_INC ${CMAKE_SOURCE_DIR}/${DATAPROC_INC} ABSOLUTE)
	endif()
endif()

if(EXISTS ${DATAPROC_INC}/${PROCESSES_INCLUDE_H})
	file(GLOB_RECURSE SOURCES_PROC_APPLY "${DATAPROC_INC}/*.cpp")
	list(APPEND SOURCES_PROC_APPLY "${CMAKE_CURRENT_SOURCE_DIR}/src/proc.cpp")
	message(STATUS ${SOURCES_PROC_APPLY})

	add_library(symphas_proc_o OBJECT ${SOURCES_PROC} ${SOURCES_PROC_APPLY})

	target_link_libraries(symphas_proc_o PUBLIC symphas_sol symphas_sym symphas_datatypes)
	target_compile_definitions(symphas_proc_o PUBLIC USING_PROC PRIVATE PROC_EXPORTS)
	target_include_directories(symphas_proc_o PUBLIC
		$<BUILD_INTERFACE:${DATAPROC_INC}> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc> $<BUILD_INTERFACE:${MODEL_INCLUDE_HEADER_DIR}>
		$<INSTALL_INTERFACE:include>)

	add_library(symphas_proc SHARED)
	target_link_libraries(symphas_proc PUBLIC symphas_proc_o)

else()
	message(STATUS "no processes header ${PROCESSES_INCLUDE_H} found in ${DATAPROC_INC}")
endif()


