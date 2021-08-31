
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

file(GLOB_RECURSE SOURCES_SYM "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

set(PRINTABLE_EQUATIONS YES CACHE BOOL "if yes, all expressions will have a printable form")

if(PRINTABLE_EQUATIONS)
	list(APPEND ADD_DEFINITIONS PRINTABLE_EQUATIONS)
endif()


add_library(symphas_sym_o OBJECT ${SOURCES_SYM})
target_include_directories(symphas_sym_o PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
	$<INSTALL_INTERFACE:include>)
target_compile_definitions(symphas_sym_o PUBLIC USING_SYM PRIVATE EXPR_EXPORTS)
target_link_libraries(symphas_sym_o PUBLIC symphas_datatypes)

add_library(symphas_sym SHARED)
target_link_libraries(symphas_sym PUBLIC symphas_sym_o)



