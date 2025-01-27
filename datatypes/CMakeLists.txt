
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

file(GLOB_RECURSE SOURCES_DATATYPES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu")


add_library(symphas_datatypes_o OBJECT ${SOURCES_DATATYPES})
target_include_directories(symphas_datatypes_o PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
	$<INSTALL_INTERFACE:include>)
target_link_libraries(symphas_datatypes_o PUBLIC symphas_lib)
target_compile_definitions(symphas_datatypes_o PRIVATE DATATYPES_EXPORTS)


add_library(symphas_datatypes SHARED)
target_link_libraries(symphas_datatypes PUBLIC symphas_datatypes_o)


