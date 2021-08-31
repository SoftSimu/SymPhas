
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

set(SOURCES_CONF "src/conf.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/confsystem.cpp")


if(PROC_ENABLED)
	set(ADD_PROC_LIB symphas_proc)
	list(APPEND SOURCES_CONF "${CMAKE_CURRENT_SOURCE_DIR}/src/confdata.cpp")
	message(STATUS "including proc lib")
endif()

add_library(symphas_conf_o OBJECT ${SOURCES_CONF})
target_include_directories(symphas_conf_o PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
	$<INSTALL_INTERFACE:include>)
target_link_libraries(symphas_conf_o PUBLIC symphas_sol ${ADD_PROC_LIB})
target_compile_definitions(symphas_conf_o PUBLIC USING_CONF)
	
add_library(symphas_conf SHARED)
target_link_libraries(symphas_conf PUBLIC symphas_conf_o)
