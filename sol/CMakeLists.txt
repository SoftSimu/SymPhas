
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

file(GLOB_RECURSE SOURCES_SOL "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SOURCES_SOL "${CMAKE_CURRENT_SOURCE_DIR}/src/colormap.cpp")

option(USE_VTK "Use the VTK graphics library to display phase field" OFF)
set(CUSTOM_SYS_IMPLEMENTATION_HEADER "" CACHE STRING "the header for additional implementations of solver and provisional systems")
set(MODEL_CALL_FUNCTION "" CACHE STRING "a function templated on the specialized model type")


if(USE_VTK)
	find_package(VTK COMPONENTS
		CommonColor
		CommonComputationalGeometry
		CommonCore
		CommonDataModel
		CommonExecutionModel
		CommonMath
		CommonMisc
		CommonSystem
		CommonTransforms
		FiltersCore
		FiltersGeneral
		FiltersSources
		ImagingCore
		InteractionStyle
		RenderingCore
		kwiml
	)
	
	if(${VTK_FOUND})
		list(APPEND ADD_DEFINITIONS VTK_ON)
		list(APPEND ADD_LIBRARIES ${VTK_LIBRARIES})
		list(APPEND SOURCES_SOL "${CMAKE_CURRENT_SOURCE_DIR}/src/colormap.cpp")

	else()
		message(STATUS "VTK package not found")
	endif()
endif()


if(CUSTOM_SYS_IMPLEMENTATION_HEADER)
	if(EXISTS ${CMAKE_SOURCE_DIR}/${CUSTOM_SYS_IMPLEMENTATION_HEADER})
		list(APPEND ADD_DEFINITIONS CUSTOM_SYS_HEADER="${CMAKE_SOURCE_DIR}/${CUSTOM_SYS_IMPLEMENTATION_HEADER}")
	elseif(SOLVER_INCLUDE_HEADER_DIR AND EXISTS ${SOLVER_INCLUDE_HEADER_DIR}/${CUSTOM_SYS_IMPLEMENTATION_HEADER})
		list(APPEND ADD_DEFINITIONS CUSTOM_SYS_HEADER="${SOLVER_INCLUDE_HEADER_DIR}/${CUSTOM_SYS_IMPLEMENTATION_HEADER}")
	else()
		message(STATUS "custom system implementation header not found (place header in cmake source directory)")
	endif()
endif()

if(MODEL_CALL_FUNCTION)
	list(APPEND ADD_DEFINITIONS MODEL_APPLY_CALL=${MODEL_CALL_FUNCTION})
	message(STATUS "model selection list is defined using definition MODEL_APPLY_CALL set to '${MODEL_CALL_FUNCTION}'")

	set(BASIC_STENCILS ON CACHE BOOL "compile model selection list with only the basic stencils")
	if(NOT BASIC_STENCILS)
		list(APPEND ADD_DEFINITIONS ALL_STENCILS)
		message(STATUS "building with all stencils available")
	else()
		message(STATUS "building with basic stencils")
	endif()
else()
	message(STATUS "model selection tree is not defined because no definition was given to MODEL_APPLY_CALL")
endif()

message(STATUS "definitions added to module 'sol' are ${ADD_DEFINITIONS}")

add_library(symphas_sol_o OBJECT ${SOURCES_SOL})
target_include_directories(symphas_sol_o PUBLIC 
	${VTK_INCLUDE_DIRS}
	${CUSTOM_SYS_INCLUDE}
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
	$<INSTALL_INTERFACE:include>)
target_link_libraries(symphas_sol_o PUBLIC symphas_sym PRIVATE ${ADD_LIBRARIES})
target_compile_definitions(symphas_sol_o PUBLIC ${ADD_DEFINITIONS})

if(USE_IO)
	target_link_libraries(symphas_sol_o PUBLIC symphas_io)
endif()

if(${VTK_FOUND})
	vtk_module_autoinit(
		TARGETS symphas_sol_o
		MODULES ${VTK_LIBRARIES}
    )
endif()


add_library(symphas_sol SHARED)
target_link_libraries(symphas_sol PUBLIC symphas_sol_o)



