
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
		InteractionStyle
		FiltersCore
		FiltersGeneral
		FiltersSources
		FiltersProgrammable
		ImagingCore
		RenderingCore
		kwiml
		RenderingGL2PSOpenGL2
		RenderingOpenGL2
		IOImage
		QUIET
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


	set(AVAILABLE_DIMENSIONS "2" CACHE STRING "comma separated list of dimensions available to simulate (may dramatically increase compile time and resources)")
	set(AVAILABLE_STENCILS_BASIC_ONLY ON CACHE BOOL "when true, compile only the basic stencils, otherwise all stencils are compiled including fully auto-generated ones")
	set(AVAILABLE_STENCILS_AUTOGENERATION OFF CACHE BOOL "generates stencil coefficients during compile time as needed (may dramatically increase compile time and resources)")

	string(REPLACE "," ";" DIMENSION_LIST ${AVAILABLE_DIMENSIONS})

	foreach(DIMENSION ${DIMENSION_LIST})
		set(AVAILABLE_STENCIL_ACCURACY_${DIMENSION}D "2" CACHE STRING "defines the list of accuracies available for the finite difference stencil in ${DIMENSION}D")
	endforeach()


	if(NOT AVAILABLE_STENCILS_BASIC_ONLY)
		list(APPEND ADD_DEFINITIONS ALL_STENCILS)
		message(STATUS "stencils: building with all stencils available")
	else()
		message(STATUS "stencils: building with basic stencils")
	endif()
	
	if(AVAILABLE_STENCILS_AUTOGENERATION)
		list(APPEND ADD_DEFINITIONS GENERATE_UNDEFINED_STENCILS_ON)
		message(STATUS "stencils: autogenerating required stencils")
	endif()

	
	foreach(DIMENSION ${DIMENSION_LIST})
		string(REPLACE "," ";" ACCURACY_LIST ${AVAILABLE_STENCIL_ACCURACY_${DIMENSION}D})
		list(APPEND ADD_DEFINITIONS ORDER_LIST_${DIMENSION}D=${AVAILABLE_STENCIL_ACCURACY_${DIMENSION}D})
		message(STATUS "stencils: available ${DIMENSION}D accuracy is {${ACCURACY_LIST}}")
		foreach(ACCURACY ${ACCURACY_LIST})
			list(APPEND ADD_DEFINITIONS ORDER_LIST_${DIMENSION}D_HAS_${ACCURACY}H)
		endforeach()
	endforeach()

	list(APPEND ADD_DEFINITIONS AVAILABLE_DIMENSIONS=${AVAILABLE_DIMENSIONS})

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
target_link_libraries(symphas_sol_o PUBLIC symphas_sym PUBLIC ${ADD_LIBRARIES})
target_compile_definitions(symphas_sol_o PUBLIC ${ADD_DEFINITIONS})
target_compile_definitions(symphas_sol_o PRIVATE SOL_EXPORTS)

if(USE_IO)
	target_link_libraries(symphas_sol_o PUBLIC symphas_io)
endif()



add_library(symphas_sol SHARED)
target_link_libraries(symphas_sol PUBLIC symphas_sol_o)

if(${VTK_FOUND})
	vtk_module_autoinit(
		TARGETS symphas_sol symphas_sol_o
		MODULES ${VTK_LIBRARIES}
    )
endif()


