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
		set(VTK_ON "True")
		list(APPEND ADD_LIBRARIES ${VTK_LIBRARIES})
		list(APPEND SOURCES_SOL "${CMAKE_CURRENT_SOURCE_DIR}/src/colormap.cpp")
	else()
		set(VTK_ON "False")
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

set(AVAILABLE_DIMENSIONS "2" CACHE STRING "comma separated list of dimensions available to simulate (may dramatically increase compile time and resources)")
set(AVAILABLE_STENCILS_BASIC_ONLY ON CACHE BOOL "when true, compile only the basic stencils, otherwise all stencils are compiled including fully auto-generated ones")
set(AVAILABLE_STENCILS_AUTOGENERATION OFF CACHE BOOL "generates stencil coefficients during compile time as needed (may dramatically increase compile time and resources)")

# Multi-stencil support: fixed point combinations for each stencil type per dimension
# These are internal defaults - users should use ENABLE_STENCIL_OVERRIDE for configuration
# 1D stencils
set(AVAILABLE_STENCIL_POINTS_LAPLACIAN_1D "3")
set(AVAILABLE_STENCIL_POINTS_GRADLAPLACIAN_1D "4") 
set(AVAILABLE_STENCIL_POINTS_BILAPLACIAN_1D "5")

# 2D stencils
set(AVAILABLE_STENCIL_POINTS_LAPLACIAN_2D "5,9")
set(AVAILABLE_STENCIL_POINTS_GRADLAPLACIAN_2D "6,8,12,16") 
set(AVAILABLE_STENCIL_POINTS_BILAPLACIAN_2D "13,17,21")

# 3D stencils  
set(AVAILABLE_STENCIL_POINTS_LAPLACIAN_3D "7,15,19,21,27")
set(AVAILABLE_STENCIL_POINTS_GRADLAPLACIAN_3D "10,12,28,36,40") 
set(AVAILABLE_STENCIL_POINTS_BILAPLACIAN_3D "21,25,41,52,57")

string(REPLACE "," ";" DIMENSION_LIST ${AVAILABLE_DIMENSIONS})

foreach(DIMENSION ${DIMENSION_LIST})
	set(AVAILABLE_STENCIL_ACCURACY_${DIMENSION}D "2" CACHE STRING "defines the list of accuracies available for the finite difference stencil in ${DIMENSION}D")
endforeach()


if(NOT AVAILABLE_STENCILS_BASIC_ONLY)
	set(ALL_STENCILS "True")
	message(STATUS "stencils: building with all stencils available")
else()
	set(ALL_STENCILS "False")
	message(STATUS "stencils: building with basic stencils")
endif()
	
if(AVAILABLE_STENCILS_AUTOGENERATION)
	set(GENERATE_UNDEFINED_STENCILS_ON "True")
	message(STATUS "stencils: autogenerating required stencils")
else()
	set(GENERATE_UNDEFINED_STENCILS_ON "False")
endif()

# Process multi-stencil point configurations per dimension
foreach(DIMENSION ${DIMENSION_LIST})
	string(REPLACE "," ";" LAPLACIAN_POINTS_LIST_${DIMENSION}D ${AVAILABLE_STENCIL_POINTS_LAPLACIAN_${DIMENSION}D})
	string(REPLACE "," ";" GRADLAPLACIAN_POINTS_LIST_${DIMENSION}D ${AVAILABLE_STENCIL_POINTS_GRADLAPLACIAN_${DIMENSION}D})
	string(REPLACE "," ";" BILAPLACIAN_POINTS_LIST_${DIMENSION}D ${AVAILABLE_STENCIL_POINTS_BILAPLACIAN_${DIMENSION}D})
endforeach()

# Generate stencil point definitions for the header
set(WRITE_STENCIL_POINT_DEFINITIONS "")
foreach(DIMENSION ${DIMENSION_LIST})
	set(WRITE_STENCIL_POINT_DEFINITIONS "${WRITE_STENCIL_POINT_DEFINITIONS}\n#ifndef AVAILABLE_LAPLACIAN_POINTS_${DIMENSION}D\n#define AVAILABLE_LAPLACIAN_POINTS_${DIMENSION}D ${AVAILABLE_STENCIL_POINTS_LAPLACIAN_${DIMENSION}D}\n#endif\n")
	set(WRITE_STENCIL_POINT_DEFINITIONS "${WRITE_STENCIL_POINT_DEFINITIONS}\n#ifndef AVAILABLE_GRADLAPLACIAN_POINTS_${DIMENSION}D\n#define AVAILABLE_GRADLAPLACIAN_POINTS_${DIMENSION}D ${AVAILABLE_STENCIL_POINTS_GRADLAPLACIAN_${DIMENSION}D}\n#endif\n")
	set(WRITE_STENCIL_POINT_DEFINITIONS "${WRITE_STENCIL_POINT_DEFINITIONS}\n#ifndef AVAILABLE_BILAPLACIAN_POINTS_${DIMENSION}D\n#define AVAILABLE_BILAPLACIAN_POINTS_${DIMENSION}D ${AVAILABLE_STENCIL_POINTS_BILAPLACIAN_${DIMENSION}D}\n#endif\n")
endforeach()

	
foreach(DIMENSION ${DIMENSION_LIST})
	string(REPLACE "," ";" ACCURACY_LIST ${AVAILABLE_STENCIL_ACCURACY_${DIMENSION}D})
	message(STATUS "stencils: available ${DIMENSION}D accuracy is {${ACCURACY_LIST}}")
	set(WRITE_ORDER_DEFINITIONS 
		"${WRITE_ORDER_DEFINITIONS}\n#ifndef ORDER_LIST_${DIMENSION}D\n#define ORDER_LIST_${DIMENSION}D ${AVAILABLE_STENCIL_ACCURACY_${DIMENSION}D}\n#endif\n")
	foreach(ACCURACY ${ACCURACY_LIST})
		set(WRITE_ORDER_ACCURACY_DEFINITIONS
			"${WRITE_ORDER_ACCURACY_DEFINITIONS}\n#ifndef ORDER_LIST_${DIMENSION}D_HAS_${ACCURACY}H\n#define ORDER_LIST_${DIMENSION}D_HAS_${ACCURACY}H\n#endif\n")
	endforeach()
endforeach()

message(STATUS "definitions added to module 'sol' are ${ADD_DEFINITIONS}")


# Stencil Discovery and Override System
# =====================================

# Discover available stencils from implementation files
# Based on actual implementations in stencilh2.h and stencilh4.h

# 1D Second-order stencils (from stencilh2.h)
set(IMPLEMENTED_LAPLACIAN_1D_2H "3")
set(IMPLEMENTED_BILAPLACIAN_1D_2H "5") 
set(IMPLEMENTED_GRADLAPLACIAN_1D_2H "4")

# 2D Second-order stencils (from stencilh2.h)
set(IMPLEMENTED_LAPLACIAN_2D_2H "5;9")
set(IMPLEMENTED_BILAPLACIAN_2D_2H "13;17;21")
set(IMPLEMENTED_GRADLAPLACIAN_2D_2H "6;8;12;16")

# 2D Fourth-order stencils (from stencilh4.h)
set(IMPLEMENTED_LAPLACIAN_2D_4H "9;17;21")
set(IMPLEMENTED_BILAPLACIAN_2D_4H "21;25;33;37")
set(IMPLEMENTED_GRADLAPLACIAN_2D_4H "14;18;26;30")

# 3D Second-order stencils (from stencilh2.h)
set(IMPLEMENTED_LAPLACIAN_3D_2H "7;15;19;21;27")
set(IMPLEMENTED_BILAPLACIAN_3D_2H "21;25;41;52;57")
set(IMPLEMENTED_GRADLAPLACIAN_3D_2H "10;12;28;36;40")

# Option to enable stencil override system
option(ENABLE_STENCIL_OVERRIDE "Enable user-configurable stencil point selection via CMake cache variables" OFF)

if(ENABLE_STENCIL_OVERRIDE)
    message(STATUS "stencils: enabling user-configurable stencil override system")
    
    # Create boolean cache variables for each IMPLEMENTED stencil, organized by dimension and accuracy
    foreach(DIMENSION ${DIMENSION_LIST})
        string(REPLACE "," ";" ACCURACY_LIST ${AVAILABLE_STENCIL_ACCURACY_${DIMENSION}D})
        foreach(ACCURACY ${ACCURACY_LIST})
            # Create cache variables for each stencil type and point count
            foreach(STENCIL_TYPE laplacian bilaplacian gradlaplacian)
                string(TOUPPER ${STENCIL_TYPE} STENCIL_TYPE_UPPER)
                if(DEFINED IMPLEMENTED_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H)
                    set(FIRST_STENCIL TRUE)
                    foreach(POINTS ${IMPLEMENTED_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H})
                        # Only enable the first (smallest) stencil by default
                        if(FIRST_STENCIL)
                            set(DEFAULT_VALUE ON)
                            set(FIRST_STENCIL FALSE)
                        else()
                            set(DEFAULT_VALUE OFF)
                        endif()
                        
                        option(STENCIL_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H_${POINTS}PTS 
                               "Enable ${DIMENSION}D ${STENCIL_TYPE} stencil with ${POINTS} points (${ACCURACY}th order accuracy)" 
                               ${DEFAULT_VALUE})
                    endforeach()
                endif()
            endforeach()
        endforeach()
    endforeach()
    
    # Generate override lists from selected boolean options
    set(WRITE_STENCIL_OVERRIDE_DEFINITIONS "")
    foreach(DIMENSION ${DIMENSION_LIST})
        string(REPLACE "," ";" ACCURACY_LIST ${AVAILABLE_STENCIL_ACCURACY_${DIMENSION}D})
        foreach(ACCURACY ${ACCURACY_LIST})
            foreach(STENCIL_TYPE laplacian bilaplacian gradlaplacian)
                string(TOUPPER ${STENCIL_TYPE} STENCIL_TYPE_UPPER)
                set(ENABLED_POINTS_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H "")
                
                if(DEFINED IMPLEMENTED_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H)
                    foreach(POINTS ${IMPLEMENTED_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H})
                        if(STENCIL_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H_${POINTS}PTS)
                            if(ENABLED_POINTS_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H)
                                set(ENABLED_POINTS_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H 
                                    "${ENABLED_POINTS_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H},${POINTS}")
                            else()
                                set(ENABLED_POINTS_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H "${POINTS}")
                            endif()
                        endif()
                    endforeach()
                    
                    # Generate override definitions for stencilincludes.h
                    if(ENABLED_POINTS_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H)
                        string(TOUPPER ${STENCIL_TYPE} STENCIL_MACRO_NAME)
                        if(STENCIL_MACRO_NAME STREQUAL "LAPLACIAN")
                            set(STENCIL_MACRO_NAME "LAP")
                        elseif(STENCIL_MACRO_NAME STREQUAL "BILAPLACIAN")
                            set(STENCIL_MACRO_NAME "BILAP")
                        elseif(STENCIL_MACRO_NAME STREQUAL "GRADLAPLACIAN")
                            set(STENCIL_MACRO_NAME "GRADLAP")
                        endif()
                        
                        set(WRITE_STENCIL_OVERRIDE_DEFINITIONS 
                            "${WRITE_STENCIL_OVERRIDE_DEFINITIONS}\n#define OVERRIDE_${STENCIL_MACRO_NAME}_${DIMENSION}D_${ACCURACY}H_POINTS (${ENABLED_POINTS_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H})\n")
                    endif()
                endif()
            endforeach()
        endforeach()
    endforeach()
    
    # Add enable flag for stencilincludes.h
    set(WRITE_STENCIL_OVERRIDE_DEFINITIONS "#define STENCIL_OVERRIDE_ENABLED\n${WRITE_STENCIL_OVERRIDE_DEFINITIONS}")
    
    message(STATUS "stencils: generated override definitions for user-selected stencils")
else()
    foreach(DIMENSION ${DIMENSION_LIST})
        string(REPLACE "," ";" ACCURACY_LIST ${AVAILABLE_STENCIL_ACCURACY_${DIMENSION}D})
        foreach(ACCURACY ${ACCURACY_LIST})
            # Unset cache variables for each stencil type and point count
            foreach(STENCIL_TYPE laplacian bilaplacian gradlaplacian)
                string(TOUPPER ${STENCIL_TYPE} STENCIL_TYPE_UPPER)
                if(DEFINED IMPLEMENTED_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H)
                    foreach(POINTS ${IMPLEMENTED_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H})
                        unset(STENCIL_${STENCIL_TYPE_UPPER}_${DIMENSION}D_${ACCURACY}H_${POINTS}PTS)
                    endforeach()
                endif()
            endforeach()
        endforeach()
    endforeach()

    set(WRITE_STENCIL_OVERRIDE_DEFINITIONS "")
    message(STATUS "stencils: using default stencil configuration")
endif()



configure_file("inc/stencildefs.h.in" "${CMAKE_SOURCE_DIR}/inc/stencildefs.h")


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
