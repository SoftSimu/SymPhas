
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

find_package(xdrfile)




# link the sources from lib

set(USE_LATEX_FORMAT OFF CACHE BOOL "set the gnuplot plot scripts to output to a latex format")

list(APPEND SOURCES_IO 
	"${CMAKE_CURRENT_SOURCE_DIR}/src/io.cpp" 
	"${CMAKE_CURRENT_SOURCE_DIR}/src/savedefines.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/writedefines.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/readgp.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/readgpc.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/writegp.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/writegpc.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/src/writecsv.cpp")

set(XDR_LIBS "")
if(xdrfile_FOUND)
	set(USING_XDR ON)
	list(APPEND XDR_LIBS ${xdrfile_LIBRARIES})
	list(APPEND XDR_INC ${xdrfile_INCLUDE_DIRS})

	list(APPEND SOURCES_IO 
		"${CMAKE_CURRENT_SOURCE_DIR}/src/readxdr.cpp"
		"${CMAKE_CURRENT_SOURCE_DIR}/src/writexdr.cpp")
else()
	set(USING_XDR OFF)
	message(STATUS "xdrfile package was not found!")
endif()

if(MSVC)
	set(GNU_PLOTTER "windows")
else()
	set(GNU_PLOTTER "qt")
endif()

configure_file("inc/plottername.h.in" "${CMAKE_SOURCE_DIR}/inc/plottername.h")
configure_file("inc/readincludes.h.in" "${CMAKE_SOURCE_DIR}/inc/readincludes.h")
configure_file("inc/writeincludes.h.in" "${CMAKE_SOURCE_DIR}/inc/writeincludes.h")


add_library(symphas_io_o OBJECT ${SOURCES_IO})
target_include_directories(symphas_io_o PUBLIC 
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
	$<INSTALL_INTERFACE:include>)
target_include_directories(symphas_io_o PRIVATE ${XDR_INC})
target_link_libraries(symphas_io_o PUBLIC symphas_datatypes PRIVATE ${XDR_LIBS})
target_compile_definitions(symphas_io_o PRIVATE IO_EXPORTS ${XDR_DEFINITIONS})

if(USE_LATEX_FORMAT)
	set(LATEX_PLOT ON)
endif()

set(USING_IO ON)
configure_file("io-defs.h.in" "${CMAKE_SOURCE_DIR}/inc/io-defs.h")

message(STATUS "LaTeX formats is ${USE_LATEX_FORMAT}")

add_library(symphas_io SHARED)
target_link_libraries(symphas_io PUBLIC symphas_io_o)



	