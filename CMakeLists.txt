
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

cmake_minimum_required(VERSION 3.14)

project(SymPhas
	VERSION 1.0
	DESCRIPTION "Symbolic Phase Field Solver"
	LANGUAGES CXX C)
	
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(BUILD_SHARED_LIBS TRUE)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE) 
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

set(COMBINE_SHARED ON CACHE BOOL "create a single dll to include in the compilation of an executable")

set(ADDITIONAL_COMPILE_OPTIONS "")

if(MSVC)
	list(APPEND ADDITIONAL_COMPILE_OPTIONS "/D_CRT_SECURE_NO_WARNINGS /bigobj ")
endif()


set(LIBDIR "lib")
set(INCLUDEDIR "include")
set(BINDIR "bin")
set(INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/${LIBDIR})
set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/${INCLUDEDIR})

set(USE_IO ON CACHE BOOL "include the module for writing and reading data files")
set(ENABLE_TIMINGS ON CACHE BOOL "enable collecting the elapsed time information of the simulation")


set(SOLVER_INCLUDE_HEADER_DIR "." CACHE PATH "the location of the header file defining the solvers to be used")
set(SOLVER_INCLUDE_HEADER_NAME "" CACHE STRING "the name of the header file within SOLVER_INCLUDE_HEADER_DIR")
set(MODEL_INCLUDE_HEADER_DIR "." CACHE PATH "the location of the header file defining the model definitions")
set(MODEL_INCLUDE_HEADER_NAME "" CACHE STRING "the name of the header file within MODEL_INCLUDE_HEADER_DIR")


if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	list(APPEND ADDITIONAL_COMPILE_OPTIONS "/F85894592")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	list(APPEND ADDITIONAL_COMPILE_OPTIONS "-Wno-unused-function")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftime-trace")
endif()

add_compile_options(${ADDITIONAL_COMPILE_OPTIONS})
message(STATUS ${ADDITIONAL_COMPILE_OPTIONS})


#
# construct a list of the base group names

set(GROUP_NAMES_0 "")
list(APPEND GROUP_NAMES_0 lib datatypes sym sol)


if(USE_IO)
	set(USE_PROC ON CACHE BOOL "include the module for data post processing")
	set(USE_CONF ON CACHE BOOL "include the module for managing configurations")
	list(APPEND GROUP_NAMES_0 io)
else()
	set(USE_PROC OFF)
	set(USE_CONF OFF)
endif()

if(USE_CONF)
	list(APPEND GROUP_NAMES_0 conf)
endif()

if(USE_PROC)
	list(APPEND GROUP_NAMES_0 proc)
endif()

foreach(GROUP_NAME ${GROUP_NAMES_0})
	
	if(${GROUP_NAME} STREQUAL "proc")
		if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${GROUP_NAME})
			
			add_subdirectory(${GROUP_NAME})
			list(APPEND GROUP_NAMES ${GROUP_NAME})

			if(TARGET symphas_proc_o)
				set(PROC_ENABLED ON)
			else()
				set(PROC_ENABLED OFF)
			endif()
		endif()
	else()
		add_subdirectory(${GROUP_NAME})
		list(APPEND GROUP_NAMES ${GROUP_NAME})
	endif()

	if(TARGET symphas_${GROUP_NAME})
		list(APPEND SYMPHAS_DEP_LIBS symphas_${GROUP_NAME})
		list(APPEND SYMPHAS_OBJ_LIBS symphas_${GROUP_NAME}_o)
		list(APPEND SYMPHAS_OBJ_SRCS $<TARGET_OBJECTS:symphas_${GROUP_NAME}_o>)
	endif()

endforeach()


set(MODULES_SOURCE "prereq.cpp")
if(USE_IO)
	list(APPEND MODULES_SOURCE "modules-io.cpp")
endif()

add_library(${PROJECT_NAME} SHARED ${MODULES_SOURCE})
target_compile_features(${PROJECT_NAME}
    PUBLIC 
		cxx_std_17
		cxx_auto_type 
		cxx_variadic_templates 
		cxx_constexpr 
		cxx_decltype_auto 
		cxx_attributes 
		cxx_inheriting_constructors)

target_link_libraries(${PROJECT_NAME} PUBLIC ${SYMPHAS_DEP_LIBS})
set_target_properties(${SYMPHAS_DEP_LIBS} ${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
target_compile_options(${PROJECT_NAME} PUBLIC ${ADDITIONAL_COMPILE_OPTIONS})


if(ENABLE_TIMINGS)
	list(APPEND ADD_DEFINITIONS PRINT_TIMINGS)
endif()

if(MODEL_INCLUDE_HEADER_NAME)
	list(APPEND ADD_DEFINITIONS MODEL_INCLUDE_HEADER="${MODEL_INCLUDE_HEADER_NAME}")
	
	if(MODEL_INCLUDE_HEADER_DIR)
		if(NOT IS_ABSOLUTE ${MODEL_INCLUDE_HEADER_DIR})
			get_filename_component(MODEL_INCLUDE_HEADER_DIR ${CMAKE_SOURCE_DIR}/${MODEL_INCLUDE_HEADER_DIR} ABSOLUTE)
		endif()
	else()
		set(MODEL_INCLUDE_HEADER_DIR ${CMAKE_SOURCE_DIR})
	endif()
	

		
	if(EXISTS ${MODEL_INCLUDE_HEADER_DIR}/${MODEL_INCLUDE_HEADER_NAME})
		message(STATUS "using model include source at '${MODEL_INCLUDE_HEADER_DIR}/${MODEL_INCLUDE_HEADER_NAME}'")
		target_include_directories(${PROJECT_NAME} PUBLIC
			$<BUILD_INTERFACE:${MODEL_INCLUDE_HEADER_DIR}>
			$<INSTALL_INTERFACE:include>)
	else()
		message(WARNING "the folder ${MODEL_INCLUDE_HEADER_DIR} does not contain ${MODEL_INCLUDE_HEADER_NAME}")
		unset(MODEL_INCLUDE_HEADER_NAME)
	endif()

else()
	set(MODEL_INCLUDE_HEADER_DIR ${CMAKE_SOURCE_DIR})

endif()


if(SOLVER_INCLUDE_HEADER_NAME)
	list(APPEND ADD_DEFINITIONS SOLVER_INCLUDE_HEADER="${SOLVER_INCLUDE_HEADER_NAME}")

	if(SOLVER_INCLUDE_HEADER_DIR)
		if(NOT IS_ABSOLUTE ${SOLVER_INCLUDE_HEADER_DIR})
			get_filename_component(SOLVER_INCLUDE_HEADER_DIR ${CMAKE_SOURCE_DIR}/${SOLVER_INCLUDE_HEADER_DIR} ABSOLUTE)
		endif()
	else()
		set(SOLVER_INCLUDE_HEADER_DIR ${CMAKE_SOURCE_DIR})
	endif()
	
	if(EXISTS ${SOLVER_INCLUDE_HEADER_DIR}/${SOLVER_INCLUDE_HEADER_NAME})
		message(STATUS "using solver include source at '${SOLVER_INCLUDE_HEADER_DIR}/${SOLVER_INCLUDE_HEADER_NAME}'")
		target_include_directories(${PROJECT_NAME} PUBLIC
			$<BUILD_INTERFACE:${SOLVER_INCLUDE_HEADER_DIR}>
			$<INSTALL_INTERFACE:include>)
	else()
		message(WARNING "the folder ${SOLVER_INCLUDE_HEADER_DIR} does not contain ${SOLVER_INCLUDE_HEADER_NAME}")
		unset(SOLVER_INCLUDE_HEADER_NAME)
	endif()

else()
	set(SOLVER_INCLUDE_HEADER_DIR ${CMAKE_SOURCE_DIR})

endif()


if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
	list(APPEND ADD_DEFINITIONS DEBUG)
	message(STATUS "using Debug mode")
else()
	message(STATUS "build type = ${CMAKE_BUILD_TYPE}")
endif()
if(BASIC_MODELS)
	list(APPEND ADD_DEFINITIONS BASIC_MODELS)
	message(STATUS "defining BASIC_MODELS")
endif()

target_compile_definitions(${PROJECT_NAME} PUBLIC ${ADD_DEFINITIONS} PRIVATE MODULES_EXPORT)
target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include>)
	
# installation of symphas
# install libraries and install header files into the include directory

get_target_property(SYMPHAS_DEFINITIONS ${PROJECT_NAME} COMPILE_DEFINITIONS)
foreach(DEP_LIB ${SYMPHAS_DEP_LIBS})
	get_target_property(DEP_DEFINITION ${DEP_LIB} COMPILE_DEFINITIONS)
	if(DEP_DEFINITION)
		list(APPEND SYMPHAS_DEFINITIONS ${DEP_DEFINITION})
	endif()
endforeach()

foreach(DEF ${SYMPHAS_DEFINITIONS})
	message("\t${DEF}")
endforeach()

list(APPEND SYMPHAS_LIBS ${SYMPHAS_DEP_LIBS} ${PROJECT_NAME})

if(COMBINE_SHARED)
	add_library(symphas_all SHARED ${MODULES_SOURCE})
	target_link_libraries(symphas_all PRIVATE ${SYMPHAS_OBJ_LIBS})
	target_compile_definitions(symphas_all PRIVATE MODULES_EXPORT)

	foreach(LIB ${SYMPHAS_LIBS})
		get_target_property(LIB_INC ${LIB} INCLUDE_DIRECTORIES)
		if(LIB_INC)
			target_include_directories(symphas_all PUBLIC ${LIB_INC})
		endif()
	endforeach()
	
	install(TARGETS symphas_all
		ARCHIVE DESTINATION ${LIBDIR}
		LIBRARY DESTINATION ${LIBDIR}
		RUNTIME DESTINATION ${BINDIR})
endif()
	



configure_file("symphas.h.in" "symphas.h")
configure_file("prereq-defs.h.in" "prereq-defs.h")

configure_file("symphas.h.in" "${CMAKE_CURRENT_SOURCE_DIR}/symphas.h")
configure_file("prereq-defs.h.in" "${CMAKE_CURRENT_SOURCE_DIR}/prereq-defs.h")

install(TARGETS ${SYMPHAS_OBJ_LIBS} EXPORT SymPhasDepends)
foreach(LIB ${SYMPHAS_LIBS})
	set_target_properties (${LIB} PROPERTIES SOVERSION ${SymPhas_VERSION} VERSION ${SymPhas_VERSION_MAJOR})
	install(TARGETS ${LIB} EXPORT SymPhasDepends
		ARCHIVE DESTINATION ${LIBDIR}
		LIBRARY DESTINATION ${LIBDIR}
		RUNTIME DESTINATION ${BINDIR})
endforeach()


foreach(GROUP_NAME ${GROUP_NAMES})
	file(GLOB_RECURSE INCLUDE_FILES "${GROUP_NAME}/inc/*.h")
	install(FILES ${INCLUDE_FILES} DESTINATION ${INCLUDEDIR})
	message(STATUS "module '${GROUP_NAME}' headers will be installed to ${INCLUDEDIR}")
endforeach()

if(SOLVER_INCLUDE_HEADER_NAME 
	AND NOT ${SOLVER_INCLUDE_HEADER_DIR} STREQUAL ${CMAKE_SOURCE_DIR})

	file(GLOB_RECURSE INCLUDE_FILES "${SOLVER_INCLUDE_HEADER_DIR}/*.h")
	install(FILES ${INCLUDE_FILES} DESTINATION ${INCLUDEDIR})
	message(STATUS "installing solver definition headers.. ${INCLUDE_FILES}")
endif()

if(MODEL_INCLUDE_HEADER_NAME
	AND NOT ${MODEL_INCLUDE_HEADER_DIR} STREQUAL ${SOLVER_INCLUDE_HEADER_DIR}
	AND NOT ${MODEL_INCLUDE_HEADER_DIR} STREQUAL ${CMAKE_SOURCE_DIR})

	file(GLOB_RECURSE INCLUDE_FILES "${MODEL_INCLUDE_HEADER_DIR}/*.h")
	install(FILES ${INCLUDE_FILES} DESTINATION ${INCLUDEDIR})
	message(STATUS "installing model definition headers.. ${CMAKE_SOURCE_DIR}")
endif()

file(GLOB_RECURSE INSTALL_HEADERS "*.h")
install(FILES ${INSTALL_HEADERS} DESTINATION ${INCLUDEDIR})


# configure the pkgconfig file

set(VERSION ${SymPhas_VERSION})
set(URL "https://github.com/MetricManifold/SymPhas")


configure_file(
   "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}.pc.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
   @ONLY)
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
	DESTINATION "${LIBDIR}/pkgconfig")


#
#
# set up the cmake configurations so this can be imported by other projects

include(CMakePackageConfigHelpers)
string(TOLOWER "${PROJECT_NAME}-${VERSION}" VERSION_PROJECT)
set(EXPORT_INSTALL_DIR "${INSTALL_LIBDIR}/cmake/${VERSION_PROJECT}")


write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Configversion.cmake"
	VERSION ${VERSION}
	COMPATIBILITY SameMajorVersion)

configure_package_config_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
	INSTALL_DESTINATION "${EXPORT_INSTALL_DIR}"
	PATH_VARS EXPORT_INSTALL_DIR)

# installation of the configured files
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Configversion.cmake"
	DESTINATION "${EXPORT_INSTALL_DIR}")

export(TARGETS ${SYMPHAS_OBJ_LIBS} ${SYMPHAS_LIBS} ${SYMPHAS_INT_LIBS} 
	NAMESPACE SYMPHAS:: 
	FILE ${PROJECT_BINARY_DIR}/SymPhasDepends.cmake)

install(EXPORT SymPhasDepends 
	NAMESPACE SYMPHAS::
	DESTINATION "${EXPORT_INSTALL_DIR}")
