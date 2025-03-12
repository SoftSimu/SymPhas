
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

cmake_minimum_required(VERSION 3.22)


project(SymPhas
	VERSION 1.0
	DESCRIPTION "Symbolic Phase Field Solver"
	LANGUAGES CXX C)
	

set(USING_CUDA OFF CACHE BOOL "use CUDA (must be installed)")
set(USING_MPI OFF CACHE BOOL "use MPI (must be installed)")

if(USING_CUDA)
	message(STATUS "CUDA is enabled")
	set(SYMPHAS_CUDA_ARCHITECTURES "" CACHE STRING "CUDA architectures to target")
	if (NOT SYMPHAS_CUDA_ARCHITECTURES)
		# exit with error code
		message(FATAL_ERROR "CUDA architectures must be set")
		#set(USING_CUDA OFF)
	endif()
	set(CMAKE_CUDA_ARCHITECTURES "${SYMPHAS_CUDA_ARCHITECTURES}")
	set(CMAKE_CUDA_STANDARD 17)
	set(CMAKE_CUDA_STANDARD_REQUIRED True)
	enable_language(CUDA)	
else()
	message(STATUS "CUDA is not enabled ${CMAKE_CUDA_COMPILER}")
	
endif()


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(BUILD_SHARED_LIBS TRUE)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE) 
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

set(COMBINE_SHARED ON CACHE BOOL "create a single library to include in the compilation of an executable")


set(LIBDIR "lib")
set(INCLUDEDIR "include")
set(BINDIR "bin")
set(INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/${LIBDIR})
set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/${INCLUDEDIR})

set(USE_IO ON CACHE BOOL "include the module for writing and reading data files")


set(SOLVER_INCLUDE_HEADER_DIR "." CACHE PATH "the location of the header file defining the solvers to be used")
set(SOLVER_INCLUDE_HEADER_NAME "" CACHE STRING "the name of the header file within SOLVER_INCLUDE_HEADER_DIR")
set(MODEL_INCLUDE_HEADER_DIR "." CACHE PATH "the location of the header file defining the model definitions")
set(MODEL_INCLUDE_HEADER_NAME "" CACHE STRING "the name of the header file within MODEL_INCLUDE_HEADER_DIR")
set(USE_TESTS OFF CACHE BOOL "include compilation of the tests")

set(MAIN_FILE "" CACHE STRING "The main file that will be compiled into an executable")

#
# construct a list of the base group names

set(GROUP_NAMES_0 "")
list(APPEND GROUP_NAMES_0 lib datatypes sym sol)


if(USE_IO)
	set(USE_PROC OFF CACHE BOOL "include the module for data post processing")
	set(USE_CONF ON CACHE BOOL "include the module for managing configurations")
	list(APPEND GROUP_NAMES_0 io)
else()
	set(USE_PROC OFF)
	set(USE_CONF OFF)
endif()

if(USE_PROC)
	list(APPEND GROUP_NAMES_0 proc)
endif()

if(USE_CONF)
	message(STATUS "getting 'nlohmann_json' from Github with FetchContent")
	include(FetchContent)
	FetchContent_Declare(nlohmann_json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
	FetchContent_MakeAvailable(nlohmann_json)

	list(APPEND GROUP_NAMES_0 conf)
endif()


set(MODULES_SOURCE "modules-models.cpp")
if(USE_IO)
	list(APPEND MODULES_SOURCE "modules-io.cpp")
endif()

if (USE_TESTS)
	file(GLOB_RECURSE SOURCES_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")

	set(CHOSEN_TESTS "" CACHE STRING "list of tests")
	set_property(CACHE CHOSEN_TESTS PROPERTY STRINGS ${SOURCES_TESTS})
	
	if(NOT "${CHOSEN_TESTS}" STREQUAL "")
		string(REPLACE " " ";" CHOSEN_TEST_LIST ${CHOSEN_TESTS})
	endif()
	
	if("${CHOSEN_TESTS}" STREQUAL "all")
		foreach(SOURCE_TESTS ${SOURCES_TESTS})
			get_filename_component(CHOSEN_TEST_NAME ${SOURCE_TESTS} NAME_WE)
			set(INCLUDE_TEST_FILES "${INCLUDE_TEST_FILES}\n#include \"${CHOSEN_TEST_NAME}.h\"")
			set(RUN_TEST_FILES "${RUN_TEST_FILES}\n\t${CHOSEN_TEST_NAME}();")
		endforeach()
	endif()

	foreach(CHOSEN_TEST ${CHOSEN_TEST_LIST})

		get_filename_component(CHOSEN_TEST_NAME ${CHOSEN_TEST} NAME_WE)
		set(CHOSEN_TEST_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/tests/${CHOSEN_TEST_NAME}.cpp")
	
		if("${CHOSEN_TEST_SOURCE}" IN_LIST SOURCES_TESTS)
			list(APPEND SOURCES_TESTS ${CHOSEN_TEST_SOURCE})
			set(INCLUDE_TEST_FILES "#include \"${CHOSEN_TEST_NAME}.h\"")
			set(RUN_TEST_FILES "\t${CHOSEN_TEST_NAME}();")
		else()
			message(STATUS "test '${CHOSEN_TEST}' does not exist")
			set(SOURCES_TESTS "")
		endif()
	endforeach()

endif()


add_library(${PROJECT_NAME} SHARED ${MODULES_SOURCE} ${SOURCES_TESTS})
target_compile_features(${PROJECT_NAME}
    PUBLIC 
		cxx_std_17
		cxx_auto_type 
		cxx_variadic_templates 
		cxx_constexpr 
		cxx_decltype_auto 
		cxx_attributes 
		cxx_inheriting_constructors)


if(MODEL_INCLUDE_HEADER_NAME)
	list(APPEND SYMPHAS_DEFINITIONS MODEL_INCLUDE_HEADER="${MODEL_INCLUDE_HEADER_NAME}")
	
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
		set(MODEL_INCLUDE_HEADER_NAME "modelmacros.h")
	endif()
else()
	if(MODEL_INCLUDE_HEADER_DIR)
		message(STATUS ${MODEL_INCLUDE_HEADER_DIR})
		if(NOT IS_ABSOLUTE ${MODEL_INCLUDE_HEADER_DIR})
			get_filename_component(MODEL_INCLUDE_HEADER_DIR ${CMAKE_SOURCE_DIR}/${MODEL_INCLUDE_HEADER_DIR} ABSOLUTE)
		endif()
		
		target_include_directories(${PROJECT_NAME} PUBLIC
			$<BUILD_INTERFACE:${MODEL_INCLUDE_HEADER_DIR}>
			$<INSTALL_INTERFACE:include>)
		message(STATUS "the folder ${MODEL_INCLUDE_HEADER_DIR} is set as the model definition folder and will be added to the includes")
	else()
		set(MODEL_INCLUDE_HEADER_DIR ${CMAKE_SOURCE_DIR})
	endif()
	set(MODEL_INCLUDE_HEADER_NAME "modelmacros.h")
endif()


if(SOLVER_INCLUDE_HEADER_NAME)
	list(APPEND SYMPHAS_DEFINITIONS SOLVER_INCLUDE_HEADER="${SOLVER_INCLUDE_HEADER_NAME}")

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
		set(SOLVER_INCLUDE_HEADER_NAME "solver.h")
	endif()
else()
	if(SOLVER_INCLUDE_HEADER_DIR)
		if(NOT IS_ABSOLUTE ${SOLVER_INCLUDE_HEADER_DIR})
			get_filename_component(SOLVER_INCLUDE_HEADER_DIR ${CMAKE_SOURCE_DIR}/${SOLVER_INCLUDE_HEADER_DIR} ABSOLUTE)
		endif()
		
		target_include_directories(${PROJECT_NAME} PUBLIC
			$<BUILD_INTERFACE:${SOLVER_INCLUDE_HEADER_DIR}>
			$<INSTALL_INTERFACE:include>)
		message(STATUS "the folder ${SOLVER_INCLUDE_HEADER_DIR} is set as the solvers folder and will be added to the includes")
	else()
		set(SOLVER_INCLUDE_HEADER_DIR ${CMAKE_SOURCE_DIR})
	endif()
	set(SOLVER_INCLUDE_HEADER_NAME "solver.h")
endif()

get_target_property(EXTRA_INCLUDE_DIRECTORIES ${PROJECT_NAME} INCLUDE_DIRECTORIES)
	
foreach(GROUP_NAME ${GROUP_NAMES_0})

	if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${GROUP_NAME}/${GROUP_NAME}-defs.h.in)
		set(CONFIGURED_DEF_FILES "${CONFIGURED_DEF_FILES}\n#include \"${GROUP_NAME}-defs.h\"")
	endif()

	if(${GROUP_NAME} STREQUAL "proc")
		if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${GROUP_NAME})
			
			add_subdirectory(${GROUP_NAME})
			#message(STATUS ${EXTRA_INCLUDE_DIRECTORIES})
			#target_compile_definitions(symphas_${GROUP_NAME} PRIVATE ${SYMPHAS_DEFINITIONS})

			list(APPEND GROUP_NAMES ${GROUP_NAME})


			if(TARGET symphas_proc_o)
				set(USING_PROC ON)
				
				target_include_directories(symphas_${GROUP_NAME}_o PRIVATE
					$<BUILD_INTERFACE:${MODEL_INCLUDE_HEADER_DIR}>
					$<INSTALL_INTERFACE:include>)
				target_include_directories(symphas_${GROUP_NAME}_o PRIVATE
					$<BUILD_INTERFACE:${SOLVER_INCLUDE_HEADER_DIR}>
					$<INSTALL_INTERFACE:include>)
			
			else()
				set(USING_PROC OFF)
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

set_target_properties(${SYMPHAS_DEP_LIBS} ${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
message(STATUS "the following def files are configured: ${CONFIGURED_DEF_FILES}")
configure_file("configured-defs.h.in" "${CMAKE_SOURCE_DIR}/inc/configured-defs.h")

if (USE_TESTS)
	target_include_directories(${PROJECT_NAME} PUBLIC 
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tests/>)
	set(USING_TESTS ON)
else()
	set(USING_TESTS OFF)
endif()

target_compile_definitions(${PROJECT_NAME} PUBLIC ${SYMPHAS_DEFINITIONS} PRIVATE MODULES_EXPORT)
target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include>)
	

target_link_libraries(${PROJECT_NAME} PUBLIC ${SYMPHAS_DEP_LIBS})
target_compile_options(${PROJECT_NAME} PUBLIC ${ADDITIONAL_COMPILE_OPTIONS})



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


if (USE_TESTS)

	configure_file("tests.cpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/tests.cpp")
	add_library(symphas_tests_o OBJECT "tests.cpp")
	target_include_directories(symphas_tests_o PUBLIC 
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
		$<INSTALL_INTERFACE:include>)
	target_link_libraries(symphas_tests_o PUBLIC ${PROJECT_NAME})

	add_library(symphas_tests SHARED)
	target_link_libraries(symphas_tests PUBLIC symphas_tests_o)

	if("${CHOSEN_TESTS}" STREQUAL "all")
		message(STATUS "using all tests in ${CMAKE_CURRENT_SOURCE_DIR}/tests/")
	else()
		message(STATUS "using test(s) '${CHOSEN_TESTS}' in ${CMAKE_CURRENT_SOURCE_DIR}/tests/")
	endif()

	add_executable(symphas_test "tests.cpp")
	target_link_libraries(symphas_test PUBLIC symphas_tests)

endif()

configure_file("prereq-defs.h.in" "${CMAKE_SOURCE_DIR}/inc/prereq-defs.h")

if (USE_CONF)
	message(STATUS "exporting 'nlohmann_json' library")
	list(APPEND SYMPHAS_LIBS nlohmann_json)
endif()

if(COMBINE_SHARED)
	add_library(symphas_all STATIC ${MODULES_SOURCE})
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
	

install(TARGETS ${SYMPHAS_OBJ_LIBS} EXPORT SymPhasDepends)
foreach(LIB ${SYMPHAS_LIBS})
	set_target_properties (${LIB} PROPERTIES SOVERSION ${SymPhas_VERSION} VERSION ${SymPhas_VERSION_MAJOR})
	install(TARGETS ${LIB} EXPORT SymPhasDepends
		ARCHIVE DESTINATION ${LIBDIR}
		LIBRARY DESTINATION ${LIBDIR}
		RUNTIME DESTINATION ${BINDIR})
endforeach()


foreach(GROUP_NAME IN LISTS GROUP_NAMES ITEMS "${CMAKE_SOURCE_DIR}")
	file(GLOB_RECURSE INCLUDE_FILES "${GROUP_NAME}/inc/*.h")
	install(FILES ${INCLUDE_FILES} DESTINATION ${INCLUDEDIR})
	message(STATUS "module '${GROUP_NAME}' headers will be installed to ${INCLUDEDIR}")
endforeach()

if(SOLVER_INCLUDE_HEADER_NAME 
	AND NOT ${SOLVER_INCLUDE_HEADER_DIR} STREQUAL ${CMAKE_SOURCE_DIR})

	file(GLOB_RECURSE INCLUDE_FILES "${SOLVER_INCLUDE_HEADER_DIR}/*.h")
	install(FILES ${INCLUDE_FILES} DESTINATION ${INCLUDEDIR})
	message(STATUS "installed solver definition headers.. ${INCLUDE_FILES}")
endif()

if(MODEL_INCLUDE_HEADER_NAME
	AND NOT ${MODEL_INCLUDE_HEADER_DIR} STREQUAL ${SOLVER_INCLUDE_HEADER_DIR}
	AND NOT ${MODEL_INCLUDE_HEADER_DIR} STREQUAL ${CMAKE_SOURCE_DIR})

	file(GLOB_RECURSE INCLUDE_FILES "${MODEL_INCLUDE_HEADER_DIR}/*.h")
	install(FILES ${INCLUDE_FILES} DESTINATION ${INCLUDEDIR})
	message(STATUS "installed model definition headers.. ${INCLUDE_FILES}")
endif()


message(STATUS "installing rest of headers..")
file(GLOB INSTALL_HEADERS "*.h")
install(FILES ${INSTALL_HEADERS} DESTINATION ${INCLUDEDIR})
foreach(GROUP_NAME ${GROUP_NAMES_0})
	file(GLOB INSTALL_HEADERS "${GROUP_NAME}/inc/*.h" "${GROUP_NAME}/inc/*.cuh")
	install(FILES ${INSTALL_HEADERS} DESTINATION ${INCLUDEDIR})
endforeach()
message(STATUS "finished installing all headers to ${INCLUDEDIR}")


# check if MAIN_FILE is set and create an executable
if (MAIN_FILE)
	add_executable(symphas_main ${MAIN_FILE})
	target_link_libraries(symphas_main ${SYMPHAS_LIBS})
	
	#if (WIN32)
	#	target_compile_options(symphas_main PRIVATE -v -fmacro-backtrace-limit=0)
	#else()
	#	target_compile_options(symphas_main PRIVATE -v -ftemplate-backtrace-limit=0)
	#endif()
	#target_compile_definitions(symphas_impl PUBLIC SOLVER_TYPE=${SOLVER_TYPE})
endif()


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
