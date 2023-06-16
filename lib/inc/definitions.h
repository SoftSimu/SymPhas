
/* ***************************************************************************
 * This file is part of the SymPhas library, a framework for implementing
 * solvers for phase-field problems with compile-time symbolic algebra.
 *
 * Copyright (c) 2018-2021 by Steven A. Silber and Mikko Karttunen
 *
 * SymPhas is free software, which can be redistributed or modified under
 * the terms of the GNU Lesser General Public License (LGPL) as published
 * by the Free Software Foundation; LGPL version 3, or later versions at
 * your choice.
 *
 * SymPhas is distributed with the faith that it will be helpful and
 * practical but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * ***************************************************************************
 *
 * MODULE:  lib
 * PURPOSE: Contains definitions of types and macros used throughout 
 * SymPhas. Certain types are standardized here in order to provide
 * portability as well as clarity in usage of various functionality.
 *
 * ***************************************************************************
 */



#pragma once

#include <complex>
#include <utility>
#include <array>

#include "macros.h"
#include "cell.h"
#include "names.h"

//! \cond

#if __GNUC__ >= 8 || _MSC_VER >= 1914
#	if defined(__has_include)
#		if __has_include(<filesystem>)
#			define FILESYSTEM_HEADER_AVAILABLE
#		endif
#	endif
#endif


#if __GNUC__ >= 9 || _MSC_VER >= 1914
#	if defined(__has_include)
#		if __has_include(<execution>)
#			define EXECUTION_HEADER_AVAILABLE
#		endif
#	endif
#endif

#ifdef _MSC_VER
#define DLLEXPORT __declspec(dllexport)
#define DLLIMPORT __declspec(dllimport)
#else
#define DLLEXPORT
#define DLLIMPORT
#endif

#ifdef LIB_EXPORTS
#define DLLLIB DLLEXPORT
#else
#define DLLLIB DLLIMPORT
#endif

#if THREADS > 1
#define MULTITHREAD
#endif
#define MULTITHREAD_TRIGGER_COUNT 128

//! \endcond



/*! 
 * \namespace symphas
 * \brief Primary namespace for the SymPhas library.
 * 
 * Primary namespace under which various elements used throughout SymPhas are defined,
 * mainly to avoid naming conflicts, and for these commonly used elements to reside under
 * the same name.
 * 
 */



namespace symphas
{
	//! Mathematical constant \f$\pi\f$.
	/*!
	 * Constant value equal to the mathematical constant \f$\pi\f$, used throughout
	 * SymPhas.
	 */
	constexpr double PI = 3.141592653589793238463;

	//! Constant \f$\epsilon\f$, near zero.
	/*!
	 * Constant which represents a value which is almost but not equal to zero.
	 */
	constexpr double EPS = 1e-8;



	//! Alias used for representing reference wrapper, for convenience purposes. 
	/*!
	 * Alias type to rename a reference wrapper from the standard library. The
	 * reference type is commonly used throughout the symbolic algebra library,
	 * primarily to allow variables to incorporate an object in an expression
	 * without copying it or using managed memory.
	 */
	template<typename T>
	using ref = std::reference_wrapper<T>;

}




/* **************************************************************************
 * SymPhas aliased types for portability and convenience.
 * **************************************************************************/


//! Alias type for an iteration variable.
/*!
 * SymPhas using declaration for a type for an iteration or
 * incrementable variable. This type is aliased for portability to 
 * choose the iterable type, if the usual type int isn't appopriate.
 */
using iter_type = int;

//! Alias type representing length.
/*! 
 * SymPhas using declaration for type representing length of a sequence
 * or array. This type is aliased for portability to 
 * choose the iterable type, if the usual type int isn't appopriate.
 */
using len_type = int;




/* **************************************************************************
 * SymPhas aliased types to standardize types in the context of phase
 * field problems.
 * **************************************************************************/


//! Alias type representing a real-valued quantity. 
/*!
 * The alias for a double to represent the type of a phase field order 
 * parameter which is real-valued. Used in SymPhas when the type of the 
 * variable or context refers to the value of an order parameter.
 */
using scalar_t = double;


//! Alias type representing a complex-valued quantity. 
/*!
 * Same as ::scalar_t, except for complex types for order parameters.
 */
using complex_t = std::complex<scalar_t>;

template<typename T, size_t D>
using any_vector_t = VECTOR_TYPE_NAME(T, D);

//! Alias type representing a vector-valued quantity. 
/*!
 * Same as ::scalar_t, except for vector types for order parameters. The vector
 * has the underlying datatype ::scalar_t.
 */
template<size_t D> 
using vector_t = VECTOR_TYPE_NAME(scalar_t, D);

//! Alias type representing a vector-valued quantity of complex type. 
/*!
 * Same as ::scalar_t, except for vector types for order parameters. The vector
 * has the underlying datatype ::scalar_t.
 */
template<size_t D>
using cvector_t = VECTOR_TYPE_NAME(complex_t, D);

//! Get the element type of a vector using type traits.
template<typename T>
struct vector_element_type;

//! Get the element type of a vector using type traits.
template<typename T, size_t D>
struct vector_element_type<any_vector_t<T, D>>
{
	using type = T;
};

//! The type of a coordinate in an axis.
/*!
 * The scalar type which is used to represent a position on one axis.
 */
using axis_coord_t = double;
#define COORD_FMT "%.4lf"


//! An nested alias for the type of a system axis.
/*!
 * An alias type using a struct templated on the dimension to represent the type
 * for the point in the axis of a system in that dimension, which can be specialized 
 * by the dimension, i.e. see axis_nd<1>. 
 * 
 * In particular, a value of this type represents a coordinate point in a system; 
 * hence the axis type in the context of systems of 2 or more dimensions has 
 * the underlying type that is equivalent to a one dimensional array of length `D`. 
 * Each coordinate value of a point is a simple scalar type, typically `double` type, 
 * so the underlying array is of this type.
 * 
 * In usage, axis coordinates are considered to be listed in order starting with 
 * the horizontal component, \f$x\f$.
 */
template<size_t D>
struct axis_nd
{
	using type = std::array<axis_coord_t, D>;
};

//! Specialization of the 1-dimensional axis. 
/*!
 * A 1-dimensional system has a single axis, namely, the coordinate system is fully
 * specified by points which are just one value. This specialization, therefore, 
 * designates the 1-dimensional axis type as the scalar type (as opposed to the  
 * usual array type in higher dimensions).
 */
template<>
struct axis_nd<1>
{
	using type = axis_coord_t;
};
//
//template<size_t D>
//struct axis_nd : std::array<axis_coord_t, D>
//{
//	using std::array<axis_coord_t, D>::array;
//	template<typename... Ts>
//	axis_nd(Ts&&... data) : std::array<axis_coord_t, D>({ std::forward<Ts>(data)... }) {}
//};
//
////! Specialization of the 1-dimensional axis. 
///*!
// * A 1-dimensional system has a single axis, namely, the coordinate system is fully
// * specified by points which are just one value. This specialization, therefore,
// * designates the 1-dimensional axis type as the scalar type (as opposed to the
// * usual array type in higher dimensions).
// */
//template<>
//struct axis_nd<1>
//{
//	axis_nd(axis_coord_t value = 0) : value{ value } {}
//
//	operator axis_coord_t const& () const
//	{
//		return value;
//	}
//
//	operator axis_coord_t& ()
//	{
//		return value;
//	}
//
//protected:
//
//	axis_coord_t value;
//};

//! Alias to nested type axis_nd::type;
/*!
 * Convenience alias to avoid specifying the nested type name within axis_nd.
 */
template<size_t D>
using axis_nd_t = typename axis_nd<D>::type;

//! Additional name for 1d axis type.
/*!
 * An alias for the 1-dimensional specialization of the axis type. Primarily used
 * for convenience purposes in referring to the 1-dimensional axis type.
 */
using axis_1d_type = axis_nd_t<1>;

//! Additional name for 2d axis type.
/*!
 * An alias for the 2-dimensional axis type. Refer to ::axis_1d_type.
 */
using axis_2d_type = axis_nd_t<2>;

//! Additional name for 3d axis type.
/*!
 * An alias for the 3-dimensional axis type. Refer to ::axis_1d_type.
 */
using axis_3d_type = axis_nd_t<3>;




/*
 * Type traits which return the type of the given arithmetic result between
 * values of generic types.
 */

//! Derive the result type of multiplication between the given types.
template<typename... Ts>
using mul_result_t = decltype((std::declval<Ts>() * ...));

//! Derive the result type of division between the given types.
template<typename... Ts>
using div_result_t = decltype((std::declval<Ts>() / ...));

//! Derive the result type of addition between the given types.
template<typename... Ts>
using add_result_t = decltype((std::declval<Ts>() + ...));

//! Derive the result type of addition between the given types.
template<typename... Ts>
using sub_result_t = decltype((std::declval<Ts>() - ...));


//! Compile time constant template variable to return maximum.
template<auto A, auto... Rest>
constexpr auto fixed_max = fixed_max<A, fixed_max<Rest...>>;
template<auto A, auto B>
constexpr auto fixed_max<A, B> = (A > B) ? A : B;
template<auto A>
constexpr auto fixed_max<A> = A;


//! Compile time constant template variable to return minimum.
template<auto A, auto... Rest>
constexpr auto fixed_min = fixed_min<A, fixed_min<Rest...>>;
template<auto A, auto B>
constexpr auto fixed_min<A, B> = (A < B) ? A : B;
template<auto A>
constexpr auto fixed_min<A> = A;
template<int I>
constexpr int fixed_abs = (I < 0) ? -I : I;
template<int I>
constexpr int fixed_sign = (I < 0) ? -1 : 1;



template<size_t N, size_t D, size_t... Rest>
constexpr size_t GCD_of = GCD_of<GCD_of<N, D>, Rest...>;
template<size_t N, size_t D>
constexpr size_t GCD_of<N, D> = (N > D) ? GCD_of<D, N - (N / D) * D> : GCD_of<N, D - (D / N) * N>;
template<size_t N>
constexpr size_t GCD_of<N, 0> = N;
template<size_t D>
constexpr size_t GCD_of<0, D> = 1;


//! Default names of order parameters that are used when unspecified.
/*!
 * An array of names that assigns names to order parameters when not explicitly
 * named by the user. The number of order parameters that have a default name
 * is unlimited, but once all the names in this pool are used, then a specific
 * format is used to assign names.
 */
inline const char* ORDER_PARAMETER_NAMES[] = {
#ifdef LATEX_PLOT
"\\psi", "\\rho"
#else
"psi", "rho"
#endif
};

//! Default name format of order parameter when other names have been assigned.
/*!
 * A string format for putting a number value into a name to allow an
 * unlimited number of order parameters to be named. This format is used once
 * all the other names are used.
 */
inline const char* ORDER_PARAMETER_NAME_EXTRA_FMT =
#ifdef LATEX_PLOT
"\\phi_%d";
#else
"op%d";
#endif


//! Default names of variables that are used when unspecified.
/*!
 * Assigns names to variables in order to be neatly displayed. Any variable
 * can be named.
 */
inline const char* VARIABLE_NAMES[] = { "Q", "W", "U", "V" };

//! Default format of variable name when other names have been assigned.
/*!
 * Uses a format to name a variable when all the other names have already been
 * assigned.
 */
inline const char* VARIABLE_NAME_EXTRA_FMT = "v%d";



/* **************************************************************************
 * Enumerations used in SymPhas in a general context
 * **************************************************************************/

//! Represents system geometry.
/*!
 * An enumeration with values which represent the symmetry of a system in the
 * context of a phase field problem. This is used in a general context in SymPhas.
 * They are tags that are associated with the coordinate axes of a given system.
 */
enum class Geometry 
{
	CARTESIAN,	//!< Cartesian geometry for both two and three dimensional systems
	POLAR,		//!< Polar geometry, applies only to two dimensional systems
	SPHERICAL,	//!< Spherical coordinates, applies only to three dimensions
	CYLINDRICAL	//!< Cylindrical coordinates, applies only to three dimensions
};



//! Values for labeling the axes of a grid.
/*!
 * Global names that are used to refer to the axes of a grid, in use cases such
 * as labeling the intervals of a given axis. They are also listed in the order
 * that the coordinate point would be written.
 *
 * The axes listed in this enumeration apply for grids up to 3 dimensions.
 */
enum class Axis
{
	X,		//!< The horizontal component of a grid.
	Y,		//!< The vertical component of a grid.
	Z,		//!< The depth component of a grid.
	T,		//!< The angle of polar coordinates.
	S,		//!< The angle of spherical coordinates.
	R,		//!< The radius of polar or spherical coordinates.
	NONE	//!< No axis is specified.
};


/* **************************************************************************
 * Macro definitions used throughout SymPhas
 * **************************************************************************/

//! Filepointer to SymPhas log.
/*! 
 * Name of the filepointer to where the regular log messages are sent.
 */
#define SYMPHAS_LOG stdout

#define SYMPHAS_DEBUG SYMPHAS_LOG
#define SYMPHAS_INFO SYMPHAS_LOG
#define SYMPHAS_WARN SYMPHAS_LOG
#define SYMPHAS_ERR stderr


//! Character log banner.
/*!
 * Definition for the string of characters representing an ASCII banner used
 * in log output.
 */
#define OUTPUT_BANNER "-------------------------------------------------------------\n"






 //! Start index of a simulation.
 /*!
  * The index defined to be the initial index of a simulation will start.
  */
#define INDEX_INIT 0

 //! Start time of a simulation.
 /*!
  * The time defined to be the start of the simulation.
  */
#define TIME_INIT 0


#define NUM_INIT_CONSTANTS 6			//!< Maximum number of parameters for the initial conditions.
#define DEFAULT_COEFF_VALUE 1.0			//!< The default value of a coefficient in the equations of motion


//! Format string for printing the timestamp.
/*!
 * The timestamp format string takes the date time information in order from
 * the largest time interval (year) to the smallest:
 * 
 * Year . Month . Day . Hour . Minute . Second
 */
#define TIMESTAMP_FMT "%d.%02d.%02d.%02d.%02d.%02d"

//! Format string of additional timestamp uniqueness identifier.
/*!
 * If the timestamp string must be unique, i.e. when creating a system file
 * name which must not already exist, then an additional string of this 
 * format can be appended to the timestamp.
 */
#define TIMESTAMP_ID_APPEND "_%d"



// ************************************************************
// The following definitions are removed for the API and are internal only


#define ERR_CODE_FILE_OPEN 1000
#define SYMPHAS_MSG_ERR_FILE_OPEN "error opening file '%s'\n"

//! \cond

#define BUFFER_LENGTH 256
#define BUFFER_LENGTH_R(N) (BUFFER_LENGTH >> N)
#define BUFFER_LENGTH_L(N) (BUFFER_LENGTH << N

#define BUFFER_LENGTH_R4 16
#define BUFFER_LENGTH_R2 64
#define BUFFER_LENGTH_R1 128
#define BUFFER_LENGTH_L2 1024
#define BUFFER_LENGTH_L3 2056
#define BUFFER_LENGTH_L4 (BUFFER_LENGTH << 4)
#define BUFFER_LENGTH_L5 (BUFFER_LENGTH << 5)


#define LINE_READ_BUFFER BUFFER_LENGTH_L3			//!< Buffer length from reading an entire line from input.
#define FILECHUNK_READ_BUFFER BUFFER_LENGTH_L(8)	//!< Buffer length for reading a large chunk of a file.

//! \endcond

