
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
 * PURPOSE: Defines functions for basic functionality throughout SymPhas.
 * This includes functions which perform common tasks.
 *
 * ***************************************************************************
 */

#pragma once


#ifdef _MSC_VER
#include <windows.h>
#undef max
#undef min
#else
#include <errno.h>
#include <sys/stat.h>
#include <libgen.h>
#endif

#ifdef _MSC_VER 
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif


#include <vector>
#include <iostream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <numeric>


//#include "symphasthread.h"
#include "definitions.h"


#ifdef FILESYSTEM_HEADER_AVAILABLE
#include <filesystem>
#endif

namespace symphas::internal
{
	template<size_t I, typename... Ss>
	struct index_of_type_match;
}




namespace symphas
{
	//! Defines math functions for arbitrary types.
	/*!
	 * Additional types are introduced by SymPhas. Therefore, mathematical 
	 * functions which take arbitrary types as the argument are implemented for 
	 * common functions including sine and cosine.
	 */
	namespace math {}

	//! Functionality internal to SymPhas.
	/*!
	 * Namespace defining all functionality which is used only internally to
	 * SymPhas.
	 */
	namespace internal {}
}

namespace symphas::internal
{

	template<size_t D, typename T>
	void sort_reorganize(axis_nd_t<D>* data_x, T* data_y, iter_type* sort_array, len_type len)
	{
		for (iter_type i = 0; i < len; ++i)
		{
			axis_nd_t<D> x;
			x = data_x[i];
			T y = data_y[i];

			iter_type j = i;
			while (true)
			{
				iter_type k = sort_array[j];
				sort_array[j] = j;

				if (k == i) break;

				data_x[j] = data_x[k];
				data_y[j] = data_y[k];

				j = k;
			}

			data_x[j] = x;
			data_y[j] = y;
		}
	}


	struct any_case_comparator
	{
		bool operator() (std::string const& s1, std::string const& s2) const
		{
			return strcasecmp(s1.c_str(), s2.c_str()) < 0;
		}
		bool operator() (std::string const& s1, const char* s2) const
		{
			return strcasecmp(s1.c_str(), s2) < 0;
		}
		bool operator() (const char* s1, std::string const& s2) const
		{
			return strcasecmp(s1, s2.c_str()) < 0;
		}
		bool operator() (const char* s1, const char* s2) const
		{
			return strcasecmp(s1, s2) < 0;
		}
	};




	template<size_t D>
	auto dims_as_tuple(const len_type* dims);

	template<>
	inline auto dims_as_tuple<1>(const len_type* dims)
	{
		return std::make_tuple(dims[0]);
	}

	template<>
	inline auto dims_as_tuple<2>(const len_type* dims)
	{
		return std::make_tuple(dims[0], dims[1]);
	}

	template<>
	inline auto dims_as_tuple<3>(const len_type* dims)
	{
		return std::make_tuple(dims[0], dims[1], dims[2]);
	}
}


/*! \addtogroup grid
 * @{
 */

//! Defines functions providing information about a grid.
namespace grid
{
	//! Return the length of a grid with the prescribed dimensions.
	/*!
	 * The number of elements in a grid with the prescribed dimensions is
	 * computed and returned.
	 * 
	 * \param dimensions The number of elements along each axis.
	 * 
	 * \param D The dimension of the grid.
	 */
	template<size_t D>
	len_type length(len_type const* dimensions);

	//! Specialization of grid::length(len_type const*).
	template<>
	inline len_type length<1>(len_type const* dimensions)
	{
		return ((dimensions != nullptr) ? dimensions[0] : 0);
	}

	template<size_t D>
	len_type length(len_type const* dimensions)
	{
		return ((dimensions != nullptr) ? dimensions[D - 1] : 0) * length<D - 1>(dimensions);
	}

	inline len_type length(len_type const* dimensions, size_t dimension)
	{
		len_type len = 1;
		for (iter_type i = 0; i < dimension; ++i)
		{
			len *= dimensions[i];
		}
		return len;
	}

}

//! @}


//! Generic template for returning the result of the math function.
/*!
 * For a value of generic type which is not accepted by any overloads,
 * it is explicitly cast to the type `complex_t` to compute the result.
 */
#define MATH_FUNCTION_OVERLOADS(NAME, FUNC) \
inline scalar_t _ ## NAME(int v) { return FUNC(static_cast<scalar_t>(v)); } \
inline scalar_t _ ## NAME(scalar_t v) { return FUNC(v); } \
inline complex_t _ ## NAME(complex_t v) { return FUNC(v); } \
template<typename T> \
auto _ ## NAME(T v) { return _ ## NAME(static_cast<complex_t>(v)); } \
template<typename T> \
auto NAME(T v) { return _ ## NAME(v); }

namespace symphas::math
{
	namespace
	{



	}
	
	//! Returns conjugate complex number for a generic type.
	/*!
	 * Since additional types are introduced in the SymPhas library, an
	 * extension of the standard library `conj` function is implemented to
	 * apply to additional types.
	 * 
	 * This implementation simply forwards its argument to the standard library
	 * `conj` function.
	 */
	template<typename T>
	complex_t conj(T const& v)
	{
		return std::conj(v);
	}

	//! Returns modulus of a value of a generic type.
	/*!
	 * The generic type has some relationship with a complex number, to allow
	 * taking the modulus to be a logically accepted result.
	 * 
	 * Since additional types are introduced in the SymPhas library, an
	 * extension of the standard library `abs` function is implemented to
	 * apply to additional types.
	 *
	 * This implementation simply forwards its argument to the standard library
	 * `abs` function.
	 */
	template<typename T>
	scalar_t modulus(T const& v)
	{
		return std::abs(v);
	}
	
	//! Returns the real part of a value of a generic type.
	/*!
	 * The generic type has some relationship with a complex number, to allow
	 * taking the real part to be a logically accepted result.
	 * 
	 * Since additional types are introduced in the SymPhas library, an
	 * extension of the standard library `real` function is implemented to
	 * apply to additional types.
	 *
	 * This implementation simply forwards its argument to the standard library
	 * `real` function.
	 */
	template<typename T>
	scalar_t real(T const& v)
	{
		using namespace std; // standard namespace brought into scope
		return real(v);
	}

	//! Returns the imaginary part of a value of a generic type.
	/*!
	 * The generic type has some relationship with a complex number, to allow
	 * taking the imaginary part to be a logically accepted result.
	 *
	 * Since additional types are introduced in the SymPhas library, an
	 * extension of the standard library `imag` function is implemented to
	 * apply to additional types.
	 *
	 * This implementation simply forwards its argument to the standard library
	 * `imag` function.
	 */
	template<typename T>
	scalar_t imag(T const& v)
	{
		using namespace std; // standard namespace brought into scope
		return imag(v);
	}


	MATH_FUNCTION_OVERLOADS(cos, std::cos);
	MATH_FUNCTION_OVERLOADS(sin, std::sin);
	MATH_FUNCTION_OVERLOADS(tan, std::tan);
	MATH_FUNCTION_OVERLOADS(cosh, std::cosh);
	MATH_FUNCTION_OVERLOADS(sinh, std::sinh);
	MATH_FUNCTION_OVERLOADS(tanh, std::tanh);
	MATH_FUNCTION_OVERLOADS(acos, std::acos);
	MATH_FUNCTION_OVERLOADS(asin, std::asin);
	MATH_FUNCTION_OVERLOADS(atan, std::atan);
	MATH_FUNCTION_OVERLOADS(acosh, std::acosh);
	MATH_FUNCTION_OVERLOADS(asinh, std::asinh);
	MATH_FUNCTION_OVERLOADS(atanh, std::atanh);
	MATH_FUNCTION_OVERLOADS(sqrt, std::sqrt);

}




// *************************************************************************************************

namespace symphas
{
	//! General support functions.
	/*!
	 * The namespace which contains all the general support functions used 
	 * throughout SymPhas.
	 */
	namespace lib {}
}

namespace symphas::lib
{
	//! Determine the length of the coordinate relative to the origin. 
	/*!
	 * Determine the length of the coordinate relative to the origin. 
	 * 
	 * \param p The axis point in a 1-dimensional grid.
	 */
	inline axis_coord_t length(axis_1d_type p)
	{
		return std::abs(p);
	}

	//! Determine the length of the coordinate relative to the origin. 
	/*!
	 * Determine the length of the coordinate relative to the origin.
	 *
	 * \param p The axis point in a 2-dimensional grid.
	 */
	inline axis_coord_t length(axis_2d_type p)
	{
		return std::sqrt(p[0] * p[0] + p[1] * p[1]);
	}

	//! Determine the length of the coordinate relative to the origin. 
	/*!
	 * Determine the length of the coordinate relative to the origin.
	 *
	 * \param p The axis point in a 3-dimensional grid.
	 */
	inline axis_coord_t length(axis_3d_type p)
	{
		return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
	}

	//! Determine the distance between two axis points.
	/*!
	 * Determine the distance between two axis points.
	 *
	 * \param p0 The first point in a 1-dimensional grid.
	 * \param p1 The second point in a 1-dimensional grid.
	 */
	inline axis_coord_t distance(axis_1d_type p0, axis_1d_type p1)
	{
		return length(axis_1d_type{ p0 - p1 });
	}

	//! Determine the distance between two axis points.
	/*!
	 * Determine the distance between two axis points.
	 *
	 * \param p0 The first point in a 2-dimensional grid.
	 * \param p1 The second point in a 2-dimensional grid.
	 */
	inline axis_coord_t distance(axis_2d_type p0, axis_2d_type p1)
	{
		return length(axis_2d_type{ p0[0] - p1[0], p0[1] - p1[1] });
	}

	//! Determine the distance between two axis points.
	/*!
	 * Determine the distance between two axis points.
	 *
	 * \param p0 The first point in a 3-dimensional grid.
	 * \param p1 The second point in a 3-dimensional grid.
	 */
	inline axis_coord_t distance(axis_3d_type p0, axis_3d_type p1)
	{
		return length(axis_3d_type{ p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2] });
	}

	//! Obtain the number of elements in a list of unique distances.
	/*!
	 * For a system of the prescribed dimensions, return the number of
	 * unique distances from the origin as computed by the Euclidean norm.
	 * 
	 * \param dims The number of elements along each axis.
	 */
	len_type radial_length(const len_type(&dims)[1]);

	//! See symphas::lib::radial_length().
	len_type radial_length(const len_type(&dims)[2]);

	//! See symphas::lib::radial_length().
	len_type radial_length(const len_type(&dims)[3]);


	//! Construct two lists representing distances and multiplicity.
	/*!
	 * For a system of the prescribed dimensions, construct two lists. The
	 * first list is a sorted list of all the unique distances of points
	 * measured from the origin by the Euclidean norm. The second list is
	 * the number of times that unique distance appears in the regular grid.
	 * The second list is typically used for renormalization when counting
	 * points of the radial average.
	 * 
	 * To obtain the mapping between points in a regular grid and
	 * their position in the distance list, use the function
	 * symphas::lib::make_radial_index_map().
	 *
	 * \param dims The number of elements along each axis.
	 * \param dx The spatial separation between points in the axis.
	 */
	std::tuple<double*, size_t*> make_radial_arrays(const len_type(&dims)[1],
		double dx = 1.0);
	//! See symphas::lib::make_radial_arrays(const len_type(&)[1]).
	std::tuple<double*, size_t*> make_radial_arrays(const len_type(&dims)[2],
		double dx = 1.0, double dy = 1.0);
	//! See symphas::lib::make_radial_arrays(const len_type(&)[1]).
	std::tuple<double*, size_t*> make_radial_arrays(const len_type(&dims)[3],
		double dx = 1.0, double dy = 1.0, double dz = 1.0);



	//! A list of values that can be indexed to get positions in radial arrays.
	/*!
	 * An array is returned, which has length equal to the length of the
	 * system of the prescribed dimensions, and each element in this array
	 * corresponds to elements of that system. The values in the array are 
	 * the positions (indices) in the lists that are returned by 
	 * symphas::lib::make_radial_arrays().
	 * Thus, one can obtain the index in the radial array of any point in a
	 * regular array.
	 * 
	 * In other words, the list maps points from the given system into the
	 * list of distances from the origin.
	 * In this way, an index \f$n\f$ for a point in a
	 * regular grid maps to an index in the distance list, representing
	 * the distance that point \f$n\f$ from the origin (where \f$n=0\f$
	 * for the origin).
	 * 
	 * \param dims The number of elements along each axis.
	 */
	iter_type* make_radial_index_map(const len_type(&dims)[1]);
	iter_type* make_radial_index_map(const len_type(&dims)[2]);
	iter_type* make_radial_index_map(const len_type(&dims)[3]);


	//! Reimplementation of `make_tuple` which differs in return type.
	/*!
	 * A reimplementation of the standard library function `make_tuple` which
	 * accepts arguments by copy instead. This is sometimes the desired method
	 * of creating a tuple because standard library `make_tuple` will convert
	 * all given types `Ti` using `std::decay<Ti>::type` to produce type `Vi`
	 * constituting the tuple types. 
	 * 
	 * This function will instead concatenate tuples produced from the copies 
	 * of the arguments using the type of the argument directly, meaning that 
	 * the resulting tuple type is exactly of the types provided.
	 * 
	 * \param[in] ts The values which are packaged into a tuple.
	 * 
	 * \tparam Ts The types of the arguments.
	 */
	template<typename... Ts>
	auto make_tuple(Ts... ts);

	//! Unpack a tuple of pairs into a pair of tuples.
	/*!
	 * Given a tuple of pairs, two new tuples are created and returned in a
	 * pair, such that the first tuple consists of elements that are the first
	 * element of the pairs in the given tuple, and the second tuple consists of
	 * elements that are the second element of the pairs in the given tuple.
	 * 
	 * \param[in] ts The list of pairs which will be split into two tuples.
	 * 
	 * \tparam As The list of types corresponding to the list of first elements
	 * in the given pair list.
	 * \tparam Bs The list of types corresponding to the list of second elements
	 * in the given pair list.
	 */
	template<typename... As, typename... Bs>
	auto unfurl_tuple(std::tuple<std::pair<As, Bs>...> const& ts);

	//! Return all elements at and beyond the given index in the tuple.
	/*!
	 * A new tuple is constructed containing the elements from the given tuple
	 * such that the index of the element is equal to and greater the given
	 * template parameter `I`.
	 * 
	 * \param[in] t The tuple from where the output elements are taken.
	 * 
	 * \tparam I The index in the given tuple from where elements are taken.
	 * \tparam Ts The list of types of elements in the tuple.
	 */
	template<size_t I, typename... Ts>
	auto get_tuple_ge(std::tuple<Ts...> const& t);

	//! Return all elements before the given index in the tuple.
	/*!
	 * A new tuple is constructed containing the elements from the given tuple
	 * such that the index of the element is (strictly) less than the given
	 * template parameter `I`.
	 *
	 * \param[in] t The tuple from where the output elements are taken.
	 *
	 * \tparam I The index in the given tuple before which elements are taken.
	 * \tparam Ts The list of types of elements in the tuple.
	 */
	template<size_t I, typename... Ts>
	auto get_tuple_lt(std::tuple<Ts...> const& t);

	//! Return all elements at and beyond the given index in the sequence.
	/*!
	 * A new sequence is constructed containing the elements from the given sequence
	 * such that the index of the element is equal to and greater the given
	 * template parameter `I`.
	 *
	 * \param[in] seq The tuple from where the output elements are taken.
	 *
	 * \tparam I The index in the given tuple from where elements are taken.
	 * \tparam Ns... The list of types of elements in the tuple.
	 */
	template<size_t I, size_t... Ns>
	auto get_seq_ge(std::index_sequence<Ns...> const& seq);

	//! Return all elements before the given index in the sequence.
	/*!
	 * A new sequence is constructed containing the elements from the given sequence
	 * such that the index of the element is (strictly) less than the given
	 * template parameter `I`.
	 *
	 * \param[in] seq The tuple from where the output elements are taken.
	 *
	 * \tparam I The index in the given tuple before which elements are taken.
	 * \tparam Ns... The list of types of elements in the tuple.
	 */
	template<size_t I, size_t... Ns>
	auto get_seq_lt(std::index_sequence<Ns...> const& seq);


	//! Constructs a list based on element types in the provided lists.
	/*!
	 * A list is constructed from two inputs, consisting of keys and values. The 
	 * content of the keys is derived from the first elements of the inputs, and
	 * the value is a pair consisting of elements from the two inputs which 
	 * correspond to the type of the key. 
	 * 
	 * In particular, a list of pairs is produced such that each pair has a key
	 * called the matched type. The types of the first elements from each of the 
	 * pairs in the input lists are considered the matched types. The value is a
	 * pair of elements that is associated with the key, and is constructed 
	 * using the second elements from the input lists. The construction is such 
	 * that the first element of the value is always associated with the first 
	 * input list, and the second element the second input list. The two 
	 * elements constituting the value must satisfy the condition that the first 
	 * element from the pair from which they are derived (i.e. the pair the
	 * elements are sourced from in the input lists) have the type same as the
	 * key type. 
	 * 
	 * The matched type is determined by searching through the first list for
	 * the pair where its first element matches the first element of the leading
	 * pair in the second list. If the elements aren't of the same type, the
	 * next pair in the first list is checked against the leading element of the
	 * second list, until either the list is exhausted or a match is found. If a
	 * match is found, a new pair is created as the key/value combination. If a
	 * match is not found, the key/value combination is created with the second
	 * element of the leading pair in the second list and a default
	 * initialized instance of the given type `empty_type`. This is also the
	 * case when the first list is fully matched and there are elements
	 * remaining in the second list. If the second list is empty and items
	 * remain in the first list, then additional pairs will be appended for
	 * which the value pair has first element sourced from the first list pairs
	 * and the second element the default initialized `empty_type` type.
	 * 
	 * The construction of each pair in the output tuple can be summarized by
	 * the following recipe:
	 * 
	 * Return a pair consisting of:
	 * 1. The key called the matched type; the first elements of each pair in
	 * both lists are the set of all keys.
	 * 2. A pair consisting of one of the following:
	 *     1. Values from both the first and second list.
	 *     2. A value from the first list and a default initialized
	 * `empty_type`.
	 *     3. A default initialized `empty_type` and a value from the second
	 * list.
	 *
	 * \param[in] as The first list of pairs, which have their first elements
	 * matched against those in the second list.
	 * \param[in] bs The second list of pairs, which have their first elements
	 * matched against those in the first list.
	 *
	 * \tparam empty_list The object type which is default initialized when
	 * there are no matches in the opposing list.
	 * \tparam match_A_types The list of all key types from the first list.
	 * \tparam As The list of all values from the first list.
	 * \tparam match_B_types The list of all key types from the second list.
	 * \tparam Bs The list of all values from the second list.
	 */
	template<typename empty_type, typename... match_A_types, typename... As, typename... match_B_types, typename... Bs>
	auto unzip_pot(
		std::tuple<std::pair<match_A_types, As>...> const& as, 
		std::tuple<std::pair<match_B_types, Bs>...> const& bs);



	//! Return the index of the given type in the type list.
	/*!
	 * Given the type list, identify the index of the `I`-th chosen type. If
	 * the type is not identified, then `-1` will be returned.
	 *
	 * \tparam Type The type to be found in the list.
	 * \tparam I The current running index of the type to find.
	 * \tparam S0 The next type in the type list to compare.
	 */
	template<typename Type, size_t I, size_t Sn, typename S0,
		typename std::enable_if_t<!std::is_same<Type, S0>::value, int> = 0>
		constexpr int index_of_type()
	{
		return -1;
	}

	//! Return the index of the given type in the type list.
	/*!
	 * Given the type list, identify the index of the `I`-th chosen type. If
	 * the type is not identified, then `-1` will be returned.
	 *
	 * \tparam Type The type to be found in the list.
	 * \tparam I The current running index of the type to find.
	 * \tparam S0 The next type in the type list to compare.
	 * \tparam Ss... The remaining types in the list.
	 */
	template<typename Type, size_t I, size_t Sn, typename S0, typename... Ss,
		typename std::enable_if_t<std::is_same<Type, S0>::value, int> = 0>
		constexpr int index_of_type()
	{
		return symphas::internal::index_of_type_match<I, S0, Ss...>::template value<Type, Sn>;
	}

	//! Return the index of the given type in the type list.
	/*!
	 * Given the type list, identify the index of the `I`-th chosen type. If
	 * the type is not identified, then `-1` will be returned.
	 *
	 * \tparam Type The type to be found in the list.
	 * \tparam I The current running index of the type to find.
	 * \tparam S0 The next type in the type list to compare.
	 * \tparam S0 The second next type in the type list to compare.
	 * \tparam Ss... The remaining types in the list.
	 */
	template<typename Type, size_t I, size_t Sn, typename S0, typename S1, typename ...Ss,
		typename std::enable_if_t<!std::is_same<Type, S0>::value, int> = 0>
		constexpr int index_of_type()
	{
		return index_of_type<Type, I, Sn, S1, Ss...>();
	}


	template<size_t I, typename... Ts>
	struct type_of_index;


	template<typename T0, typename... Ts>
	struct type_of_index<0, T0, Ts...>
	{
		using type = T0;
	};

	template<size_t I, typename T0, typename... Ts>
	struct type_of_index<I, T0, Ts...>
	{
		using type = typename type_of_index<I - 1, Ts...>::type;
	};

	template<size_t I, typename... Ts>
	struct type_of_index<I, std::tuple<Ts...>>
	{
		using type = typename type_of_index<I, Ts...>::type;
	};

	// **************************************************************************************

	//! Joining a single sequence simply returns it.
	template<typename T, size_t... Ys>
	constexpr auto seq_join(std::integer_sequence<T, Ys...>)
	{
		return std::integer_sequence<T, Ys...>{};
	}

	//! Joins two index sequences.
	template<typename T, size_t... Ys, size_t... Qs>
	constexpr auto seq_join(std::integer_sequence<T, Ys...>, std::integer_sequence<T, Qs...>)
	{
		return std::integer_sequence<T, Ys..., Qs...>{};
	}

	//! Joins two index sequences.
	template<typename T, size_t... Ys, size_t... Qs, typename... Seqs>
	constexpr auto seq_join(std::integer_sequence<T, Ys...>, std::integer_sequence<T, Qs...>, Seqs... seqs)
	{
		return seq_join(std::integer_sequence<T, Ys..., Qs...>{}, seqs...);
	}


	template<typename T, size_t... Ys>
	constexpr auto seq_neg(std::integer_sequence<int, Ys...>)
	{
		return std::integer_sequence<int, -Ys...>{};
	}

	template<size_t... Ys>
	constexpr auto seq_neg(std::index_sequence<Ys...>)
	{
		return std::integer_sequence<int, -int(Ys)...>{};
	}


	//! Adding a single sequence simply returns it.
	template<typename T, size_t... Ys>
	constexpr auto seq_add(std::integer_sequence<T, Ys...>)
	{
		return std::integer_sequence<T, Ys...>{};
	}

	//! Adds two index sequences.
	/*!
	 * The values are added pointwise, between two index sequences. When
	 * sequences are not of equal size, the shortest one is considered
	 * to have 0s in the remaining entries, so that a sequence equal in length
	 * to the longest sequence is always returned.
	 */
	template<typename T, size_t... Ys, size_t... Qs>
	constexpr auto seq_add(std::integer_sequence<T, Ys...>, std::integer_sequence<T, Qs...>)
	{
		if constexpr (sizeof...(Ys) == sizeof...(Qs))
		{
			return std::integer_sequence<T, (Qs + Ys)...>{};
		}
		else if constexpr (sizeof...(Ys) > sizeof...(Qs))
		{
			return seq_join(
				seq_add(
					symphas::lib::get_seq_lt<sizeof...(Qs)>(std::integer_sequence<T, Ys...>{}),
					std::integer_sequence<T, Qs...>{}),
				symphas::lib::get_seq_ge<sizeof...(Qs)>(std::integer_sequence<T, Ys...>{}));
		}
		else
		{
			return seq_join(
				seq_add(
					std::integer_sequence<T, Ys...>{},
					symphas::lib::get_seq_lt<sizeof...(Ys)>(std::integer_sequence<T, Qs...>{})),
				symphas::lib::get_seq_ge<sizeof...(Ys)>(std::integer_sequence<T, Qs...>{}));
		}
	}

	//! Adds multiple index sequences.
	/*!
	 * The values are added pointwise, between two index sequences. When
	 * sequences are not of equal size, the shortest one is considered
	 * to have 0s in the remaining entries, so that a sequence equal in length
	 * to the longest sequence is always returned.
	 */
	template<typename T, size_t... Ys, size_t... Qs, typename... Seqs>
	constexpr auto seq_add(std::integer_sequence<T, Ys...>, std::integer_sequence<T, Qs...>, Seqs... seqs)
	{
		return seq_add(seq_add(std::integer_sequence<T, Ys...>{}, std::integer_sequence<T, Qs...>{}), seqs...);
	}


	//! Subtracting a single sequence simply returns it.
	template<typename T, size_t... Ys>
	constexpr auto seq_sub(std::integer_sequence<T, Ys...>)
	{
		return std::integer_sequence<T, Ys...>{};
	}

	//! Subtracts multiple index sequences from the first one.
	/*!
	 * The values are subtracted pointwise, between two index sequences. When
	 * sequences are not of equal size, the shortest one is considered
	 * to have 0s in the remaining entries, so that a sequence equal in length
	 * to the longest sequence is always returned. The sequences are
	 * all subtracted from the first one.
	 */
	template<typename T, size_t... Ys, size_t... Qs, typename... Seqs>
	constexpr auto seq_sub(std::integer_sequence<T, Ys...>, std::integer_sequence<T, Qs...>, Seqs... seqs)
	{
		return seq_sub(seq_add(std::integer_sequence<int, int(Ys)...>{}, seq_neg(std::integer_sequence<T, Qs...>{})), seqs...);
	}


	//! The index sequence result type of joining multiple sequences.
	template<typename Seq, typename... Seqs>
	struct seq_join_result
	{
		using type = decltype(seq_join(std::declval<Seq>(), std::declval<Seqs>()...));
	};

	//! The index sequence result type of adding multiple sequences.
	template<typename Seq, typename... Seqs>
	struct seq_add_result
	{
		using type = decltype(seq_add(std::declval<Seq>(), std::declval<Seqs>()...));
	};

	//! The index sequence result type of adding multiple sequences.
	template<typename Seq, typename... Seqs>
	struct seq_sub_result
	{
		using type = decltype(seq_sub(std::declval<Seq>(), std::declval<Seqs>()...));
	};

	//! Alias for the join result of multiple sequences.
	template<typename... Seqs>
	using seq_join_t = typename seq_join_result<Seqs...>::type;

	//! Alias for the add result of multiple sequences.
	template<typename... Seqs>
	using seq_add_t = typename seq_add_result<Seqs...>::type;

	//! Alias for the add result of multiple sequences.
	template<typename... Seqs>
	using seq_sub_t = typename seq_sub_result<Seqs...>::type;


	template<size_t N, typename T>
	struct value_in_seq;

	template<size_t N>
	struct value_in_seq<N, std::index_sequence<>>
	{
		static const bool value = false;
	};

	template<size_t N, size_t I0, size_t... Is>
	struct value_in_seq<N, std::index_sequence<I0, Is...>>
	{
		static const bool value = (I0 == N) || (value_in_seq<N, std::index_sequence<Is...>>::value);
	};




	namespace internal
	{


		template<typename Seq, typename... Seqs>
		static size_t constexpr seq_len_product()
		{
			if constexpr (sizeof...(Seqs) == 0)
			{
				return Seq::size();
			}
			else
			{
				return Seq::size() * seq_len_product<Seqs...>();
			}
		}


		template<typename T, size_t V>
		static size_t constexpr get_value_from_seq(std::integer_sequence<T, V>)
		{
			return V;
		}

		template<size_t N, typename T, T... Es, typename = std::enable_if_t<(N < sizeof...(Es)), int>>
		static size_t constexpr seq_value(std::integer_sequence<T, Es...>)
		{
			return get_value_from_seq(std::tuple_element_t<N, std::tuple<std::integer_sequence<T, Es>...>>{});
		}


		template<typename T>
		struct CrossProductFunctions
		{
			template<T... Es>
			using seq_t = std::integer_sequence<T, Es...>;

			// **********************************************************
			// Expand the cross list into a full tuple of sequences representing the row of combined values.

			template<T E1>
			static auto constexpr expand2(seq_t<E1>, seq_t<>)
			{
				return std::make_tuple();
			}

			template<T E1, T E2, T... E2s>
			static auto constexpr expand2(seq_t<E1>, seq_t<E2, E2s...>)
			{
				return std::tuple_cat(
					std::make_tuple(seq_t<E1, E2>{}),
					expand2(seq_t<E1>{}, seq_t<E2s...>{}));
			}

			template<T E1, T... E1s>
			static auto constexpr expand1(seq_t<E1, E1s...>, seq_t<>)
			{
				return std::tuple<>{};
			}

			template<T E2, T... E2s>
			static auto constexpr expand1(seq_t<>, seq_t<E2, E2s...>)
			{
				return std::tuple<>{};
			}

			static auto constexpr expand1(seq_t<>, seq_t<>)
			{
				return std::tuple<>{};
			}

			template<T E1, T... E1s, T E2, T... E2s>
			static auto constexpr expand1(seq_t<E1, E1s...>, seq_t<E2, E2s...>)
			{
				return std::tuple_cat(
					expand2(seq_t<E1>{}, seq_t<E2, E2s...>{}),
					expand1(seq_t<E1s...>{}, seq_t<E2, E2s...>{}));
			}

			template<T... E1s>
			static auto constexpr expand33(seq_t<E1s...>, seq_t<>)
			{
				return std::make_tuple();
			}


			template<T... E1s, T E2, T... E2s>
			static auto constexpr expand33(seq_t<E1s...>, seq_t<E2, E2s...>)
			{
				return std::tuple_cat(
					std::make_tuple(symphas::lib::seq_join(seq_t<E1s...>{}, seq_t<E2>{})),
					expand33(seq_t<E1s...>{}, seq_t<E2s...>{}));
			}

			template<T E, T... Es>
			static auto constexpr expand22(std::tuple<>, seq_t<E, Es...>)
			{
				return std::tuple<>{};
			}

			template<typename Row, typename... Rows, T E, T... Es>
			static auto constexpr expand22(std::tuple<Row, Rows...>, seq_t<E, Es...>)
			{
				return std::tuple_cat(
					expand33(Row{}, seq_t<E, Es...>{}),
					expand22(std::tuple<Rows...>{}, seq_t<E, Es...>{}));
			}

			template<typename Row, typename... Rows, T E, T... Es, typename List0, typename... Lists>
			static auto constexpr expand22(std::tuple<Row, Rows...>, seq_t<E, Es...>, List0, Lists...)
			{
				return expand22(
					std::tuple_cat(
						expand33(Row{}, seq_t<E, Es...>{}),
						expand22(std::tuple<Rows...>{}, seq_t<E, Es...>{})),
					List0{}, Lists{}...);
			}

			template<T E, T... Es, typename List0, typename... Lists>
			static auto constexpr expand11(seq_t<E, Es...>, List0, Lists...)
			{
				return expand22(expand1(seq_t<E, Es...>{}, List0{}), Lists{}...);
			}


			// **********************************************************
			// Selects only a single row without constructing the whole cross list.

			template<size_t N>
			static auto constexpr select(seq_t<>)
			{
				return std::integer_sequence<T>{};
			}

			template<size_t N, T E, T... Es>
			static auto constexpr select(seq_t<E, Es...>)
			{
				return std::integer_sequence < T, seq_value<N>(seq_t<E, Es...>{}) > {};
			}

			template<size_t N, T E, T... Es, typename Seq, typename... Seqs, size_t L = seq_len_product<Seq, Seqs...>(), size_t N0 = N / L, size_t N1 = N - N0 * L>
			static auto constexpr select(seq_t<E, Es...>, Seq, Seqs...)
			{
				return symphas::lib::seq_join(std::integer_sequence < T, seq_value<N0>(seq_t<E, Es...>{}) > {}, select<N1>(Seq{}, Seqs{}...));
			}
		};
	}

	/*!
	 * Generate the cross product or cross join of all the numeric elements in the
	 * provided std::integer_sequence types.
	 *
	 * \tparam Lists The std::integer_sequence containing a list of values which are
	 * cross joined.
	 */
	template<typename... Lists>
	struct CrossProductList;

	template<>
	struct CrossProductList<>
	{
		//using type = std::tuple<>;
		static const size_t count = 0;
		static const size_t rank = 0;
	};

	template<typename T, T... Es>
	struct CrossProductList<std::integer_sequence<T, Es...>>
	{
		//using type = std::tuple<std::integer_sequence<T, Es>...>;

		static const size_t count = sizeof...(Es);
		static const size_t rank = 1;

		template<size_t N>
		static const size_t size = std::integer_sequence<T, Es...>::size();

		template<size_t N, typename = std::enable_if_t<(N < count), int>>
		using row = std::tuple_element_t<N, std::tuple<std::integer_sequence<T, Es>...>>;

	};


	template<typename T, T... E1s, T... E2s>
	struct CrossProductList<std::integer_sequence<T, E1s...>, std::integer_sequence<T, E2s...>>
	{
		//using type = decltype(internal::CrossProductFunctions<T>::expand1(
		//	std::declval<std::integer_sequence<T, E1s...>>(),
		//	std::declval<std::integer_sequence<T, E2s...>>()));

		static const size_t count = (sizeof...(E1s) * sizeof...(E2s));
		static const size_t rank = 2;

		template<size_t N>
		static const size_t size = std::tuple_element_t<N, std::tuple<std::integer_sequence<T, E1s...>, std::integer_sequence<T, E2s...>>>::size();

		template<size_t N, typename = std::enable_if_t<(N < count), int>>
		using row = decltype(internal::CrossProductFunctions<T>::select<N>(
			std::declval<std::integer_sequence<T, E1s...>>(),
			std::declval<std::integer_sequence<T, E2s...>>()));
	};

	template<typename T, T... Es, typename List1, typename List2, typename... Lists>
	struct CrossProductList<std::integer_sequence<T, Es...>, List1, List2, Lists...>
	{
		//using type = decltype(internal::CrossProductFunctions<T>::expand11(
		//	std::declval<std::integer_sequence<T, Es...>>(),
		//	std::declval<List1>(),
		//	std::declval<List2>(),
		//	std::declval<Lists>()...));

		static const size_t count = internal::seq_len_product<std::integer_sequence<T, Es...>, List1, List2, Lists...>();
		static const size_t rank = 3 + sizeof...(Lists);

		template<size_t N>
		static const size_t size = std::tuple_element_t<N, std::tuple<std::integer_sequence<T, Es...>, List1, List2, Lists...>>::size();

		template<size_t N, typename = std::enable_if_t<(N < count), int>>
		using row = decltype(internal::CrossProductFunctions<T>::template select<N>(
			std::declval<std::integer_sequence<T, Es...>>(),
			std::declval<List1>(),
			std::declval<List2>(),
			std::declval<Lists>()...));

	};



	//! Concatenates two tuple lists together into one tuple type.
	/*
	 * In general, the type is not defined; tuples need to be concatenated.
	 */
	template<typename A, typename B>
	struct combine_types
	{
		using type = std::tuple<A, B>;
	};

	//! Concatenates two tuple lists together into one tuple type.
	/*
	 * The underlying types of two tuples are joined together in one tuple.
	 */
	template<typename... G1s, typename... G2s>
	struct combine_types<std::tuple<G1s...>, std::tuple<G2s...>>
	{
		using type = std::tuple<G1s..., G2s...>;
	};

	//! Collect and count like types from a tuple list.
	/*!
	 * Give the list of types, the list is made unique and number of times
	 * each type appears in the list is repeated.
	 */
	template<typename... Gs>
	struct cc_like_types;

	namespace
	{

		template<typename G>
		struct cc_like_types_apply
		{
			using type = G;
			using count = std::index_sequence<1>;
			using discard = std::tuple<>;
		};

		template<typename G0>
		struct cc_like_types_apply<std::tuple<G0>>
		{
			using type = std::tuple<G0>;
			using count = std::index_sequence<1>;
			using discard = std::tuple<>;
		};

		template<typename G0, typename... Gs>
		struct cc_like_types_apply<std::tuple<G0, G0, Gs...>>
		{
		protected:
			using cft0 = cc_like_types_apply<std::tuple<G0, Gs...>>;

		public:
			using type = typename cft0::type;
			using count = seq_add_t<std::index_sequence<1>, typename cft0::count>;
			using discard = typename cft0::discard;
		};

		template<typename G0, typename G1, typename... Gs>
		struct cc_like_types_apply<std::tuple<G0, G1, Gs...>>
		{
		protected:
			using cft0 = cc_like_types_apply<std::tuple<G0, Gs...>>;

		public:
			using type = typename cft0::type;
			using count = typename cft0::count;
			using discard = typename combine_types<std::tuple<G1>, typename cft0::discard>::type;
		};
	}

	template<>
	struct cc_like_types<std::tuple<>>
	{
		using type = std::tuple<>;
		using count = std::index_sequence<>;
	};

	template<typename G0, typename... Gs>
	struct cc_like_types<std::tuple<G0, Gs...>>
	{
	protected:

		using cc0 = cc_like_types_apply<std::tuple<G0, Gs...>>;
		using cc1 = cc_like_types<typename cc0::discard>;

	public:

		using type = typename combine_types<
			typename cc0::type,
			typename cc1::type
		>::type;

		using count = seq_join_t<
			typename cc0::count,
			typename cc1::count>;
	};

	template<typename... Gs>
	struct cc_like_types
	{
		using type = typename cc_like_types<std::tuple<Gs...>>::type;
		using count = typename cc_like_types<std::tuple<Gs...>>::count;
	};


	//! Returns the unique list of types from the list.
	/*!
	 * Returns the unique list of types is returned from joining together
	 * the list of types.
	 */
	template<typename A, typename B>
	struct combine_types_unique
	{
		using type = typename cc_like_types<typename combine_types<A, B>::type>::type;
	};

	namespace
	{
		// helper function that iterates over the tuple
		template<typename... As, typename... Bs, size_t... Is>
		auto unfurl_tuple(std::tuple<std::pair<As, Bs>...> const& ts, std::index_sequence<Is...>)
		{
			auto first_elements = symphas::lib::make_tuple(std::get<Is>(ts).first...);
			auto second_elements = symphas::lib::make_tuple(std::get<Is>(ts).second...);
			return std::make_pair(first_elements, second_elements);
		}

		// helper function that takes index sequence and rebuilds the tuple from the 
		// given index onwards
		template<size_t I, typename... Ts, size_t... Is>
		auto get_tuple_ge(std::tuple<Ts...> const& ts, std::index_sequence<Is...>)
		{
			return symphas::lib::make_tuple(std::get<I + Is>(ts)...);
		}

		// helper function that takes index sequence and rebuilds the tuple up to and
		// not including the given index
		template<typename... Ts, size_t... Is>
		auto get_tuple_lt(std::tuple<Ts...> const& ts, std::index_sequence<Is...>)
		{
			return symphas::lib::make_tuple(std::get<Is>(ts)...);
		}


		template<typename T, size_t... Ns>
		auto get_seq_from_tuple(std::tuple<std::integer_sequence<T, Ns>...> const&)
		{
			return std::integer_sequence<T, Ns...>{};
		}


		template<size_t I, size_t... Ns, size_t... Is>
		auto get_seq_ge(std::tuple<std::index_sequence<Ns>...> const& ts, std::index_sequence<Is...>)
		{
			return get_seq_from_tuple(get_tuple_ge<I>(ts, std::index_sequence<Is...>{}));
		}

		template<size_t... Ns, size_t... Is>
		auto get_seq_lt(std::tuple<std::index_sequence<Ns>...> const& ts, std::index_sequence<Is...>)
		{
			return get_seq_from_tuple(get_tuple_lt(ts, std::index_sequence<Is...>{}));
		}


		// **************************************************************************************


		//! End the recursion when there is nothing remaining to match against.
		template<typename empty_type, typename match_type, typename A, typename... match_A_types, typename... As, typename B>
		auto unzip_pot_one(
			std::tuple<std::pair<match_type, A>, std::pair<match_A_types, As>...> const& as, 
			std::pair<match_type, B> const& b0)
		{
			auto matched = std::make_tuple(std::make_pair(b0.first, std::make_pair(std::get<0>(as).second, b0.second)));
			auto rest = symphas::lib::get_tuple_ge<1>(as);

			return std::make_pair(matched, rest);
		}

		//! End the recursion when there is nothing remaining to match against.
		template<typename empty_type, typename match_B_type, typename B>
		auto unzip_pot_one(
			std::tuple<> const&, 
			std::pair<match_B_type, B> const& b0)
		{
			auto matched = std::make_tuple(std::make_pair(b0.first, std::make_pair(empty_type{}, b0.second)));
			return std::make_pair(matched, std::make_tuple());
		}

		//! End the recursion when there is nothing remaining to match against.
		template<typename empty_type, typename match_A_type, typename A, typename... match_A_types, typename... As,
			typename match_B_type, typename B, typename std::enable_if<!std::is_same<match_A_type, match_B_type>::value, int>::type = 0>
		auto unzip_pot_one(
			std::tuple<std::pair<match_A_type, A>, std::pair<match_A_types, As>...> const& as, 
			std::pair<match_B_type, B> const& b0)
		{
			auto [matched, unmatched] = unzip_pot_one<empty_type>(symphas::lib::get_tuple_ge<1>(as), b0);
			auto rest = std::tuple_cat(std::make_tuple(std::get<0>(as)), unmatched);

			return std::make_pair(matched, rest);
		}


		//! Construct the `b` data element.
		/*!
		 * Each pair of the tuple has the match type as the first member and
		 * a pair as the second member, which contains the content from the
		 * given `as` as the first element and an default constructed
		 * instance of the type `empty_type` as the second element. The 
		 * resulting list of pairs is thus the same length as the input list. 
		 */
		template<typename empty_type, typename... match_B_types, typename... Bs, size_t... Is>
		auto unzip_pot_b_apply(
			std::tuple<std::pair<match_B_types, Bs>...> const& bs, 
			std::index_sequence<Is...>)
		{
			return std::make_tuple(
				std::make_pair(
					std::get<Is>(bs).first, 
					std::make_pair(empty_type{}, std::get<Is>(bs).second)
				)...);
		}

		//! Construct the `a` data element.
		/*!
		 * Each pair of the tuple has the match type as the first member and
		 * a pair as the second member, which contains the content from the
		 * given `bs` as the second element and an default constructed
		 * instance of the type `empty_type` as the first element. The resulting
		 * list of pairs is thus the same length as the input list. 
		 */
		template<typename empty_type, typename... match_A_types, typename... As, size_t... Is>
		auto unzip_pot_a_apply(
			std::tuple<std::pair<match_A_types, As>...> const& as, 
			std::index_sequence<Is...>)
		{
			return std::make_tuple(
				std::make_pair(
					std::get<Is>(as).first, 
					std::make_pair(std::get<Is>(as).second, empty_type{})
				)...);
		}

		//! Construct the `b` data element.
		template<typename empty_type, typename... match_B_types, typename... Bs>
		auto unzip_pot_b(std::tuple<std::pair<match_B_types, Bs>...> const& bs)
		{
			using bs_type = std::tuple<std::pair<match_B_types, Bs>...>;
			constexpr size_t bs_len = std::tuple_size<bs_type>::value;

			return unzip_pot_b_apply<empty_type>(bs, std::make_index_sequence<bs_len>{});
		}

		//! Construct the `a` data element.
		template<typename empty_type, typename... match_A_types, typename... As>
		auto unzip_pot_a(std::tuple<std::pair<match_A_types, As>...> const& as)
		{
			using as_type = std::tuple<std::pair<match_A_types, As>...>;
			constexpr size_t as_len = std::tuple_size<as_type>::value;

			return unzip_pot_a_apply<empty_type>(as, std::make_index_sequence<as_len>{});
		}




		//! Apply the algorithm to unzip a list based on the types.
		/*!
		 * Entry point for the recursive unzip_pot function, carries through 
		 * the tuple as it is being unzipped, 
		 * and the final return operation unfurls the tuple.
		 */
		template<typename empty_type, typename... Rs>
		auto unzip_pot_carry_apply(
			std::tuple<Rs...> const& matched, 
			std::tuple<> const&, 
			std::tuple<> const&)
		{
			return matched;
		}

		//! Apply the algorithm to unzip a list based on the types.
		/*!
		 * See symphas::lib::unzip_pot_carry_apply.
		 */
		template<typename empty_type, typename... Rs, typename match_A_type, typename A, typename... match_A_types, typename... As>
		auto unzip_pot_carry_apply(
			std::tuple<Rs...> const& matched,
			std::tuple<std::pair<match_A_type, A>, std::pair<match_A_types, As>...> const& as, 
			std::tuple<> const&)
		{
			auto pot_as = unzip_pot_a<empty_type>(as);
			auto matched_list = std::tuple_cat(matched, pot_as);

			return matched_list;
		}


		//! Apply the algorithm to unzip a list based on the types.
		/*!
		 * See symphas::lib::unzip_pot_carry_apply.
		 */
		template<typename empty_type, typename... Rs, typename match_B_type, typename B, typename... match_B_types, typename... Bs>
		auto unzip_pot_carry_apply(
			std::tuple<Rs...> const& matched, 
			std::tuple<> const&, 
			std::tuple<std::pair<match_B_type, B>, std::pair<match_B_types, Bs>...> const& bs)
		{
			auto pot_bs = unzip_pot_b<empty_type>(bs);
			auto matched_list = std::tuple_cat(matched, pot_bs);

			return matched_list;
		}


		//! Apply the algorithm to unzip a list based on the types.
		/*!
		 * See symphas::lib::unzip_pot_carry_apply.
		 */
		template<typename empty_type, typename... Rs, typename match_A_type, typename A, typename... match_A_types, typename... As, typename match_B_type, typename B, typename... match_B_types, typename... Bs>
		auto unzip_pot_carry_apply(
			std::tuple<Rs...> const& matched,
			std::tuple<std::pair<match_A_type, A>, std::pair<match_A_types, As>...> const& as, 
			std::tuple<std::pair<match_B_type, B>, std::pair<match_B_types, Bs>...> const& bs)
		{
			auto [one_matched, as_rest] = unzip_pot_one<empty_type>(as, std::get<0>(bs));
			auto match_append = std::tuple_cat(matched, one_matched);
			auto bs_rest = symphas::lib::get_tuple_ge<1>(bs);
			
			return unzip_pot_carry_apply<empty_type>(
				match_append, 
				as_rest,
				bs_rest);
		}

	}


	//! Implements function to return identity.
	/*!
	 * Implements a `static constexpr` function in order to return the known
	 * identity at compile time. Each specialization of this struct
	 * implements the identity for its respective type.
	 */
	template<typename T>
	struct identity;

	//! Identity for the scalar type, always equal to `1.0`.
	template<>
	struct identity<scalar_t>
	{
		constexpr scalar_t operator()()
		{
			return 1.0;
		}
	};

	//! Identity for the complex type, always equal to `1.0 + 0.0*i`.
	template<>
	struct identity<complex_t>
	{
		constexpr complex_t operator()()
		{
			return { 1.0, 0.0 };
		}
	};

	//! The identity specialization for VectorValue.
	/*!
	 * An identity for a vector is not defined. However, it is naively
	 * assumed that an identity element is considered in the context of
	 * point-wise scalar multiplication. There is no
	 * identity for the dot product of two vectors. Instead, the identity
	 * is defined such that the dot product
	 * of two vector identities results in
	 * the scalar identity of the underlying type, multiplied by `D`.
	 * 
	 * \tparam T The underlying type of the vector.
	 * \tparam D The dimension of the vector.
	 */
	template<typename T, size_t D>
	struct identity<VectorValue<T, D>>
	{
		constexpr VectorValue<T, D> operator()()
		{
			VectorValue<T, D> idty;
			std::fill(idty.v, idty.v + D, identity<T>{}());
			return idty;
		}
	};

	//! The identity for the alias type ::vector_t.
	/*!
	 * An identity for a vector is not well defined. See
	 * identity<VectorValue>.
	 * 
	 * The identity is defined for the vector type given by the alias
	 * ::vector_t. All that is assumed about the underlying type of
	 * ::vector_t is that there is an access operator `[]` and that it
	 * can be default constructed.
	 * 
	 * Using that,
	 * each element of the ::vector_t data is set to the identity given by
	 * the underlying data type of the vector, `scalar_t`.
	 * 
	 * \tparam D The dimension of the vector.
	 */
	template<size_t D>
	struct identity<vector_t<D>>
	{
		constexpr vector_t<D> operator()()
		{
			vector_t<D> idty;
			for (iter_type i = 0; i < D; ++i)
			{
				idty[i] = identity<scalar_t>{}();
			}
			return idty;
		}
	};

	//! Helper function to return the identity.
	/*!
	 * Determines the identity of the given type by using the specialized
	 * structure. An identity is provided for commonly used SymPhas types.
	 */
	template<typename T>
	inline auto constexpr get_identity()
	{
		return identity<T>{}();
	}




	//! Generalized function to expand the given list.
	/*!
	 * All the elements from the input are added to the out array, which is
	 * appropriately resized in order to accommodate the appended elements.
	 * Since the array is resized, the array needs to be taken by reference to
	 * be properly reassigned.
	 *
	 * \param[in] in The buffer list which will be copied to out.
	 * \param[out] out The elements from the buffer parameter `in` are added to
	 * the end of this list and its memory is reallocated to fit the new size.
	 * \param in_len The total length of the input, which is how many elements
	 * will be appended to the output.
	 * \param out_len The current length of out, so that the memory can be
	 * appropriately reallocated with a larger length.
	 */
	template<typename T>
	void expand_append_list(const T* in, T*& out, size_t in_len, size_t out_len)
	{
		size_t total_len = out_len + in_len;
		T* out_buffer = new T[total_len];
		std::copy(out, out + out_len, out_buffer);
		std::copy(in, in + in_len, out_buffer + out_len);

		delete[] out;
		out = out_buffer;

	}


	//! Generalized function to expand the given list.
	/*!
	 * A given number of copies of the provided input value are added to the out 
	 * array, which is appropriately resized in order to accommodate the 
	 * appended elements. Since the array is resized, the array needs to be 
	 * taken by reference to be properly reassigned.
	 *
	 * \param[in] in The element of which copies will be appended to the output.
	 * \param[out] out The elements from the buffer parameter `in` are added to
	 * the end of this list and its memory is reallocated to fit the new size.
	 * \param in_len The total length of the input, which is how many elements
	 * will be appended to the output.
	 * \param out_len The current length of out, so that the memory can be
	 * appropriately reallocated with a larger length.
	 */
	template<typename T>
	void expand_append_list(T in, T*& out, size_t in_len, size_t out_len = 0)
	{
		size_t total_len = out_len + in_len;
		T* out_buffer = new T[total_len];
		if (out_len > 0)
		{
			std::copy(out, out + out_len, out_buffer);
		}
		std::fill(out_buffer + out_len, out_buffer + total_len, in);
		
		delete[] out;
		out = out_buffer;
	}



	//! Determine the maximum value from a list of values.
	/*!
	 * Provided a list of values of generic types, the maximum value is
	 * returned by using the comparator operator `>`. Thus, for custom defined
	 * types, this operator needs to be defined.
	 * 
	 * The algorithm for multiple values terminates at a single value.
	 * 
	 * \param val0 Terminating value which is simply returned.
	 * 
	 * \tparam T0 Type of terminating value.
	 */
	template<typename T0>
	constexpr auto max_value(T0 val0)
	{
		return val0;
	}

	//! Determine the maximum value from a list of values.
	/*!
	 * Provided a list of values of generic types, the maximum value is
	 * returned by using the comparator operator `>`. Thus, for custom defined
	 * types, this operator needs to be defined.
	 * 
	 * The return type is derived by determining the result type of addition
	 * between instances of all types.
	 * 
	 * \param val0 One of the values to compare.
	 * \param val1 A second value to compare.
	 * \param vals Further list of values which will be compared to determine
	 * the maximum value.
	 * 
	 * \tparam T0 Type of first value.
	 * \tparam T1 Type of second value.
	 * \tparam Ts List of types of the additional compared value list.
	 */
	template<typename T0, typename T1, typename... Ts>
	constexpr auto max_value(T0 val0, T1 val1, Ts ...vals) -> add_result_t<T0, T1, Ts...>
	{
		auto amin = max_value(val1, vals...);
		return (val0 > amin) ? val0 : amin;
	}

	//! Determine the minimum value from a list of values.
	/*!
	 * Provided a list of values of generic types, the minimum value is
	 * returned by using the comparator operator `>`. Thus, for custom defined
	 * types, this operator needs to be defined.
	 *
	 * The algorithm for multiple values terminates at a single value.
	 *
	 * \param val0 Terminating value which is simply returned.
	 *
	 * \tparam T0 Type of terminating value.
	 */
	template<typename T0>
	constexpr auto min_value(T0 val0)
	{
		return val0;
	}

	//! Determine the minimum value from a list of values.
	/*!
	 * Provided a list of values of generic types, the minimum value is
	 * returned by using the comparator operator `>`. Thus, for custom defined
	 * types, this operator needs to be defined.
	 *
	 * The return type is derived by determining the result type of addition
	 * between instances of all types.
	 *
	 * \param val0 One of the values to compare.
	 * \param val1 A second value to compare.
	 * \param vals Further list of values which will be compared to determine
	 * the minimum value.
	 *
	 * \tparam T0 Type of first value.
	 * \tparam T1 Type of second value.
	 * \tparam Ts List of types of the additional compared value list.
	 */
	template<typename T0, typename T1, typename... Ts>
	constexpr auto min_value(T0 val0, T1 val1, Ts ...vals) -> add_result_t<T0, T1, Ts...>
	{
		auto amin = min_value(val1, vals...);
		return (val0 < amin) ? val0 : amin;
	}


	//! Transform a string into a file name format.
	/* 
	 * Transform a string to lowercase and remove spaces. Additional
	 * transformations are also made to convert any string to an appropriate
	 * file name. The string is not changed in place, but is copied
	 * into the given output, which is expected to be at least as long as
	 * input.
	 * 
	 * \param[in] in The input string which is transformed.
	 * \param[out] out The result of converting the input to a file name format.
	 * The function expects that this string is already allocated to a length
	 * at least as long as the input.
	 */
	void to_file_name(char const* in, char* out, size_t count = 0);

	//! Transform the string in place to lowercase.
	/*!
	 * The given string is manipulated so that all characters are converted
	 * into lowercase.
	 * 
	 * \param[inout] str The input string which is manipulated in place.
	 */
	void to_lower(char* str);

	//! Transform the string in place to lowercase.
	/*!
	 * The given string is transformed into the output string so that all 
	 * characters are converted into lowercase.
	 * 
	 * \param[in] str The input string which is converted to lowercase.
	 * \param[out] out The result of converting the input to lowercase.
	 */
	void to_lower(const char* str, char* out);

	//! Transform the string in place to uppercase.
	/*!
	 * The given string is manipulated so that all characters are converted
	 * into uppercase.
	 *
	 * \param[inout] str The input string which is manipulated in place.
	 */
	void to_upper(char* str);

	//! Transform the string in place to uppercase.
	/*!
	 * The given string is transformed into the output string so that all
	 * characters are converted into uppercase.
	 *
	 * \param[in] str The input string which is converted to uppercase.
	 * \param[out] out The result of converting the input to uppercase.
	 */
	void to_upper(const char* str, char* out);

	//! Trim whitespace from the beginning and end.
	/*!
	 * The given string is manipulated so that all whitespace is removed from
	 * the beginning and end. The whitespace characters includes `\n` `\t` and
	 * `\r` in addition to space, ` `.
	 *
	 * \param[inout] str The input string which is manipulated in place.
	 */
	void str_trim(char* str);

	//! Trim whitespace from the beginning and end.
	/*!
	 * The given string is manipulated so that all whitespace is removed from
	 * the beginning and end. The whitespace characters includes `\n` `\t` and
	 * `\r` in addition to space, ` `.
	 *
	 * \param[in] str The input string which is trimmed.
	 * \param[out] out The result of trimming the input.
	 */
	void str_trim(const char* str, char* out);


	//! Returns offset to next non-digit character, which may the first one.
	/*!
	 * Returns offset to a position in the given name string to the next
	 * character which is not a digit, which can also be the end of the string.
	 * In the case that the offset is the end of the string, then the offset
	 * will be equal to the result of `std::strlen`. So this case, the result
	 * of dereferencing the input by the result would be the terminating
	 * character.
	 * 
	 * \param content The pointer to the start of the string that will be iterated
	 * over until the current index no longer refers to a digit.
	 */
	iter_type pos_after_digit(const char* content);

	//! Returns offset to next character after matching the start of the string.
	/*!
	 * Returns offset to a position in the given string to the next
	 * character that follows the provided token, if the token exists. This
	 * can also be the end of the string, see pos_after_digit. In the case the
	 * input does not start with the given token, the function will return
	 * 0. Likewise if there is no complete tokens, an offset to a partial match 
	 * (where the token no longer matches) is still considered 0. The match
	 * is performed index by index.
	 * 
	 * \param content the pointer to the start of the string that will be compared with
	 * the token index by index to either return the position or return 0 if
	 * the token is non-matching.
	 * \param token The token which is matched against the start of the given
	 * string.
	 */
	iter_type pos_after_token(const char* content, const char* token);

	//! Returns offset to next character after a token in the set if matching.
	/*!
	 * Returns offset to a position in the given name string to the next
	 * character which follows the end of the token string for one of the tokens
	 * in the provided set. Regarding return value and end of string behaviour, 
	 * see pos_after_digit.
	 *
	 * \param content the pointer to the start of the string that will be compared with
	 * the token index by index to either return the position or return 0 if
	 * the token is non-matching.
	 * \param tokens The tokens which are matched against the start of the given
	 * string.
	 */
	iter_type pos_after_token(const char* content, std::initializer_list<const char*> tokens);



	//! Return the number of digits in a positive integer.
	/*!
	 * Count the number of digits in a positive integer. Returns the
	 * value as a compile time constant.
	 */
	template<size_t N>
	size_t constexpr num_digits()
	{
		if (N <= 9)
		{
			return 1;
		}

		return 1 + num_digits<N / 10>();
	}

	//! Return the number of digits in a positive integer.
	/*!
	 * Count the number of digits in a positive integer.
	 */
	size_t num_digits(int n);
	size_t num_digits(size_t n);


	//! Compute the power between two positive integers.
	/*!
	 * Returns a compile time constant result for taking the power of one
	 * positive integer to the other.
	 * 
	 * \param b Base positive integer.
	 * \param p Power positive integer.
	 */
	constexpr size_t power_ull(size_t b, size_t p)
	{
		size_t e = 1;
		for (iter_type i = 0; i < p; ++i)
		{
			e *= b;
		}
		return e;
	}



	inline iter_type next_block_i(iter_type i, len_type full_len, len_type block_count)
	{
		iter_type
			count = full_len / block_count,
			rem = full_len - (count * block_count);

		return count * i + std::min(i, rem);
	}

	// ****************************************************************************************




	//! Applies binary sorting to axis values centered at the origin.
	bool origin_sort(axis_coord_t a, axis_coord_t b, len_type axis_len);
	bool regular_sort(axis_coord_t a, axis_coord_t b, len_type);


	//! Sort the given data series.
	/*!
	 * A given list of data is sorted in place. In particular, the order of
	 * the given data is modified.
	 *
	 * \param data The data which is modified by sorting.
	 * \param L The length of the data series.
	 * \param f A binary function used to order the data.
	 */
	template<typename T, typename F>
	void sort_data(std::vector<std::pair<axis_1d_type, T>>& data, len_type L, F f)
	{
		std::sort(
			data.begin(), data.end(),
			[&](auto a, auto b) {
				return f(a.first, b.first, L); }
		);
	}


	template<typename T, typename F>
	void sort_data(std::vector<std::pair<axis_2d_type, T>>& data, len_type L, len_type M, F f)
	{
		std::sort(
			data.begin(), data.end(),
			[&](auto a, auto b) {
				return (a.first[1] == b.first[1])
					? f(a.first[0], b.first[0], L)
					: f(a.first[1], b.first[1], M); }
		);
	}

	template<typename T, typename F>
	void sort_data(std::vector<std::pair<axis_3d_type, T>>& data, len_type L, len_type M, len_type N, F f)
	{
		std::sort(
			data.begin(), data.end(),
			[&](auto a, auto b) {
				return (a.first[2] == b.first[2] && a.first[1] == b.first[1])
					? f(a.first[0], b.first[0], L)
					: (a.first[2] == b.first[2])
					? f(a.first[1], b.first[1], M)
					: f(a.first[2], b.first[2], N); }
		);
	}

	template<typename T>
	void sort_data(std::vector<std::pair<axis_1d_type, T>>& data, len_type L)
	{
		sort_data(data, L, regular_sort);
	}
	template<typename T>
	void sort_data(std::vector<std::pair<axis_1d_type, T>>& data, len_type L, len_type M)
	{
		sort_data(data, L, M, regular_sort);
	}
	template<typename T>
	void sort_data(std::vector<std::pair<axis_1d_type, T>>& data, len_type L, len_type M, len_type N)
	{
		sort_data(data, L, M, N, regular_sort);
	}


	template<typename T, typename F>
	void sort_data(axis_1d_type* data_x, T* data_y, len_type L, F f)
	{
		iter_type* sort_array = new iter_type[L];
		std::iota(sort_array, sort_array + L, 0);

		std::sort(
			sort_array, sort_array + L,
			[&](auto a, auto b) {
				return f(data_x[a], data_x[b], L); }
		);

		symphas::internal::sort_reorganize<1>(data_x, data_y, sort_array, L);
		delete[] sort_array;
	}

	template<typename T, typename F>
	void sort_data(axis_2d_type* data_x, T* data_y, len_type L, len_type M, F f)
	{
		len_type len = L * M;
		iter_type* sort_array = new iter_type[len];
		std::iota(sort_array, sort_array + len, 0);

		std::sort(
			sort_array, sort_array + len,
			[&](auto a, auto b) {
				return (data_x[a][1] == data_x[b][1])
					? f(data_x[a][0], data_x[b][0], L)
					: f(data_x[a][1], data_x[b][1], M); }
		);

		symphas::internal::sort_reorganize<2>(data_x, data_y, sort_array, len);
		delete[] sort_array;
	}

	template<typename T, typename F>
	void sort_data(axis_3d_type* data_x, T* data_y, len_type L, len_type M, len_type N, F f)
	{
		len_type len = L * M * N;
		iter_type* sort_array = new iter_type[len];
		std::iota(sort_array, sort_array + len, 0);

		std::sort(
			sort_array, sort_array + len,
			[&](auto a, auto b) {
				return (data_x[a][2] == data_x[b][2] && data_x[a][1] == data_x[b][1])
					? f(data_x[a][0], data_x[b][0], L)
					: (data_x[a][2] == data_x[b][2])
					? f(data_x[a][1], data_x[b][1], M)
					: f(data_x[a][2], data_x[b][2], N); }
		);

		symphas::internal::sort_reorganize<3>(data_x, data_y, sort_array, len);
		delete[] sort_array;
	}

	template<typename X, typename T>
	void sort_data(X&& data_x, T&& data_y, len_type L)
	{
		sort_data(data_x, data_y, L, regular_sort);
	}
	template<typename X, typename T>
	void sort_data(X&& data_x, T&& data_y, len_type L, len_type M)
	{
		sort_data(data_x, data_y, L, M, regular_sort);
	}
	template<typename X, typename T>
	void sort_data(X&& data_x, T&& data_y, len_type L, len_type M, len_type N)
	{
		sort_data(data_x, data_y, L, M, N, regular_sort);
	}


	//! The given data is put into a vector of data.
	/*!
	 * The values of the data along the individual axes are provided and then
	 * assorted into a list where the values of the two axes are matched
	 * point wise. This assumes that the two data lists already correspond,
	 * and does not perform any special consideration of dimension.
	 * 
	 * \param data_x The \f$x\f$ axis values of the data.
	 * \param data_y The \f$y\f$ axis values of the data.
	 * \param len The number of values in the data set.
	 */
	template<typename X, typename Y>
	std::vector<std::pair<X, Y>> combine_data(const X* data_x, const Y* data_y, len_type len)
	{
		std::vector<std::pair<X, Y>> out;
		out.reserve(len);

		for (iter_type i = 0; i < len; ++i)
		{
			out.emplace_back(data_x[i], data_y[i]);
		}
		return out;
	}

	template<typename X, typename Y, size_t N>
	std::vector<std::pair<X, Y[N]>> combine_data(const X* data_x, const Y(*data_y)[N], len_type len)
	{
		std::vector<std::pair<X, Y[N]>> out(len);

		for (iter_type i = 0; i < len; ++i)
		{
			out[i].first = data_x[i];
			for (iter_type n = 0; n < N; ++n)
			{
				out[i].second[n] = data_y[i][n];
			}
		}
		return out;
	}





	//! For an list of data, return the dimensions.
	/*!
	 * For an list of data, return the dimensions. The data is assumed to be
	 * sorted.
	 *
	 * \param sorted_data The sorted list of data from which the dimensions
	 * are inferred.
	 */
	template<typename T>
	inline std::tuple<len_type> get_sorted_dimensions(std::vector<std::pair<axis_1d_type, T>> sorted_data)
	{
		return { static_cast<iter_type>(sorted_data.size()) };
	}


	template<typename T>
	inline std::tuple<len_type, len_type> get_sorted_dimensions(std::vector<std::pair<axis_2d_type, T>> sorted_data)
	{
		len_type L = 1, M = 0;
		double y0 = sorted_data.front().first[1];
		for (auto it = sorted_data.begin() + 1; it < sorted_data.end() && y0 == it->first[1]; ++it)
		{
			++L;
		}
		M = static_cast<len_type>(sorted_data.size()) / L;

		return { L, M };

	}



	template<typename T>
	std::tuple<len_type, len_type, len_type> get_sorted_dimensions(std::vector<std::pair<axis_3d_type, T>> sorted_data)
	{
		len_type L = 0, M = 0, N = 0;
		double
			z0 = sorted_data.front().first[2],
			y0 = sorted_data.front().first[1];
		for (auto it = sorted_data.begin() + 1; it < sorted_data.end(); ++it)
		{
			if (y0 != it->first[1])
			{
				L = static_cast<iter_type>(it - sorted_data.begin());
				break;
			}
		}
		for (auto it = sorted_data.begin() + L; it < sorted_data.end(); it += L)
		{
			if (z0 != it->first[2])
			{
				M = static_cast<iter_type>(it - sorted_data.begin()) / L;
				break;
			}
		}
		N = static_cast<len_type>(sorted_data.size()) / (M * L);

		return { L, M, N };

	}


	inline std::tuple<len_type, len_type, len_type> get_sorted_dimensions(const axis_3d_type* sorted_data, len_type len)
	{
		len_type L = 0, M = 0, N = 0;
		double
			z0 = (*sorted_data)[2],
			y0 = (*sorted_data)[1];
		for (auto* it = sorted_data + 1; it < sorted_data + len; ++it)
		{
			if (y0 != (*it)[1])
			{
				L = static_cast<iter_type>(it - sorted_data);
				break;
			}
		}
		for (auto *it = sorted_data + L; it < sorted_data + len; it += L)
		{
			if (z0 != (*it)[2])
			{
				M = static_cast<iter_type>(it - sorted_data) / L;
				break;
			}
		}
		N = len / (M * L);

		return { L, M, N };

	}

	inline std::tuple<len_type, len_type> get_sorted_dimensions(const axis_2d_type* sorted_data, len_type len)
	{
		len_type L = 0, M = 0;
		double y0 = (*sorted_data)[1];
		while (y0 == (*sorted_data++)[1])
		{
			++L;
		}
		M = len / L;

		return { L, M };

	}

	inline std::tuple<len_type> get_sorted_dimensions(const axis_1d_type*, len_type len)
	{
		return { len };
	}

	//! Return the min, max and separation distance of the given axis.
	template<size_t D>
	auto get_axis_mms(const axis_nd_t<D>* axis_values, len_type len)
	{
		std::array<axis_coord_t, D> min{}, max{}, sep{};
		std::fill(min.begin(), min.end(), std::numeric_limits<axis_coord_t>::max());
		std::fill(max.begin(), max.end(), std::numeric_limits<axis_coord_t>::min());
		std::fill(sep.begin(), sep.end(), std::numeric_limits<axis_coord_t>::max());

		for (iter_type i = 0; i < len; ++i)
		{
			for (iter_type n = 0; n < D; ++n)
			{
				double dd = std::abs(min[n] - axis_values[i][n]);
				sep[n] = (dd > 0 && dd < sep[n]) ? dd : sep[n];

				min[n] = (axis_values[i][n] < min[n]) ? axis_values[i][n] : min[n];
				max[n] = (axis_values[i][n] > max[n]) ? axis_values[i][n] : max[n];
			}
		}

		return std::make_tuple(min, max, sep);
	}

	template<>
	inline auto get_axis_mms<1>(const axis_nd_t<1>* axis_values, len_type len)
	{
		std::array<axis_coord_t, 1> min = { std::numeric_limits<axis_coord_t>::max() };
		std::array<axis_coord_t, 1> max = { std::numeric_limits<axis_coord_t>::min() };
		std::array<axis_coord_t, 1> sep = { std::numeric_limits<axis_coord_t>::max() };

		for (iter_type i = 0; i < len; ++i)
		{
			double dd = std::abs(min[0] - axis_values[i]);
			sep[0] = (dd > 0 && dd < sep[0]) ? dd : sep[0];

			min[0] = (axis_values[i] < min[0]) ? axis_values[i] : min[0];
			max[0] = (axis_values[i] > max[0]) ? axis_values[i] : max[0];
		}

		return std::make_tuple(min, max, sep);
	}



	template<size_t D>
	auto get_dimensions(const axis_nd_t<D>* axis_values, len_type len)
	{
		auto [min, max, sep] = get_axis_mms<D>(axis_values, len);
		len_type dim[D];
		for (iter_type n = 0; n < D; ++n)
		{
			dim[n] = static_cast<len_type>(std::round((max[n] - min[n]) / sep[n])) + 1;
		}
		return symphas::internal::dims_as_tuple<D>(dim);
	}


	template<size_t D, typename T>
	auto get_dimensions(std::vector<std::pair<axis_nd_t<D>, T>> const& data)
	{
		axis_nd_t<D>* axis_values = new axis_nd_t<D>[data.size()]{};
		for (iter_type i = 0; i < data.size(); ++i)
		{
			axis_values[i] = data[i].first;
		}

		auto&& dims = get_dimensions<D>(axis_values, static_cast<len_type>(data.size()));
		delete[] axis_values;
		return dims;
	}



	//! Return information on the interval of the given axis data.
	/*!
	 * The given data is expected to be sorted.
	 */
	template<typename T>
	void fill_sorted_ranges(std::vector<std::pair<axis_3d_type, T>> const& data, axis_coord_t(&ranges)[6])
	{
		if (data.size() > 0)
		{
			ranges[0] = data.front().first[0];
			ranges[1] = data.back().first[0];
			ranges[2] = data.front().first[1];
			ranges[3] = data.back().first[1];
			ranges[4] = data.front().first[2];
			ranges[5] = data.back().first[2];
		}
	}

	template<typename T>
	void fill_sorted_ranges(std::vector<std::pair<axis_2d_type, T>> const& data, axis_coord_t(&ranges)[4])
	{
		if (data.size() > 0)
		{
			ranges[0] = data.front().first[0];
			ranges[1] = data.back().first[0];
			ranges[2] = data.front().first[1];
			ranges[3] = data.back().first[1];
		}
	}

	template<typename T>
	void fill_sorted_ranges(std::vector<std::pair<axis_1d_type, T>> const& data, axis_coord_t(&ranges)[2])
	{
		if (data.size() > 0)
		{
			ranges[0] = data.front().first;
			ranges[1] = data.back().first;
		}
	}

	template<size_t D, typename T>
	void fill_sorted_ranges(std::vector<std::pair<axis_nd_t<D / 2>, T>> const& data, axis_coord_t(&ranges)[D], T(&extrema)[2])
	{
		if (data.size() > 0)
		{
			fill_sorted_ranges(data, ranges);
			extrema[0] = std::max_element(data.begin(), data.end(), [&](auto a, auto b) { return a.second > b.second; })->second;
			extrema[1] = std::max_element(data.begin(), data.end(), [&](auto a, auto b) { return a.second < b.second; })->second;
		}	
	}

	


	inline void fill_sorted_ranges(const axis_3d_type* data_x, len_type len, axis_coord_t(&ranges)[6])
	{
		if (len > 0)
		{
			auto [L, M, N] = get_sorted_dimensions(data_x, len);
			ranges[0] = (data_x[0])[0];
			ranges[1] = (data_x[L - 1])[0];
			ranges[2] = (data_x[0])[1];
			ranges[3] = (data_x[L * (M - 1)])[1];
			ranges[4] = (data_x[0])[2];
			ranges[5] = (data_x[L * M * (N - 1)])[2];
		}
	}

	inline void fill_sorted_ranges(const axis_2d_type* data_x, len_type len, axis_coord_t(&ranges)[4])
	{
		if (len > 0)
		{
			auto [L, M] = get_sorted_dimensions(data_x, len);
			ranges[0] = (data_x[0])[0];
			ranges[1] = (data_x[L - 1])[0];
			ranges[2] = (data_x[0])[1];
			ranges[3] = (data_x[L * (M - 1)])[1];
		}
	}

	inline void fill_sorted_ranges(const axis_1d_type* data_x, len_type len, axis_coord_t(&ranges)[2])
	{
		if (len > 0)
		{
			ranges[0] = *(data_x);
			ranges[1] = *(data_x + len - 1);
		}
	}

	template<typename T>
	void fill_sorted_ranges(const T* data, len_type len, T(&ranges)[2])
	{
		if (len > 0)
		{
			auto max = std::max_element(data, data + len, [&](auto a, auto b) { return a < b; });
			auto min = std::max_element(data, data + len, [&](auto a, auto b) { return a > b; });
			ranges[0] = min;
			ranges[1] = max;
		}
	}
	



	template<typename T>
	void fill_ranges(std::vector<std::tuple<iter_type, axis_1d_type, T>> const& data, T(&ranges)[4])
	{
		if (data.size() > 0)
		{
			auto max_it_y = std::max_element(data.begin(), data.end(), [&](auto a, auto b) { return std::get<2>(a) < std::get<2>(b); });
			auto min_it_y = std::max_element(data.begin(), data.end(), [&](auto a, auto b) { return std::get<2>(a) > std::get<2>(b); });
			auto max_it_x = std::max_element(data.begin(), data.end(), [&](auto a, auto b) { return std::get<1>(a) < std::get<1>(b); });
			auto min_it_x = std::max_element(data.begin(), data.end(), [&](auto a, auto b) { return std::get<1>(a) > std::get<1>(b); });

			auto max_y = std::get<2>(*max_it_y);
			auto min_y = std::get<2>(*min_it_y);
			auto max_x = std::get<1>(*max_it_x);
			auto min_x = std::get<1>(*min_it_x);

			ranges[0] = static_cast<T>(max_x);
			ranges[1] = static_cast<T>(min_x);
			ranges[2] = max_y;
			ranges[3] = min_y;
		}
	}

	inline void fill_ranges(axis_3d_type const* data_x, len_type len, axis_coord_t (&ranges)[6])
	{
		// initialize infinite of opposite sign on ranges 
		// so that the checks with data_x can work.

		ranges[0] = +INFINITY;
		ranges[2] = +INFINITY;
		ranges[4] = +INFINITY;

		ranges[1] = -INFINITY;
		ranges[3] = -INFINITY;
		ranges[5] = -INFINITY;

		for (iter_type i = 0; i < len; ++i)
		{
			if (data_x[i][0] < ranges[0])
			{
				ranges[0] = data_x[i][0];
			}
			if (data_x[i][0] > ranges[1])
			{
				ranges[1] = data_x[i][0];
			}
			if (data_x[i][1] < ranges[2])
			{
				ranges[2] = data_x[i][1];
			}
			if (data_x[i][1] > ranges[3])
			{
				ranges[3] = data_x[i][1];
			}
			if (data_x[i][2] < ranges[4])
			{
				ranges[4] = data_x[i][2];
			}
			if (data_x[i][2] > ranges[5])
			{
				ranges[5] = data_x[i][2];
			}
		}
	}


	inline void fill_ranges(axis_2d_type const* data_x, len_type len, axis_coord_t(&ranges)[4])
	{
		// initialize infinite of opposite sign on ranges 
		// so that the checks with data_x can work.

		ranges[0] = +INFINITY;
		ranges[2] = +INFINITY;

		ranges[1] = -INFINITY;
		ranges[3] = -INFINITY;

		for (iter_type i = 0; i < len; ++i)
		{
			if (data_x[i][0] < ranges[0])
			{
				ranges[0] = data_x[i][0];
			}
			if (data_x[i][0] > ranges[1])
			{
				ranges[1] = data_x[i][0];
			}
			if (data_x[i][1] < ranges[2])
			{
				ranges[2] = data_x[i][1];
			}
			if (data_x[i][1] > ranges[3])
			{
				ranges[3] = data_x[i][1];
			}
		}
	}



	inline void fill_ranges(axis_1d_type const* data_x, len_type len, axis_coord_t(&ranges)[2])
	{
		// initialize infinite of opposite sign on ranges 
		// so that the checks with data_x can work.

		ranges[0] = +INFINITY;
		ranges[1] = -INFINITY;

		for (iter_type i = 0; i < len; ++i)
		{
			if (data_x[i] < ranges[0])
			{
				ranges[0] = data_x[i];
			}
			if (data_x[i] > ranges[1])
			{
				ranges[1] = data_x[i];
			}
		}
	}




	/*
	 *
	 *
	 * Implementations of tuple functions.
	 */


	template<typename... As, typename... Bs>
	auto unfurl_tuple(std::tuple<std::pair<As, Bs>...> const& ts)
	{
		return unfurl_tuple(ts, std::make_index_sequence<std::tuple_size<std::tuple<std::pair<As, Bs>...>>::value>{});
	}

	template<size_t I, typename... Ts>
	auto get_tuple_ge(std::tuple<Ts...> const& ts)
	{
		return get_tuple_ge<I>(ts, std::make_index_sequence<sizeof...(Ts) - I>{});
	}

	template<size_t I, typename... Ts>
	auto get_tuple_lt(std::tuple<Ts...> const& ts)
	{
		return get_tuple_lt(ts, std::make_index_sequence<I>{});
	}

	template<size_t I, size_t... Ns>
	auto get_seq_ge(std::index_sequence<Ns...> const& seq)
	{
		if constexpr (sizeof...(Ns) < I)
		{
			return std::index_sequence<>{};
		}
		else
		{
			return get_seq_ge<I>(std::tuple<std::index_sequence<Ns>...>{}, std::make_index_sequence<sizeof...(Ns) - I>{});
		}
	}

	template<size_t I, size_t... Ns>
	auto get_seq_lt(std::index_sequence<Ns...> const& seq)
	{
		return get_seq_lt(std::tuple<std::index_sequence<Ns>...>{}, std::make_index_sequence<I>{});
	}

	template<typename... Ts>
	auto make_tuple(Ts... ts)
	{
		return std::tuple_cat(std::tuple<Ts>(ts)...);
	}

	template<typename empty_type, typename... match_A_types, typename... As, typename... match_B_types, typename... Bs>
	auto unzip_pot(std::tuple<std::pair<match_A_types, As>...> const& a, std::tuple<std::pair<match_B_types, Bs>...> const& b)
	{
		return unzip_pot_carry_apply<empty_type>(std::make_tuple(), a, b);
	}


	//! Print the timestamp into the given buffer string.
	/*! 
	 * Print the timestamp into the given buffer string.
	 */
	void write_ts_str(char* buffer);

#ifdef FILESYSTEM_HEADER_AVAILABLE
	void make_directory(std::filesystem::path dir, int err_no = ERR_CODE_FILE_OPEN);
#endif
	
	//! Create a directory given by the string argument
	/*!
	 * Create a directory given by the string argument. If any parent path
	 * components of the directory do not exist, these will also be generated.
	 */
	void make_directory(const char* dir, int err_no = ERR_CODE_FILE_OPEN);

	//! Create a directory given by the string argument
	/*!
	 * Create a directory given by the string argument. If any parent path
	 * components of the directory do not exist, these will also be generated.
	 */
	void make_directory_for_file(const char* dir, int err_no = ERR_CODE_FILE_OPEN);

	//! Return the directory component of the given path. 
	/*!
	 * Standardized method to get the directory component of the
	 * given path.
	 */
	void get_parent_directory(char* path, char* &basepath);

#ifdef FILESYSTEM_HEADER_AVAILABLE
	std::filesystem::path get_parent_directory(std::filesystem::path dir);
#endif


	//! Based on the condition, return the value.
	/*!
	 * Used for initializing constant variables by evaluating the constant
	 * expression condition and returning the first given parameter
	 * if true.
	 * 
	 * \param v1 The value which is returned if the condition is true.
	 * \param v2 The value return if the condition is false.
	 * 
	 * \tparam condition The condition that is evaluated to
	 * determine which value is returned.
	 */
	template<bool condition, typename C0, typename C1>
	constexpr decltype(auto) conditional_return_value(C0&& v1, C1&& v2)
	{
		if constexpr (condition) return v1;
		else return v2;
	}

}

template<typename S0, typename ...Ss>
struct symphas::internal::index_of_type_match<0, S0, Ss...>
{
	template<typename Type, size_t Sn>
	static const int value = Sn - sizeof...(Ss) - 1;
};

template<size_t I, typename S0, typename... Ss>
struct symphas::internal::index_of_type_match<I, S0, Ss...>
{
	template<typename Type, size_t Sn>
	static const int value = symphas::lib::index_of_type<Type, I - 1, Sn, Ss...>();
};

template<size_t I, typename... Ss>
struct symphas::internal::index_of_type_match
{
	template<typename Type, size_t Sn>
	static const int value = -1;
};



// **************************************************************************************


