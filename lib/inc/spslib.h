
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
#include "timer.h"

#ifdef FILESYSTEM_HEADER_AVAILABLE
#include <filesystem>
#endif




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


}


/*! \addtogroup grid
 * @{
 */

//! Defines functions providing information about a grid.
namespace grid
{
	template<size_t D>
	void set_dimensions(const len_type* dimensions, len_type* dims)
	{
		if (dimensions == nullptr)
		{
			std::fill(dims, dims + D, 0);
		}
		else
		{
			for (iter_type i = 0; i < D; ++i)
			{
				dims[i] = dimensions[i];
			}
		}
	}

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
		return ((dimensions != nullptr) ? std::max(0, dimensions[0]) : 0);
	}

	template<size_t D>
	len_type length(len_type const* dimensions)
	{
		return ((dimensions != nullptr) ? std::max(0, dimensions[D - 1]) : 0) * length<D - 1>(dimensions);
	}

	inline len_type length(len_type const* dimensions, size_t dimension)
	{
		len_type len = 1;
		for (iter_type i = 0; i < dimension; ++i)
		{
			len *= std::max(0, dimensions[i]);
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
#define MATH_FUNCTION_OVERLOADS_IMPL(NAME, RETURN) \
inline scalar_t _ ## NAME(int vv) { auto v = static_cast<scalar_t>(vv); RETURN; } \
inline scalar_t _ ## NAME(scalar_t v) { RETURN; } \
inline complex_t _ ## NAME(complex_t v) { RETURN; } \
template<typename T> \
auto _ ## NAME(T v) { return _ ## NAME(static_cast<complex_t>(v)); } \
template<typename T> \
T NAME(T v) { return _ ## NAME(v); }

#define MATH_FUNCTION_OVERLOADS(NAME, FUNC) MATH_FUNCTION_OVERLOADS_IMPL(NAME, return FUNC(v))

namespace symphas::math
{
	
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

	//! Returns conjugate complex number for a generic type.
	/*!
	 * Since additional types are introduced in the SymPhas library, an
	 * extension of the standard library `conj` function is implemented to
	 * apply to additional types.
	 *
	 * This implementation simply forwards its argument to the standard library
	 * `conj` function.
	 */
	template<typename T, size_t D>
	auto conj(any_vector_t<T, D> const& v)
	{
		using std::conj;

		any_vector_t<T, D> vc;
		for (iter_type i = 0; i < D; ++i)
		{
			vc[i] = conj(v[i]);
		}
		return vc;
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

	template<typename T, size_t D>
	T abs(const T(&value)[D])
	{
		T result{};
		for (iter_type i = 0; i < D; ++i)
		{
			result += value[i] * value[i];
		}

		using std::sqrt;
		return sqrt(result);
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
	MATH_FUNCTION_OVERLOADS(log, std::log);

	MATH_FUNCTION_OVERLOADS_IMPL(sec, return 1. / std::cos(v));
	MATH_FUNCTION_OVERLOADS_IMPL(csc, return 1. / std::sin(v));
	MATH_FUNCTION_OVERLOADS_IMPL(cot, return 1. / std::tan(v));



	template<size_t N, size_t E>
	constexpr size_t fixed_pow = N * fixed_pow<N, E - 1>;
	template<size_t N>
	constexpr size_t fixed_pow<N, 0> = 1;
	template<size_t N>
	constexpr size_t fixed_pow<N, 1> = N;



	template<typename T, typename E>
	auto pow(T&& base, E&& exponent)
	{
		using std::pow;
		return pow(std::forward<T>(base), std::forward<E>(exponent));
	}


	template<typename T, size_t... Is>
	auto pow(T const& base, std::index_sequence<Is...>)
	{
		auto f = [&] (size_t) { return base; };
		return (f(Is) * ...);
	}

	template<size_t N, typename T>
	auto pow(T const& base)
	{
		if constexpr (N <= 6)
		{
			return pow(base, std::make_index_sequence<N>{});
		}
		else
		{
			//using symphas::math::pow;
			return pow(base, N);
		}
	}

	template<size_t O, typename T, size_t D>
	auto pow(any_vector_t<T, D> const& v)
	{
		using std::pow;
		using symphas::math::pow;

		if constexpr (O == 1)
		{
			return v;
		}
		else if constexpr (O == 2)
		{
			return dot(v, v);
		}
		else if constexpr (O % 2 == 0)
		{
			return pow<O / 2>(dot(v, v));
		}
		else
		{
			return pow<O / 2>(dot(v, v)) * v;
		}
	}

	template<typename T, size_t D>
	auto pow(any_vector_t<T, D> const& v, double exponent)
	{
		throw;
	}

	//! Cross product function for vectors.
	/*!
	 * Compound assignment operator to calculate the conventional cross product and assign it
	 * to the calling instance.
	 *
	 * Assignment operator to compute the cross product, which is only applicable to the
	 * 3-dimensional VectorValue template specialization. Operator is called by the left hand
	 * instance, taking data from the right hand VectorValue to compute the cross product and
	 * assigning the result to the left hand side.
	 *
	 * \param rhs The VectorValue instance where data is taken from to compute the cross product
	 * the data from the left hand instance.
	 */
	template<typename T>
	any_vector_t<T, 3> cross(any_vector_t<T, 3> const& lhs, any_vector_t<T, 3> const& rhs)
	{
		any_vector_t<T, 3> result;

		result[0] = (lhs[1] * rhs[2] - lhs[2] * rhs[1]);
		result[1] = (lhs[2] * rhs[0] - lhs[0] * rhs[2]);
		result[2] = (lhs[0] * rhs[1] - lhs[1] * rhs[0]);

		return result;
	}


	//template<size_t D>
	//auto dot(any_vector_t<scalar_t, D> const& lhs, any_vector_t<complex_t, D> const& rhs)
	//{
	//	complex_t sum = 0;
	//	for (iter_type i = 0; i < D; ++i)
	//	{
	//		sum += lhs[i] * rhs[i];
	//	}
	//	return sum;
	//}

	//template<size_t D>
	//auto dot(any_vector_t<complex_t, D> const& lhs, any_vector_t<scalar_t, D> const& rhs)
	//{
	//	complex_t sum = 0;
	//	for (iter_type i = 0; i < D; ++i)
	//	{
	//		sum += lhs[i] * rhs[i];
	//	}
	//	return sum;
	//}


	template<size_t... Os>
	constexpr auto sum()
	{
		return (Os + ...);
	}

	template<size_t... Os>
	constexpr auto product()
	{
		return (Os * ...);
	}
}


namespace symphas::internal
{
	template<size_t, typename T>
	decltype(auto) repeat_value(T&& value)
	{
		return std::forward<T>(value);
	}

	template<typename, typename T>
	decltype(auto) repeat_value(T&& value)
	{
		return std::forward<T>(value);
	}


	template<typename T>
	struct check_is_simple_data
	{
	protected:

		constexpr static std::true_type _get_value(complex_t)
		{
			return {};
		}

		constexpr static std::true_type _get_value(scalar_t)
		{
			return {};
		}

		constexpr static std::true_type _get_value(int)
		{
			return {};
		}

		template<typename T0>
		constexpr static std::false_type _get_value(T0)
		{
			return {};
		}

		constexpr static auto get_value(T ptr)
		{
			return _get_value(ptr);
		}


	public:

		static const bool value = std::invoke_result_t<decltype(&check_is_simple_data<T>::get_value), T>::value;
	};

	template<typename T>
	struct check_is_simple_data<T&> : check_is_simple_data<T> {};
	template<typename T>
	struct check_is_simple_data<T const> : check_is_simple_data<T> {};
}


template<typename T>
constexpr bool is_simple_data = symphas::internal::check_is_simple_data<T>::value;

template<typename T>
constexpr bool is_non_vector = true;
template<typename T, size_t D>
constexpr bool is_non_vector<any_vector_t<T, D>> = false;

template<typename T, size_t N, size_t M>
using any_matrix_t = any_vector_t<any_vector_t<T, M>, N>;
template<typename T, size_t N>
using any_row_vector_t = any_matrix_t<T, 1, N>;


template<typename T, size_t D, typename std::enable_if_t<is_non_vector<T>, int> = 0>
auto operator*(any_vector_t<T, D> const& lhs, any_vector_t<T, D> const& rhs)
{
	using namespace std;
	using namespace symphas::math;
	return dot(lhs, rhs);
}

template<typename T, typename S, size_t D, typename std::enable_if_t<is_non_vector<T> && is_non_vector<S>, int> = 0>
auto operator*(any_vector_t<T, D> const& lhs, any_vector_t<S, D> const& rhs)
{
	using namespace std;
	using namespace symphas::math;
	return dot(lhs, rhs);
}

template<typename T, typename S, size_t D, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator*(any_row_vector_t<T, D> const& lhs, any_vector_t<S, D> const& rhs)
{
	any_vector_t<T, D> lhs0;
	for (iter_type i = 0; i < D; ++i)
	{
		lhs0[i] = lhs[0][i];
	}

	using namespace std;
	using namespace symphas::math;
	return dot(lhs0, rhs);
}

template<typename T, typename S, size_t N, typename std::enable_if_t<(is_non_vector<T> && is_simple_data<S>), int> = 0>
auto operator*(any_vector_t<T, N> const& lhs, S const& rhs)
{
	any_vector_t<mul_result_t<T, S>, N> result;
	for (int i = 0; i < N; ++i)
	{
		result[i] = lhs[i] * rhs;
	}
	return result;
}

template<typename T, typename S, size_t N, typename std::enable_if_t<(is_simple_data<T> && is_non_vector<S>), int> = 0>
auto operator*(T const& lhs, any_vector_t<S, N> const& rhs)
{
	any_vector_t<mul_result_t<T, S>, N> result;
	for (int i = 0; i < N; ++i)
	{
		result[i] = lhs * rhs[i];
	}
	return result;
}

template<typename T, typename S, size_t N, typename std::enable_if_t<(is_non_vector<T> && is_simple_data<S>), int> = 0>
auto operator*(any_row_vector_t<T, N> const& lhs, S const& rhs)
{
	any_row_vector_t<mul_result_t<T, S>, N> result;
	for (int i = 0; i < N; ++i)
	{
		result[0][i] = lhs[0][i] * rhs;
	}
	return result;
}

template<typename T, typename S, size_t N, typename std::enable_if_t<(is_simple_data<T> && is_non_vector<S>), int> = 0>
auto operator*(T const& lhs, any_row_vector_t<S, N> const& rhs)
{
	any_row_vector_t<mul_result_t<T, S>, N> result;
	for (int i = 0; i < N; ++i)
	{
		result[0][i] = lhs * rhs[0][i];
	}
	return result;
}

template<typename T, typename S, size_t N, size_t M, typename std::enable_if_t<(is_non_vector<T> && is_simple_data<S>), int> = 0>
auto operator*(any_matrix_t<T, N, M> const& lhs, S const& rhs)
{
	any_matrix_t<mul_result_t<T, S>, N, M> result;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			result[i][i] = lhs[i][j] * rhs;
		}
	}
	return result;
}

template<typename T, typename S, size_t N, size_t M, typename std::enable_if_t<(is_simple_data<T> && is_non_vector<S>), int> = 0>
auto operator*(T const& lhs, any_matrix_t<S, N, M> const& rhs)
{
	any_matrix_t<mul_result_t<T, S>, N, M> result;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			result[i][i] = lhs * rhs[i][j];
		}
	}
	return result;
}

template<typename T, typename S, size_t L, size_t M, size_t N, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator*(any_matrix_t<T, L, M> const& lhs, any_matrix_t<S, M, N> const& rhs)
{
	any_matrix_t<mul_result_t<T, S>, L, N> result;
	for (int i = 0; i < L; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			mul_result_t<T, S> sum{};
			for (int k = 0; k < M; ++k)
			{
				sum += lhs[i][k] * rhs[k][j];
			}
			result[i][j] = sum;
		}
	}
	return result;
}

template<typename T, typename S, size_t L, size_t M, typename std::enable_if_t<(is_non_vector<T>&& is_non_vector<S>), int> = 0>
auto operator*(any_matrix_t<T, L, M> const& lhs, any_vector_t<S, M> const& rhs)
{
	any_vector_t<mul_result_t<T, S>, M> result;
	for (int i = 0; i < L; ++i)
	{
		mul_result_t<T, S> sum{};
		for (int k = 0; k < M; ++k)
		{
			sum += lhs[i][k] * rhs[k];
		}
		result[i] = sum;
	}
	return result;
}

template<typename T, typename S, size_t L, size_t M, typename std::enable_if_t<(is_non_vector<T>&& is_non_vector<S>), int> = 0>
auto operator*(any_row_vector_t<T, L> const& lhs, any_matrix_t<S, L, M> const& rhs)
{
	any_row_vector_t<mul_result_t<T, S>, M> result;
	for (int i = 0; i < M; ++i)
	{
		mul_result_t<T, S> sum{};
		for (int k = 0; k < L; ++k)
		{
			sum += lhs[k] * rhs[k][i];
		}
		result[i] = sum;
	}
	return result;
}

template<typename T, typename S, size_t N, typename std::enable_if_t<(is_non_vector<T> && is_simple_data<S>), int> = 0>
auto operator/(any_vector_t<T, N> const& lhs, S const& rhs)
{
	any_vector_t<div_result_t<T, S>, N> result;
	for (int i = 0; i < N; ++i)
	{
		result[i] = lhs[i] / rhs;
	}
	return result;
}

template<typename T, typename S, size_t N, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator*(any_vector_t<T, 1> const& lhs, any_vector_t<S, N> const& rhs)
{
	return lhs.v[0] * rhs;
}

template<typename T, typename S, size_t N, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator*(any_vector_t<T, N> const& lhs, any_vector_t<S, 1> const& rhs)
{
	return lhs * rhs.v[0];
}

template<typename T, typename S, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator*(any_vector_t<T, 1> const& lhs, any_vector_t<S, 1> const& rhs)
{
	return lhs.v[0] * rhs.v[0];
}

template<typename T, typename std::enable_if_t<(is_non_vector<T>), int> = 0>
auto operator*(any_vector_t<T, 1> const& lhs, any_vector_t<T, 1> const& rhs)
{
	return lhs.v[0] * rhs.v[0];
}

template<typename T, typename S, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator+(any_vector_t<T, 1> const& lhs, S const& rhs)
{
	any_vector_t<add_result_t<T, S>, 1> result;
	result[0] = lhs[0] + rhs;
	return result;
}

template<typename T, typename S, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator+(T const& lhs, any_vector_t<S, 1> const& rhs)
{
	any_vector_t<add_result_t<T, S>, 1> result;
	result[0] = lhs + rhs[0];
	return result;
}

template<typename T, typename S, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator+(any_row_vector_t<T, 1> const& lhs, S const& rhs)
{
	any_row_vector_t<add_result_t<T, S>, 1> result;
	result[0][0] = lhs[0][0] + rhs;
	return result;
}

template<typename T, typename S, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator+(T const& lhs, any_row_vector_t<S, 1> const& rhs)
{
	any_row_vector_t<add_result_t<T, S>, 1> result;
	result[0][0] = lhs + rhs[0][0];
	return result;
}

template<typename T, typename S, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator-(any_vector_t<T, 1> const& lhs, S const& rhs)
{
	any_vector_t<add_result_t<T, S>, 1> result;
	result[0] = lhs[0] - rhs;
	return result;
}

template<typename T, typename S, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator-(T const& lhs, any_vector_t<S, 1> const& rhs)
{
	any_vector_t<add_result_t<T, S>, 1> result;
	result[0] = lhs - rhs[0];
	return result;
}

template<typename T, typename S, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator-(any_row_vector_t<T, 1> const& lhs, S const& rhs)
{
	any_row_vector_t<add_result_t<T, S>, 1> result;
	result[0][0] = lhs[0][0] - rhs;
	return result;
}

template<typename T, typename S, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator-(T const& lhs, any_row_vector_t<S, 1> const& rhs)
{
	any_row_vector_t<add_result_t<T, S>, 1> result;
	result[0][0] = lhs - rhs[0][0];
	return result;
}


template<typename T, typename S, size_t N, size_t M, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator+(any_matrix_t<T, N, M> const& lhs, any_matrix_t<S, N, M> const& rhs)
{
	any_matrix_t<add_result_t<T, S>, N, M> result;
	for (int i = 0; i < N; ++i)
	{
		result[0][i] = lhs[0][i] + rhs[0][i];
	}
	return result;
}

template<typename T, typename S, size_t N, typename std::enable_if_t<(is_non_vector<T> && is_non_vector<S>), int> = 0>
auto operator+(any_vector_t<T, N> const& lhs, any_vector_t<S, N> const& rhs)
{
	any_vector_t<add_result_t<T, S>, N> result;
	for (int i = 0; i < N; ++i)
	{
		result[i] = lhs[i] + rhs[i];
	}
	return result;
}

template<typename T, typename S, size_t N, size_t M, typename std::enable_if_t<(is_non_vector<T>&& is_non_vector<S>), int> = 0>
auto operator-(any_matrix_t<T, N, M> const& lhs, any_matrix_t<S, N, M> const& rhs)
{
	any_matrix_t<sub_result_t<T, S>, N, M> result;
	for (int i = 0; i < N; ++i)
	{
		result[0][i] = lhs[0][i] - rhs[0][i];
	}
	return result;
}

template<typename T, typename S, size_t N, typename std::enable_if_t<(is_non_vector<T>&& is_non_vector<S>), int> = 0>
auto operator-(any_vector_t<T, N> const& lhs, any_vector_t<S, N> const& rhs)
{
	any_vector_t<sub_result_t<T, S>, N> result;
	for (int i = 0; i < N; ++i)
	{
		result[i] = lhs[i] - rhs[i];
	}
	return result;
}

inline axis_nd_t<2> operator+(axis_nd_t<2> const& a, axis_nd_t<2> const& b)
{
	return { a[0] + b[0], a[1] + b[1] };
}

inline axis_nd_t<2> operator-(axis_nd_t<2> const& a, axis_nd_t<2> const& b)
{
	return { a[0] - b[0], a[1] - b[1] };
}

inline axis_nd_t<3> operator+(axis_nd_t<3> const& a, axis_nd_t<3> const& b)
{
	return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

inline axis_nd_t<3> operator-(axis_nd_t<3> const& a, axis_nd_t<3> const& b)
{
	return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
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


	template<typename... Ts>
	struct types_list {};

	template<typename... Ts>
	auto to_types_list(std::tuple<Ts...> const&)
	{
		return types_list<Ts...>{};
	}

	template<typename T0, typename T1, typename... Ts>
	auto to_types_list(T0 const&, T1 const&, Ts const&...)
	{
		return types_list<T0, T1, Ts...>{};
	}

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
	template<size_t I, typename T, T... Ns>
	auto get_seq_ge(std::integer_sequence<T, Ns...> const& seq);

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
	template<size_t I, typename T, T... Ns>
	auto get_seq_lt(std::integer_sequence<T, Ns...> const& seq);


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


	//! Joining a single sequence simply returns it.
	template<typename T, T... Ys>
	constexpr auto seq_join(std::integer_sequence<T, Ys...>)
	{
		return std::integer_sequence<T, Ys...>{};
	}

	//! Joins two index sequences.
	template<typename T, T... Ys, T... Qs>
	constexpr auto seq_join(std::integer_sequence<T, Ys...>, std::integer_sequence<T, Qs...>)
	{
		return std::integer_sequence<T, Ys..., Qs...>{};
	}

	//! Joins two index sequences.
	template<typename T, T... Ys, T... Qs, typename... Seqs>
	constexpr auto seq_join(std::integer_sequence<T, Ys...>, std::integer_sequence<T, Qs...>, Seqs... seqs)
	{
		return seq_join(std::integer_sequence<T, Ys..., Qs...>{}, seqs...);
	}

	template<typename T, T... Ys>
	constexpr auto seq_neg(std::integer_sequence<T, Ys...>)
	{
		return std::integer_sequence<int, -Ys...>{};
	}

	template<size_t... Ys>
	constexpr auto seq_neg(std::index_sequence<Ys...>)
	{
		return std::integer_sequence<int, -int(Ys)...>{};
	}


	//! Adding a single sequence simply returns it.
	template<typename T, T... Ys>
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
	template<typename T0, T0... Ys, typename T1, T1... Qs, typename T = std::conditional_t<std::is_same<T0, T1>::value, T0, int>>
	constexpr auto seq_add(std::integer_sequence<T0, Ys...>, std::integer_sequence<T1, Qs...>)
	{
		if constexpr (sizeof...(Ys) == sizeof...(Qs))
		{
			return std::integer_sequence<T, (T)(Qs + Ys)...>{};
		}
		else if constexpr (sizeof...(Ys) > sizeof...(Qs))
		{
			return seq_join(
				seq_add(
					symphas::lib::get_seq_lt<sizeof...(Qs)>(std::integer_sequence<T0, Ys...>{}),
					std::integer_sequence<T1, Qs...>{}),
				symphas::lib::get_seq_ge<sizeof...(Qs)>(std::integer_sequence<T0, Ys...>{}));
		}
		else
		{
			return seq_join(
				seq_add(
					std::integer_sequence<T0, Ys...>{},
					symphas::lib::get_seq_lt<sizeof...(Ys)>(std::integer_sequence<T1, Qs...>{})),
				symphas::lib::get_seq_ge<sizeof...(Ys)>(std::integer_sequence<T1, Qs...>{}));
		}
	}

	//! Adds multiple index sequences.
	/*!
	 * The values are added pointwise, between two index sequences. When
	 * sequences are not of equal size, the shortest one is considered
	 * to have 0s in the remaining entries, so that a sequence equal in length
	 * to the longest sequence is always returned.
	 */
	template<typename T0, T0... Ys, typename T1, T1... Qs, typename... Seqs>
	constexpr auto seq_add(std::integer_sequence<T0, Ys...>, std::integer_sequence<T1, Qs...>, Seqs... seqs)
	{
		return seq_add(seq_add(std::integer_sequence<T0, Ys...>{}, std::integer_sequence<T1, Qs...>{}), seqs...);
	}


	//! Subtracting a single sequence simply returns it.
	template<typename T, T... Ys>
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
	template<typename T0, T0... Ys, typename T1, T1... Qs, typename... Seqs>
	constexpr auto seq_sub(std::integer_sequence<T0, Ys...>, std::integer_sequence<T1, Qs...>, Seqs... seqs)
	{
		return seq_sub(seq_add(std::integer_sequence<int, int(Ys)...>{}, seq_neg(std::integer_sequence<T1, Qs...>{})), seqs...);
	}

	template<typename T, T... Ys>
	constexpr auto filter_seq(std::integer_sequence<T, Ys...>)
	{
		return std::integer_sequence<T, Ys...>{};
	}

	namespace
	{
		template<typename T, T... Ys1, T Y0, T... Ys2, T Q0, T... Qs>
		constexpr auto _filter_seq(std::integer_sequence<T, Ys1...>, std::integer_sequence<T, Y0, Ys2...>, std::integer_sequence<T, Q0>, std::integer_sequence<T, Qs...>);

		template<typename T, T... Ys1, T Q0, T... Qs>
		constexpr auto _filter_seq(std::integer_sequence<T>, std::integer_sequence<T, Ys1...>, std::integer_sequence<T>, std::integer_sequence<T, Q0, Qs...>)
		{
			return _filter_seq(std::integer_sequence<T>{}, std::integer_sequence<T, Ys1...>{}, std::integer_sequence<T, Q0>{}, std::integer_sequence<T, Qs...>{});
		}

		template<typename T, T... Ys1>
		constexpr auto _filter_seq(std::integer_sequence<T>, std::integer_sequence<T, Ys1...>, std::integer_sequence<T>, std::integer_sequence<T>)
		{
			return std::integer_sequence<T, Ys1...>{};
		}

		template<typename T, T Y0, T... Ys1>
		constexpr auto _filter_seq(std::integer_sequence<T, Y0, Ys1...>, std::integer_sequence<T>, std::integer_sequence<T>, std::integer_sequence<T>)
		{
			return std::integer_sequence<T, Y0, Ys1...>{};
		}

		template<typename T, T... Ys1, T Q0>
		constexpr auto _filter_seq(std::integer_sequence<T, Ys1...>, std::integer_sequence<T>, std::integer_sequence<T, Q0>, std::integer_sequence<T>)
		{
			return std::integer_sequence<T, Ys1...>{};
		}

		template<typename T, T... Ys1, T Q0, T Q1, T... Qs>
		constexpr auto _filter_seq(std::integer_sequence<T, Ys1...>, std::integer_sequence<T>, std::integer_sequence<T, Q0>, std::integer_sequence<T, Q1, Qs...>)
		{
			return _filter_seq(std::integer_sequence<T>{}, std::integer_sequence<T, Ys1...>{}, std::integer_sequence<T, Q1>{}, std::integer_sequence<T, Qs...>{});
		}

		template<typename T, T... Ys1, T Y0, T... Ys2, T Q0, T... Qs>
		constexpr auto _filter_seq(std::integer_sequence<T, Ys1...>, std::integer_sequence<T, Y0, Ys2...>, std::integer_sequence<T, Q0>, std::integer_sequence<T, Qs...>)
		{
			if constexpr (Y0 == Q0)
			{
				return _filter_seq(std::integer_sequence<T, Ys1...>{}, std::integer_sequence<T, Ys2...>{}, std::integer_sequence<T, Q0>{}, std::integer_sequence<T, Qs...>{});
			}
			else
			{
				return _filter_seq(std::integer_sequence<T, Ys1..., Y0>{}, std::integer_sequence<T, Ys2...>{}, std::integer_sequence<T, Q0>{}, std::integer_sequence<T, Qs...>{});
			}
		}
	}

	//! Filters from the first sequence, all values in the subsequent sequences.
	/*!
	 * All values that appear in the second and subsequence sequences are removed from the first sequence. If
	 * there are any values in the second sequence that do not appear in the first, there is no
	 * effect. If the intersection of elements shared between the sets is empty, then it
	 * will return the first sequence.
	 * 
	 * Note: Each element is filtered only once. This will not filter all repeating elements.
	 */
	template<typename T, T... Ys, T... Qs, typename... Seqs>
	constexpr auto filter_seq(std::integer_sequence<T, Ys...>, std::integer_sequence<T, Qs...>, Seqs... seqs)
	{
		return filter_seq(_filter_seq(std::integer_sequence<T>{}, std::integer_sequence<T, Ys...>{}, std::integer_sequence<T>{}, std::integer_sequence<T, Qs...>{}), seqs...);
	}


	//! Makes a new sequence with values shared between all the sequences.
	/*!
	 * All values that are shared by all provided sequences will be put into a new
	 * sequence.
	 */
	template<typename T, T... Ys>
	constexpr auto intersect_seq(std::integer_sequence<T, Ys...>)
	{
		return std::integer_sequence<T, Ys...>{};
	}

	//! Makes a new sequence with values shared between all the sequences.
	/*!
	 * All values that are shared by all provided sequences will be put into a new
	 * sequence.
	 */
	template<typename T, T... Ys, T... Qs, typename... Seqs>
	constexpr auto intersect_seq(std::integer_sequence<T, Ys...>, std::integer_sequence<T, Qs...>, Seqs... seqs);



	//! The index sequence result type of joining multiple sequences.
	template<typename Seq, typename... Seqs>
	struct seq_join_result
	{
		using type = decltype(seq_join(std::declval<Seq>(), std::declval<Seqs>()...));
	};

	//! The index sequence result type of joining multiple sequences.
	template<typename Seq, typename... Seqs>
	struct seq_join_result<types_list<Seq, Seqs...>>
	{
		using type = decltype(seq_join(std::declval<Seq>(), std::declval<Seqs>()...));
	};

	//! The index sequence result type of joining multiple sequences.
	template<>
	struct seq_join_result<types_list<>>
	{
		using type = std::index_sequence<>;
	};

	//! The index sequence result type of adding multiple sequences.
	template<typename Seq, typename... Seqs>
	struct seq_add_result
	{
		using type = decltype(seq_add(std::declval<Seq>(), std::declval<Seqs>()...));
	};

	//! The index sequence result type of adding a value to each element of the sequence.
	template<size_t N, typename Seq>
	struct seq_offset_result;

	//! The index sequence result type of adding a value to each element of the sequence.
	template<size_t N, size_t... Is>
	struct seq_offset_result<N, std::index_sequence<Is...>>
	{
		using type = std::index_sequence<(N + Is)...>;
	};

	//! The index sequence result type of adding multiple sequences.
	template<typename Seq>
	struct seq_from_to_result;

	//! The index sequence result type of adding multiple sequences.
	template<size_t I0, size_t J0>
	struct seq_from_to_result<std::index_sequence<I0, J0>>
	{
		using type = typename seq_offset_result<I0, std::make_index_sequence<J0 - I0>>::type;
	};


	//! The index sequence result type of adding multiple sequences.
	template<typename Seq, typename... Seqs>
	struct seq_sub_result
	{
		using type = decltype(seq_sub(std::declval<Seq>(), std::declval<Seqs>()...));
	};

	//! The index sequence result type of filtering multiple sequences.
	/*!
	 * The first sequence that is given has the values removed from it that are found in the
	 * other sequences.
	 */
	template<typename... Seqs>
	struct filter_seq_result;

	//! The index sequence result type of filtering multiple sequences.
	template<typename T, T... Is>
	struct filter_seq_result<std::integer_sequence<T, Is...>>
	{
		using type = std::integer_sequence<T, Is...>;
	};

	//! The index sequence result type of filtering multiple sequences.
	template<typename T>
	struct filter_seq_result<std::integer_sequence<T>, std::integer_sequence<T>>
	{
		using type = std::integer_sequence<T>;
	};

	template<typename T, T I0, T... Fs>
	struct test_ne
	{
		static const bool value = ((I0 != Fs) && ...);
	};

	//! The index sequence result type of filtering multiple sequences.
	template<typename T, T I0, T... Is, T... Fs>
	struct filter_seq_result<std::integer_sequence<T, I0, Is...>, std::integer_sequence<T, Fs...>>
	{
		using type = typename seq_join_result<
			std::conditional_t<test_ne<T, I0, Fs...>::value, std::integer_sequence<T, I0>, std::integer_sequence<T>>,
			std::conditional_t<test_ne<T, Is, Fs...>::value, std::integer_sequence<T, Is>, std::integer_sequence<T>>...>::type;
	};

	//! The index sequence result type of filtering multiple sequences.
	template<typename T, T I0, T... Is, bool F0, bool... Fs>
	struct filter_seq_result<std::integer_sequence<T, I0, Is...>, std::integer_sequence<bool, F0, Fs...>>
	{
		using type = typename seq_join_result<
			std::conditional_t<F0, std::integer_sequence<T, I0>, std::integer_sequence<T>>,
			std::conditional_t<Fs, std::integer_sequence<T, Is>, std::integer_sequence<T>>...>::type;
	};

	//! The index sequence result type of filtering multiple sequences.
	template<typename Seq0, typename Seq1, typename Seq2, typename... Seqs>
	struct filter_seq_result<Seq0, Seq1, Seq2, Seqs...>
	{
		using type = typename filter_seq_result<typename filter_seq_result<Seq0, Seq1>::type, Seq2, Seqs...>::type;
	};

	//! The index sequence result type of intersecting multiple sequences.
	template<typename... Seqs>
	struct intersect_seq_result
	{
		using type = decltype(intersect_seq(std::declval<Seqs>()...));
	};

	//! Alias for the join result of multiple sequences.
	template<typename Seq, typename... Seqs>
	using seq_join_t = typename seq_join_result<Seq, Seqs...>::type;

	//! Alias for the addition result of multiple sequences.
	template<typename... Seqs>
	using seq_add_t = typename seq_add_result<Seqs...>::type;

	//! The index sequence result type of adding a value to each element of the sequence.
	template<size_t N, typename Seq>
	using seq_offset_t = typename seq_offset_result<N, Seq>::type;

	template<typename Seq>
	using seq_from_to_t = typename seq_from_to_result<Seq>::type;

	//! Alias for the subtraction result of multiple sequences.
	template<typename... Seqs>
	using seq_sub_t = typename seq_sub_result<Seqs...>::type;

	//! Alias for the filter result of multiple sequences.
	/*!
	 * Filter from the first sequence everything that's in the subsequent sequences. That is,
	 * any element in the sequences after the first will be removed from the first.
	 */
	template<typename... Seqs>
	using filter_seq_t = typename filter_seq_result<Seqs...>::type;

	//! Alias for the intersection result of multiple sequences.
	template<typename... Seqs>
	using intersect_seq_t = typename intersect_seq_result<Seqs...>::type;


	template<typename T, typename Seq>
	struct change_seq_type_impl;

	template<typename T, typename T0, T0... Is>
	struct change_seq_type_impl<T, std::integer_sequence<T0, Is...>>
	{
		using type = std::integer_sequence<T, T(Is)...>;
	};

	template<typename T, typename Seq>
	using change_seq_type_t = typename change_seq_type_impl<T, Seq>::type;

	//! Returns whether there is the given value in the sequence.
	/*!
	 * Returns true if the value `N` is in the given sequence.
	 * 
	 * \tparam N The existence of this value in the sequence is checked.
	 * \tparam T The sequence to check.
	 */
	template<typename T, T N, typename Seq>
	struct is_value_in_seq;

	template<typename T, T N>
	struct is_value_in_seq<T, N, std::integer_sequence<T>>
	{
		static const bool value = false;
	};

	template<typename T, T N, T I0, T... Is>
	struct is_value_in_seq<T, N, std::integer_sequence<T, I0, Is...>>
	{
		static const bool value = (I0 == N) || (is_value_in_seq<T, N, std::integer_sequence<T, Is...>>::value);
	};


	//! Get the value at the specified index in the sequence.
	/*!
	 * Provides the value at the given index in the sequence.
	 * 
	 * \tparam N The index to get the value from.
	 * \tparam T The sequence.
	 */
	template<size_t N, typename T>
	struct seq_index_value_impl;

	template<size_t N, typename T>
	struct seq_index_value_impl<N, std::integer_sequence<T>>;

	template<typename T, T I0, T... Is>
	struct seq_index_value_impl<0, std::integer_sequence<T, I0, Is...>>
	{
		static const T value = I0;
	};

	template<size_t N, typename T, T I0, T... Is>
	struct seq_index_value_impl<N, std::integer_sequence<T, I0, Is...>>
	{
		static const T value = seq_index_value_impl<N - 1, std::integer_sequence<T, Is...>>::value;
	};

	template<int N, typename T>
	using seq_index_value = std::conditional_t<(N < 0),
		seq_index_value_impl<size_t(T::size() + N), T>,
		seq_index_value_impl<size_t(N), T>>;

	template<typename Pick, typename From>
	struct select_seq_impl;

	template<typename T, size_t... Is, T... Ns>
	struct select_seq_impl<std::index_sequence<Is...>, std::integer_sequence<T, Ns...>>
	{
		using type = std::integer_sequence<T, seq_index_value<Is, std::integer_sequence<T, Ns...>>::value...>;
	};

	template<typename Pick, typename From>
	using select_seq = typename select_seq_impl<Pick, From>::type;


	//! Create a sequence made of a single value.
	/*!
	 * The sequence is equal to repeating the given parameter `I` a 
	 * total of `N` times.
	 * 
	 * \tparam T The type of the repeated element.
	 * \tparam N The number of repeated elements.
	 * \tparam I The value of the repeated element.
	 */
	template<size_t N, typename T, T I>
	struct seq_repeating_value
	{
	protected:

		template<size_t N0>
		static constexpr T V = I;

		template<size_t... Ns>
		static constexpr auto build_sequence(std::index_sequence<Ns...>)
		{
			return std::integer_sequence<T, V<Ns>...>{};
		}

	public:

		using type = decltype(build_sequence(std::make_index_sequence<N>{}));
	};

	template<size_t N, typename T, T I>
	using seq_repeating_value_t = typename seq_repeating_value<N, T, I>::type;


	//! Makes a new sequence with values shared between all the sequences.
	/*!
	 * All values that are shared by all provided sequences will be put into a new
	 * sequence.
	 */
	template<typename T, T... Ys, T... Qs, typename... Seqs>
	constexpr auto intersect_seq(std::integer_sequence<T, Ys...>, std::integer_sequence<T, Qs...>, Seqs... seqs)
	{
		using filtered_t = seq_join_t<std::integer_sequence<T>, 
			std::conditional_t<
				is_value_in_seq<T, Ys, std::integer_sequence<T, Qs...>>::value,
				std::integer_sequence<T, Ys>,
				std::integer_sequence<T>>...
			>;
		return intersect_seq(filtered_t{}, seqs...);
	}

	namespace internal
	{

		template<typename Seq, typename... Seqs>
		size_t constexpr seq_len_product()
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
		size_t constexpr get_value_from_seq(std::integer_sequence<T, V>)
		{
			return V;
		}

		template<size_t N, typename T, T... Es, typename std::enable_if_t<(N < sizeof...(Es)), int> = 0>
		size_t constexpr seq_value(std::integer_sequence<T, Es...>);


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
				return std::integer_sequence<T, seq_value<N>(seq_t<E, Es...>{})>{};
			}

			template<size_t N, T E, T... Es, typename Seq, typename... Seqs, size_t L = seq_len_product<Seq, Seqs...>()>
			static auto constexpr select(seq_t<E, Es...>, Seq, Seqs...)
			{
				if constexpr (L == 0)
				{
					return select_non_empty_seq<N>(symphas::lib::types_list<>{}, symphas::lib::types_list<seq_t<E, Es...>, Seq, Seqs...>{});
				}
				else if constexpr (N < L)
				{
					constexpr size_t N0 = N / L;
					constexpr size_t N1 = N - N0 * L;
					return symphas::lib::seq_join(std::integer_sequence<T, seq_value<N0>(seq_t<E, Es...>{})>{}, select<N1>(Seq{}, Seqs{}...));
				}
				else
				{
					return select<N % L>(seq_t<E, Es...>{}, Seq{}, Seqs{}...);
				}
			}

			template<size_t N, typename... Seqs0, T I0, T... Is, typename... Seqs>
			static auto constexpr select_non_empty_seq(symphas::lib::types_list<Seqs0...>, symphas::lib::types_list<std::integer_sequence<T, I0, Is...>, Seqs...>)
			{
				return select_non_empty_seq<N>(symphas::lib::types_list<Seqs0..., std::integer_sequence<T, I0, Is...>>{}, symphas::lib::types_list<Seqs...>{});
			}

			template<size_t N, typename... Seqs0, typename... Seqs>
			static auto constexpr select_non_empty_seq(symphas::lib::types_list<Seqs0...>, symphas::lib::types_list<std::integer_sequence<T>, Seqs...>)
			{
				return select_non_empty_seq<N>(symphas::lib::types_list<Seqs0...>{}, symphas::lib::types_list<Seqs...>{});
			}

			template<size_t N, typename... Seqs0, typename... Seqs>
			static auto constexpr select_non_empty_seq(symphas::lib::types_list<Seqs0...>, symphas::lib::types_list<>)
			{
				return select<N>(Seqs0{}...);
			}
		};
	}


	template<typename T>
	struct types_list_size;

	template<typename... Ts>
	struct types_list_size<types_list<Ts...>>
	{
		static const size_t value = sizeof...(Ts);
	};

	template<typename... Ts>
	struct expand_types_list_impl;

	template<>
	struct expand_types_list_impl<>
	{
		using type = types_list<>;
	};

	template<typename T>
	struct expand_types_list_impl<T>
	{
		using type = types_list<T>;
	};

	template<>
	struct expand_types_list_impl<types_list<>>
	{
		using type = types_list<>;
	};

	template<typename T>
	struct expand_types_list_impl<types_list<T>>
	{
		using type = typename expand_types_list_impl<T>::type;
	};

	template<typename... T1s, typename... T2s>
	struct expand_types_list_impl<types_list<T1s...>, types_list<T2s...>>
	{
		using type = types_list<T1s..., T2s...>;
	};

	template<typename T1, typename... T2s>
	struct expand_types_list_impl<T1, types_list<T2s...>>
	{
		using type = types_list<T1, T2s...>;
	};

	template<typename... T1s, typename T2>
	struct expand_types_list_impl<types_list<T1s...>, T2>
	{
		using type = types_list<T1s..., T2>;
	};

	template<typename T1, typename T2>
	struct expand_types_list_impl<T1, T2>
	{
		using type = types_list<T1, T2>;
	};

	template<typename T0, typename T1, typename... Ts>
	struct expand_types_list_impl<T0, T1, Ts...>
	{
		using type = typename expand_types_list_impl<typename expand_types_list_impl<T0, T1>::type, Ts...>::type;
	};

	template<typename... Ts>
	using expand_types_list = typename expand_types_list_impl<Ts...>::type;

	template<typename Seq0, typename Seq1>
	struct make_seq_pairs;

	template<size_t... Is, size_t... Js>
	struct make_seq_pairs<std::index_sequence<Is...>, std::index_sequence<Js...>>
	{
		using type = expand_types_list<std::conditional_t<(Is < Js), std::index_sequence<Is, Js>, types_list<>>...>;
	};

	template<size_t N, typename... Seqs>
	struct seq_skip_indices_impl;

	template<>
	struct seq_skip_indices_impl<0, std::index_sequence<0>>
	{
		using type = std::index_sequence<>;
	};

	template<size_t N>
	struct seq_skip_indices_impl<N, std::index_sequence<N - 1>>
	{
		using type = std::make_index_sequence<N - 1>;
	};

	template<size_t N>
	struct seq_skip_indices_impl<N, std::index_sequence<>>
	{
		using type = std::make_index_sequence<N>;
	};

	template<size_t N>
	struct seq_skip_indices_impl<N, std::index_sequence<0>>
	{
		using type = seq_offset_t<1, std::make_index_sequence<N - 1>>;
	};

	template<size_t N, size_t I0>
	struct seq_skip_indices_impl<N, std::index_sequence<I0>>
	{
		using type = seq_join_t<std::make_index_sequence<I0>, seq_offset_t<I0 + 1, std::make_index_sequence<N - I0 - 1>>>;
	};

	template<size_t N, size_t I0, size_t J0>
	struct seq_skip_indices_impl<N, types_list<std::index_sequence<I0, J0>>>
	{
		using type = seq_join_t<seq_from_to_t<std::index_sequence<I0, J0>>>;
	};

	template<size_t N, typename Seq0, typename Seq1, typename... Seqs>
	struct seq_skip_indices_impl<N, Seq0, Seq1, Seqs...>
	{
		using type = seq_join_t<seq_from_to_t<Seq0>, seq_from_to_t<Seq1>, seq_from_to_t<Seqs>...>;
	};

	template<size_t N, typename Seq0, typename Seq1, typename... Seqs>
	struct seq_skip_indices_impl<N, types_list<Seq0, Seq1, Seqs...>>
	{
		using type = typename seq_skip_indices_impl<N, Seq0, Seq1, Seqs...>::type;
	};

	template<size_t N>
	struct seq_skip_indices_impl<N, types_list<>>
	{
		using type = std::index_sequence<>;
	};


	template<size_t N, size_t I0, size_t I1, size_t... Is>
	struct seq_skip_indices_impl<N, std::index_sequence<I0, I1, Is...>>
	{
        using type = std::conditional_t<
            std::is_same<std::make_index_sequence<N>, std::index_sequence<I0, I1, Is...>>::value,
            std::index_sequence<>, 
            typename seq_skip_indices_impl<N, 
                typename make_seq_pairs<
                    std::index_sequence<0, I0 + 1, I1 + 1, (Is + 1)...>, 
                    std::index_sequence<I0, I1, Is..., N>>::type>::type>;
	};

	template<size_t N, typename Seq>
	using seq_skip_indices = typename seq_skip_indices_impl<N, Seq>::type;

	template<size_t N, typename T>
	struct list_repeating_type
	{
	protected:

		template<size_t N0>
		using repeat_type = T;

		template<size_t... Ns>
		static constexpr auto build_list(std::index_sequence<Ns...>)
		{
			return types_list<repeat_type<Ns>...>{};
		}

	public:

		using type = decltype(build_list(std::make_index_sequence<N>{}));
	};

	template<size_t N, typename T>
	using list_repeating_type_t = typename list_repeating_type<N, T>::type;



	template<size_t I, typename Seq0>
	struct seq_ge_impl;

	template<size_t I, typename T, T... Ns>
	struct seq_ge_impl<I, std::integer_sequence<T, Ns...>>
	{
		using type = std::conditional_t<
			(sizeof...(Ns) <= I),
			std::integer_sequence<T>,
			symphas::lib::select_seq<seq_offset_t<I, std::make_index_sequence<fixed_min<sizeof...(Ns) - I, sizeof...(Ns)>>>, std::integer_sequence<T, Ns...>>
		>;
	};

	template<size_t I, typename Seq0>
	using seq_ge_t = typename seq_ge_impl<I, Seq0>::type;

	template<size_t I, typename Seq0>
	struct seq_lt_impl;

	template<size_t I, typename T, T... Ns>
	struct seq_lt_impl<I, std::integer_sequence<T, Ns...>>
	{
		using type = symphas::lib::select_seq<std::make_index_sequence<fixed_min<I, sizeof...(Ns)>>, std::integer_sequence<T, Ns...>>;
	};

	template<size_t I, typename Seq0>
	using seq_lt_t = typename seq_lt_impl<I, Seq0>::type;


	template<typename T, T N, typename... Seqs>
	struct find_index_lt;

	// find the largest index of the sequence that is smaller than N
	template<typename T, T N, T I0>
	struct find_index_lt<T, N, std::integer_sequence<T, I0>>
	{
		static const size_t value = (I0 < N) ? 0 : -1;
	};

	// find the largest index of the sequence that is smaller than N
	template<typename T, T N, T I0, T I1, T... Is>
	struct find_index_lt<T, N, std::integer_sequence<T, I0, I1, Is...>>
	{
		static const size_t value = 
			(seq_index_value<sizeof...(Is) / 2 + 1, std::integer_sequence<T, I0, I1, Is...>>::value <= N)
				? sizeof...(Is) / 2 + 1 + find_index_lt<T, N, seq_ge_t<sizeof...(Is) / 2 + 1, std::integer_sequence<T, I0, I1, Is...>>>::value
				: find_index_lt<T, N, seq_lt_t<sizeof...(Is) / 2 + 1, std::integer_sequence<T, I0, I1, Is...>>>::value;
	};

	template<typename T, T N, typename... Seqs>
	struct find_index_ge;

	// find the largest index of the sequence that is smaller than N
	template<typename T, T N, T I0>
	struct find_index_ge<T, N, std::integer_sequence<T, I0>>
	{
		static const size_t value = (I0 >= N) ? 0 : 1;
	};

	// find the largest index of the sequence that is smaller than N
	template<typename T, T N, T I0, T I1, T... Is>
	struct find_index_ge<T, N, std::integer_sequence<T, I0, I1, Is...>>
	{
		static const size_t value =
			(seq_index_value<sizeof...(Is) / 2, std::integer_sequence<T, I0, I1, Is...>>::value < N)
			? sizeof...(Is) / 2 + 1 + find_index_ge<T, N, seq_ge_t<sizeof...(Is) / 2 + 1, std::integer_sequence<T, I0, I1, Is...>>>::value
			: find_index_ge<T, N, seq_lt_t<sizeof...(Is) / 2 + 1, std::integer_sequence<T, I0, I1, Is...>>>::value;
	};


	template<typename T, T I, typename Seq0>
	struct seq_before_impl;

	template<typename T, T I, T... Ns>
	struct seq_before_impl<T, I, std::integer_sequence<T, Ns...>>
	{
		using type = symphas::lib::select_seq<
			std::make_index_sequence<find_index_lt<T, I, std::integer_sequence<T, Ns...>>::value + 1>,
			std::integer_sequence<T, Ns...>>;
	};

	template<typename T, T I, typename Seq0>
	using seq_before_t = typename seq_before_impl<T, I, Seq0>::type;


	template<typename T, T I, typename Seq0>
	struct seq_after_at_impl;

	template<typename T, T I, T... Ns>
	struct seq_after_at_impl<T, I, std::integer_sequence<T, Ns...>>
	{
	protected:

		static const size_t N = find_index_ge<T, I, std::integer_sequence<T, Ns...>>::value;

	public:

		using type = symphas::lib::select_seq<
			seq_offset_t<N, std::make_index_sequence<sizeof...(Ns) - N>>, 
			std::integer_sequence<T, Ns...>>;
	};

	template<typename T, T I, typename Seq0>
	using seq_after_at_t = typename seq_after_at_impl<T, I, Seq0>::type;



	template<bool, typename Seq0, typename Seq1>
	struct select_merge_sort_seq_impl;

	template<typename Seq0, typename Seq1>
	struct merge_sort_seq_impl;

	template<typename T, T N0, T... Ns>
	struct merge_sort_seq_impl<std::integer_sequence<T, N0, Ns...>, std::integer_sequence<T>>
	{
		using type = std::integer_sequence<T, N0, Ns...>;
	};

	template<typename T, T M0, T... Ms>
	struct merge_sort_seq_impl<std::integer_sequence<T>, std::integer_sequence<T, M0, Ms...>>
	{
		using type = std::integer_sequence<T, M0, Ms...>;
	};

	template<typename T, T N0, T... Ns, T M0, T... Ms>
	struct merge_sort_seq_impl<std::integer_sequence<T, N0, Ns...>, std::integer_sequence<T, M0, Ms...>>
	{
		using type = typename select_merge_sort_seq_impl<(N0 < M0), std::integer_sequence<T, N0, Ns...>, std::integer_sequence<T, M0, Ms...>>::type;
	};


	template<typename T, T N0, T... Ns, T M0, T... Ms>
	struct select_merge_sort_seq_impl<true, std::integer_sequence<T, N0, Ns...>, std::integer_sequence<T, M0, Ms...>>
	{
		using type = typename merge_sort_seq_impl<
			seq_join_t<seq_before_t<T, M0, std::integer_sequence<T, N0, Ns...>>, std::integer_sequence<T, M0>, seq_after_at_t<T, M0, std::integer_sequence<T, N0, Ns...>>>,
			std::integer_sequence<T, Ms...>>::type;
	};

	template<typename T, T N0, T... Ns, T M0, T... Ms>
	struct select_merge_sort_seq_impl<false, std::integer_sequence<T, N0, Ns...>, std::integer_sequence<T, M0, Ms...>>
	{
		using type = typename merge_sort_seq_impl<
			seq_join_t<seq_before_t<T, N0, std::integer_sequence<T, M0, Ms...>>, std::integer_sequence<T, N0>, seq_after_at_t<T, N0, std::integer_sequence<T, M0, Ms...>>>,
			std::integer_sequence<T, Ns...>>::type;
	};

	template<typename T, T N0, T... Ns, T... Ms>
	struct select_merge_sort_seq_impl<false, std::integer_sequence<T, N0, Ns...>, std::integer_sequence<T, N0, Ms...>>
	{
		using type = typename merge_sort_seq_impl<
			std::integer_sequence<T, N0, N0, Ns...>,
			std::integer_sequence<T, Ms...>>::type;
	};

	template<typename Seq0, typename Seq1>
	using merge_sort_seq = typename merge_sort_seq_impl<Seq0, Seq1>::type;

	template<typename Seq>
	struct sorted_seq_impl;

	template<typename T>
	struct sorted_seq_impl<std::integer_sequence<T>>
	{
		using type = std::integer_sequence<T>;
	};

	template<typename T, T N0>
	struct sorted_seq_impl<std::integer_sequence<T, N0>>
	{
		using type = std::integer_sequence<T, N0>;
	};

	template<typename T, T N0, T N1, T... Ns>
	struct sorted_seq_impl<std::integer_sequence<T, N0, N1, Ns...>>
	{
		using type = merge_sort_seq<
			typename sorted_seq_impl<seq_ge_t<sizeof...(Ns) / 2 + 1, std::integer_sequence<T, N0, N1, Ns...>>>::type,
			typename sorted_seq_impl<seq_lt_t<sizeof...(Ns) / 2 + 1, std::integer_sequence<T, N0, N1, Ns...>>>::type>;
	};

	template<typename Seq0>
	using sorted_seq = typename sorted_seq_impl<Seq0>::type;


	//! Puts the nth value of all sequences into the result.
	template<size_t N, typename... Ts>
	struct nth_value_of_seqs
	{
		using type = std::integer_sequence<bool, symphas::lib::seq_index_value<N, Ts>::value...>;
	};

	template<size_t N, typename... Ts>
	using nth_value_of_seqs_t = typename nth_value_of_seqs<N, Ts...>::type;

	//! Collect and count like types from a tuple list.
	/*!
	 * Give the list of types, the list is made unique and number of times
	 * each type appears in the list is repeated.
	 */
	template<typename... Gs>
	struct filter_types_impl;

	template<typename... G0s>
	struct filter_types_impl<types_list<G0s...>>
	{
		using type = types_list<G0s...>;
	};

	template<typename... G0s>
	struct filter_types_impl<types_list<G0s...>, types_list<>>
	{
		using type = types_list<G0s...>;
	};

	template<bool... flags, typename... G0s>
	struct filter_types_impl<types_list<G0s...>, std::integer_sequence<bool, flags...>>
	{
		using type = expand_types_list<std::conditional_t<flags, types_list<G0s>, types_list<>>...>;
	};

	template<typename... G0s, typename G1, typename... G1s>
	struct filter_types_impl<types_list<G0s...>, types_list<G1, G1s...>>
	{
	protected:

		using filter_mask = std::integer_sequence<bool, !std::is_same<G1, G0s>::value...>;
		using filtered_1 = typename filter_types_impl<types_list<G0s...>, filter_mask>::type;

	public:

		using type = typename filter_types_impl<filtered_1, types_list<G1s...>>::type;

	};

	template<typename... G0s, typename... G1s, typename... Rest>
	struct filter_types_impl<types_list<G0s...>, types_list<G1s...>, Rest...>
	{
		using type = typename filter_types_impl<
			typename filter_types_impl<types_list<G0s...>, types_list<G1s...>>::type,
			Rest...>::type;

	};

	//! Filter types from the first type list.
	/*!
	 * Filter from the first list everything that's in the subsequent lists. That is,
	 * any element in the lists after the first will be removed from the first.
	 * 
	 * If an std::integer_sequence<bool> type is passed as the second argument, it will act
	 * as a mask.
	 */
	template<typename... Gs>
	using filter_types = typename filter_types_impl<Gs...>::type;

	template<typename... Gs>
	struct cc_like_types;

	template<>
	struct cc_like_types<>
	{
		using type = types_list<>;
		using count_seq = std::index_sequence<>;
	};

	template<typename G0, typename... Gs>
	struct cc_like_types<G0, Gs...>
	{
	protected:

		using filtered = symphas::lib::filter_types<types_list<G0, Gs...>, types_list<G0>>;
		using like_types_rest = cc_like_types<filtered>;
		static const size_t count = (((std::is_same<G0, Gs>::value) ? 1 : 0) + ... + 1);

	public:

		using type = expand_types_list<G0, typename like_types_rest::type>;
		using count_seq = seq_join_t<std::index_sequence<count>, typename like_types_rest::count_seq>;
	};


	template<typename... Gs>
	struct cc_like_types<types_list<Gs...>>
	{
		using type = typename cc_like_types<Gs...>::type;
		using count_seq = typename cc_like_types<Gs...>::count_seq;
	};

	//! Returns the unique list of types from the list.
	/*!
	 * Returns the unique list of types is returned from joining together
	 * the list of types.
	 */
	template<typename... Gs>
	struct combine_types_unique
	{
		using type = typename cc_like_types<expand_types_list<Gs...>>::type;
	};

	template<bool flag, size_t N, size_t D, size_t I, typename T0, typename... Ts>
	struct nth_periodic_shift_impl;

	template<size_t N, size_t D, typename T0, typename... Ts>
	struct nth_periodic_shift_impl<true, N, D, D, T0, Ts...>
	{
		using type = typename nth_value_of_seqs<N, T0, Ts...>::type;
	};

	template<size_t N, size_t D, size_t I, typename T0, typename... Ts>
	struct nth_periodic_shift_impl<false, N, D, I, T0, Ts...>
	{
		using type = typename nth_periodic_shift_impl<D == I + 1, N, D, I + 1,
			symphas::lib::seq_join_t<T0, T0>,
			symphas::lib::seq_join_t<Ts, Ts>...,
			symphas::lib::seq_join_t<symphas::lib::seq_repeating_value_t<T0::size(), bool, 0>, symphas::lib::seq_repeating_value_t<T0::size(), bool, 1>>
		>::type;
	};

	template<size_t N, size_t D>
	struct nth_periodic_shift
	{
		using type = typename nth_periodic_shift_impl<D == 1, N, D, 1, std::integer_sequence<bool, 0, 1>>::type;
	};

	template<size_t N, size_t D>
	using nth_periodic_shift_t = typename nth_periodic_shift<N, D>::type;


	namespace
	{
		// helper function that iterates over the tuple
		template<typename... As, typename... Bs, typename T, T... Is>
		auto unfurl_tuple(std::tuple<std::pair<As, Bs>...> const& ts, std::integer_sequence<T, Is...>)
		{
			auto first_elements = symphas::lib::make_tuple(std::get<Is>(ts).first...);
			auto second_elements = symphas::lib::make_tuple(std::get<Is>(ts).second...);
			return std::make_pair(first_elements, second_elements);
		}

		// helper function that takes index sequence and rebuilds the tuple from the 
		// given index onwards
		template<size_t I, typename... Ts, size_t... Is>
		auto get_tuple_ge(std::tuple<Ts...> const& ts, std::index_sequence<Is... >)
		{
			return symphas::lib::make_tuple(std::get<I + Is>(ts)...);
		}

		// helper function that takes index sequence and rebuilds the tuple up to and
		// not including the given index
		template<typename... Ts, typename T, T... Is>
		auto get_tuple_lt(std::tuple<Ts...> const& ts, std::integer_sequence<T, Is...>)
		{
			return symphas::lib::make_tuple(std::get<Is>(ts)...);
		}


		template<typename T, T... Ns>
		auto get_seq_from_tuple(std::tuple<std::integer_sequence<T, Ns>...> const&)
		{
			return std::integer_sequence<T, Ns...>{};
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

	template<typename data_type>
	struct basic_forward_iterator_container;


	template<typename data_type>
	struct basic_forward_iterator
	{
		using iterator_category = std::forward_iterator_tag;
		using difference_type = int;
		using value_type = basic_forward_iterator_container<data_type>;
		using pointer = value_type*;
		using reference = int;

		template<typename T>
		basic_forward_iterator(T&& data, iter_type pos = 0) :
			ptr{ std::forward<T>(data), pos } {}

		basic_forward_iterator_container<data_type> operator*() const
		{
			return ptr;
		}

		// Prefix increment
		basic_forward_iterator& operator++()
		{
			++ptr.pos;
			return *this;
		}

		// Postfix increment
		basic_forward_iterator operator++(int)
		{
			basic_forward_iterator tmp = *this;
			++(*this);
			return tmp;
		}

		friend bool operator==(basic_forward_iterator const& a, basic_forward_iterator const& b)
		{
			return a.ptr.pos == b.ptr.pos;
		}

		friend bool operator!=(basic_forward_iterator const& a, basic_forward_iterator const& b)
		{
			return !(a == b);
		}

		basic_forward_iterator_container<data_type> ptr;
	};


	template<typename... Ts>
	struct zip_iterator
	{
		using iterator_category = std::forward_iterator_tag;
		using difference_type = int;
		using value_type = std::tuple<Ts...>;
		using pointer = value_type*; 
		using reference = int;  

		using seq = std::make_index_sequence<sizeof...(Ts)>;

		zip_iterator(Ts const&... ts) :
			data{ std::begin(ts)... }, n{ 0 }//, end{ std::end(ts)... }
		{
		}

		decltype(auto) operator*() const
		{ 
			return get_data(seq{});
		}

		decltype(auto) operator[](iter_type offset) const
		{
			return get_data(seq{}, offset);
		}
		
		// Prefix increment
		zip_iterator<Ts...>& operator++() 
		{ 
			n += 1;
			return *this; 
		}

		// Postfix increment
		zip_iterator<Ts...> operator++(int)
		{ 								 
			zip_iterator<Ts...> tmp = *this;
			++(*this); 
			return tmp; 
		}

		friend bool operator==(zip_iterator<Ts...> const& a, zip_iterator<Ts...> const& b)
		{ 
			return a.check_all_same(b, seq{});
		}

		friend bool operator!=(zip_iterator<Ts...> const& a, zip_iterator<Ts...> const& b)
		{
			return a.check_different(b, seq{});
		}

		template<size_t... Is>
		bool check_all_same(zip_iterator<Ts...> const& b, std::index_sequence<Is...>) const
		{
			return (check_one_same<Is>(b) && ...);
		}

		template<size_t... Is>
		bool check_different(zip_iterator<Ts...> const& b, std::index_sequence<Is...>) const
		{
			return (!check_one_same<Is>(b) || ...);
		}

		template<size_t I>
		bool check_one_same(zip_iterator<Ts...> const& b) const
		{
			if (n == b.n)
			{
				return std::get<I>(data) == std::get<I>(b.data);
			}
			else
			{
				return std::get<I>(b.data) - std::get<I>(data) == b.n - n;
			}
		}

		template<size_t... Is>
		decltype(auto) get_data(std::index_sequence<Is...>, iter_type offset = 0) const
		{
			return std::tie(std::get<Is>(data)[n + offset]...);
		}

		template<size_t... Is>
		decltype(auto) get_data(std::index_sequence<Is...>, iter_type offset = 0)
		{
			return std::tie(std::get<Is>(data)[n + offset]...);
		}
		
	protected:
		
		template<typename T>
		static auto begin(T const& t)
		{
			return std::begin(t);
		}

	public:
		
		std::tuple<std::invoke_result_t<decltype(&zip_iterator<Ts...>::begin<Ts>), Ts>...> data;
		iter_type n;

	};


	template<typename T0, typename... Ts>
	struct zip_container
	{

		zip_container(T0 const& t0, Ts const&... ts) :
			iter{ t0, ts... }, len{ static_cast<len_type>(std::end(t0) - std::begin(t0)) } {}


		auto begin()
		{
			return iter;
		}

		auto end()
		{
			zip_iterator end(iter);
			end.n = len;
			return end;
		}

		auto begin() const
		{
			return iter;
		}

		auto end() const
		{
			zip_iterator end(iter);
			end.n = len;
			return end;
		}

		zip_iterator<T0, Ts...> iter;
		len_type len;
	};


	template<typename T0, typename... Ts>
	struct zip_container<T0&&, Ts&&...>
	{
		zip_container(T0&& t0, Ts&&... ts) :
			data{ t0, ts... },
			iter{ get_iter() },
			len{ static_cast<len_type>(std::end(std::get<0>(data)) - std::begin(std::get<0>(data))) }
		{
		}

		auto begin()
		{
			return iter;
		}

		auto end()
		{
			zip_iterator end(iter);
			end.n = len;
			return end;
		}

		auto begin() const
		{
			return iter;
		}

		auto end() const
		{
			zip_iterator end(iter);
			end.n = len;
			return end;
		}

		std::tuple<T0, Ts...> data;
		zip_iterator<T0, Ts...> iter;
		len_type len;

	protected:

		auto get_iter() const
		{
			return get_iter(std::make_index_sequence<1 + sizeof...(Ts)>{});
		}
		
		template<size_t... Is>
		auto get_iter(std::index_sequence<Is...>) const
		{
			return zip_iterator(std::get<Is>(data)...);
		}
	};


	template<typename... Ts>
	zip_container(Ts const&...) -> zip_container<Ts...>;

	template<typename... Ts>
	zip_container(Ts&&...) -> zip_container<Ts&&...>;


	template<typename... Ts>
	auto zip(Ts&&... ts)
	{
		return zip_container(std::forward<Ts>(ts)...);
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


	// **************************************************************************************

	//! Return the index of the given type in the type list.
	/*!
	 * Given the type list, identify the index of the `I`-th chosen type. If
	 * the type is not identified, then `-1` will be returned.
	 *
	 * \tparam Type The type to be found in the list.
	 * \tparam I The index of the type to find out of those that match.
	 * \tparam Ss... The types in the list.
	 */
	template<typename Type, size_t I, typename... Ss, size_t... Is>
	constexpr int index_of_type_impl(std::index_sequence<Is...>)
	{
		using mask_t = seq_join_t<std::index_sequence<>,
			std::conditional_t<
			(std::is_same<Ss, Type>::value),
			std::index_sequence<Is>,
			std::index_sequence<>>...>;

		if constexpr (I < mask_t::size())
		{
			return seq_index_value<I, mask_t>::value;
		}
		else
		{
			return -1;
		}
	}

	template<typename T>
	struct unroll_types_list {};

	//! Return the index of the given type in the type list.
	/*!
	 * Given the type list, identify the index of the `I`-th chosen type. If
	 * the type is not identified, then `-1` will be returned.
	 *
	 * \tparam Type The type to be found in the list.
	 * \tparam I The index of the type to find out of those that match.
	 * \tparam S0 The next type in the type list to compare.
	 * \tparam Ss... The remaining types in the list.
	 */
	template<typename Type, size_t N, typename ...Ss>
	constexpr int nth_index_of_type = index_of_type_impl<Type, N, Ss...>(std::make_index_sequence<sizeof...(Ss)>{});

	//! Return the index of the given type in the type list.
	/*!
	 * Given the type list, identify the index of the `I`-th chosen type. If
	 * the type is not identified, then `-1` will be returned.
	 *
	 * \tparam Type The type to be found in the list.
	 * \tparam S0 The next type in the type list to compare.
	 * \tparam Ss... The remaining types in the list.
	 */
	template<typename Type, typename ...Ss>
	constexpr int index_of_type = nth_index_of_type<Type, 0, Ss...>;

	template<typename Type>
	constexpr int index_of_type<Type> = -1;

	template<typename Type, Type N, Type... Is>
	constexpr int index_of_value = index_of_type<std::integer_sequence<Type, N>, std::integer_sequence<Type, Is>...>;

	template<typename Type, size_t N, typename S0, typename ...Ss>
	constexpr int nth_index_of_type<Type, N, types_list<S0, Ss...>> = nth_index_of_type<Type, N, S0, Ss...>;

	template<typename Type, typename S0, typename ...Ss>
	constexpr int index_of_type<Type, types_list<S0, Ss...>> = index_of_type<Type, S0, Ss...>;

	template<typename Type, typename S0, typename ...Ss>
	constexpr int index_of_type<Type, unroll_types_list<types_list<S0, Ss...>>> = index_of_type<Type, S0, Ss...>;



	template<size_t I, typename... Ts>
	struct type_at_index_impl;

	//template<>
	//struct type_at_index_impl<0, types_list<>>
	//{
	//	using type = void;
	//};

	template<typename T0>
	struct type_at_index_impl<0, T0>
	{
		using type = T0;
	};

	template<typename T0, typename T1, typename... Ts>
	struct type_at_index_impl<0, T0, T1, Ts...>
	{
		using type = T0;
	};

	template<size_t I, typename T0, typename... Ts>
	struct type_at_index_impl<I, T0, Ts...>
	{
		using type = typename type_at_index_impl<I - 1, Ts...>::type;
	};

	template<>
	struct type_at_index_impl<0, unroll_types_list<types_list<>>>
	{
		using type = void;
	};

	template<typename T0, typename... Ts>
	struct type_at_index_impl<0, unroll_types_list<types_list<T0, Ts...>>>
	{
		using type = T0;
	};

	template<size_t I, typename... Ts>
	struct type_at_index_impl<I, unroll_types_list<types_list<Ts...>>>
	{
		using type = typename type_at_index_impl<I, Ts...>::type;
	};

	template<size_t I, typename... Ts>
	using type_at_index = typename type_at_index_impl<I, Ts...>::type;

	template<size_t I, typename... Ts>
	struct types_after_at_index_impl;


	template<size_t I>
	struct types_after_at_index_impl<I>
	{
		using type = types_list<>;
	};

	template<size_t I, typename T0, typename... Ts>
	struct types_after_at_index_impl<I, T0, Ts...>
	{
		using type = std::conditional_t<
			(I == 0),
			types_list<T0, Ts...>,
			typename types_after_at_index_impl<I - 1, Ts...>::type>;
	};

	template<size_t I, typename... Ts>
	struct types_after_at_index_impl<I, unroll_types_list<types_list<Ts...>>>
	{
		using type = typename types_after_at_index_impl<I, Ts...>::type;
	};

	template<size_t I, typename... Ts>
	using types_after_at_index = typename types_after_at_index_impl<I, Ts...>::type;


	template<typename Seq, typename... Ts>
	struct select_types_impl;

	template<size_t... Is, typename... Ts>
	struct select_types_impl<std::index_sequence<Is...>, Ts...>
	{
		using type = types_list<type_at_index<Is, Ts...>...>;
	};

	template<size_t... Is, typename... Ts>
	struct select_types_impl<std::index_sequence<Is...>, unroll_types_list<types_list<Ts...>>>
	{
		using type = types_list<type_at_index<Is, Ts...>...>;
	};
	

	template<typename Seq, typename... Ts>
	using select_types = typename select_types_impl<Seq, Ts...>::type;

	template<typename Seq, typename... Ts>
	struct filter_types_on_index_impl;

	template<size_t... Is, typename... Ts>
	struct filter_types_on_index_impl<std::index_sequence<Is...>, Ts...>
	{
		using type = select_types<symphas::lib::filter_seq_t<std::make_index_sequence<sizeof...(Ts)>, std::index_sequence<Is...>>, Ts...>;
	};

	template<size_t... Is, typename... Ts>
	struct filter_types_on_index_impl<std::index_sequence<Is...>, unroll_types_list<types_list<Ts...>>>
	{
		using type = select_types<symphas::lib::filter_seq_t<std::make_index_sequence<sizeof...(Ts)>, std::index_sequence<Is...>>, Ts...>;
	};

	template<typename Seq, typename... Ts>
	using filter_types_on_index = typename filter_types_on_index_impl<Seq, Ts...>::type;


	template<size_t I, typename T, typename... Ts>
	struct types_before_index_impl;

	template<size_t I, typename... T0s>
	struct types_before_index_impl<I, types_list<T0s...>>
	{
		using type = types_list<T0s...>;
	};

	template<typename... T0s, typename T0, typename... Ts>
	struct types_before_index_impl<0, types_list<T0s...>, T0, Ts...>
	{
		using type = types_list<T0s...>;
	};

	template<size_t I, typename... T0s, typename T0, typename... Ts>
	struct types_before_index_impl<I, types_list<T0s...>, T0, Ts...>
	{
		using type = typename types_before_index_impl<I - 1, types_list<T0s..., T0>, Ts...>::type;
	};

	template<size_t I, typename... Ts>
	struct types_before_index_impl<I, types_list<>, unroll_types_list<types_list<Ts...>>>
	{
		using type = typename types_before_index_impl<I, types_list<>, Ts...>::type;
	};

	template<typename T0, typename... Ts>
	struct types_before_index_impl<0, types_list<>, unroll_types_list<types_list<T0, Ts...>>>
	{
		using type = types_list<T0>;
	};

	template<size_t I, typename... Ts>
	using types_before_index = typename types_before_index_impl<I, types_list<>, Ts...>::type;


	template<size_t I0, size_t I, typename... Ts>
	struct types_between_index_impl
	{
		using type = types_after_at_index<I0, unroll_types_list<types_before_index<I, Ts...>>>;
	};

	template<size_t I0, size_t I, typename... Ts>
	using types_between_index = typename types_between_index_impl<I0, I, Ts...>::type;



	template<typename T, typename Seq>
	struct reverse_types_list_impl;

	template<typename... Es, size_t... Is>
	struct reverse_types_list_impl<types_list<Es...>, std::index_sequence<Is...>>
	{
		using type = types_list<type_at_index<sizeof...(Is) - 1 - Is, Es...>...>;
	};

	template<typename... Es>
	struct reverse_types_list_impl<types_list<Es...>, void>
	{
		using type = typename reverse_types_list_impl<types_list<Es...>, std::make_index_sequence<sizeof...(Es)>>::type;
	};

	//! Reverse a ::types_list.
	template<typename T>
	using reverse_types_list = typename reverse_types_list_impl<T, void>::type;


	// **************************************************************************************

    template<size_t N, typename T, T... Es, typename std::enable_if_t<(N < sizeof...(Es)), int>>
    size_t constexpr internal::seq_value(std::integer_sequence<T, Es...>)
    {
        return get_value_from_seq(symphas::lib::type_at_index<N, std::integer_sequence<T, Es>...>{});
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
		static const size_t count = 0;
		static const size_t rank = 0;
	};

	template<typename T, T... Es>
	struct CrossProductList<std::integer_sequence<T, Es...>>
	{
		static const size_t count = sizeof...(Es);
		static const size_t rank = 1;

		template<size_t N>
		static const size_t size = std::integer_sequence<T, Es...>::size();

		template<size_t N, typename std::enable_if_t<(N < count), int> = 0>
		using row = type_at_index<N, unroll_types_list<types_list<std::integer_sequence<T, Es>...>>>;

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
		static const size_t size = type_at_index<N, unroll_types_list<types_list<std::integer_sequence<T, E1s...>, std::integer_sequence<T, E2s...>>>>::size();

		template<size_t N>
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
		static const size_t size = type_at_index<N, unroll_types_list<types_list<std::integer_sequence<T, Es...>, List1, List2, Lists...>>>::size();

		template<size_t N>
		using row = decltype(internal::CrossProductFunctions<T>::template select<N>(
			std::declval<std::integer_sequence<T, Es...>>(),
			std::declval<List1>(),
			std::declval<List2>(),
			std::declval<Lists>()...));
	};

	template<typename T>
	struct array_container
	{
		T* data;
		size_t n;

		array_container() : array_container(0) {}
		array_container(size_t n) : data{ (n > 0) ? new T[n]{} : nullptr }, n{ n } {}
		array_container(const T* data, size_t n) : array_container(n)
		{
			if (data)
			{
				for (iter_type i = 0; i < n; ++i)
				{
					this->data[i] = data[i];
				}
			}
		}

		array_container(array_container<T> const& other) noexcept : array_container(other.data, other.n) {}

		array_container(array_container<T>&& other) noexcept : array_container()
		{
			swap(*this, other);
		}

		array_container& operator=(array_container other)
		{
			swap(*this, other);
			return *this;
		}

		~array_container()
		{
			delete[] data;
		}

		operator const T* () const
		{
			return data;
		}

		operator T* ()
		{
			return data;
		}

		T& operator[](iter_type i)
		{
			return data[i];
		}

		const T& operator[](iter_type i) const
		{
			return data[i];
		}

		friend void swap(array_container<T>& first, array_container<T>& second)
		{
			using std::swap;
			swap(first.data, second.data);
			swap(first.n, second.n);
		}

		const T* begin() const
		{
			return data;
		}

		const T* end() const
		{
			return data + n;
		}

		T* begin()
		{
			return data;
		}

		T* end()
		{
			return data + n;
		}

	};

	using string = array_container<char>;

	// **************************************************************************************


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


	inline int exponent10(double value)
	{
		double p = std::log10(std::abs(value));
		return (p < 0) ? int(p) - 1 : int(p);
	}

	inline double base10(double value)
	{
		return value * std::pow(10, -exponent10(value));
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

	//! Get the elements of the tuple greater or equal to the given index.
	/*!
	 * Returns a new tuple consisting of all the elements with index greater than or equal to
	 * the given index. If the given index is 0, then the original tuple will be returned.
	 */
	template<size_t I, typename... Ts>
	auto get_tuple_ge(std::tuple<Ts...> const& ts)
	{
		return get_tuple_ge<I>(ts, std::make_index_sequence<sizeof...(Ts) - I>{});
	}

	//! Get the elements of the tuple less than the given index.
	/*!
	 * Returns a new tuple consisting of all the elements with index strictly less than the
	 * given index. If the index is equal to the length of the tuple, the original tuple will 
	 * be returned.
	 */
	template<size_t I, typename... Ts>
	auto get_tuple_lt(std::tuple<Ts...> const& ts)
	{
		return get_tuple_lt(ts, std::make_index_sequence<I>{});
	}
	
	//! Get the elements of the tuple in the given range, not including the last index.
	/*!
	 * Returns a new tuple consisting of all the elements between the first and last indices
	 * provided, not including the last index. That is, this function will return the same tuple
	 * the last index is equal to the size of the tuple. 
	 */
	template<size_t I0, size_t In, typename... Ts, typename std::enable_if_t<(I0 < In), int> = 0>
	auto get_tuple_bw(std::tuple<Ts...> const& ts)
	{
		return get_tuple_lt(ts, symphas::lib::seq_add(symphas::lib::seq_repeating_value_t<In - I0, size_t, I0>{}, std::make_index_sequence<In - I0>{}));
	}

	template<size_t I0, size_t In, typename... Ts, typename std::enable_if_t<(I0 >= In), int> = 0>
	auto get_tuple_bw(std::tuple<Ts...> const& ts)
	{
		return std::tuple<>{};
	}

	template<size_t I, typename T, T... Ns>
	auto get_seq_ge(std::integer_sequence<T, Ns...> const& seq)
	{
		if constexpr (sizeof...(Ns) <= I)
		{
			return std::integer_sequence<T>{};
		}
		else
		{
			return symphas::lib::select_seq<seq_offset_t<I, std::make_index_sequence<sizeof...(Ns) - I>>, std::integer_sequence<T, Ns...>>{};
		}
	}

	template<size_t I, typename T, T... Ns>
	auto get_seq_lt(std::integer_sequence<T, Ns...> const& seq)
	{
		return symphas::lib::select_seq<std::make_index_sequence<I>, std::integer_sequence<T, Ns...>>{};
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
	void get_parent_directory(const char* path, const char* &basepath);

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



	template<Axis... axs>
	struct axis_list {};

	template<Side... sides>
	struct side_list {};


	//! Put the first `D` axes in a types list.
	/*!
	 * Construct a types list with the first `D` axes in the list.
	 */
	template<size_t D>
	auto make_axis_list();


	template<>
	inline auto make_axis_list<1>()
	{
		return axis_list<Axis::X>{};
	}

	template<>
	inline auto make_axis_list<2>()
	{
		return axis_list<Axis::X, Axis::Y>{};
	}

	template<>
	inline auto make_axis_list<3>()
	{
		return axis_list<Axis::X, Axis::Y, Axis::Z>{};
	}

	template<size_t D, typename G>
	auto make_axis_list()
	{
		return symphas::lib::types_list<G, decltype(make_axis_list<D>())>{};
	}
}


// **************************************************************************************


