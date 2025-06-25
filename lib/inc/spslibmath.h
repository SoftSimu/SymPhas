
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
#include <libgen.h>
#include <sys/stat.h>
#endif

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <numeric>
#include <vector>

// #include "symphasthread.h"
#include "definitions.h"
#include "timer.h"

namespace symphas {
//! Defines math functions for arbitrary types.
/*!
 * Additional types are introduced by SymPhas. Therefore, mathematical
 * functions which take arbitrary types as the argument are implemented for
 * common functions including sine and cosine.
 */
namespace math {}
}  // namespace symphas

//! Generic template for returning the result of the math function.
/*!
 * For a value of generic type which is not accepted by any overloads,
 * it is explicitly cast to the type `complex_t` to compute the result.
 */
#define MATH_FUNCTION_OVERLOADS_IMPL(NAME, EVAL)   \
  inline auto _##NAME(int vv) {                    \
    auto v = static_cast<scalar_t>(vv);            \
    return EVAL;                                   \
  }                                                \
  inline auto _##NAME(scalar_t v) { return EVAL; } \
  inline auto _##NAME(complex_t c) {               \
    std::complex<scalar_t> v(c.real(), c.imag());  \
    return EVAL;                                   \
  }                                                \
  template <typename T>                            \
  auto _##NAME(T v) {                              \
    return _##NAME(static_cast<complex_t>(v));     \
  }                                                \
  template <typename T>                            \
  auto NAME(T v) {                                 \
    return _##NAME(v);                             \
  }

#define MATH_FUNCTION_OVERLOADS(NAME, FUNC) \
  MATH_FUNCTION_OVERLOADS_IMPL(NAME, FUNC(v))

namespace symphas::math {

//! Returns conjugate complex number for a generic type.
/*!
 * Since additional types are introduced in the SymPhas library, an
 * extension of the standard library `conj` function is implemented to
 * apply to additional types.
 *
 * This implementation simply forwards its argument to the standard library
 * `conj` function.
 */
template <typename T>
complex_t conj(T const& v) {
  std::complex<scalar_t> c = std::conj(v);
  return complex_t(c.real(), c.imag());
}

//! Returns conjugate complex number for a internal complex type.
/*!
 * Since additional types are introduced in the SymPhas library, an
 * extension of the standard library `conj` function is implemented to
 * apply to complex_t type.
 */
template<>
inline complex_t conj<>(complex_t const& v) {
  return complex_t(v.real(), -v.imag());
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
template <typename T, size_t D>
auto conj(any_vector_t<T, D> const& v) {
  using std::conj;

  any_vector_t<T, D> vc;
  for (iter_type i = 0; i < D; ++i) {
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
template <typename T>
scalar_t modulus(T const& v) {
  return std::abs(v);
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
template <>
inline scalar_t modulus<>(complex_t const& v) {
  return std::sqrt(v.real() * v.real() + v.imag() * v.imag());
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
template <typename T>
scalar_t real(T const& v) {
  using namespace std;  // standard namespace brought into scope
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
template <typename T>
scalar_t imag(T const& v) {
  using namespace std;  // standard namespace brought into scope
  return imag(v);
}

template <typename T, size_t D>
__host__ __device__ T abs(const T (&value)[D]) {
  T result{};
  for (iter_type i = 0; i < D; ++i) {
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
MATH_FUNCTION_OVERLOADS(abs, std::abs);

MATH_FUNCTION_OVERLOADS_IMPL(sec, 1. / std::cos(v));
MATH_FUNCTION_OVERLOADS_IMPL(csc, 1. / std::sin(v));
MATH_FUNCTION_OVERLOADS_IMPL(cot, 1. / std::tan(v));

template <size_t N, size_t E>
constexpr size_t fixed_pow = N * fixed_pow<N, E - 1>;
template <size_t N>
constexpr size_t fixed_pow<N, 0> = 1;
template <size_t N>
constexpr size_t fixed_pow<N, 1> = N;

template <typename T, typename E>
__host__ __device__ auto pow(T&& base, E&& exponent) {
  using std::pow;
  return pow(std::forward<T>(base), std::forward<E>(exponent));
}

template <typename T, size_t... Is>
__host__ __device__ auto pow(T const& base, std::index_sequence<Is...>) {
  auto f = [&](size_t) { return base; };
  return (f(Is) * ...);
}

template <size_t N, typename T>
__host__ __device__ auto pow(T const& base) {
  if constexpr (N <= 6) {
    return pow(base, std::make_index_sequence<N>{});
  } else {
    // using symphas::math::pow;
    return pow(base, N);
  }
}

template <size_t O, typename T, size_t D>
__host__ __device__ auto pow(any_vector_t<T, D> const& v) {
  using std::pow;
  using symphas::math::pow;

  if constexpr (O == 1) {
    return v;
  } else if constexpr (O == 2) {
    return dot(v, v);
  } else if constexpr (O % 2 == 0) {
    return pow<O / 2>(dot(v, v));
  } else {
    return pow<O / 2>(dot(v, v)) * v;
  }
}

template <typename T, size_t D>
__host__ __device__ auto pow(any_vector_t<T, D> const& v, double exponent) {
  throw;
}

//! Cross product function for vectors.
/*!
 * Compound assignment operator to calculate the conventional cross product and
 * assign it to the calling instance.
 *
 * Assignment operator to compute the cross product, which is only applicable to
 * the 3-dimensional VectorValue template specialization. Operator is called by
 * the left hand instance, taking data from the right hand VectorValue to
 * compute the cross product and assigning the result to the left hand side.
 *
 * \param rhs The VectorValue instance where data is taken from to compute the
 * cross product the data from the left hand instance.
 */
template <typename T>
any_vector_t<T, 3> cross(any_vector_t<T, 3> const& lhs,
                         any_vector_t<T, 3> const& rhs) {
  any_vector_t<T, 3> result;

  result[0] = (lhs[1] * rhs[2] - lhs[2] * rhs[1]);
  result[1] = (lhs[2] * rhs[0] - lhs[0] * rhs[2]);
  result[2] = (lhs[0] * rhs[1] - lhs[1] * rhs[0]);

  return result;
}

// template<size_t D>
// auto dot(any_vector_t<scalar_t, D> const& lhs, any_vector_t<complex_t, D>
// const& rhs)
//{
//	complex_t sum = 0;
//	for (iter_type i = 0; i < D; ++i)
//	{
//		sum += lhs[i] * rhs[i];
//	}
//	return sum;
// }

// template<size_t D>
// auto dot(any_vector_t<complex_t, D> const& lhs, any_vector_t<scalar_t, D>
// const& rhs)
//{
//	complex_t sum = 0;
//	for (iter_type i = 0; i < D; ++i)
//	{
//		sum += lhs[i] * rhs[i];
//	}
//	return sum;
// }

template <size_t... Os>
constexpr auto sum() {
  return (Os + ...);
}

template <size_t... Os>
constexpr auto product() {
  return (Os * ...);
}
}  // namespace symphas::math

namespace symphas::internal {

template <typename T>
struct check_is_simple_data {
 protected:
  constexpr static std::true_type _get_value(complex_t) { return {}; }

  constexpr static std::true_type _get_value(scalar_t) { return {}; }

  constexpr static std::true_type _get_value(int) { return {}; }

  template <typename T0>
  constexpr static std::false_type _get_value(T0) {
    return {};
  }

  constexpr static auto get_value(T ptr) { return _get_value(ptr); }

 public:
  static const bool value =
      std::invoke_result_t<decltype(&check_is_simple_data<T>::get_value),
                           T>::value;
};

template <typename T>
struct check_is_simple_data<T&> : check_is_simple_data<T> {};
template <typename T>
struct check_is_simple_data<T const> : check_is_simple_data<T> {};
}  // namespace symphas::internal

namespace symphas {
template <typename T>
constexpr bool is_simple_data =
    symphas::internal::check_is_simple_data<T>::value;

template <typename T>
constexpr bool is_non_vector = true;
template <typename T, size_t D>
constexpr bool is_non_vector<any_vector_t<T, D>> = false;
}  // namespace symphas

template <typename T, size_t N, size_t M>
using any_matrix_t = any_vector_t<any_vector_t<T, M>, N>;
template <typename T, size_t N>
using any_row_vector_t = any_matrix_t<T, 1, N>;

template <typename T, size_t D,
          typename std::enable_if_t<symphas::is_non_vector<T>, int> = 0>
__host__ __device__ auto operator*(any_vector_t<T, D> const& lhs,
                                   any_vector_t<T, D> const& rhs) {
  using namespace std;
  using namespace symphas::math;
  return dot(lhs, rhs);
}

template <typename T, typename S, size_t D,
          typename std::enable_if_t<
              symphas::is_non_vector<T> && symphas::is_non_vector<S>, int> = 0>
__host__ __device__ auto operator*(any_vector_t<T, D> const& lhs,
                                   any_vector_t<S, D> const& rhs) {
  using namespace std;
  using namespace symphas::math;
  return dot(lhs, rhs);
}

template <
    typename T, typename S, size_t D,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator*(any_row_vector_t<T, D> const& lhs,
                                   any_vector_t<S, D> const& rhs) {
  any_vector_t<mul_result_t<T, S>, D> lhs0;
  for (iter_type i = 0; i < D; ++i) {
    lhs0[i] = lhs[0][i];
  }

  using namespace std;
  using namespace symphas::math;
  return dot(lhs0, rhs);
}

template <
    typename T, typename S, size_t D,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator*(any_vector_t<T, D> const& lhs,
                                   any_row_vector_t<S, D> const& rhs) {
  any_matrix_t<mul_result_t<T, S>, D, D> lhs0;
  for (iter_type i = 0; i < D; ++i) {
    for (iter_type j = 0; j < D; ++j) {
      lhs0[i][j] = lhs[i] * rhs[0][j];
    }
  }
  return lhs0;
}

template <
    typename T, typename S, size_t N,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_simple_data<S>), int> = 0>
__host__ __device__ auto operator*(any_vector_t<T, N> const& lhs,
                                   S const& rhs) {
  any_vector_t<mul_result_t<T, S>, N> result;
  for (int i = 0; i < N; ++i) {
    result[i] = lhs[i] * rhs;
  }
  return result;
}

template <typename T, size_t N,
          typename std::enable_if_t<symphas::is_non_vector<T>, int> = 0>
__host__ __device__ auto operator*(any_vector_t<T, N> const& lhs,
                                   scalar_t const& rhs) {
  any_vector_t<mul_result_t<T, scalar_t>, N> result;
  for (int i = 0; i < N; ++i) {
    result[i] = lhs[i] * rhs;
  }
  return result;
}

template <
    typename T, typename S, size_t N,
    typename std::enable_if_t<
        (symphas::is_simple_data<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator*(T const& lhs,
                                   any_vector_t<S, N> const& rhs) {
  any_vector_t<mul_result_t<T, S>, N> result;
  for (int i = 0; i < N; ++i) {
    result[i] = lhs * rhs[i];
  }
  return result;
}

template <
    typename T, typename S, size_t N,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_simple_data<S>), int> = 0>
__host__ __device__ auto operator*(any_row_vector_t<T, N> const& lhs,
                                   S const& rhs) {
  any_row_vector_t<mul_result_t<T, S>, N> result;
  for (int i = 0; i < N; ++i) {
    result[0][i] = lhs[0][i] * rhs;
  }
  return result;
}

template <
    typename T, typename S, size_t N,
    typename std::enable_if_t<
        (symphas::is_simple_data<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator*(T const& lhs,
                                   any_row_vector_t<S, N> const& rhs) {
  any_row_vector_t<mul_result_t<T, S>, N> result;
  for (int i = 0; i < N; ++i) {
    result[0][i] = lhs * rhs[0][i];
  }
  return result;
}

template <
    typename T, typename S, size_t N, size_t M,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_simple_data<S>), int> = 0>
__host__ __device__ auto operator*(any_matrix_t<T, N, M> const& lhs,
                                   S const& rhs) {
  any_matrix_t<mul_result_t<T, S>, N, M> result;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      result[i][i] = lhs[i][j] * rhs;
    }
  }
  return result;
}

template <
    typename T, typename S, size_t N, size_t M,
    typename std::enable_if_t<
        (symphas::is_simple_data<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator*(T const& lhs,
                                   any_matrix_t<S, N, M> const& rhs) {
  any_matrix_t<mul_result_t<T, S>, N, M> result;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      result[i][i] = lhs * rhs[i][j];
    }
  }
  return result;
}

template <
    typename T, typename S, size_t L, size_t M, size_t N,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator*(any_matrix_t<T, L, M> const& lhs,
                                   any_matrix_t<S, M, N> const& rhs) {
  any_matrix_t<mul_result_t<T, S>, L, N> result;
  for (int i = 0; i < L; ++i) {
    for (int j = 0; j < N; ++j) {
      mul_result_t<T, S> sum{};
      for (int k = 0; k < M; ++k) {
        sum += lhs[i][k] * rhs[k][j];
      }
      result[i][j] = sum;
    }
  }
  return result;
}

template <
    typename T, typename S, size_t L, size_t M,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator*(any_matrix_t<T, L, M> const& lhs,
                                   any_vector_t<S, M> const& rhs) {
  any_vector_t<mul_result_t<T, S>, M> result;
  for (int i = 0; i < L; ++i) {
    mul_result_t<T, S> sum{};
    for (int k = 0; k < M; ++k) {
      sum += lhs[i][k] * rhs[k];
    }
    result[i] = sum;
  }
  return result;
}

template <
    typename T, typename S, size_t L, size_t M,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator*(any_row_vector_t<T, L> const& lhs,
                                   any_matrix_t<S, L, M> const& rhs) {
  any_row_vector_t<mul_result_t<T, S>, M> result;
  for (int i = 0; i < M; ++i) {
    mul_result_t<T, S> sum{};
    for (int k = 0; k < L; ++k) {
      sum += lhs[k] * rhs[k][i];
    }
    result[i] = sum;
  }
  return result;
}

template <
    typename T, typename S, size_t N,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_simple_data<S>), int> = 0>
__host__ __device__ auto operator/(any_vector_t<T, N> const& lhs,
                                   S const& rhs) {
  any_vector_t<div_result_t<T, S>, N> result;
  for (int i = 0; i < N; ++i) {
    result[i] = lhs[i] / rhs;
  }
  return result;
}

template <
    typename T, typename S, size_t N,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator*(any_vector_t<T, 1> const& lhs,
                                   any_vector_t<S, N> const& rhs) {
  return lhs.v[0] * rhs;
}

template <
    typename T, typename S, size_t N,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator*(any_vector_t<T, N> const& lhs,
                                   any_vector_t<S, 1> const& rhs) {
  return lhs * rhs.v[0];
}

template <
    typename T, typename S,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator*(any_vector_t<T, 1> const& lhs,
                                   any_vector_t<S, 1> const& rhs) {
  return lhs.v[0] * rhs.v[0];
}

template <typename T,
          typename std::enable_if_t<(symphas::is_non_vector<T>), int> = 0>
__host__ __device__ auto operator*(any_vector_t<T, 1> const& lhs,
                                   any_vector_t<T, 1> const& rhs) {
  return lhs.v[0] * rhs.v[0];
}

template <
    typename T, typename S,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator+(any_vector_t<T, 1> const& lhs,
                                   S const& rhs) {
  any_vector_t<add_result_t<T, S>, 1> result;
  result[0] = lhs[0] + rhs;
  return result;
}

template <
    typename T, typename S,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator+(T const& lhs,
                                   any_vector_t<S, 1> const& rhs) {
  any_vector_t<add_result_t<T, S>, 1> result;
  result[0] = lhs + rhs[0];
  return result;
}

template <
    typename T, typename S,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator+(any_row_vector_t<T, 1> const& lhs,
                                   S const& rhs) {
  any_row_vector_t<add_result_t<T, S>, 1> result;
  result[0][0] = lhs[0][0] + rhs;
  return result;
}

template <
    typename T, typename S,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator+(T const& lhs,
                                   any_row_vector_t<S, 1> const& rhs) {
  any_row_vector_t<add_result_t<T, S>, 1> result;
  result[0][0] = lhs + rhs[0][0];
  return result;
}

template <
    typename T, typename S,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator-(any_vector_t<T, 1> const& lhs,
                                   S const& rhs) {
  any_vector_t<add_result_t<T, S>, 1> result;
  result[0] = lhs[0] - rhs;
  return result;
}

template <
    typename T, typename S,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator-(T const& lhs,
                                   any_vector_t<S, 1> const& rhs) {
  any_vector_t<add_result_t<T, S>, 1> result;
  result[0] = lhs - rhs[0];
  return result;
}

template <
    typename T, typename S,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator-(any_row_vector_t<T, 1> const& lhs,
                                   S const& rhs) {
  any_row_vector_t<add_result_t<T, S>, 1> result;
  result[0][0] = lhs[0][0] - rhs;
  return result;
}

template <
    typename T, typename S,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator-(T const& lhs,
                                   any_row_vector_t<S, 1> const& rhs) {
  any_row_vector_t<add_result_t<T, S>, 1> result;
  result[0][0] = lhs - rhs[0][0];
  return result;
}

template <
    typename T, typename S, size_t N, size_t M,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator+(any_matrix_t<T, N, M> const& lhs,
                                   any_matrix_t<S, N, M> const& rhs) {
  any_matrix_t<add_result_t<T, S>, N, M> result;
  for (int i = 0; i < N; ++i) {
    result[0][i] = lhs[0][i] + rhs[0][i];
  }
  return result;
}

template <
    typename T, typename S, size_t N,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator+(any_vector_t<T, N> const& lhs,
                                   any_vector_t<S, N> const& rhs) {
  any_vector_t<add_result_t<T, S>, N> result;
  for (int i = 0; i < N; ++i) {
    result[i] = lhs[i] + rhs[i];
  }
  return result;
}

template <
    typename T, typename S, size_t N, size_t M,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator-(any_matrix_t<T, N, M> const& lhs,
                                   any_matrix_t<S, N, M> const& rhs) {
  any_matrix_t<sub_result_t<T, S>, N, M> result;
  for (int i = 0; i < N; ++i) {
    result[0][i] = lhs[0][i] - rhs[0][i];
  }
  return result;
}

template <
    typename T, typename S, size_t N,
    typename std::enable_if_t<
        (symphas::is_non_vector<T> && symphas::is_non_vector<S>), int> = 0>
__host__ __device__ auto operator-(any_vector_t<T, N> const& lhs,
                                   any_vector_t<S, N> const& rhs) {
  any_vector_t<sub_result_t<T, S>, N> result;
  for (int i = 0; i < N; ++i) {
    result[i] = lhs[i] - rhs[i];
  }
  return result;
}

inline axis_nd_t<2> operator+(axis_nd_t<2> const& a, axis_nd_t<2> const& b) {
  return {a[0] + b[0], a[1] + b[1]};
}

inline axis_nd_t<2> operator-(axis_nd_t<2> const& a, axis_nd_t<2> const& b) {
  return {a[0] - b[0], a[1] - b[1]};
}

inline axis_nd_t<3> operator+(axis_nd_t<3> const& a, axis_nd_t<3> const& b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

inline axis_nd_t<3> operator-(axis_nd_t<3> const& a, axis_nd_t<3> const& b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

namespace symphas {
namespace lib {

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
template <typename T0>
constexpr auto max_value(T0 val0) {
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
template <typename T0, typename T1, typename... Ts>
constexpr auto max_value(T0 val0, T1 val1, Ts... vals)
    -> add_result_t<T0, T1, Ts...> {
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
template <typename T0>
constexpr auto min_value(T0 val0) {
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
template <typename T0, typename T1, typename... Ts>
constexpr auto min_value(T0 val0, T1 val1, Ts... vals)
    -> add_result_t<T0, T1, Ts...> {
  auto amin = min_value(val1, vals...);
  return (val0 < amin) ? val0 : amin;
}

//! Implements function to return identity.
/*!
 * Implements a `static constexpr` function in order to return the known
 * identity at compile time. Each specialization of this struct
 * implements the identity for its respective type.
 */
template <typename T>
struct identity;

//! Identity for the scalar type, always equal to `1.0`.
template <>
struct identity<scalar_t> {
  __host__ __device__ constexpr scalar_t operator()() { return 1.0; }
};

//! Identity for the complex type, always equal to `1.0 + 0.0*i`.
template <>
struct identity<complex_t> {
  __host__ __device__ constexpr complex_t operator()() { return {1.0, 0.0}; }
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
template <typename T, size_t D>
struct identity<VectorValue<T, D>> {
  __host__ __device__ constexpr VectorValue<T, D> operator()() {
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
template <size_t D>
struct identity<vector_t<D>> {
  __host__ __device__ constexpr vector_t<D> operator()() {
    vector_t<D> idty;
    for (iter_type i = 0; i < D; ++i) {
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
template <typename T>
__host__ __device__ inline auto constexpr get_identity() {
  return identity<T>{}();
}

//! Return the number of digits in a positive integer.
/*!
 * Count the number of digits in a positive integer. Returns the
 * value as a compile time constant.
 */
template <size_t N>
size_t constexpr num_digits() {
  if (N <= 9) {
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
constexpr size_t power_ull(size_t b, size_t p) {
  size_t e = 1;
  for (iter_type i = 0; i < p; ++i) {
    e *= b;
  }
  return e;
}

inline int exponent10(double value) {
  double p = std::log10(std::abs(value));
  return (p < 0) ? int(p) - 1 : int(p);
}

inline double base10(double value) {
  return value * std::pow(10, -exponent10(value));
}

}  // namespace lib
}  // namespace symphas
