
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
 * MODULE:  io
 * PURPOSE: Contains algorithms to transform datasets, such as lists
 * of vectors.
 *
 * ***************************************************************************
 */

#pragma once

#include <algorithm>
#include <vector>

#include "definitions.h"

namespace symphas::lib {
//! Computes the absolute value of a list of values.
/*!
 * Given a list of complex values, the absolute value of each element is
 * computed and then written to the same position into the output array.
 *
 * \param[in] data The array of complex type which is the input. The
 * values of this array are transformed to scalar values and copied to out
 * using the absolute value function provided by the standard library.
 * \param[out] out The list of the absolute value of the elements from data.
 * \param len Length of the array pointed to by data.
 */
inline void make_scalar_values(const complex_t* data, scalar_t* out,
                               size_t len) {
  std::transform(data, data + len, out,
                 [](auto v) { return symphas::math::abs(v); });
}

//! Computes the absolute value of a list of values.
/*!
 * Overload of make_scalar_values for an input of an array of type ::vector_t.
 * See make_scalar_values(const complex_t*, scalar_t*, size_t)
 *
 * \param[in] data The array of vector type which is the input. The
 * values of this array are transformed to scalar values and copied to out
 * using an existing absolute value function, abs(VectorValue const&).
 * \param[out] out The list of the absolute value of the elements from data.
 * \param len Length of the array pointed to by data.
 */
template <size_t D>
inline void make_scalar_values(const vector_t<D>* data, scalar_t* out,
                               size_t len) {
  auto abs_f = [](vector_t<D> const* v) { return abs(*v); };
  std::transform(data, data + len, out, abs_f);
}

//! Computes the absolute value of a list of values.
/*!
 * Overload of make_scalar_values for an input of an array of type ::scalar_t.
 * See make_scalar_values(const complex_t*, scalar_t*, size_t)
 *
 * \param[in] data The array of scalar type which is the input. The
 * values of this array are copied to out.
 * \param[out] out The list containing a copy of elements from data.
 * \param len Length of the array pointed to by data.
 */
inline void make_scalar_values(const scalar_t* data, scalar_t* out,
                               size_t len) {
  std::copy(data, data + len, out);
}

//! Returns a new list of the real values from input.
/*!
 * Generates a new list where each value is a real value from the input
 * list. Where t input list is also real, it only copies the list.
 *
 * \param[in] data The array of scalar type which is the input. The
 * values of this array are copied.
 */
template <size_t D>
inline std::vector<std::pair<axis_nd_t<D>, scalar_t>> to_scalar_field(
    std::vector<std::pair<axis_nd_t<D>, scalar_t>> const& data) {
  return data;
}

//! Returns a new list of the real values from input.
/*!
 * Generates a new list where each value is a real value from the input
 * list. Computes the real value by taking the magnitude of each complex
 * value.
 *
 * \param[in] data The array of complex type which is the input. The
 * values of this array are made into a real-valued list by taking the
 * complex magnitude.
 */
template <size_t D>
inline std::vector<std::pair<axis_nd_t<D>, scalar_t>> to_scalar_field(
    std::vector<std::pair<axis_nd_t<D>, complex_t>> const& data) {
  std::vector<std::pair<axis_nd_t<D>, scalar_t>> out;
  out.reserve(data.size());
  for (auto const& d : data) {
    out.emplace_back(d.first, std::abs(d.second));
  }
  return out;
}

}  // namespace symphas::lib