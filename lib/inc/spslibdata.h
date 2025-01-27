
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

#include <algorithm>

#include "definitions.h"
#include "types.h"

namespace symphas {
//! Defines math functions for arbitrary types.
/*!
 * Additional types are introduced by SymPhas. Therefore, mathematical
 * functions which take arbitrary types as the argument are implemented for
 * common functions including sine and cosine.
 */

//! Functionality internal to SymPhas.
/*!
 * Namespace defining all functionality which is used only internally to
 * SymPhas.
 */
namespace internal {}
}  // namespace symphas

namespace symphas::lib {

template <size_t D, typename T>
void sort_reorganize(axis_nd_t<D>* data_x, T* data_y, iter_type* sort_array,
                     len_type len) {
  for (iter_type i = 0; i < len; ++i) {
    axis_nd_t<D> x;
    x = data_x[i];
    T y = data_y[i];

    iter_type j = i;
    while (true) {
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

}  // namespace symphas::lib

namespace symphas {
//! General support functions.
/*!
 * The namespace which contains all the general support functions used
 * throughout SymPhas.
 */
namespace lib {

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
std::tuple<double*, size_t*> make_radial_arrays(const len_type (&dims)[1],
                                                double dx = 1.0);
//! See symphas::lib::make_radial_arrays(const len_type(&)[1]).
std::tuple<double*, size_t*> make_radial_arrays(const len_type (&dims)[2],
                                                double dx = 1.0,
                                                double dy = 1.0);
//! See symphas::lib::make_radial_arrays(const len_type(&)[1]).
std::tuple<double*, size_t*> make_radial_arrays(const len_type (&dims)[3],
                                                double dx = 1.0,
                                                double dy = 1.0,
                                                double dz = 1.0);

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
iter_type* make_radial_index_map(const len_type (&dims)[1]);
iter_type* make_radial_index_map(const len_type (&dims)[2]);
iter_type* make_radial_index_map(const len_type (&dims)[3]);

// **************************************************************************************

template <typename T>
void assign(T* to, iter_type ii, const T* from, iter_type n) {
  if (to != nullptr) {
    to[ii] = from[n];
  }
}

template <typename T>
void assign(T* to, iter_type ii) {
  if (to != nullptr) {
    to[ii] = T{};
  }
}

template <typename T, size_t D>
void assign(T* (&to)[D], iter_type ii, const T* (&from)[D], iter_type n) {
  if (*to != nullptr) {
    for (iter_type i = 0; i < D; ++i) {
      to[i][ii] = from[i][n];
    }
  }
}

template <typename T, size_t D>
void assign(T* (&to)[D], iter_type ii) {
  if (*to != nullptr) {
    for (iter_type i = 0; i < D; ++i) {
      to[i][ii] = T{};
    }
  }
}

inline void assign(double_arr2* to, iter_type ii, double_arr2* from,
                   iter_type n) {
  if (*to != nullptr) {
    for (iter_type i = 0; i < 2; ++i) {
      to[i][ii] = from[i][n];
    }
  }
}

inline void assign(double_arr2* to, iter_type ii) {
  if (*to != nullptr) {
    for (iter_type i = 0; i < 2; ++i) {
      to[i][ii] = double{};
    }
  }
}

template <typename T>
bool is_null(T* to) {
  return (to == nullptr);
}

template <typename T, size_t D>
bool is_null(T* (&to)[D]) {
  return (*to == nullptr);
}

inline bool is_null(double_arr2* to) { return (*to == nullptr); }

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
template <typename T, typename F>
void sort_data(std::vector<std::pair<axis_1d_type, T>>& data, len_type L, F f) {
  std::sort(data.begin(), data.end(),
            [&](auto a, auto b) { return f(a.first, b.first, L); });
}

template <typename T, typename F>
void sort_data(std::vector<std::pair<axis_2d_type, T>>& data, len_type L,
               len_type M, F f) {
  std::sort(data.begin(), data.end(), [&](auto a, auto b) {
    return (a.first[1] == b.first[1]) ? f(a.first[0], b.first[0], L)
                                      : f(a.first[1], b.first[1], M);
  });
}

template <typename T, typename F>
void sort_data(std::vector<std::pair<axis_3d_type, T>>& data, len_type L,
               len_type M, len_type N, F f) {
  std::sort(data.begin(), data.end(), [&](auto a, auto b) {
    return (a.first[2] == b.first[2] && a.first[1] == b.first[1])
               ? f(a.first[0], b.first[0], L)
           : (a.first[2] == b.first[2]) ? f(a.first[1], b.first[1], M)
                                        : f(a.first[2], b.first[2], N);
  });
}

template <typename T>
void sort_data(std::vector<std::pair<axis_1d_type, T>>& data, len_type L) {
  sort_data(data, L, regular_sort);
}
template <typename T>
void sort_data(std::vector<std::pair<axis_1d_type, T>>& data, len_type L,
               len_type M) {
  sort_data(data, L, M, regular_sort);
}
template <typename T>
void sort_data(std::vector<std::pair<axis_1d_type, T>>& data, len_type L,
               len_type M, len_type N) {
  sort_data(data, L, M, N, regular_sort);
}

template <typename T, typename F>
void sort_data(axis_1d_type* data_x, T* data_y, len_type L, F f) {
  iter_type* sort_array = new iter_type[L];
  std::iota(sort_array, sort_array + L, 0);

  std::sort(sort_array, sort_array + L,
            [&](auto a, auto b) { return f(data_x[a], data_x[b], L); });

  symphas::lib::sort_reorganize<1>(data_x, data_y, sort_array, L);
  delete[] sort_array;
}

template <typename T, typename F>
void sort_data(axis_2d_type* data_x, T* data_y, len_type L, len_type M, F f) {
  len_type len = L * M;
  iter_type* sort_array = new iter_type[len];
  std::iota(sort_array, sort_array + len, 0);

  std::sort(sort_array, sort_array + len, [&](auto a, auto b) {
    return (data_x[a][1] == data_x[b][1]) ? f(data_x[a][0], data_x[b][0], L)
                                          : f(data_x[a][1], data_x[b][1], M);
  });

  symphas::lib::sort_reorganize<2>(data_x, data_y, sort_array, len);
  delete[] sort_array;
}

template <typename T, typename F>
void sort_data(axis_3d_type* data_x, T* data_y, len_type L, len_type M,
               len_type N, F f) {
  len_type len = L * M * N;
  iter_type* sort_array = new iter_type[len];
  std::iota(sort_array, sort_array + len, 0);

  std::sort(sort_array, sort_array + len, [&](auto a, auto b) {
    return (data_x[a][2] == data_x[b][2] && data_x[a][1] == data_x[b][1])
               ? f(data_x[a][0], data_x[b][0], L)
           : (data_x[a][2] == data_x[b][2]) ? f(data_x[a][1], data_x[b][1], M)
                                            : f(data_x[a][2], data_x[b][2], N);
  });

  symphas::lib::sort_reorganize<3>(data_x, data_y, sort_array, len);
  delete[] sort_array;
}

template <typename X, typename T>
void sort_data(X&& data_x, T&& data_y, len_type L) {
  sort_data(data_x, data_y, L, regular_sort);
}
template <typename X, typename T>
void sort_data(X&& data_x, T&& data_y, len_type L, len_type M) {
  sort_data(data_x, data_y, L, M, regular_sort);
}
template <typename X, typename T>
void sort_data(X&& data_x, T&& data_y, len_type L, len_type M, len_type N) {
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
template <typename X, typename Y>
std::vector<std::pair<X, Y>> combine_data(const X* data_x, const Y* data_y,
                                          len_type len) {
  std::vector<std::pair<X, Y>> out;
  out.reserve(len);

  for (iter_type i = 0; i < len; ++i) {
    out.emplace_back(data_x[i], data_y[i]);
  }
  return out;
}

template <typename X, typename Y, size_t N>
std::vector<std::pair<X, Y[N]>> combine_data(const X* data_x,
                                             const Y (*data_y)[N],
                                             len_type len) {
  std::vector<std::pair<X, Y[N]>> out(len);

  for (iter_type i = 0; i < len; ++i) {
    out[i].first = data_x[i];
    for (iter_type n = 0; n < N; ++n) {
      out[i].second[n] = data_y[i][n];
    }
  }
  return out;
}

}  // namespace lib
}  // namespace symphas

// **************************************************************************************
