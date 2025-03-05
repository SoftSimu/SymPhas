
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
#include <filesystem>

#include "spslibdata.h"
#include "spslibiter.h"
#include "spsliblist.h"
#include "spslibmath.h"
#include "spslibstr.h"
#include "spslibtuple.h"
#include "timer.h"
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

/*! \addtogroup grid
 * @{
 */

//! Defines functions providing information about a grid.
namespace grid {
template <size_t D>
__host__ __device__ void set_dimensions(const len_type* dimensions,
                                        len_type* dims) {
  if (dimensions == nullptr) {
    std::fill(dims, dims + D, 0);
  } else {
    for (iter_type i = 0; i < D; ++i) {
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
template <size_t D>
__host__ __device__ len_type length(len_type const* dimensions);

//! Specialization of grid::length(len_type const*).
template <>
__host__ __device__ inline len_type length<1>(len_type const* dimensions) {
  return ((dimensions != nullptr && dimensions[0] > 0) ? dimensions[0] : 0);
}

template <size_t D>
__host__ __device__ len_type length(len_type const* dimensions) {
  return ((dimensions != nullptr && dimensions[D - 1] > 0) ? dimensions[D - 1]
                                                           : 0) *
         length<D - 1>(dimensions);
}

__host__ __device__ inline len_type length(len_type const* dimensions,
                                           size_t dimension) {
  len_type len = 1;
  for (iter_type i = 0; i < dimension; ++i) {
    len *= (dimensions[i] > 0) ? dimensions[i] : 0;
  }
  return len;
}

}  // namespace grid

/*!
 * @}
 */

namespace symphas {
namespace lib {

namespace internal {

template <typename T>
struct CrossProductFunctions {
  template <T... Es>
  using seq_t = std::integer_sequence<T, Es...>;

  // **********************************************************
  // Expand the cross list into a full tuple of sequences representing the row
  // of combined values.

  template <T E1>
  static auto constexpr expand2(seq_t<E1>, seq_t<>) {
    return std::make_tuple();
  }

  template <T E1, T E2, T... E2s>
  static auto constexpr expand2(seq_t<E1>, seq_t<E2, E2s...>) {
    return std::tuple_cat(std::make_tuple(seq_t<E1, E2>{}),
                          expand2(seq_t<E1>{}, seq_t<E2s...>{}));
  }

  template <T E1, T... E1s>
  static auto constexpr expand1(seq_t<E1, E1s...>, seq_t<>) {
    return std::tuple<>{};
  }

  template <T E2, T... E2s>
  static auto constexpr expand1(seq_t<>, seq_t<E2, E2s...>) {
    return std::tuple<>{};
  }

  static auto constexpr expand1(seq_t<>, seq_t<>) { return std::tuple<>{}; }

  template <T E1, T... E1s, T E2, T... E2s>
  static auto constexpr expand1(seq_t<E1, E1s...>, seq_t<E2, E2s...>) {
    return std::tuple_cat(expand2(seq_t<E1>{}, seq_t<E2, E2s...>{}),
                          expand1(seq_t<E1s...>{}, seq_t<E2, E2s...>{}));
  }

  template <T... E1s>
  static auto constexpr expand33(seq_t<E1s...>, seq_t<>) {
    return std::make_tuple();
  }

  template <T... E1s, T E2, T... E2s>
  static auto constexpr expand33(seq_t<E1s...>, seq_t<E2, E2s...>) {
    return std::tuple_cat(
        std::make_tuple(seq_join(seq_t<E1s...>{}, seq_t<E2>{})),
        expand33(seq_t<E1s...>{}, seq_t<E2s...>{}));
  }

  template <T E, T... Es>
  static auto constexpr expand22(std::tuple<>, seq_t<E, Es...>) {
    return std::tuple<>{};
  }

  template <typename Row, typename... Rows, T E, T... Es>
  static auto constexpr expand22(std::tuple<Row, Rows...>, seq_t<E, Es...>) {
    return std::tuple_cat(expand33(Row{}, seq_t<E, Es...>{}),
                          expand22(std::tuple<Rows...>{}, seq_t<E, Es...>{}));
  }

  template <typename Row, typename... Rows, T E, T... Es, typename List0,
            typename... Lists>
  static auto constexpr expand22(std::tuple<Row, Rows...>, seq_t<E, Es...>,
                                 List0, Lists...) {
    return expand22(
        std::tuple_cat(expand33(Row{}, seq_t<E, Es...>{}),
                       expand22(std::tuple<Rows...>{}, seq_t<E, Es...>{})),
        List0{}, Lists{}...);
  }

  template <T E, T... Es, typename List0, typename... Lists>
  static auto constexpr expand11(seq_t<E, Es...>, List0, Lists...) {
    return expand22(expand1(seq_t<E, Es...>{}, List0{}), Lists{}...);
  }

  // **********************************************************
  // Selects only a single row without constructing the whole cross list.

  template <size_t N>
  static auto constexpr select(seq_t<>) {
    return std::integer_sequence<T>{};
  }

  template <size_t N, T E, T... Es>
  static auto constexpr select(seq_t<E, Es...>) {
    return std::integer_sequence<T, seq_value<N>(seq_t<E, Es...>{})>{};
  }

  template <size_t N, T E, T... Es, typename Seq, typename... Seqs,
            size_t L = seq_len_product<Seq, Seqs...>::value>
  static auto constexpr select(seq_t<E, Es...>, Seq, Seqs...) {
    if constexpr (L == 0) {
      return select_non_empty_seq<N>(
          types_list<>{}, types_list<seq_t<E, Es...>, Seq, Seqs...>{});
    } else if constexpr (N < L) {
      constexpr size_t N0 = N / L;
      constexpr size_t N1 = N - N0 * L;
      return seq_join(
          std::integer_sequence<T, seq_value<N0>(seq_t<E, Es...>{})>{},
          select<N1>(Seq{}, Seqs{}...));
    } else {
      return select<N % L>(seq_t<E, Es...>{}, Seq{}, Seqs{}...);
    }
  }

  template <size_t N, typename... Seqs0, T I0, T... Is, typename... Seqs>
  static auto constexpr select_non_empty_seq(
      types_list<Seqs0...>,
      types_list<std::integer_sequence<T, I0, Is...>, Seqs...>) {
    return select_non_empty_seq<N>(
        types_list<Seqs0..., std::integer_sequence<T, I0, Is...>>{},
        types_list<Seqs...>{});
  }

  template <size_t N, typename... Seqs0, typename... Seqs>
  static auto constexpr select_non_empty_seq(
      types_list<Seqs0...>, types_list<std::integer_sequence<T>, Seqs...>) {
    return select_non_empty_seq<N>(types_list<Seqs0...>{},
                                   types_list<Seqs...>{});
  }

  template <size_t N, typename... Seqs0>
  static auto constexpr select_non_empty_seq(types_list<Seqs0...>,
                                             types_list<>) {
    return select<N>(Seqs0{}...);
  }
};

}  // namespace internal

/*!
 * \brief Generate the cross product or cross join of all the numeric elements
 * in the provided std::integer_sequence types.
 *
 * This struct template is used to generate the cross product of multiple
 * std::integer_sequence types. It provides the count of the total number of
 * combinations, the rank (number of lists in the cross product), and a `row`
 * type alias that represents a specific combination.
 *
 * \tparam Lists The std::integer_sequence containing a list of values which are
 * cross joined.
 */
template <typename... Lists>
struct CrossProductList;

template <>
struct CrossProductList<> {
  static const size_t count = 0;
  static const size_t rank = 0;
};

template <typename T, T... Es>
struct CrossProductList<std::integer_sequence<T, Es...>> {
  static const size_t count = sizeof...(Es);
  static const size_t rank = 1;

  template <size_t N>
  static const size_t size = std::integer_sequence<T, Es...>::size();

  template <size_t N, typename std::enable_if_t<(N < count), int> = 0>
  using row = type_at_index<
      N, unroll_types_list<types_list<std::integer_sequence<T, Es>...>>>;
};

template <typename T, T... E1s, T... E2s>
struct CrossProductList<std::integer_sequence<T, E1s...>,
                        std::integer_sequence<T, E2s...>> {
  static const size_t count = (sizeof...(E1s) * sizeof...(E2s));
  static const size_t rank = 2;

  template <size_t N>
  static const size_t size = type_at_index<
      N,
      unroll_types_list<types_list<std::integer_sequence<T, E1s...>,
                                   std::integer_sequence<T, E2s...>>>>::size();

  template <size_t N>
  using row = decltype(internal::CrossProductFunctions<T>::template select<N>(
      std::declval<std::integer_sequence<T, E1s...>>(),
      std::declval<std::integer_sequence<T, E2s...>>()));
};

template <typename T, T... Es, typename List1, typename List2,
          typename... Lists>
struct CrossProductList<std::integer_sequence<T, Es...>, List1, List2,
                        Lists...> {
  static const size_t count = seq_len_product<std::integer_sequence<T, Es...>,
                                              List1, List2, Lists...>::value;
  static const size_t rank = 3 + sizeof...(Lists);

  template <size_t N>
  static const size_t size = type_at_index<
      N, unroll_types_list<types_list<std::integer_sequence<T, Es...>, List1,
                                      List2, Lists...>>>::size();

  template <size_t N>
  using row = decltype(internal::CrossProductFunctions<T>::template select<N>(
      std::declval<std::integer_sequence<T, Es...>>(), std::declval<List1>(),
      std::declval<List2>(), std::declval<Lists>()...));
};
}  // namespace lib
}  // namespace symphas
// *************************************************************************************************

namespace symphas {
//! General support functions.
/*!
 * The namespace which contains all the general support functions used
 * throughout SymPhas.
 */
namespace lib {

//! Determine the length of the coordinate relative to the origin.
/*!
 * Determine the length of the coordinate relative to the origin.
 *
 * \param p The axis point in a 1-dimensional grid.
 */
inline axis_coord_t length(axis_1d_type p) { return std::abs(p); }

//! Determine the length of the coordinate relative to the origin.
/*!
 * Determine the length of the coordinate relative to the origin.
 *
 * \param p The axis point in a 2-dimensional grid.
 */
inline axis_coord_t length(axis_2d_type p) {
  return std::sqrt(p[0] * p[0] + p[1] * p[1]);
}

//! Determine the length of the coordinate relative to the origin.
/*!
 * Determine the length of the coordinate relative to the origin.
 *
 * \param p The axis point in a 3-dimensional grid.
 */
inline axis_coord_t length(axis_3d_type p) {
  return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
}

//! Determine the distance between two axis points.
/*!
 * Determine the distance between two axis points.
 *
 * \param p0 The first point in a 1-dimensional grid.
 * \param p1 The second point in a 1-dimensional grid.
 */
inline axis_coord_t distance(axis_1d_type p0, axis_1d_type p1) {
  return length(axis_1d_type{p0 - p1});
}

//! Determine the distance between two axis points.
/*!
 * Determine the distance between two axis points.
 *
 * \param p0 The first point in a 2-dimensional grid.
 * \param p1 The second point in a 2-dimensional grid.
 */
inline axis_coord_t distance(axis_2d_type p0, axis_2d_type p1) {
  return length(axis_2d_type{p0[0] - p1[0], p0[1] - p1[1]});
}

//! Determine the distance between two axis points.
/*!
 * Determine the distance between two axis points.
 *
 * \param p0 The first point in a 3-dimensional grid.
 * \param p1 The second point in a 3-dimensional grid.
 */
inline axis_coord_t distance(axis_3d_type p0, axis_3d_type p1) {
  return length(axis_3d_type{p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2]});
}

//! Obtain the number of elements in a list of unique distances.
/*!
 * For a system of the prescribed dimensions, return the number of
 * unique distances from the origin as computed by the Euclidean norm.
 *
 * \param dims The number of elements along each axis.
 */
len_type radial_length(const len_type (&dims)[1]);

//! See symphas::lib::radial_length().
len_type radial_length(const len_type (&dims)[2]);

//! See symphas::lib::radial_length().
len_type radial_length(const len_type (&dims)[3]);

//! Print the timestamp into the given buffer string.
/*!
 * Print the timestamp into the given buffer string.
 */
void write_ts_str(char* buffer);

void make_directory(std::filesystem::path dir, int err_no = ERR_CODE_FILE_OPEN);

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
 * given path. The returned pointer must be deleted to avoid a memory leak.
 */
char* get_parent_directory(const char* path);

std::filesystem::path get_parent_directory(std::filesystem::path dir);

template <Axis... axs>
struct axis_list {};

template <Side... sides>
struct side_list {};

//! Put the first `D` axes in a types list.
/*!
 * Construct a types list with the first `D` axes in the list.
 */
template <size_t D>
auto make_axis_list();

template <>
inline auto make_axis_list<1>() {
  return axis_list<Axis::X>{};
}

template <>
inline auto make_axis_list<2>() {
  return axis_list<Axis::X, Axis::Y>{};
}

template <>
inline auto make_axis_list<3>() {
  return axis_list<Axis::X, Axis::Y, Axis::Z>{};
}

template <size_t D, typename G>
auto make_axis_list() {
  return symphas::lib::types_list<G, decltype(make_axis_list<D>())>{};
}
}  // namespace lib
}  // namespace symphas

// **************************************************************************************
