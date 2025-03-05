
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

#include <type_traits>

namespace symphas {

/*!
 * \brief A template struct that holds an array of values of type T with a
 * specified dimension.
 *
 * The `array_values_type` struct is a utility for managing arrays of a fixed
 * size. It provides constructors for initializing the array and implicit
 * conversion operators to treat the struct as a reference to the underlying
 * array.
 *
 * \tparam T The type of the elements in the array.
 * \tparam D The dimension (size) of the array.
 */
template <typename T, size_t D>
struct array_values_type {
  using arr_t = T[D];

  arr_t range;
  template <typename... value_types>
  __device__ __host__ array_values_type(value_types... values)
      : range{values...} {};

  __device__ __host__ operator const arr_t&() const { return range; }
  __device__ __host__ operator arr_t&() { return range; }

  __device__ __host__ const arr_t& as_array() const { return range; }
  __device__ __host__ arr_t& as_array() { return range; }

  __device__ __host__ const T& operator[](iter_type i) const {
    return range[i];
  }
  __device__ __host__ T& operator[](iter_type i) { return range[i]; }
};

template <typename T, size_t D0, size_t D1>
struct array_2_values_type {
  using arr2_t = T[D0][D1];

  arr2_t range;
  template <typename... value_types>
  __device__ __host__ array_2_values_type(value_types&&... values) : range{} {
    initialize(std::make_index_sequence<sizeof...(value_types)>{},
               std::forward<value_types>(values)...);
  };

  __device__ __host__ operator const arr2_t&() const { return range; }
  __device__ __host__ operator arr2_t&() { return range; }

  __device__ __host__ const arr2_t& as_array() const { return range; }
  __device__ __host__ arr2_t& as_array() { return range; }

 private:
  template <size_t I0, size_t... I0s>
  void fill(array_values_type<T, D1> const& value,
            std::index_sequence<I0s...>) {
    ((range[I0][I0s] = value[I0s]), ...);
  }

  template <size_t I0>
  void fill(array_values_type<T, D1> const& value) {
    fill<I0>(value, std::make_index_sequence<D1>{});
  }

  template <typename... value_types, size_t... I0s>
  void initialize(std::index_sequence<I0s...>, value_types&&... values) {
    ((fill<I0s>(std::forward<value_types>(values))), ...);
  }
};

/*!
 * \brief Type alias for an array_values_type struct with double elements and a
 * dimension of 2.
 *
 * This alias simplifies the usage of `array_values_type` for the common case of
 * an array of two doubles.
 */
using arr_d2_t = array_values_type<double, 2>;

/*!
 * \brief Type alias for an array_values_type struct with iter_type elements and
 * a specified dimension.
 *
 * This alias allows for creating `array_values_type` structs with `iter_type`
 * elements and a variable dimension.
 *
 * \tparam D The dimension (size) of the array.
 */
template <size_t D>
using arr_n_t = array_values_type<iter_type, D>;

/*!
 * \brief Type alias for an array_values_type struct with len_type elements and
 * a specified dimension.
 *
 * This alias allows for creating `array_values_type` structs with `len_type`
 * elements and a variable dimension.
 *
 * \tparam D The dimension (size) of the array.
 */
template <size_t D>
using arr_l_t = array_values_type<len_type, D>;

/*!
 * \brief Type alias for an array_values_type struct with len_type elements and
 * a specified dimension.
 *
 * This alias allows for creating `array_values_type` structs with `len_type`
 * elements and a variable dimension.
 *
 * \tparam D The dimension (size) of the array.
 */
template <size_t D>
using arr2_l_t = array_2_values_type<len_type, D, 2>;

template <auto... values>
constexpr bool disjunction_values = (values || ...);

template <auto... values>
constexpr bool conjunction_values = (values && ...);

template <typename test_type, typename... Ts>
constexpr bool are_all_same_v =
    std::conjunction_v<std::is_same<test_type, Ts>...>;

template <size_t I, typename T>
struct type_ignore_index_impl {
  using type = T;
};

template <typename I, typename T>
struct type_ignore_type_impl {
  using type = T;
};

/*!
 * \brief A template alias that ignores its first parameter and returns the
 * second one.
 *
 * This template alias is used in metaprogramming contexts where a template with
 * a certain number of parameters is needed, but some of those parameters are
 * not used.
 *
 * \tparam I An index parameter that is ignored.
 * \tparam T The type that is returned by this template alias.
 */
template <size_t I, typename T>
using type_ignore_index = typename type_ignore_index_impl<I, T>::type;

/*!
 * \brief A template alias that ignores its first parameter and returns the
 * second one.
 *
 * This template alias is used in metaprogramming contexts where a template with
 * a certain number of parameters is needed, but some of those parameters are
 * not used.
 *
 * \tparam I An index parameter that is ignored.
 * \tparam T The type that is returned by this template alias.
 */
template <typename I, typename T>
using type_ignore_type = typename type_ignore_type_impl<I, T>::type;

/*!
 * \brief A template alias that ignores its first parameter and returns the
 * second one.
 *
 * This template alias is used in metaprogramming contexts where a template with
 * a certain number of parameters is needed, but some of those parameters are
 * not used.
 *
 * \tparam I An index parameter that is ignored.
 * \tparam T The type that is returned by this template alias.
 */
template <size_t I, typename T>
constexpr size_t index_ignore_type = I;

}  // namespace symphas

namespace symphas {
namespace lib {

template <size_t, typename T>
decltype(auto) repeat_value(T&& value) {
  return std::forward<T>(value);
}

template <typename, typename T>
decltype(auto) repeat_value(T&& value) {
  return std::forward<T>(value);
}

template <typename... Ts>
struct types_list {};
template <typename T>
struct unroll_types_list {};

template <typename T>
struct types_list_size;

template <typename... Ts>
struct types_list_size<types_list<Ts...>> {
  static const size_t value = sizeof...(Ts);
};

template <typename... Ts>
struct expand_types_list_impl;

template <>
struct expand_types_list_impl<> {
  using type = types_list<>;
};

template <typename T>
struct expand_types_list_impl<T> {
  using type = types_list<T>;
};

template <>
struct expand_types_list_impl<types_list<>> {
  using type = types_list<>;
};

template <typename T>
struct expand_types_list_impl<types_list<T>> {
  using type = typename expand_types_list_impl<T>::type;
};

template <typename... T1s, typename... T2s>
struct expand_types_list_impl<types_list<T1s...>, types_list<T2s...>> {
  using type = types_list<T1s..., T2s...>;
};

template <typename T1, typename... T2s>
struct expand_types_list_impl<T1, types_list<T2s...>> {
  using type = types_list<T1, T2s...>;
};

template <typename... T1s, typename T2>
struct expand_types_list_impl<types_list<T1s...>, T2> {
  using type = types_list<T1s..., T2>;
};

template <typename T1, typename T2>
struct expand_types_list_impl<T1, T2> {
  using type = types_list<T1, T2>;
};

template <typename T0, typename T1, typename... Ts>
struct expand_types_list_impl<T0, T1, Ts...> {
  using type = typename expand_types_list_impl<
      typename expand_types_list_impl<T0, T1>::type, Ts...>::type;
};

template <typename... Ts>
using expand_types_list = typename expand_types_list_impl<Ts...>::type;

// **************************************************************************************

template <size_t I, typename... Ts>
struct type_at_index_impl;

// template<>
// struct type_at_index_impl<0, types_list<>>
//{
//	using type = void;
// };

template <typename T0>
struct type_at_index_impl<0, T0> {
  using type = T0;
};

template <typename T0, typename T1, typename... Ts>
struct type_at_index_impl<0, T0, T1, Ts...> {
  using type = T0;
};

template <size_t I, typename T0, typename... Ts>
struct type_at_index_impl<I, T0, Ts...> {
  using type = typename type_at_index_impl<I - 1, Ts...>::type;
};

template <>
struct type_at_index_impl<0, unroll_types_list<types_list<>>> {
  using type = void;
};

template <typename T0, typename... Ts>
struct type_at_index_impl<0, unroll_types_list<types_list<T0, Ts...>>> {
  using type = T0;
};

template <size_t I, typename... Ts>
struct type_at_index_impl<I, unroll_types_list<types_list<Ts...>>> {
  using type = typename type_at_index_impl<I, Ts...>::type;
};

template <size_t I, typename... Ts>
using type_at_index = typename type_at_index_impl<I, Ts...>::type;

// Primary template for type_at_index_impl
template <size_t I, typename... Ts>
struct direct_type_at_index_impl;

// Specialization for the base case when I is 0
template <typename T0, typename... Ts>
struct direct_type_at_index_impl<0, T0, Ts...> {
  using type = T0;
};

// Recursive case: reduce the index and continue
template <size_t I, typename T0, typename... Ts>
struct direct_type_at_index_impl<I, T0, Ts...>
    : direct_type_at_index_impl<I - 1, Ts...> {
  using parent_type = direct_type_at_index_impl<I - 1, Ts...>;
  using type = typename parent_type::type;
};

// Type alias for easier usage
template <size_t I, typename... Ts>
using direct_type_at_index = typename direct_type_at_index_impl<I, Ts...>::type;

template <size_t I, typename... Ts>
struct types_after_at_index_impl;

template <size_t I>
struct types_after_at_index_impl<I> {
  using type = types_list<>;
};

template <size_t I, typename T0, typename... Ts>
struct types_after_at_index_impl<I, T0, Ts...> {
  using type = std::conditional_t<
      (I == 0), types_list<T0, Ts...>,
      typename types_after_at_index_impl<I - 1, Ts...>::type>;
};

template <size_t I, typename... Ts>
struct types_after_at_index_impl<I, unroll_types_list<types_list<Ts...>>> {
  using type = typename types_after_at_index_impl<I, Ts...>::type;
};

template <size_t I, typename... Ts>
using types_after_at_index = typename types_after_at_index_impl<I, Ts...>::type;

template <typename Seq, typename... Ts>
struct select_types_impl;

template <size_t... Is, typename... Ts>
struct select_types_impl<std::index_sequence<Is...>, Ts...> {
  using type = types_list<type_at_index<Is, Ts...>...>;
};

template <size_t... Is, typename... Ts>
struct select_types_impl<std::index_sequence<Is...>,
                         unroll_types_list<types_list<Ts...>>> {
  using type = types_list<type_at_index<Is, Ts...>...>;
};

template <typename Seq, typename... Ts>
using select_types = typename select_types_impl<Seq, Ts...>::type;

template <size_t I, typename T, typename... Ts>
struct types_before_index_impl;

template <size_t I, typename... T0s>
struct types_before_index_impl<I, types_list<T0s...>> {
  using type = types_list<T0s...>;
};

template <typename... T0s, typename T0, typename... Ts>
struct types_before_index_impl<0, types_list<T0s...>, T0, Ts...> {
  using type = types_list<T0s...>;
};

template <size_t I, typename... T0s, typename T0, typename... Ts>
struct types_before_index_impl<I, types_list<T0s...>, T0, Ts...> {
  using type = typename types_before_index_impl<I - 1, types_list<T0s..., T0>,
                                                Ts...>::type;
};

template <size_t I, typename... Ts>
struct types_before_index_impl<I, types_list<>,
                               unroll_types_list<types_list<Ts...>>> {
  using type = typename types_before_index_impl<I, types_list<>, Ts...>::type;
};

template <typename T0, typename... Ts>
struct types_before_index_impl<0, types_list<>,
                               unroll_types_list<types_list<T0, Ts...>>> {
  using type = types_list<T0>;
};

template <size_t I, typename... Ts>
using types_before_index =
    typename types_before_index_impl<I, types_list<>, Ts...>::type;

template <size_t I0, size_t I, typename... Ts>
struct types_between_index_impl {
  using type =
      types_after_at_index<I0, unroll_types_list<types_before_index<I, Ts...>>>;
};

template <size_t I0, size_t I, typename... Ts>
using types_between_index =
    typename types_between_index_impl<I0, I, Ts...>::type;

template <typename T, typename Seq>
struct reverse_types_list_impl;

template <typename... Es, size_t... Is>
struct reverse_types_list_impl<types_list<Es...>, std::index_sequence<Is...>> {
  using type = types_list<type_at_index<sizeof...(Is) - 1 - Is, Es...>...>;
};

template <typename... Es>
struct reverse_types_list_impl<types_list<Es...>, void> {
  using type = typename reverse_types_list_impl<
      types_list<Es...>, std::make_index_sequence<sizeof...(Es)>>::type;
};

//! Reverse a ::types_list.
template <typename T>
using reverse_types_list = typename reverse_types_list_impl<T, void>::type;

}  // namespace lib
}  // namespace symphas
