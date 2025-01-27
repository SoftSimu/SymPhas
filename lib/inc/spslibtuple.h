
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

#include "definitions.h"

namespace symphas {
namespace lib {

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
template <typename... Ts>
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
template <typename... As, typename... Bs>
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
template <size_t I, typename... Ts>
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
template <size_t I, typename... Ts>
auto get_tuple_lt(std::tuple<Ts...> const& t);

namespace {
// helper function that iterates over the tuple
template <typename... As, typename... Bs, typename T, T... Is>
auto unfurl_tuple(std::tuple<std::pair<As, Bs>...> const& ts,
                  std::integer_sequence<T, Is...>) {
  auto first_elements = symphas::lib::make_tuple(std::get<Is>(ts).first...);
  auto second_elements = symphas::lib::make_tuple(std::get<Is>(ts).second...);
  return std::make_pair(first_elements, second_elements);
}

// helper function that takes index sequence and rebuilds the tuple from the
// given index onwards
template <size_t I, typename... Ts, size_t... Is>
auto get_tuple_ge(std::tuple<Ts...> const& ts, std::index_sequence<Is...>) {
  return symphas::lib::make_tuple(std::get<I + Is>(ts)...);
}

// helper function that takes index sequence and rebuilds the tuple up to and
// not including the given index
template <typename... Ts, typename T, T... Is>
auto get_tuple_lt(std::tuple<Ts...> const& ts,
                  std::integer_sequence<T, Is...>) {
  return symphas::lib::make_tuple(std::get<Is>(ts)...);
}

template <typename T, T... Ns>
auto get_seq_from_tuple(std::tuple<std::integer_sequence<T, Ns>...> const&) {
  return std::integer_sequence<T, Ns...>{};
}

//! End the recursion when there is nothing remaining to match against.
template <typename empty_type, typename match_type, typename A,
          typename... match_A_types, typename... As, typename B>
auto unzip_pot_one(std::tuple<std::pair<match_type, A>,
                              std::pair<match_A_types, As>...> const& as,
                   std::pair<match_type, B> const& b0) {
  auto matched = std::make_tuple(std::make_pair(
      b0.first, std::make_pair(std::get<0>(as).second, b0.second)));
  auto rest = symphas::lib::get_tuple_ge<1>(as);

  return std::make_pair(matched, rest);
}

//! End the recursion when there is nothing remaining to match against.
template <typename empty_type, typename match_B_type, typename B>
auto unzip_pot_one(std::tuple<> const&, std::pair<match_B_type, B> const& b0) {
  auto matched = std::make_tuple(
      std::make_pair(b0.first, std::make_pair(empty_type{}, b0.second)));
  return std::make_pair(matched, std::make_tuple());
}

//! End the recursion when there is nothing remaining to match against.
template <typename empty_type, typename match_A_type, typename A,
          typename... match_A_types, typename... As, typename match_B_type,
          typename B,
          typename std::enable_if<
              !std::is_same<match_A_type, match_B_type>::value, int>::type = 0>
auto unzip_pot_one(std::tuple<std::pair<match_A_type, A>,
                              std::pair<match_A_types, As>...> const& as,
                   std::pair<match_B_type, B> const& b0) {
  auto [matched, unmatched] =
      unzip_pot_one<empty_type>(symphas::lib::get_tuple_ge<1>(as), b0);
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
template <typename empty_type, typename... match_B_types, typename... Bs,
          size_t... Is>
auto unzip_pot_b_apply(std::tuple<std::pair<match_B_types, Bs>...> const& bs,
                       std::index_sequence<Is...>) {
  return std::make_tuple(
      std::make_pair(std::get<Is>(bs).first,
                     std::make_pair(empty_type{}, std::get<Is>(bs).second))...);
}

//! Construct the `a` data element.
/*!
 * Each pair of the tuple has the match type as the first member and
 * a pair as the second member, which contains the content from the
 * given `bs` as the second element and an default constructed
 * instance of the type `empty_type` as the first element. The resulting
 * list of pairs is thus the same length as the input list.
 */
template <typename empty_type, typename... match_A_types, typename... As,
          size_t... Is>
auto unzip_pot_a_apply(std::tuple<std::pair<match_A_types, As>...> const& as,
                       std::index_sequence<Is...>) {
  return std::make_tuple(
      std::make_pair(std::get<Is>(as).first,
                     std::make_pair(std::get<Is>(as).second, empty_type{}))...);
}

//! Construct the `b` data element.
template <typename empty_type, typename... match_B_types, typename... Bs>
auto unzip_pot_b(std::tuple<std::pair<match_B_types, Bs>...> const& bs) {
  using bs_type = std::tuple<std::pair<match_B_types, Bs>...>;
  constexpr size_t bs_len = std::tuple_size<bs_type>::value;

  return unzip_pot_b_apply<empty_type>(bs, std::make_index_sequence<bs_len>{});
}

//! Construct the `a` data element.
template <typename empty_type, typename... match_A_types, typename... As>
auto unzip_pot_a(std::tuple<std::pair<match_A_types, As>...> const& as) {
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
template <typename empty_type, typename... Rs>
auto unzip_pot_carry_apply(std::tuple<Rs...> const& matched,
                           std::tuple<> const&, std::tuple<> const&) {
  return matched;
}

//! Apply the algorithm to unzip a list based on the types.
/*!
 * See symphas::lib::unzip_pot_carry_apply.
 */
template <typename empty_type, typename... Rs, typename match_A_type,
          typename A, typename... match_A_types, typename... As>
auto unzip_pot_carry_apply(
    std::tuple<Rs...> const& matched,
    std::tuple<std::pair<match_A_type, A>,
               std::pair<match_A_types, As>...> const& as,
    std::tuple<> const&) {
  auto pot_as = unzip_pot_a<empty_type>(as);
  auto matched_list = std::tuple_cat(matched, pot_as);

  return matched_list;
}

//! Apply the algorithm to unzip a list based on the types.
/*!
 * See symphas::lib::unzip_pot_carry_apply.
 */
template <typename empty_type, typename... Rs, typename match_B_type,
          typename B, typename... match_B_types, typename... Bs>
auto unzip_pot_carry_apply(
    std::tuple<Rs...> const& matched, std::tuple<> const&,
    std::tuple<std::pair<match_B_type, B>,
               std::pair<match_B_types, Bs>...> const& bs) {
  auto pot_bs = unzip_pot_b<empty_type>(bs);
  auto matched_list = std::tuple_cat(matched, pot_bs);

  return matched_list;
}

//! Apply the algorithm to unzip a list based on the types.
/*!
 * See symphas::lib::unzip_pot_carry_apply.
 */
template <typename empty_type, typename... Rs, typename match_A_type,
          typename A, typename... match_A_types, typename... As,
          typename match_B_type, typename B, typename... match_B_types,
          typename... Bs>
auto unzip_pot_carry_apply(
    std::tuple<Rs...> const& matched,
    std::tuple<std::pair<match_A_type, A>,
               std::pair<match_A_types, As>...> const& as,
    std::tuple<std::pair<match_B_type, B>,
               std::pair<match_B_types, Bs>...> const& bs) {
  auto [one_matched, as_rest] = unzip_pot_one<empty_type>(as, std::get<0>(bs));
  auto match_append = std::tuple_cat(matched, one_matched);
  auto bs_rest = symphas::lib::get_tuple_ge<1>(bs);

  return unzip_pot_carry_apply<empty_type>(match_append, as_rest, bs_rest);
}

}  // namespace

template <typename... As, typename... Bs>
auto unfurl_tuple(std::tuple<std::pair<As, Bs>...> const& ts) {
  return unfurl_tuple(
      ts, std::make_index_sequence<
              std::tuple_size<std::tuple<std::pair<As, Bs>...>>::value>{});
}

//! Get the elements of the tuple greater or equal to the given index.
/*!
 * Returns a new tuple consisting of all the elements with index greater than or
 * equal to the given index. If the given index is 0, then the original tuple
 * will be returned.
 */
template <size_t I, typename... Ts>
auto get_tuple_ge(std::tuple<Ts...> const& ts) {
  return get_tuple_ge<I>(ts, std::make_index_sequence<sizeof...(Ts) - I>{});
}

//! Get the elements of the tuple less than the given index.
/*!
 * Returns a new tuple consisting of all the elements with index strictly less
 * than the given index. If the index is equal to the length of the tuple, the
 * original tuple will be returned.
 */
template <size_t I, typename... Ts>
auto get_tuple_lt(std::tuple<Ts...> const& ts) {
  return get_tuple_lt(ts, std::make_index_sequence<I>{});
}

//! Get the elements of the tuple in the given range, not including the last
//! index.
/*!
 * Returns a new tuple consisting of all the elements between the first and last
 * indices provided, not including the last index. That is, this function will
 * return the same tuple the last index is equal to the size of the tuple.
 */
template <size_t I0, size_t In, typename... Ts,
          typename std::enable_if_t<(I0 < In), int> = 0>
auto get_tuple_bw(std::tuple<Ts...> const& ts) {
  return get_tuple_lt(
      ts, symphas::lib::seq_add(
              symphas::lib::seq_repeating_value_t<In - I0, size_t, I0>{},
              std::make_index_sequence<In - I0>{}));
}

template <size_t I0, size_t In, typename... Ts,
          typename std::enable_if_t<(I0 >= In), int> = 0>
auto get_tuple_bw(std::tuple<Ts...> const& ts) {
  return std::tuple<>{};
}

template <size_t I, typename T, T... Ns>
auto get_seq_ge(std::integer_sequence<T, Ns...> const& seq) {
  if constexpr (sizeof...(Ns) <= I) {
    return std::integer_sequence<T>{};
  } else {
    return symphas::lib::select_seq<
        seq_offset_t<I, std::make_index_sequence<sizeof...(Ns) - I>>,
        std::integer_sequence<T, Ns...>>{};
  }
}

template <size_t I, typename T, T... Ns>
auto get_seq_lt(std::integer_sequence<T, Ns...> const& seq) {
  return symphas::lib::select_seq<std::make_index_sequence<I>,
                                  std::integer_sequence<T, Ns...>>{};
}

template <typename... Ts>
auto make_tuple(Ts... ts) {
  return std::tuple_cat(std::tuple<Ts>(ts)...);
}

template <typename empty_type, typename... match_A_types, typename... As,
          typename... match_B_types, typename... Bs>
auto unzip_pot(std::tuple<std::pair<match_A_types, As>...> const& a,
               std::tuple<std::pair<match_B_types, Bs>...> const& b) {
  return unzip_pot_carry_apply<empty_type>(std::make_tuple(), a, b);
}

}  // namespace lib
}  // namespace symphas
