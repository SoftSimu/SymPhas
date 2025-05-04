
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

#include "spslibmetatypes.h"

namespace symphas {
namespace lib {

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
template <size_t I, typename T, T... Ns>
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
template <size_t I, typename T, T... Ns>
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
template <typename empty_type, typename... match_A_types, typename... As,
          typename... match_B_types, typename... Bs>
auto unzip_pot(std::tuple<std::pair<match_A_types, As>...> const& as,
               std::tuple<std::pair<match_B_types, Bs>...> const& bs);

//! Joining a single sequence simply returns it.
template <typename T, T... Ys>
constexpr auto seq_join(std::integer_sequence<T, Ys...>) {
  return std::integer_sequence<T, Ys...>{};
}

//! Joins two index sequences.
template <typename T, T... Ys, T... Qs>
constexpr auto seq_join(std::integer_sequence<T, Ys...>,
                        std::integer_sequence<T, Qs...>) {
  return std::integer_sequence<T, Ys..., Qs...>{};
}

//! Joins two index sequences.
template <typename T, T... Ys, T... Qs, typename... Seqs>
constexpr auto seq_join(std::integer_sequence<T, Ys...>,
                        std::integer_sequence<T, Qs...>, Seqs... seqs) {
  return seq_join(std::integer_sequence<T, Ys..., Qs...>{}, seqs...);
}

template <typename T, T... Ys>
constexpr auto seq_neg(std::integer_sequence<T, Ys...>) {
  return std::integer_sequence<int, -Ys...>{};
}

template <size_t... Ys>
constexpr auto seq_neg(std::index_sequence<Ys...>) {
  return std::integer_sequence<int, -int(Ys)...>{};
}

//! Adding a single sequence simply returns it.
template <typename T, T... Ys>
constexpr auto seq_add(std::integer_sequence<T, Ys...>) {
  return std::integer_sequence<T, Ys...>{};
}

template <typename T>
struct sequence_type_impl;

template <typename T, T... Ns>
struct sequence_type_impl<std::integer_sequence<T, Ns...>> {
  using type = T;
};

template <typename T>
using sequence_type = typename sequence_type_impl<T>::type;

//! Adds two index sequences.
/*!
 * The values are added pointwise, between two index sequences. When
 * sequences are not of equal size, the shortest one is considered
 * to have 0s in the remaining entries, so that a sequence equal in length
 * to the longest sequence is always returned.
 */
template <typename T0, T0... Ys, typename T1, T1... Qs,
          typename T = std::common_type_t<T0, T1>>
constexpr auto seq_add(std::integer_sequence<T0, Ys...>,
                       std::integer_sequence<T1, Qs...>) {
  if constexpr (sizeof...(Ys) == sizeof...(Qs)) {
    return std::integer_sequence<T, (T)(Qs + Ys)...>{};
  } else if constexpr (sizeof...(Ys) > sizeof...(Qs)) {
    return seq_join(seq_add(symphas::lib::get_seq_lt<sizeof...(Qs)>(
                                std::integer_sequence<T0, Ys...>{}),
                            std::integer_sequence<T1, Qs...>{}),
                    symphas::lib::get_seq_ge<sizeof...(Qs)>(
                        std::integer_sequence<T0, Ys...>{}));
  } else {
    return seq_join(seq_add(std::integer_sequence<T0, Ys...>{},
                            symphas::lib::get_seq_lt<sizeof...(Ys)>(
                                std::integer_sequence<T1, Qs...>{})),
                    symphas::lib::get_seq_ge<sizeof...(Ys)>(
                        std::integer_sequence<T1, Qs...>{}));
  }
}

//! Adds multiple index sequences.
/*!
 * The values are added pointwise, between two index sequences. When
 * sequences are not of equal size, the shortest one is considered
 * to have 0s in the remaining entries, so that a sequence equal in length
 * to the longest sequence is always returned.
 */
template <typename T0, T0... Ys, typename T1, T1... Qs, typename... Seqs>
constexpr auto seq_add(std::integer_sequence<T0, Ys...>,
                       std::integer_sequence<T1, Qs...>, Seqs... seqs) {
  return seq_add(seq_add(std::integer_sequence<T0, Ys...>{},
                         std::integer_sequence<T1, Qs...>{}),
                 seqs...);
}

//! Subtracting a single sequence simply returns it.
template <typename T, T... Ys>
constexpr auto seq_sub(std::integer_sequence<T, Ys...>) {
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
template <typename T0, T0... Ys, typename T1, T1... Qs, typename... Seqs>
constexpr auto seq_sub(std::integer_sequence<T0, Ys...>,
                       std::integer_sequence<T1, Qs...>, Seqs... seqs) {
  return seq_sub(seq_add(std::integer_sequence<int, int(Ys)...>{},
                         seq_neg(std::integer_sequence<T1, Qs...>{})),
                 seqs...);
}

template <typename T, T... Ys>
constexpr auto filter_seq(std::integer_sequence<T, Ys...>) {
  return std::integer_sequence<T, Ys...>{};
}

namespace {
template <typename T, T... Ys1, T Y0, T... Ys2, T Q0, T... Qs>
constexpr auto _filter_seq(std::integer_sequence<T, Ys1...>,
                           std::integer_sequence<T, Y0, Ys2...>,
                           std::integer_sequence<T, Q0>,
                           std::integer_sequence<T, Qs...>);

template <typename T, T... Ys1, T Q0, T... Qs>
constexpr auto _filter_seq(std::integer_sequence<T>,
                           std::integer_sequence<T, Ys1...>,
                           std::integer_sequence<T>,
                           std::integer_sequence<T, Q0, Qs...>) {
  return _filter_seq(
      std::integer_sequence<T>{}, std::integer_sequence<T, Ys1...>{},
      std::integer_sequence<T, Q0>{}, std::integer_sequence<T, Qs...>{});
}

template <typename T, T... Ys1>
constexpr auto _filter_seq(std::integer_sequence<T>,
                           std::integer_sequence<T, Ys1...>,
                           std::integer_sequence<T>, std::integer_sequence<T>) {
  return std::integer_sequence<T, Ys1...>{};
}

template <typename T, T Y0, T... Ys1>
constexpr auto _filter_seq(std::integer_sequence<T, Y0, Ys1...>,
                           std::integer_sequence<T>, std::integer_sequence<T>,
                           std::integer_sequence<T>) {
  return std::integer_sequence<T, Y0, Ys1...>{};
}

template <typename T, T... Ys1, T Q0>
constexpr auto _filter_seq(std::integer_sequence<T, Ys1...>,
                           std::integer_sequence<T>,
                           std::integer_sequence<T, Q0>,
                           std::integer_sequence<T>) {
  return std::integer_sequence<T, Ys1...>{};
}

template <typename T, T... Ys1, T Q0, T Q1, T... Qs>
constexpr auto _filter_seq(std::integer_sequence<T, Ys1...>,
                           std::integer_sequence<T>,
                           std::integer_sequence<T, Q0>,
                           std::integer_sequence<T, Q1, Qs...>) {
  return _filter_seq(
      std::integer_sequence<T>{}, std::integer_sequence<T, Ys1...>{},
      std::integer_sequence<T, Q1>{}, std::integer_sequence<T, Qs...>{});
}

template <typename T, T... Ys1, T Y0, T... Ys2, T Q0, T... Qs>
constexpr auto _filter_seq(std::integer_sequence<T, Ys1...>,
                           std::integer_sequence<T, Y0, Ys2...>,
                           std::integer_sequence<T, Q0>,
                           std::integer_sequence<T, Qs...>) {
  if constexpr (Y0 == Q0) {
    return _filter_seq(
        std::integer_sequence<T, Ys1...>{}, std::integer_sequence<T, Ys2...>{},
        std::integer_sequence<T, Q0>{}, std::integer_sequence<T, Qs...>{});
  } else {
    return _filter_seq(std::integer_sequence<T, Ys1..., Y0>{},
                       std::integer_sequence<T, Ys2...>{},
                       std::integer_sequence<T, Q0>{},
                       std::integer_sequence<T, Qs...>{});
  }
}
}  // namespace

//! Filters from the first sequence, all values in the subsequent sequences.
/*!
 * All values that appear in the second and subsequence sequences are removed
 * from the first sequence. If there are any values in the second sequence that
 * do not appear in the first, there is no effect. If the intersection of
 * elements shared between the sets is empty, then it will return the first
 * sequence.
 *
 * Note: Each element is filtered only once. This will not filter all repeating
 * elements.
 */
template <typename T, T... Ys, T... Qs, typename... Seqs>
constexpr auto filter_seq(std::integer_sequence<T, Ys...>,
                          std::integer_sequence<T, Qs...>, Seqs... seqs) {
  return filter_seq(
      _filter_seq(std::integer_sequence<T>{}, std::integer_sequence<T, Ys...>{},
                  std::integer_sequence<T>{},
                  std::integer_sequence<T, Qs...>{}),
      seqs...);
}

//! Makes a new sequence with values shared between all the sequences.
/*!
 * All values that are shared by all provided sequences will be put into a new
 * sequence.
 */
template <typename T, T... Ys>
constexpr auto intersect_seq(std::integer_sequence<T, Ys...>) {
  return std::integer_sequence<T, Ys...>{};
}

//! Makes a new sequence with values shared between all the sequences.
/*!
 * All values that are shared by all provided sequences will be put into a new
 * sequence.
 */
template <typename T, T... Ys, T... Qs, typename... Seqs>
constexpr auto intersect_seq(std::integer_sequence<T, Ys...>,
                             std::integer_sequence<T, Qs...>, Seqs... seqs);

//! The index sequence result type of joining multiple sequences.
template <typename Seq, typename... Seqs>
struct seq_join_result {
  using type = decltype(seq_join(std::declval<Seq>(), std::declval<Seqs>()...));
};

//! The index sequence result type of joining multiple sequences.
template <typename Seq, typename... Seqs>
struct seq_join_result<types_list<Seq, Seqs...>> {
  using type = decltype(seq_join(std::declval<Seq>(), std::declval<Seqs>()...));
};

//! The index sequence result type of joining multiple sequences.
template <>
struct seq_join_result<types_list<>> {
  using type = std::index_sequence<>;
};

//! The index sequence result type of adding multiple sequences.
template <typename Seq, typename... Seqs>
struct seq_add_result {
  using type = decltype(seq_add(std::declval<Seq>(), std::declval<Seqs>()...));
};

//! The index sequence result type of adding a value to each element of the
//! sequence.
template <size_t N, typename Seq>
struct seq_offset_result;

//! The index sequence result type of adding a value to each element of the
//! sequence.
template <size_t N, size_t... Is>
struct seq_offset_result<N, std::index_sequence<Is...>> {
  using type = std::index_sequence<(N + Is)...>;
};

//! The index sequence result type of adding multiple sequences.
template <typename Seq>
struct seq_from_to_result;

//! The index sequence result type of adding multiple sequences.
template <size_t I0, size_t J0>
struct seq_from_to_result<std::index_sequence<I0, J0>> {
  using type =
      typename seq_offset_result<I0, std::make_index_sequence<J0 - I0>>::type;
};

//! The index sequence result type of adding multiple sequences.
template <typename Seq, typename... Seqs>
struct seq_sub_result {
  using type = decltype(seq_sub(std::declval<Seq>(), std::declval<Seqs>()...));
};

//! The index sequence result type of filtering multiple sequences.
/*!
 * The first sequence that is given has the values removed from it that are
 * found in the other sequences.
 */
template <typename... Seqs>
struct filter_seq_result;

//! The index sequence result type of filtering multiple sequences.
template <typename T, T... Is>
struct filter_seq_result<std::integer_sequence<T, Is...>> {
  using type = std::integer_sequence<T, Is...>;
};

//! The index sequence result type of filtering multiple sequences.
template <typename T>
struct filter_seq_result<std::integer_sequence<T>, std::integer_sequence<T>> {
  using type = std::integer_sequence<T>;
};

template <typename T, T I0, T... Fs>
struct test_ne {
  static const bool value = ((I0 != Fs) && ...);
};

//! The index sequence result type of filtering multiple sequences.
template <typename T, T I0, T... Is, T... Fs>
struct filter_seq_result<std::integer_sequence<T, I0, Is...>,
                         std::integer_sequence<T, Fs...>> {
  using type = typename seq_join_result<
      std::conditional_t<test_ne<T, I0, Fs...>::value,
                         std::integer_sequence<T, I0>,
                         std::integer_sequence<T>>,
      std::conditional_t<test_ne<T, Is, Fs...>::value,
                         std::integer_sequence<T, Is>,
                         std::integer_sequence<T>>...>::type;
};

//! The index sequence result type of filtering multiple sequences.
template <typename T, T I0, T... Is, bool F0, bool... Fs>
struct filter_seq_result<std::integer_sequence<T, I0, Is...>,
                         std::integer_sequence<bool, F0, Fs...>> {
  using type = typename seq_join_result<
      std::conditional_t<F0, std::integer_sequence<T, I0>,
                         std::integer_sequence<T>>,
      std::conditional_t<Fs, std::integer_sequence<T, Is>,
                         std::integer_sequence<T>>...>::type;
};

//! The index sequence result type of filtering multiple sequences.
template <typename Seq0, typename Seq1, typename Seq2, typename... Seqs>
struct filter_seq_result<Seq0, Seq1, Seq2, Seqs...> {
  using type =
      typename filter_seq_result<typename filter_seq_result<Seq0, Seq1>::type,
                                 Seq2, Seqs...>::type;
};

//! The index sequence result type of intersecting multiple sequences.
template <typename... Seqs>
struct intersect_seq_result {
  using type = decltype(intersect_seq(std::declval<Seqs>()...));
};

//! Alias for the join result of multiple sequences.
template <typename Seq, typename... Seqs>
using seq_join_t = typename seq_join_result<Seq, Seqs...>::type;

//! Alias for the addition result of multiple sequences.
template <typename... Seqs>
using seq_add_t = typename seq_add_result<Seqs...>::type;

//! The index sequence result type of adding a value to each element of the
//! sequence.
template <size_t N, typename Seq>
using seq_offset_t = typename seq_offset_result<N, Seq>::type;

template <typename Seq>
using seq_from_to_t = typename seq_from_to_result<Seq>::type;

//! Alias for the subtraction result of multiple sequences.
template <typename... Seqs>
using seq_sub_t = typename seq_sub_result<Seqs...>::type;

//! Alias for the filter result of multiple sequences.
/*!
 * Filter from the first sequence everything that's in the subsequent sequences.
 * That is, any element in the sequences after the first will be removed from
 * the first.
 */
template <typename... Seqs>
using filter_seq_t = typename filter_seq_result<Seqs...>::type;

//! Alias for the intersection result of multiple sequences.
template <typename... Seqs>
using intersect_seq_t = typename intersect_seq_result<Seqs...>::type;

template <typename T, typename Seq>
struct change_seq_type_impl;

template <typename T, typename T0, T0... Is>
struct change_seq_type_impl<T, std::integer_sequence<T0, Is...>> {
  using type = std::integer_sequence<T, T(Is)...>;
};

template <typename T, typename Seq>
using change_seq_type_t = typename change_seq_type_impl<T, Seq>::type;

//! Returns whether there is the given value in the sequence.
/*!
 * Returns true if the value `N` is in the given sequence.
 *
 * \tparam N The existence of this value in the sequence is checked.
 * \tparam T The sequence to check.
 */
template <typename T, T N, typename Seq>
struct is_value_in_seq;

template <typename T, T N>
struct is_value_in_seq<T, N, std::integer_sequence<T>> {
  static const bool value = false;
};

template <typename T, T N, T I0, T... Is>
struct is_value_in_seq<T, N, std::integer_sequence<T, I0, Is...>> {
  static const bool value =
      (I0 == N) ||
      (is_value_in_seq<T, N, std::integer_sequence<T, Is...>>::value);
};

template <size_t N, typename T, typename Default = void>
struct seq_index_value_impl;

template <size_t N, typename T>
struct seq_index_value_impl<N, std::integer_sequence<T>, void>;

template <size_t N, typename T, T Default>
struct seq_index_value_impl<N, std::integer_sequence<T>,
                            std::integer_sequence<T, Default>> {
  static const T value = Default;
};

template <typename T, T I0, T... Is, typename Default>
struct seq_index_value_impl<0, std::integer_sequence<T, I0, Is...>, Default> {
  static const T value = I0;
};

template <size_t N, typename T, T I0, T... Is, typename Default>
struct seq_index_value_impl<N, std::integer_sequence<T, I0, Is...>, Default> {
  static const auto value =
      seq_index_value_impl<N - 1, std::integer_sequence<T, Is...>,
                           Default>::value;
};

template <int N, typename Seq, typename Default>
using seq_index_value_any_default = std::conditional_t<
    (N < 0), seq_index_value_impl<size_t(int(Seq::size()) + N), Seq, Default>,
    seq_index_value_impl<size_t(N), Seq, Default>>;

/*!
 * \brief Get the N-th value of an integer sequence, with a default value.
 *
 * This type alias is for a struct that retrieves the N-th value from an integer
 * sequence. If N is negative, it counts from the end of the sequence. If the
 * index is out of bounds, it returns the specified default value.
 *
 * \tparam N The index in the sequence.
 * \tparam Seq The type of the integer sequence.
 * \tparam Default The default value to return if the index is out of bounds.
 */
template <int N, typename Seq, sequence_type<Seq> Default>
using seq_index_value_with_default = seq_index_value_any_default<
    N, Seq, std::integer_sequence<sequence_type<Seq>, Default>>;

/*!
 * \brief Get the N-th value of an integer sequence.
 *
 * This type alias is for a struct that retrieves the N-th value from an integer
 * sequence. If N is negative, it counts from the end of the sequence. If the
 * index is out of bounds, it results in undefined behavior.
 *
 * \tparam N The index in the sequence.
 * \tparam Seq The type of the integer sequence.
 */
template <int N, typename Seq>
using seq_index_value = seq_index_value_any_default<N, Seq, void>;

template <typename Pick, typename From>
struct select_seq_impl;

template <typename T, size_t... Is, T... Ns>
struct select_seq_impl<std::index_sequence<Is...>,
                       std::integer_sequence<T, Ns...>> {
  using type = std::integer_sequence<
      T, seq_index_value<Is, std::integer_sequence<T, Ns...>>::value...>;
};

/*!
 * \brief Select elements from an integer sequence based on another sequence.
 *
 * This type alias is for a struct that takes two integer sequences, `Pick` and
 * `From`. It generates a new sequence where each element is selected from
 * `From` at the index specified in `Pick`.
 *
 * \tparam Pick The sequence of indices to select from `From`.
 * \tparam From The sequence from which elements are selected.
 */
template <typename Pick, typename From>
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
template <size_t N, typename T, T I>
struct seq_repeating_value {
 protected:
  template <size_t N0>
  static constexpr T V = I;

  template <size_t... Ns>
  static constexpr auto build_sequence(std::index_sequence<Ns...>) {
    return std::integer_sequence<T, V<Ns>...>{};
  }

 public:
  using type = decltype(build_sequence(std::make_index_sequence<N>{}));
};

/*!
 * \brief Generate a sequence with a repeating value.
 *
 * This type alias is for a struct that generates an integer sequence of length
 * `N`, where each element is the value `I`.
 *
 * \tparam N The length of the sequence.
 * \tparam T The type of the integer sequence.
 * \tparam I The value to be repeated in the sequence.
 */
template <size_t N, typename T, T I>
using seq_repeating_value_t = typename seq_repeating_value<N, T, I>::type;

//! Makes a new sequence with values shared between all the sequences.
/*!
 * All values that are shared by all provided sequences will be put into a new
 * sequence.
 */
template <typename T, T... Ys, T... Qs, typename... Seqs>
constexpr auto intersect_seq(std::integer_sequence<T, Ys...>,
                             std::integer_sequence<T, Qs...>, Seqs... seqs) {
  using filtered_t = seq_join_t<
      std::integer_sequence<T>,
      std::conditional_t<
          is_value_in_seq<T, Ys, std::integer_sequence<T, Qs...>>::value,
          std::integer_sequence<T, Ys>, std::integer_sequence<T>>...>;
  return intersect_seq(filtered_t{}, seqs...);
}

/*!
 * \brief Calculates the product of the sizes of a number of sequences.
 *
 * This function uses template recursion to calculate the product of the sizes
 * of one or more sequences. It takes one or more template parameters, each
 * representing a sequence. The `Seq` parameter represents the first sequence,
 * and `Seqs...` represents the rest of the sequences.
 *
 * \tparam Seq The first sequence.
 * \tparam Seqs The rest of the sequences.
 * \return The product of the sizes of the sequences.
 */
template <typename Seq, typename... Seqs>
struct seq_len_product;

template <typename Seq>
struct seq_len_product<Seq> {
  static constexpr size_t value = Seq::size();
};

template <typename Seq, typename Seq2, typename... Seqs>
struct seq_len_product<Seq, Seq2, Seqs...> {
  static constexpr size_t value =
      Seq::size() * seq_len_product<Seq2, Seqs...>::value;
};

template <typename T, size_t V>
size_t constexpr get_value_from_seq(std::integer_sequence<T, V>) {
  return V;
}

template <size_t N, typename T, T... Es,
          typename std::enable_if_t<(N < sizeof...(Es)), int> = 0>
size_t constexpr seq_value(std::integer_sequence<T, Es...>);

template <size_t N, typename T, T... Es,
          typename std::enable_if_t<(N < sizeof...(Es)), int>>
size_t constexpr seq_value(std::integer_sequence<T, Es...>) {
  return get_value_from_seq(
      symphas::lib::type_at_index<N, std::integer_sequence<T, Es>...>{});
}

template <typename Seq0, typename Seq1>
struct make_seq_pairs;

template <size_t... Is, size_t... Js>
struct make_seq_pairs<std::index_sequence<Is...>, std::index_sequence<Js...>> {
  using type = expand_types_list<std::conditional_t<
      (Is < Js), std::index_sequence<Is, Js>, types_list<>>...>;
};

template <size_t N, typename... Seqs>
struct seq_skip_indices_impl;

template <>
struct seq_skip_indices_impl<0, std::index_sequence<0>> {
  using type = std::index_sequence<>;
};

template <size_t N>
struct seq_skip_indices_impl<N, std::index_sequence<N - 1>> {
  using type = std::make_index_sequence<N - 1>;
};

template <size_t N>
struct seq_skip_indices_impl<N, std::index_sequence<>> {
  using type = std::make_index_sequence<N>;
};

template <size_t N>
struct seq_skip_indices_impl<N, std::index_sequence<0>> {
  using type = seq_offset_t<1, std::make_index_sequence<N - 1>>;
};

template <size_t N, size_t I0>
struct seq_skip_indices_impl<N, std::index_sequence<I0>> {
  using type =
      seq_join_t<std::make_index_sequence<I0>,
                 seq_offset_t<I0 + 1, std::make_index_sequence<N - I0 - 1>>>;
};

template <size_t N, size_t I0, size_t J0>
struct seq_skip_indices_impl<N, types_list<std::index_sequence<I0, J0>>> {
  using type = seq_join_t<seq_from_to_t<std::index_sequence<I0, J0>>>;
};

template <size_t N, typename Seq0, typename Seq1, typename... Seqs>
struct seq_skip_indices_impl<N, Seq0, Seq1, Seqs...> {
  using type = seq_join_t<seq_from_to_t<Seq0>, seq_from_to_t<Seq1>,
                          seq_from_to_t<Seqs>...>;
};

template <size_t N, typename Seq0, typename Seq1, typename... Seqs>
struct seq_skip_indices_impl<N, types_list<Seq0, Seq1, Seqs...>> {
  using type = typename seq_skip_indices_impl<N, Seq0, Seq1, Seqs...>::type;
};

template <size_t N>
struct seq_skip_indices_impl<N, types_list<>> {
  using type = std::index_sequence<>;
};

template <size_t N, size_t I0, size_t I1, size_t... Is>
struct seq_skip_indices_impl<N, std::index_sequence<I0, I1, Is...>> {
  using type = std::conditional_t<
      std::is_same<std::make_index_sequence<N>,
                   std::index_sequence<I0, I1, Is...>>::value,
      std::index_sequence<>,
      typename seq_skip_indices_impl<
          N, typename make_seq_pairs<
                 std::index_sequence<0, I0 + 1, I1 + 1, (Is + 1)...>,
                 std::index_sequence<I0, I1, Is..., N>>::type>::type>;
};

template <size_t N, typename Seq>
using seq_skip_indices = typename seq_skip_indices_impl<N, Seq>::type;

template <typename Seq, typename T>
struct list_repeating_type_seq;

template <size_t... Is, typename T>
struct list_repeating_type_seq<std::index_sequence<Is...>, T> {
  using type = types_list<type_ignore_index<Is, T>...>;
};

template <size_t N, typename T>
struct list_repeating_type {
  using type =
      typename list_repeating_type_seq<std::make_index_sequence<N>, T>::type;
};

template <size_t N, typename T>
using list_repeating_type_t = typename list_repeating_type<N, T>::type;

template <size_t N, typename T>
struct repeat_each_element;

template <size_t N, typename T, T... Is>
struct repeat_each_element<N, std::integer_sequence<T, Is...>> {
  using type =
      seq_join_t<std::integer_sequence<T>, seq_repeating_value_t<N, T, Is>...>;
};

/*!
 * \brief Repeat each element in an integer sequence.
 *
 * This type alias is for a struct that takes an integer sequence and a number
 * `N`, and generates a new sequence where each element from the original
 * sequence is repeated `N` times in a row.
 *
 * \tparam N The number of times each element is repeated.
 * \tparam T The type of the integer sequence.
 */

template <size_t N, typename T>
using repeat_each_element_t = typename repeat_each_element<N, T>::type;

template <size_t N, typename T>
struct repeat_all_elements;

template <size_t N, typename T, T... Is>
struct repeat_all_elements<N, std::integer_sequence<T, Is...>> {
  template <size_t... Ns>
  static constexpr auto build_sequence(std::index_sequence<Ns...>) {
    return seq_join_t<
        std::integer_sequence<T>,
        symphas::type_ignore_index<Ns, std::integer_sequence<T, Is...>>...>{};
  }

 public:
  using type = decltype(build_sequence(std::make_index_sequence<N>{}));
};

/*!
 * \brief Repeat the entire integer sequence.
 *
 * This type alias is for a struct that takes an integer sequence and a number
 * `N`, and generates a new sequence where the entire original sequence is
 * repeated `N` times.
 *
 * \tparam N The number of times the entire sequence is repeated.
 * \tparam T The type of the integer sequence.
 */
template <size_t N, typename T>
using repeat_all_elements_t = typename repeat_all_elements<N, T>::type;

template <size_t I, typename Seq0>
struct seq_ge_impl;

template <size_t I, typename T, T... Ns>
struct seq_ge_impl<I, std::integer_sequence<T, Ns...>> {
  using type = std::conditional_t<
      (sizeof...(Ns) <= I), std::integer_sequence<T>,
      symphas::lib::select_seq<
          seq_offset_t<I, std::make_index_sequence<
                              fixed_min<sizeof...(Ns) - I, sizeof...(Ns)>>>,
          std::integer_sequence<T, Ns...>>>;
};

template <size_t I, typename Seq0>
using seq_ge_t = typename seq_ge_impl<I, Seq0>::type;

template <size_t I, typename Seq0>
struct seq_lt_impl;

template <size_t I, typename T, T... Ns>
struct seq_lt_impl<I, std::integer_sequence<T, Ns...>> {
  using type = symphas::lib::select_seq<
      std::make_index_sequence<fixed_min<I, sizeof...(Ns)>>,
      std::integer_sequence<T, Ns...>>;
};

template <size_t I, typename Seq0>
using seq_lt_t = typename seq_lt_impl<I, Seq0>::type;

template <typename T, T N, typename... Seqs>
struct find_index_lt;

// find the largest index of the sequence that is smaller than N
template <typename T, T N, T I0>
struct find_index_lt<T, N, std::integer_sequence<T, I0>> {
  static const size_t value = (I0 < N) ? 0 : ~0ull;
};

// find the largest index of the sequence that is smaller than N
template <typename T, T N, T I0, T I1, T... Is>
struct find_index_lt<T, N, std::integer_sequence<T, I0, I1, Is...>> {
  static const size_t value =
      (seq_index_value<sizeof...(Is) / 2 + 1,
                       std::integer_sequence<T, I0, I1, Is...>>::value <= N)
          ? sizeof...(Is) / 2 + 1 +
                find_index_lt<
                    T, N,
                    seq_ge_t<sizeof...(Is) / 2 + 1,
                             std::integer_sequence<T, I0, I1, Is...>>>::value
          : find_index_lt<
                T, N,
                seq_lt_t<sizeof...(Is) / 2 + 1,
                         std::integer_sequence<T, I0, I1, Is...>>>::value;
};

template <typename T, T N, typename... Seqs>
struct find_index_ge;

// find the largest index of the sequence that is smaller than N
template <typename T, T N, T I0>
struct find_index_ge<T, N, std::integer_sequence<T, I0>> {
  static const size_t value = (I0 >= N) ? 0 : 1;
};

// find the largest index of the sequence that is smaller than N
template <typename T, T N, T I0, T I1, T... Is>
struct find_index_ge<T, N, std::integer_sequence<T, I0, I1, Is...>> {
  static const size_t value =
      (seq_index_value<sizeof...(Is) / 2,
                       std::integer_sequence<T, I0, I1, Is...>>::value < N)
          ? sizeof...(Is) / 2 + 1 +
                find_index_ge<
                    T, N,
                    seq_ge_t<sizeof...(Is) / 2 + 1,
                             std::integer_sequence<T, I0, I1, Is...>>>::value
          : find_index_ge<
                T, N,
                seq_lt_t<sizeof...(Is) / 2 + 1,
                         std::integer_sequence<T, I0, I1, Is...>>>::value;
};

template <typename T, T I, typename Seq0>
struct seq_before_impl;

template <typename T, T I, T... Ns>
struct seq_before_impl<T, I, std::integer_sequence<T, Ns...>> {
  using type = symphas::lib::select_seq<
      std::make_index_sequence<
          find_index_lt<T, I, std::integer_sequence<T, Ns...>>::value + 1>,
      std::integer_sequence<T, Ns...>>;
};

template <typename T, T I, typename Seq0>
using seq_before_t = typename seq_before_impl<T, I, Seq0>::type;

template <typename T, T I, typename Seq0>
struct seq_after_at_impl;

template <typename T, T I, T... Ns>
struct seq_after_at_impl<T, I, std::integer_sequence<T, Ns...>> {
 protected:
  static const size_t N =
      find_index_ge<T, I, std::integer_sequence<T, Ns...>>::value;

 public:
  using type = symphas::lib::select_seq<
      seq_offset_t<N, std::make_index_sequence<sizeof...(Ns) - N>>,
      std::integer_sequence<T, Ns...>>;
};

template <typename T, T I, typename Seq0>
using seq_after_at_t = typename seq_after_at_impl<T, I, Seq0>::type;

template <bool, typename Seq0, typename Seq1>
struct select_merge_sort_seq_impl;

template <typename Seq0, typename Seq1>
struct merge_sort_seq_impl;

template <typename T, T N0, T... Ns>
struct merge_sort_seq_impl<std::integer_sequence<T, N0, Ns...>,
                           std::integer_sequence<T>> {
  using type = std::integer_sequence<T, N0, Ns...>;
};

template <typename T, T M0, T... Ms>
struct merge_sort_seq_impl<std::integer_sequence<T>,
                           std::integer_sequence<T, M0, Ms...>> {
  using type = std::integer_sequence<T, M0, Ms...>;
};

template <typename T, T N0, T... Ns, T M0, T... Ms>
struct merge_sort_seq_impl<std::integer_sequence<T, N0, Ns...>,
                           std::integer_sequence<T, M0, Ms...>> {
  using type = typename select_merge_sort_seq_impl<
      (N0 < M0), std::integer_sequence<T, N0, Ns...>,
      std::integer_sequence<T, M0, Ms...>>::type;
};

template <typename T, T N0, T... Ns, T M0, T... Ms>
struct select_merge_sort_seq_impl<true, std::integer_sequence<T, N0, Ns...>,
                                  std::integer_sequence<T, M0, Ms...>> {
  using type = typename merge_sort_seq_impl<
      seq_join_t<seq_before_t<T, M0, std::integer_sequence<T, N0, Ns...>>,
                 std::integer_sequence<T, M0>,
                 seq_after_at_t<T, M0, std::integer_sequence<T, N0, Ns...>>>,
      std::integer_sequence<T, Ms...>>::type;
};

template <typename T, T N0, T... Ns, T M0, T... Ms>
struct select_merge_sort_seq_impl<false, std::integer_sequence<T, N0, Ns...>,
                                  std::integer_sequence<T, M0, Ms...>> {
  using type = typename merge_sort_seq_impl<
      seq_join_t<seq_before_t<T, N0, std::integer_sequence<T, M0, Ms...>>,
                 std::integer_sequence<T, N0>,
                 seq_after_at_t<T, N0, std::integer_sequence<T, M0, Ms...>>>,
      std::integer_sequence<T, Ns...>>::type;
};

template <typename T, T N0, T... Ns, T... Ms>
struct select_merge_sort_seq_impl<false, std::integer_sequence<T, N0, Ns...>,
                                  std::integer_sequence<T, N0, Ms...>> {
  using type =
      typename merge_sort_seq_impl<std::integer_sequence<T, N0, N0, Ns...>,
                                   std::integer_sequence<T, Ms...>>::type;
};

template <typename Seq0, typename Seq1>
using merge_sort_seq = typename merge_sort_seq_impl<Seq0, Seq1>::type;

template <typename Seq>
struct sorted_seq_impl;

template <typename T>
struct sorted_seq_impl<std::integer_sequence<T>> {
  using type = std::integer_sequence<T>;
};

template <typename T, T N0>
struct sorted_seq_impl<std::integer_sequence<T, N0>> {
  using type = std::integer_sequence<T, N0>;
};

template <typename T, T N0, T N1, T... Ns>
struct sorted_seq_impl<std::integer_sequence<T, N0, N1, Ns...>> {
  using type = merge_sort_seq<
      typename sorted_seq_impl<
          seq_ge_t<sizeof...(Ns) / 2 + 1,
                   std::integer_sequence<T, N0, N1, Ns...>>>::type,
      typename sorted_seq_impl<
          seq_lt_t<sizeof...(Ns) / 2 + 1,
                   std::integer_sequence<T, N0, N1, Ns...>>>::type>;
};

template <typename Seq0>
using sorted_seq = typename sorted_seq_impl<Seq0>::type;

//! Puts the nth value of all sequences into the result.
template <size_t N, typename... Ts>
struct nth_value_of_seqs {
  using type =
      std::integer_sequence<bool,
                            symphas::lib::seq_index_value<N, Ts>::value...>;
};

template <size_t N, typename... Ts>
using nth_value_of_seqs_t = typename nth_value_of_seqs<N, Ts...>::type;

template <bool flag, size_t N, size_t D, size_t I, typename T0, typename... Ts>
struct nth_periodic_shift_impl;

template <size_t N, size_t D, typename T0, typename... Ts>
struct nth_periodic_shift_impl<true, N, D, D, T0, Ts...> {
  using type = typename nth_value_of_seqs<N, T0, Ts...>::type;
};

template <size_t N, size_t D, size_t I, typename T0, typename... Ts>
struct nth_periodic_shift_impl<false, N, D, I, T0, Ts...> {
  using type = typename nth_periodic_shift_impl<
      D == I + 1, N, D, I + 1, symphas::lib::seq_join_t<T0, T0>,
      symphas::lib::seq_join_t<Ts, Ts>...,
      symphas::lib::seq_join_t<
          symphas::lib::seq_repeating_value_t<T0::size(), bool, 0>,
          symphas::lib::seq_repeating_value_t<T0::size(), bool, 1>>>::type;
};

template <size_t N, size_t D>
struct nth_periodic_shift {
  using type =
      typename nth_periodic_shift_impl<D == 1, N, D, 1,
                                       std::integer_sequence<bool, 0, 1>>::type;
};

template <size_t N, size_t D>
using nth_periodic_shift_t = typename nth_periodic_shift<N, D>::type;

}  // namespace lib
}  // namespace symphas
