
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
 * MODULE:  expr
 * PURPOSE: Defines ways to transform expressions, usually from one
 * expression to the other based on some rules or the type of expression.
 *
 * ***************************************************************************
 */

#pragma once

#include <fstream>
#include <iostream>
#include <type_traits>

#include "expressionapply.h"

// ******************************************************************************************

namespace expr::split {

namespace {

/* given the derivative index, retrieve the lowest order
 */

template <size_t I>
struct min_order_from_index {
 protected:
  static constexpr size_t get_value() {
    iter_type OD = 0;
    while (true) {
      if (I & (1ull << OD)) {
        return OD;
      }
      ++OD;
    }
  }

 public:
  static const size_t value = get_value();
};

/* convenience functions to pack an expression into either side of a pair
 */

template <typename E>
auto pack_left(OpExpression<E> const& e) {
  return std::make_pair(*static_cast<const E*>(&e), OpVoid{});
}

template <typename E>
auto pack_right(OpExpression<E> const& e) {
  return std::make_pair(OpVoid{}, *static_cast<const E*>(&e));
}

template <typename E>
auto pack_left(OpOperator<E> const& e) {
  return std::make_pair(*static_cast<const E*>(&e), OpVoid{});
}

template <typename E>
auto pack_right(OpOperator<E> const& e) {
  return std::make_pair(OpVoid{}, *static_cast<const E*>(&e));
}

template <typename... E0s, typename... E1s>
auto adds_expand_pair(std::pair<E0s, E1s> const&... pairs) {
  return std::make_pair((pairs.first + ...), (pairs.second + ...));
}

template <typename E0, typename E1, typename... E0s, typename... E1s>
auto adds_expand_pair_no_first(std::pair<E0, E1> const& pair0,
                               std::pair<E0s, E1s> const&... pairs) {
  return std::make_pair(pair0.first, (pair0.second + ... + pairs.second));
}
}  // namespace

//! Determine the smallest derivative order of the expression.
/*!
 * Determine the smallest derivative order of the expression.
 *
 * \tparam E The expression type from which to get the minimum derivative
 * order.
 */
template <typename E>
struct min_derivative_order {
  static const size_t value =
      min_order_from_index<expr::derivative_index<0, E>::value>::value;
};

// **************************************************************************************

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename Dd, typename V, typename E, typename Sp>
auto by_linear(OpDerivative<Dd, V, E, Sp> const& e);

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename A1, typename A2, typename E>
auto by_linear(OpCombination<A1, A2, E> const& e);

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename A1, typename A2, typename E>
auto by_linear(OpChain<A1, A2, E> const& e);

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename T, typename G>
auto by_linear(OpTerm<T, G> const& e);

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename V, typename E1, typename E2>
auto by_linear(OpConvolution<V, E1, E2> const& e);

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename E>
auto by_linear(OpExpression<E> const& e) {
  return pack_right(*static_cast<E const*>(&e));
}

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename E>
auto by_linear(OpOperator<E> const& e) {
  return pack_right(*static_cast<E const*>(&e));
}

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename... Es>
auto by_linear(OpAdd<Es...> const& e);

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename A, typename B>
auto by_linear(OpBinaryDiv<A, B> const& e);

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename G, typename V, typename E>
auto by_linear(OpMap<G, V, E> const& e);

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename V, typename sub_t, typename E, typename... Ts>
auto by_linear(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e);

//! Split all linear variables from nonlinear variables.
/*!
 * Using the available predicates to determine linearity/nonlinearity of
 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
 * for linearity. Similarly, the `NL` part is nonlinear.
 *
 * \param e The expression which is split.
 */
template <typename V, typename E, typename... Ts, int... I0s, int... P0s,
          typename A, typename... As, typename B, typename C>
auto by_linear(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>,
                     symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
                     SymbolicFunction<A, As...>, B, C> const& series);

// derivative
namespace {

template <typename Dd, typename V, typename E, typename T,
          typename std::enable_if_t<(expression_predicate<E>::linear), int> = 0>
auto by_linear_derivative(OpDerivative<Dd, V, E, T> const& e) {
  return pack_left(e);
}

template <typename Dd, typename V, typename E, typename T,
          typename std::enable_if_t<(!expression_predicate<E>::linear &&
                                     expression_predicate<E>::combination),
                                    int> = 0>
auto by_linear_derivative(OpDerivative<Dd, V, E, T> const& e) {
  return expr::split::by_linear(expr::apply_operators(e));
}

template <typename Dd, typename V, typename E, typename T,
          typename std::enable_if_t<(!expression_predicate<E>::linear &&
                                     !expression_predicate<E>::combination),
                                    int> = 0>
auto by_linear_derivative(OpDerivative<Dd, V, E, T> const& e) {
  return pack_right(e);
}
}  // namespace

template <typename Dd, typename V, typename E, typename Sp>
auto by_linear(OpDerivative<Dd, V, E, Sp> const& e) {
  return by_linear_derivative(e);
}

// op combination

namespace {
template <
    typename A1, typename A2, typename E,
    typename std::enable_if_t<
        (expression_predicate<OpCombination<A1, A2, E>>::linear), int> = 0>
auto by_linear_combination(OpCombination<A1, A2, E> const& e) {
  return pack_left(e);
}

template <typename A1, typename A2, typename E,
          typename std::enable_if_t<
              (!expression_predicate<OpCombination<A1, A2, E>>::linear &&
               (!expression_predicate<E>::nonlinear &&
                expression_predicate<E>::combination)),
              int> = 0>
auto by_linear_combination(OpCombination<A1, A2, E> const& e) {
  if constexpr (expr::is_linear<A1>::value && expr::is_linear<A2>::value) {
    auto [l, r] = expr::split::by_linear(expr::get_enclosed_expression(e));
    return std::make_pair(e.combination * l, e.combination * r);
  } else {
    return pack_right(e);
  }
}

template <
    typename A1, typename A2, typename E,
    typename std::enable_if_t<expression_predicate<E>::nonlinear, int> = 0>
auto by_linear_combination(OpCombination<A1, A2, E> const& e) {
  return pack_right(e);
}
}  // namespace

template <typename A1, typename A2, typename E>
auto by_linear(OpCombination<A1, A2, E> const& e) {
  return by_linear_combination(e);
}

// op chain

namespace {
template <typename A1, typename A2, typename E,
          typename std::enable_if_t<
              (expression_predicate<OpChain<A1, A2, E>>::linear), int> = 0>
auto by_linear_chain(OpChain<A1, A2, E> const& e) {
  return pack_left(e);
}

template <typename A1, typename A2, typename E,
          typename std::enable_if_t<
              (!expression_predicate<OpChain<A1, A2, E>>::linear &&
               expression_predicate<E>::combination),
              int> = 0>
auto by_linear_chain(OpChain<A1, A2, E> const& e) {
  if constexpr (expr::is_linear<A1>::value && expr::is_linear<A2>::value) {
    auto [l, r] = expr::split::by_linear(expr::get_enclosed_expression(e));
    return std::make_pair(e.combination * l, e.combination * r);
  } else {
    return pack_right(e);
  }
}

template <typename A1, typename A2, typename E,
          typename std::enable_if_t<
              (expression_predicate<OpChain<A1, A2, E>>::nonlinear), int> = 0>
auto by_linear_chain(OpChain<A1, A2, E> const& e) {
  return pack_right(e);
}
}  // namespace

template <typename A1, typename A2, typename E>
auto by_linear(OpChain<A1, A2, E> const& e) {
  return by_linear_chain(e);
}

// convolution

namespace {

template <typename V, typename E1, typename E2,
          typename std::enable_if_t<
              expression_predicate<OpConvolution<V, E1, E2>>::linear, int> = 0>
auto by_linear_convolution(OpConvolution<V, E1, E2> const& e) {
  return pack_left(e);
}

template <typename V, typename E1, typename E2,
          typename std::enable_if_t<
              !expression_predicate<OpConvolution<V, E1, E2>>::linear, int> = 0>
auto by_linear_convolution(OpConvolution<V, E1, E2> const& e) {
  return pack_right(e);
}

}  // namespace

template <typename V, typename E1, typename E2>
auto by_linear(OpConvolution<V, E1, E2> const& e) {
  return by_linear_convolution(e);
}

// handling general expressions (recursion termination) and binary operations
// (recursive)

template <typename V, typename G>
auto by_linear(OpTerm<V, G> const& e) {
  return pack_left(e);
}

template <typename V, typename... Gs, expr::exp_key_t... Xs>
auto by_linear(OpTerms<V, Term<Gs, Xs>...> const& e) {
  using index_list_t = symphas::internal::select_unique_i_<
      expr::op_types_t<OpTerms<V, Term<Gs, Xs>...>>>;
  using term_list_t =
      symphas::lib::filter_types<symphas::lib::types_list<Gs...>, index_list_t>;

  constexpr int N =
      sum_factors<OpTerms<V, Term<Gs, Xs>...>, term_list_t>::value;
  if constexpr (N == 1) {
    return pack_left(e);
  } else {
    return pack_right(e);
  }
}

namespace {
template <typename... Es, size_t... Is>
auto by_linear_adds(OpAdd<Es...> const& e, std::index_sequence<Is...>) {
  return adds_expand_pair(by_linear(expr::get<Is>(e))...);
  // return std::make_pair((std::get<0>(std::get<Is>(a)) + ...),
  // (std::get<1>(std::get<Is>(a)) + ...));
}

// tries to find if there's an expression in the numerator that is linear
// through the denominator.
// **************************************************************************************

template <typename... E0s, typename... Ts>
auto filter_linear(std::tuple<E0s...> const& terms0, Ts const&... ts);

template <typename... E0s, size_t... Is>
auto filter_linear(std::tuple<E0s...> const& terms,
                   std::index_sequence<Is...>) {
  return std::make_tuple(std::get<Is>(terms)...);
}

template <typename... E0s, typename... E1s, typename... Ts, size_t... Is>
auto filter_linear(std::tuple<E0s...> const& terms0, std::index_sequence<Is...>,
                   std::tuple<E1s...> const& terms1, Ts const&... ts) {
  using filtered_indices_t = symphas::lib::seq_join_t<
      std::index_sequence<>,
      std::conditional_t<(symphas::lib::index_of_type<E0s, E1s...> >= 0 &&
                          !std::is_same<E0s, OpVoid>::value),
                         std::index_sequence<Is>, std::index_sequence<>>...>;

  return filter_linear(filter_linear(terms0, filtered_indices_t{}), ts...);
}

template <typename... E0s, typename... Ts>
auto filter_linear(std::tuple<E0s...> const& terms0, Ts const&... ts) {
  return filter_linear(terms0, std::make_index_sequence<sizeof...(E0s)>{},
                       ts...);
}

template <typename E>
auto get_linear_term(OpExpression<E> const& e) {
  auto&& [l, _] = expr::split::by_linear(*static_cast<E const*>(&e));
  return std::make_tuple(l);
}

template <typename E>
auto get_linear_term(OpOperator<E> const& e) {
  auto&& [l, _] = expr::split::by_linear(*static_cast<E const*>(&e));
  return std::make_tuple(l);
}

template <typename A, typename B>
auto get_linear_term(OpBinaryDiv<A, B> const& e) {
  return std::make_tuple(OpVoid{});
}

// \p a are the terms from the numerator, and \p e is the term from the
// denominator which is being checked that it can divide at least one of them.
template <typename E, typename A>
auto check_term_in_num(OpExpression<E> const& e, OpExpression<A> const& a,
                       std::index_sequence<>) {
  auto div = (*static_cast<A const*>(&a)) / (*static_cast<E const*>(&e));
  return get_linear_term(div);
}

// \p a are the terms from the numerator, and \p e is the term from the
// denominator which is being checked that it can divide at least one of them.
template <typename E, typename A>
auto check_term_in_num(OpOperator<E> const& e, OpExpression<A> const& a,
                       std::index_sequence<>) {
  return std::make_tuple(OpVoid{});
}

// \p a are the terms from the numerator, and \p e is the term from the
// denominator which is being checked that it can divide at least one of them.
template <typename E, typename A>
auto check_term_in_num(OpExpression<E> const& e, OpOperator<A> const& a,
                       std::index_sequence<>) {
  return std::make_tuple(OpVoid{});
}

// \p a are the terms from the numerator, and \p e is the term from the
// denominator which is being checked that it can divide at least one of them.
template <typename E, typename A>
auto check_term_in_num(OpOperator<E> const& e, OpOperator<A> const& a,
                       std::index_sequence<>) {
  return std::make_tuple(OpVoid{});
}

// \p a are the terms from the numerator, and \p e is the term from the
// denominator which is being checked that it can divide at least one of them.
template <typename E, typename... As, size_t... Is>
auto check_term_in_num(OpExpression<E> const& e, OpAdd<As...> const& a,
                       std::index_sequence<Is...>) {
  auto coeff =
      check_term_in_num(*static_cast<E const*>(&e), expr::get<sizeof...(Is)>(a),
                        std::index_sequence<>{});
  auto rest =
      check_term_in_num(*static_cast<E const*>(&e), (expr::get<Is>(a) + ...),
                        std::make_index_sequence<sizeof...(Is) - 1>{});

  if constexpr (std::is_same<std::tuple<OpVoid>, decltype(coeff)>::value) {
    return rest;
  } else {
    return std::tuple_cat(coeff, rest);
  }
}

// \p a are the terms from the numerator, and \p e is the term from the
// denominator which is being checked that it can divide at least one of them.
template <typename E, typename... As, size_t... Is>
auto check_term_in_num(OpOperator<E> const& e, OpAdd<As...> const& a,
                       std::index_sequence<Is...>) {
  return std::make_tuple(OpVoid{});
}

// \p a are the terms from the numerator, and \p e is the term from the
// denominator which is being checked that it can divide at least one of them.
template <typename... Es, typename... As, size_t... Js>
auto check_term_in_num(OpAdd<Es...> const& e, OpAdd<As...> const& a,
                       std::index_sequence<Js...>) {
  auto coeffs = filter_linear(check_term_in_num(
      expr::get<Js>(e), a, std::make_index_sequence<sizeof...(As) - 1>{})...);
  return coeffs;
}

template <typename A, typename... Bs>
auto by_linear_divs(OpBinaryDiv<A, OpAdd<Bs...>> const& e) {
  return std::tuple<OpVoid>{};
}

template <typename T0, typename... As, typename B>
auto by_linear_divs(OpBinaryDiv<OpAdd<As...>, B> const& e) {
  return check_term_in_num(e.b, e.a, std::make_index_sequence<sizeof...(As)>{});
}

template <typename... As, typename... Bs>
auto by_linear_divs(OpBinaryDiv<OpAdd<As...>, OpAdd<Bs...>> const& e) {
  return check_term_in_num(e.b, e.a, std::make_index_sequence<sizeof...(Bs)>{});
}

template <typename... Ls, size_t... Is>
auto add_linear_terms(std::tuple<Ls...> const& linear_terms,
                      std::index_sequence<Is...>) {
  return (std::get<Is>(linear_terms) + ... + OpVoid{});
}
}  // namespace

template <typename... Es>
auto by_linear(OpAdd<Es...> const& e) {
  return by_linear_adds(e, std::make_index_sequence<sizeof...(Es)>{});
}

template <typename A, typename B>
auto by_linear(OpBinaryDiv<A, B> const& e) {
  auto linear_terms = by_linear_divs(e);
  using seq_t =
      std::make_index_sequence<std::tuple_size<decltype(linear_terms)>::value>;

  auto l = add_linear_terms(linear_terms, seq_t{});
  return std::make_pair(l, e - l);
}

template <typename G, typename V, typename E>
auto by_linear(OpMap<G, V, E> const& e) {
  if constexpr (expr::is_combination<OpMap<G, V, E>>::value) {
    auto [l, r] = by_linear(expr::get_enclosed_expression(e));
    auto lm = expr::make_map<G>(expr::coeff(e) * l);
    auto rm = expr::make_map<G>(expr::coeff(e) * r);
    return std::make_pair(lm, rm);
  } else {
    return pack_right(e);
  }
}

template <typename V, typename sub_t, typename E, typename... Ts>
auto by_linear(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e) {
  auto&& [l, nl] = by_linear(e.f.e);
  auto fl = (expr::function_of(Ts{}...) = l);
  auto fnl = (expr::function_of(Ts{}...) = nl);
  return std::make_pair(
      symphas::internal::make_symbolic_eval(expr::coeff(e), e.data, fl),
      symphas::internal::make_symbolic_eval(expr::coeff(e), e.data, fnl));
}

template <typename V, typename E, typename... Ts, int... I0s, int... P0s,
          typename A, typename... As, typename B, typename C>
auto by_linear(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>,
                     symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
                     SymbolicFunction<A, As...>, B, C> const& series) {
  auto&& [l, nl] = by_linear(series.data.e);
  return std::make_pair(
      expr::coeff(series) * expr::recreate_series(l, series.data),
      expr::coeff(series) * expr::recreate_series(nl, series.data));
}

// **************************************************************************************

// there are several different cases a compound operator might happen:
// 1. E is a function of only V(Z)
// 2. E is a function of V(Z) and others, and contains a term that's only of
// V(Z)
// 3. E is a function of V(Z) and others, but contains no terms only of V(Z)
// 4. E is not a function of V(Z)

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename... Es>
auto separate_var(OpAdd<Es...> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename E1, typename E2>
auto separate_var(OpBinaryMul<E1, E2> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename E1, typename E2>
auto separate_var(OpBinaryDiv<E1, E2> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename Dd, typename V, typename E, typename Sp>
auto separate_var(OpDerivative<Dd, V, E, Sp> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename Dd, typename V, typename G, typename Sp>
auto separate_var(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename V, typename E, typename T>
auto separate_var(OpIntegral<V, E, T> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename A1, typename A2, typename E>
auto separate_var(OpChain<A1, A2, E> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename A1, typename A2, typename E>
auto separate_var(OpCombination<A1, A2, E> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename V, typename E1, typename E2>
auto separate_var(OpConvolution<V, E1, E2> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename V, size_t D, template <typename, size_t> typename grid_type, typename E>
auto separate_var(OpConvolution<V, GaussianSmoothing<D, grid_type>, E> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename G, typename V, typename E>
auto separate_var(OpMap<G, V, E> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, auto f, typename V, typename E>
auto separate_var(OpFunctionApply<f, V, E> const& e);

//! Separates the expressions based on existence of the variable index.
/*!
 * Separates terms of the expression based on whether the term contains
 * only variable terms of the prescribed index. It will form two
 * expressions, `A` and `B`, such that `A` is an expression containing
 * only terms of variable index `Z`, and `B` is everything else, which
 * may also contain variables of index `Z`. If a term can't be separated
 * by subtracting or adding, then it will not be separated and
 * remain in the expression `B`. An example would be multiplying
 * the variable index `Z` by another variable with a different index.
 *
 * \param e The expression which is split.
 *
 * \tparam Z The index of the variable to separate.
 */
template <size_t Z, typename V, typename E, typename F, typename Arg0,
          typename... Args>
auto separate_var(OpFunction<V, E, F, Arg0, Args...> const& e);

//

template <size_t Z, typename E>
auto separate_var(E const& e) {
  return pack_right(e);
}

template <size_t Z, typename T, typename G>
auto separate_var(OpTerm<T, Variable<Z, G>> const& e) {
  return pack_left(e);
}

namespace {
/*
 * Used in avoiding "if constexpr" constructs.
 */

template <size_t Z, typename E1, typename E2>
constexpr bool svc_pred = (expr::vars<E1>::template only_id<Z>() &&
                           expr::vars<E2>::template only_id<Z>());
template <size_t Z, typename E>
constexpr bool svcg_pred = (expr::vars<E>::template only_id<Z>());
template <size_t Z, typename E1, typename E2>
constexpr bool svm_pred = svc_pred<Z, E1, E2>;
template <size_t Z, typename E1, typename E2>
constexpr bool svd_pred =
    (!expr::vars<E1>::template only_id<Z>() &&
     expr::vars<E2>::template only_id<Z>() && !svm_pred<Z, E1, E2>);
template <size_t Z, typename E>
constexpr bool svd_pred_1 = (expr::vars<E>::template only_id<Z>());
template <size_t Z, typename E>
constexpr bool svd_pred_2 = (expression_predicate<E>::combination &&
                             expr::vars<E>::template has_id<Z>());
template <size_t Z, typename E>
constexpr bool svcc_pred_1 = (expr::vars<E>::template only_id<Z>());
template <size_t Z, typename E>
constexpr bool svcc_pred_2 = (expression_predicate<E>::combination &&
                              expr::vars<E>::template has_id<Z>());

}  // namespace

namespace {

template <size_t Z, typename V, typename E1, typename E2,
          typename std::enable_if_t<svc_pred<Z, E1, E2>, int> = 0>
auto separate_var_convolution(OpConvolution<V, E1, E2> const& e) {
  return pack_left(e);
}

template <size_t Z, typename V, typename E1, typename E2,
          typename std::enable_if_t<!svc_pred<Z, E1, E2>, int> = 0>
auto separate_var_convolution(OpConvolution<V, E1, E2> const& e) {
  return pack_right(e);
}
}  // namespace

template <size_t Z, typename V, typename E1, typename E2>
auto separate_var(OpConvolution<V, E1, E2> const& e) {
  return separate_var_convolution<Z>(e);
}

namespace {
template <size_t Z, typename V, size_t D, template <typename, size_t> typename grid_type, typename E,
          typename std::enable_if_t<svcg_pred<Z, E>, int> = 0>
auto separate_var_convolution_g(
    OpConvolution<V, GaussianSmoothing<D, grid_type>, E> const& e) {
  return pack_left(e);
}

template <size_t Z, typename V, size_t D, template <typename, size_t> typename grid_type, typename E,
          typename std::enable_if_t<!svcg_pred<Z, E>, int> = 0>
auto separate_var_convolution_g(
    OpConvolution<V, GaussianSmoothing<D, grid_type>, E> const& e) {
  return pack_right(e);
}
}  // namespace

template <size_t Z, typename V, size_t D, template <typename, size_t> typename grid_type, typename E>
auto separate_var(OpConvolution<V, GaussianSmoothing<D, grid_type>, E> const& e) {
  return separate_var_convolution_g<Z>(e);
}

namespace {
template <size_t Z, typename... Es, size_t... Is>
auto separate_var_adds(OpAdd<Es...> const& e, std::index_sequence<Is...>) {
  return adds_expand_pair(separate_var<Z>(expr::get<Is>(e))...);
  // return std::make_pair((std::get<0>(std::get<Is>(a)) + ...),
  // (std::get<1>(std::get<Is>(a)) + ...));
}
}  // namespace

template <size_t Z, typename... Es>
auto separate_var(OpAdd<Es...> const& e) {
  return separate_var_adds<Z>(e, std::make_index_sequence<sizeof...(Es)>{});
}

namespace {
template <size_t Z, typename E1, typename E2,
          typename std::enable_if_t<svm_pred<Z, E1, E2>, int> = 0>
auto separate_var_mul(OpBinaryMul<E1, E2> const& e) {
  return pack_left(e);
}

template <size_t Z, typename E1, typename E2,
          typename std::enable_if_t<!svm_pred<Z, E1, E2>, int> = 0>
auto separate_var_mul(OpBinaryMul<E1, E2> const& e) {
  return pack_right(e);
}
}  // namespace

template <size_t Z, typename E1, typename E2>
auto separate_var(OpBinaryMul<E1, E2> const& e) {
  return separate_var_mul<Z>(e);
}

namespace {
template <size_t Z, typename E1, typename E2,
          typename std::enable_if_t<svd_pred<Z, E1, E2>, int> = 0>
auto separate_var_div(OpBinaryDiv<E1, E2> const& e) {
  auto [a, b] = separate_var<Z>(e.a);
  return std::make_pair(a / e.b, b / e.b);
}

template <size_t Z, typename E1, typename E2,
          typename std::enable_if_t<!svm_pred<Z, E1, E2>, int> = 0>
auto separate_var_div(OpBinaryDiv<E1, E2> const& e) {
  return pack_right(e);
}

template <size_t Z, typename E1, typename E2,
          typename std::enable_if_t<svm_pred<Z, E1, E2>, int> = 0>
auto separate_var_div(OpBinaryDiv<E1, E2> const& e) {
  return pack_left(e);
}
}  // namespace

template <size_t Z, typename E1, typename E2>
auto separate_var(OpBinaryDiv<E1, E2> const& e) {
  return separate_var_div<Z>(e);
}

namespace {
template <size_t Z, typename Dd, typename V, typename E, typename T,
          typename std::enable_if_t<svd_pred_1<Z, E>, int> = 0>
auto separate_var_derivative(OpDerivative<Dd, V, E, T> const& e) {
  return pack_left(e);
}

template <
    size_t Z, typename Dd, typename V, typename E, typename T,
    typename std::enable_if_t<(!svd_pred_1<Z, E> && svd_pred_2<Z, E>), int> = 0>
auto separate_var_derivative(OpDerivative<Dd, V, E, T> const& e) {
  return separate_var<Z>(expr::apply_operators(e));
}

template <
    size_t Z, typename Dd, typename V, typename E, typename T,
    typename std::enable_if_t<!(svd_pred_1<Z, E> || svd_pred_2<Z, E>), int> = 0>
auto separate_var_derivative(OpDerivative<Dd, V, E, T> const& e) {
  return pack_right(e);
}
}  // namespace

template <size_t Z, typename Dd, typename V, typename E, typename Sp>
auto separate_var(OpDerivative<Dd, V, E, Sp> const& e) {
  return separate_var_derivative<Z>(e);
}

namespace {
template <size_t Z, typename Dd, typename V, typename G, typename T,
          typename std::enable_if_t<svd_pred_1<Z, G>, int> = 0>
auto separate_var_derivative_lop(
    OpDerivative<Dd, V, OpTerm<OpIdentity, G>, T> const& e) {
  return pack_left(e);
}

template <size_t Z, typename Dd, typename V, typename G, typename T,
          typename std::enable_if_t<!svd_pred_1<Z, G>, int> = 0>
auto separate_var_derivative_lop(
    OpDerivative<Dd, V, OpTerm<OpIdentity, G>, T> const& e) {
  return pack_right(e);
}
}  // namespace

template <size_t Z, typename Dd, typename V, typename G, typename Sp>
auto separate_var(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e) {
  return separate_var_derivative_lop<Z>(e);
}

namespace {
template <size_t Z, typename V, typename E, typename T,
          typename std::enable_if_t<svd_pred_1<Z, E>, int> = 0>
auto separate_var_integral(OpIntegral<V, E, T> const& e) {
  return pack_left(e);
}

template <size_t Z, typename V, typename E, typename T,
          typename std::enable_if_t<!svd_pred_1<Z, E>, int> = 0>
auto separate_var_integral(OpIntegral<V, E, T> const& e) {
  return pack_right(e);
}
}  // namespace

template <size_t Z, typename V, typename E, typename T>
auto separate_var(OpIntegral<V, E, T> const& e) {
  return separate_var_integral<Z>(e);
}

namespace {
template <size_t Z, typename A1, typename A2, typename E,
          typename std::enable_if_t<svcc_pred_1<Z, E>, int> = 0>
auto separate_var_chain(OpChain<A1, A2, E> const& e) {
  return pack_left(e);
}

template <size_t Z, typename A1, typename A2, typename E,
          typename std::enable_if_t<(!svcc_pred_1<Z, E> && svcc_pred_2<Z, E>),
                                    int> = 0>
auto separate_var_chain(OpChain<A1, A2, E> const& e) {
  return separate_var<Z>(expr::apply_operators(e));
}
}  // namespace

template <size_t Z, typename A1, typename A2, typename E,
          typename std::enable_if_t<!(svcc_pred_1<Z, E> || svcc_pred_2<Z, E>),
                                    int> = 0>
auto separate_var_chain(OpChain<A1, A2, E> const& e) {
  return pack_right(e);
}

template <size_t Z, typename A1, typename A2, typename E>
auto separate_var(OpChain<A1, A2, E> const& e) {
  return separate_var_chain<Z>(e);
}

namespace {
template <size_t Z, typename A1, typename A2, typename E,
          typename std::enable_if_t<svcc_pred_1<Z, E>, int> = 0>
auto separate_var_combination(OpCombination<A1, A2, E> const& e) {
  return pack_left(e);
}

template <size_t Z, typename A1, typename A2, typename E,
          typename std::enable_if_t<(!svcc_pred_1<Z, E> && svcc_pred_2<Z, E>),
                                    int> = 0>
auto separate_var_combination(OpCombination<A1, A2, E> const& e) {
  return separate_var<Z>(expr::apply_operators(e));
}

template <size_t Z, typename A1, typename A2, typename E,
          typename std::enable_if_t<!(svcc_pred_1<Z, E> || svcc_pred_2<Z, E>),
                                    int> = 0>
auto separate_var_combination(OpCombination<A1, A2, E> const& e) {
  return pack_right(e);
}
}  // namespace

template <size_t Z, typename A1, typename A2, typename E>
auto separate_var(OpCombination<A1, A2, E> const& e) {
  return separate_var_combination<Z>(e);
}

template <size_t Z, typename G, typename V, typename E>
auto separate_var(OpMap<G, V, E> const& e) {
  if constexpr (expr::is_combination<OpMap<G, V, E>>::value) {
    auto [l, r] = separate_var<Z>(expr::get_enclosed_expression(e));
    auto lm = expr::make_map<G>(expr::coeff(e) * l);
    auto rm = expr::make_map<G>(expr::coeff(e) * r);
    return std::make_pair(lm, rm);
  } else {
    if constexpr (svcg_pred<Z, E>) {
      return pack_left(e);
    } else {
      return pack_right(e);
    }
  }
}

template <size_t Z, expr::exp_key_t X, typename V, typename E>
auto separate_var(OpPow<X, V, E> const& e) {
  if constexpr (svcg_pred<Z, E>) {
    return pack_left(e);
  } else {
    return pack_right(e);
  }
}

template <size_t Z, auto f, typename V, typename E>
auto separate_var(OpFunctionApply<f, V, E> const& e) {
  if constexpr (svcg_pred<Z, E>) {
    return pack_left(e);
  } else {
    return pack_right(e);
  }
}

template <size_t Z, typename V, typename E, typename F, typename Arg0,
          typename... Args>
auto separate_var(OpFunction<V, E, F, Arg0, Args...> const& e) {
  if constexpr (svcg_pred<Z, E>) {
    return pack_left(e);
  } else {
    return pack_right(e);
  }
}

// **************************************************************************************

//! Separates an operator from the expression, if there is one applied to it.
/*!
 * The operator is returned in the first entry of the pair, and the expression
 * it is applied to is returned in the second. The original expression that is
 * passed to this function is recovered by applying the operator to the enclosed
 * expression. For example, given an expression U which is obtained by applying
 * an operator D to an enclosed expression E, then U = D(E). Note that in
 * general, it is NOT the case that U = D * U, since this is taking the dot
 * product.
 *
 * The return value of this function will be (D, E), so that U is recovered by
 * applying D to E.
 */
template <typename E>
auto separate_operator(OpExpression<E> const& e) {
  return pack_right(*static_cast<E const*>(&e));
}

template <typename E>
auto separate_operator(OpOperator<E> const& e) {
  return pack_left(*static_cast<E const*>(&e));
}

namespace {
template <typename Sp, size_t... Os, Axis... axs>
auto separate_mixed(Sp const& solver,
                    typename Sp::template mixed_derivative<Os...>,
                    symphas::internal::axis_list<axs...>) {
  return (OpIdentity{} * ... *
          expr::make_operator_directional_derivative<axs, Os>(solver));
}

template <typename Sp, size_t... Os>
auto separate_mixed(Sp const& solver,
                    typename Sp::template mixed_derivative<Os...>) {
  return separate_mixed<Sp>(solver,
                            typename Sp::template mixed_derivative<Os...>{},
                            symphas::lib::make_axis_list<sizeof...(Os)>());
}
}  // namespace

//! Separates an operator from the expression, if there is one applied to it.
/*!
 * The operator is returned in the first entry of the pair, and the expression
 * it is applied to is returned in the second. The original expression that is
 * passed to this function is recovered by applying the operator to the enclosed
 * expression. For example, given an expression U which is obtained by applying
 * an operator D to an enclosed expression E, then U = D(E). Note that in
 * general, it is NOT the case that U = D * U, since this is taking the dot
 * product.
 *
 * The return value of this function will be (D, E), so that U is recovered by
 * applying D to E.
 */
template <typename Dd, typename V, typename E, typename Sp,
          size_t R = expr::eval_type<E>::rank>
auto separate_operator(OpDerivative<Dd, V, E, Sp> const& e) {
  constexpr size_t order = OpDerivative<Dd, V, E, Sp>::order;
  constexpr Axis axis = OpDerivative<Dd, V, E, Sp>::axis;

  if constexpr (Dd::is_directional) {
    if constexpr (Dd::is_mixed) {
      auto op = separate_mixed<Sp>(e.solver, Dd{});
      return std::make_pair(op,
                            expr::coeff(e) * expr::get_enclosed_expression(e));
    } else {
      auto op = expr::make_operator_directional_derivative<Dd::axis, Dd::order>(
          e.solver);
      return std::make_pair(op,
                            expr::coeff(e) * expr::get_enclosed_expression(e));
    }
  } else {
    if constexpr (order % 2 == 1) {
      auto op_ax =
          expr::make_operator_directional_derivative<axis, 1>(e.solver);
      if constexpr (order == 1) {
        return std::make_pair(
            op_ax, expr::coeff(e) * expr::get_enclosed_expression(e));
      } else {
        auto op = op_ax * expr::make_operator_derivative<order - 1>(e.solver);
        return std::make_pair(
            op, expr::coeff(e) * expr::get_enclosed_expression(e));
      }
    } else {
      return std::make_pair(expr::make_operator_derivative<order>(e.solver),
                            expr::coeff(e) * expr::get_enclosed_expression(e));
    }
  }
}

//! Separates an operator from the expression, if there is one applied to it.
/*!
 * The operator is returned in the first entry of the pair, and the expression
 * it is applied to is returned in the second. The original expression that is
 * passed to this function is recovered by applying the operator to the enclosed
 * expression. For example, given an expression U which is obtained by applying
 * an operator D to an enclosed expression E, then U = D(E). Note that in
 * general, it is NOT the case that U = D * U, since this is taking the dot
 * product.
 *
 * The return value of this function will be (D, E), so that U is recovered by
 * applying D to E.
 */
template <Axis ax, size_t O, typename A1, typename A2, typename E,
          typename std::enable_if<
              min_derivative_order<OpOperatorCombination<A1, A2>>::value == O,
              int>::type = 0>
auto separate_operator(OpCombination<A1, A2, E> const& e) {
  return std::make_pair(e.combination, e.e);
}

//! Separates an operator from the expression, if there is one applied to it.
/*!
 * The operator is returned in the first entry of the pair, and the expression
 * it is applied to is returned in the second. The original expression that is
 * passed to this function is recovered by applying the operator to the enclosed
 * expression. For example, given an expression U which is obtained by applying
 * an operator D to an enclosed expression E, then U = D(E). Note that in
 * general, it is NOT the case that U = D * U, since this is taking the dot
 * product.
 *
 * The return value of this function will be (D, E), so that U is recovered by
 * applying D to E.
 */
template <typename A1, typename A2, typename E>
auto separate_operator(OpCombination<A1, A2, E> const& e) {
  return separate_operator<Axis::X, 0>(e);
}

//! Separates an operator from the expression, if there is one applied to it.
/*!
 * The operator is returned in the first entry of the pair, and the expression
 * it is applied to is returned in the second. The original expression that is
 * passed to this function is recovered by applying the operator to the enclosed
 * expression. For example, given an expression U which is obtained by applying
 * an operator D to an enclosed expression E, then U = D(E). Note that in
 * general, it is NOT the case that U = D * U, since this is taking the dot
 * product.
 *
 * The return value of this function will be (D, E), so that U is recovered by
 * applying D to E.
 */
template <Axis ax, size_t O, typename A1, typename A2, typename E,
          typename std::enable_if<
              min_derivative_order<OpOperatorChain<A1, A2>>::value == O,
              int>::type = 0>
auto separate_operator(OpChain<A1, A2, E> const& e) {
  return std::make_pair(e.combination, e.e);
}

//! Separates an operator from the expression, if there is one applied to it.
/*!
 * The operator is returned in the first entry of the pair, and the expression
 * it is applied to is returned in the second. The original expression that is
 * passed to this function is recovered by applying the operator to the enclosed
 * expression. For example, given an expression U which is obtained by applying
 * an operator D to an enclosed expression E, then U = D(E). Note that in
 * general, it is NOT the case that U = D * U, since this is taking the dot
 * product.
 *
 * The return value of this function will be (D, E), so that U is recovered by
 * applying D to E.
 */
template <typename A1, typename A2, typename E>
auto separate_operator(OpChain<A1, A2, E> const& e) {
  return separate_operator<Axis::X, 0>(e);
}

// **************************************************************************************

namespace {

// factors the expression by the given variable
// the OpTerms and opnlvariables with matching type C are factored
// N is the number of types of C to factor out

template <typename C>
struct divide_variable {
  template <typename T, typename G,
            typename std::enable_if<(expr::factor_count<C, G>::value > 0),
                                    int>::type = 0>
  auto operator()(OpTerm<T, G> const& e) {
    return std::make_pair(OpTerm<OpIdentity, G>(e.data), expr::coeff(e));
  }

  template <typename T, typename G,
            typename std::enable_if<(expr::factor_count<C, G>::value == 0),
                                    int>::type = 0>
  auto operator()(OpTerm<T, G> const& e) {
    return std::make_pair(OpIdentity{}, e);
  }
};

template <size_t Z>
struct divide_variable<Variable<Z>> {
  template <typename T, typename G>
  auto operator()(OpTerm<T, Variable<Z, G>> const& e) {
    return std::make_pair(OpTerm<OpIdentity, Variable<Z, G>>(e.data),
                          expr::coeff(e));
  }

  template <typename T, typename G>
  auto operator()(OpTerm<T, G> const& e) {
    return std::make_pair(OpIdentity{}, e);
  }
};

template <size_t N, typename C, typename... Es,
          typename Enable = typename std::enable_if<
              (expr::factor_count<C, OpAdd<Es...>>::value >= N)>::type>
auto _factor(OpAdd<Es...> const& e);
template <size_t N, typename C, typename E1, typename E2,
          typename Enable = typename std::enable_if<
              (expr::factor_count<C, OpBinaryMul<E1, E2>>::value >= N)>::type>
auto _factor(OpBinaryMul<E1, E2> const& e);
template <size_t N, typename C, typename E1, typename E2,
          typename Enable = typename std::enable_if<
              (expr::factor_count<C, OpBinaryDiv<E1, E2>>::value >= N)>::type>
auto _factor(OpBinaryDiv<E1, E2> const& e);
template <size_t N, typename C, typename V, typename... Gs, exp_key_t... Xs>
auto _factor(OpTerms<V, Term<Gs, Xs>...> const& e);

template <size_t N, typename C, typename E,
          typename std::enable_if<
              (N == 0 || expr::factor_count<C, E>::value < N), int>::type = 0>
auto _factor(OpExpression<E> const& e) {
  return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
}

template <size_t N, typename C, typename E>
auto _factor(OpOperator<E> const& e) {
  return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
}

// Helper functions for factoring.
// **************************************************************************************

template <size_t N, typename C, typename... Es, size_t... Is>
auto _factor_adds(OpAdd<Es...> const& e, std::index_sequence<Is...>) {
  return adds_expand_pair_no_first(_factor<N, C>(expr::get<Is>(e))...);
  // return std::make_pair(std::get<0>(a).first,
  // expr::add_all(std::get<Is>(a).second...));
}

template <typename V, typename G0, exp_key_t X0, typename... Gs,
          exp_key_t... Xs>
auto _make_terms(V value, Term<G0, X0> const& term0,
                 Term<Gs, Xs> const&... rest) {
  return OpTerms(value, term0, rest...);
}

template <typename V, typename G0, typename... Gs, exp_key_t... Xs>
auto _make_terms(V value, Term<G0, 0> const& term0,
                 Term<Gs, Xs> const&... rest) {
  return expr::transform::sift_term(value, rest...);
}

template <size_t N, typename C, typename V, typename... Gs, exp_key_t... Xs,
          size_t... Is>
auto _select_terms(OpTerms<V, Term<Gs, Xs>...> const& e,
                   std::index_sequence<Is...>, std::index_sequence<>) {
  return std::pair(OpIdentity{}, e);
}

template <size_t N, typename C, typename V, typename... Gs, exp_key_t... Xs,
          size_t... Is, size_t I>
auto _select_terms(OpTerms<V, Term<Gs, Xs>...> const& e,
                   std::index_sequence<Is...>, std::index_sequence<I>) {
  // If the type can factor another type more than once (e.g. |k|^2 factors
  // |k|^4 twice), then if we want to factor |k|^2 from |k|^4, we need to divide
  // |k|^4 by |k|^2.
  constexpr size_t N0 =
      expr::factor_count<
          C, symphas::lib::direct_type_at_index<I - 1, Gs...>>::value -
      1;

  auto factor_data = expr::get<I>(e).data();
  auto factor_term = Term(factor_data);
  auto factor_eval =
      (factor_term * ~(factor_term.template pow<N0>())).template pow<N>();
  auto nonfactor_term = expr::get<I>(e) * (~factor_eval);

  return std::make_pair(
      _make_terms(OpIdentity{}, factor_term),
      _make_terms(e.term, nonfactor_term, expr::get<Is>(e)...));
}

template <size_t N, typename C, typename V, typename... Gs, exp_key_t... Xs,
          size_t... Is, bool... fs>
auto _factor(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<Is...>,
             std::integer_sequence<bool, fs...>) {
  using symphas::lib::seq_join_t;
  using pick_factor_t =
      seq_join_t<std::index_sequence<>,
                 std::conditional_t<fs, std::index_sequence<Is + 1>,
                                    std::index_sequence<>>...>;
  using pick_nonfactor_t =
      seq_join_t<std::index_sequence<>,
                 std::conditional_t<!fs, std::index_sequence<Is + 1>,
                                    std::index_sequence<>>...>;
  return _select_terms<N, C>(e, pick_nonfactor_t{}, pick_factor_t{});
}

template <size_t N, typename C, typename V, typename... Gs, exp_key_t... Xs>
auto _factor(OpTerms<V, Term<Gs, Xs>...> const& e) {
  using seq_t = std::make_index_sequence<sizeof...(Gs)>;
  using mask_t =
      std::integer_sequence<bool, (factor_count<C, Term<Gs, Xs>>::value >= N &&
                                   expr::is_combinable<Gs>)...>;
  return _factor<N, C>(e, seq_t{}, mask_t{});
}

template <size_t N, typename C, typename... Es, typename Enable>
auto _factor(OpAdd<Es...> const& e) {
  return _factor_adds<N, C>(e, std::make_index_sequence<sizeof...(Es)>{});
}

template <size_t N, typename C, typename E1, typename E2,
          size_t AN = expr::factor_count<C, E1>::value,
          typename std::enable_if_t<(N > AN), int> = 0>
auto _factor_sift(OpBinaryMul<E1, E2> const& e) {
  auto a = _factor<AN, C>(e.a);
  auto b = _factor<N - AN, C>(e.b);
  return std::make_pair(a.first * b.first, a.second * b.second);
}

template <size_t N, typename C, typename E1, typename E2,
          size_t AN = expr::factor_count<C, E1>::value,
          typename std::enable_if_t<(N <= AN), int> = 0>
auto _factor_sift(OpBinaryMul<E1, E2> const& e) {
  auto a = _factor<N, C>(e.a);
  return std::make_pair(a.first, a.second * e.b);
}

template <size_t N, typename C, typename E1, typename E2, typename Enable>
auto _factor(OpBinaryMul<E1, E2> const& e) {
  /* this algorithm takes into account the edge case that the same type of
   * variable is separated by a strict multiplication, although in general, a
   * multiplication rule will put them together into an nlvariable
   */
  return _factor_sift<N, C>(e);
}

template <size_t N, typename C, typename E1, typename E2, typename Enable>
auto _factor(OpBinaryDiv<E1, E2> const& e) {
  auto a = _factor<N, C>(e.a);
  return std::make_pair(a.first, a.second / e.b);
}

template <typename C>
struct is_variable_data_factor {
  static const bool value = false;
};

template <size_t Z>
struct is_variable_data_factor<Variable<Z>> {
  static const bool value = true;
};

template <typename C0, typename... Cs>
struct factor_pred {
  static const bool value = ((std::is_same<C0, Cs>::value || ...));
};

template <typename C0>
struct factor_pred<C0> {
  static const bool value = false;
};

// automatically selects the largest possible number of datas to factor
template <typename C0, typename E,
          typename std::enable_if<(expr::grid_can_combine<C0>::value ||
                                   is_variable_data_factor<C0>::value),
                                  int>::type = 0>
auto _factor(OpExpression<E> const& e) {
  constexpr size_t N = expr::factor_count<C0, E>::value;
  if constexpr (N > 0) {
    return _factor<N, C0>(*static_cast<E const*>(&e));
  } else {
    return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
  }
}

template <typename C0, typename E,
          typename std::enable_if<(!expr::grid_can_combine<C0>::value &&
                                   !is_variable_data_factor<C0>::value),
                                  int>::type = 0>
auto _factor(OpExpression<E> const& e) {
  return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
}

template <typename C0, typename E>
auto _factor(OpOperator<E> const& e) {
  return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
}
}  // namespace

//! Factor an expression by the given terms.
/*!
 * Factor an expression by the given terms which are represented by
 * the given types, not passed explicitly.
 * Returns the result of factoring the given term in the first entry of the
 * pair, and the product of all the factored terms in the second entry of the
 * pair.
 *
 * The types to factor by must all be unique. This produces a maximally
 * factored expression where all the types have been potentially factored.
 *
 * \param e The expression to factor.
 *
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <typename C0, typename... Cs, typename E,
          typename std::enable_if_t<
              (sizeof...(Cs) == 0 && !factor_pred<C0>::value), int> = 0>
auto factor(OpExpression<E> const& e) {
  auto a = _factor<C0>(*static_cast<const E*>(&e));
  return std::make_pair(a.first, a.second);
}

//! Factor an expression by the given terms.
/*!
 * Factor an expression by the given terms which are represented by
 * the given types, not passed explicitly.
 * Returns the result of factoring the given term in the first entry of the
 * pair, and the product of all the factored terms in the second entry of the
 * pair.
 *
 * The types to factor by must all be unique. This produces a maximally
 * factored expression where all the types have been potentially factored.
 *
 * \param e The expression to factor.
 *
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <typename C0, typename... Cs, typename E,
          typename std::enable_if_t<
              (sizeof...(Cs) > 0 && !factor_pred<C0, Cs...>::value), int> = 0>
auto factor(OpExpression<E> const& e) {
  auto a = _factor<C0>(*static_cast<const E*>(&e));
  auto b = factor<Cs...>(a.second);
  return std::make_pair(a.first * b.first, b.second);
}

//! Factor an expression by the given terms.
/*!
 * Factor an expression by the given terms which are represented by
 * the given types, not passed explicitly.
 * Returns the result of factoring the given term in the first entry of the
 * pair, and the product of all the factored terms in the second entry of the
 * pair.
 *
 * The types to factor by must all be unique. This produces a maximally
 * factored expression where all the types have been potentially factored.
 *
 * \param e The expression to factor.
 *
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <typename C0, typename... Cs, typename E>
auto factor(OpOperator<E> const& e) {
  return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
}

//! Factor an expression by the given terms the given number of times.
/*!
 * Factor an expression by the given terms which are represented by
 * the given types, not passed explicitly. Each term is factored out
 * the given number of times.
 * Returns the result of factoring the given term in the first entry of the
 * pair, and the product of all the factored terms in the second entry of the
 * pair.
 *
 * The types to factor by must all be unique.
 *
 * \param e The expression to factor.
 *
 * \tparam N The number of times to factor each type out.
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <size_t N, typename C0, typename E,
          typename std::enable_if_t<(!factor_pred<C0>::value), int> = 0>
auto factor(OpExpression<E> const& e) {
  constexpr size_t min_order = fixed_min<expr::factor_count<C0, E>::value, N>;
  if constexpr (min_order > 0) {
    return _factor<min_order, C0>(*static_cast<const E*>(&e));
  } else {
    return std::make_pair(OpIdentity{}, *static_cast<const E*>(&e));
  }
}

//! Factor an expression by the given terms the given number of times.
/*!
 * Factor an expression by the given terms which are represented by
 * the given types, not passed explicitly. Each term is factored out
 * the given number of times.
 * Returns the result of factoring the given term in the first entry of the
 * pair, and the product of all the factored terms in the second entry of the
 * pair.
 *
 * The types to factor by must all be unique.
 *
 * \param e The expression to factor.
 *
 * \tparam N The number of times to factor each type out.
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <
    size_t N, typename C0, typename C1, typename... Cs, typename E,
    typename std::enable_if_t<(!factor_pred<C0, C1, Cs...>::value), int> = 0>
auto factor(OpExpression<E> const& e) {
  constexpr size_t min_order = fixed_min<expr::factor_count<C0, E>::value, N>;
  if constexpr (min_order > 0) {
    auto a = _factor<min_order, C0>(*static_cast<const E*>(&e));
    auto b = factor<N, C1, Cs...>(a.second);
    return std::make_pair(a.first * b.first, b.second);
  } else {
    return factor<N, C1, Cs...>(*static_cast<const E*>(&e));
  }
}

//! Factor an expression by the given terms the given number of times.
/*!
 * Factor an expression by the given terms which are represented by
 * the given types, not passed explicitly. Each term is factored out
 * the given number of times.
 * Returns the result of factoring the given term in the first entry of the
 * pair, and the product of all the factored terms in the second entry of the
 * pair.
 *
 * The types to factor by must all be unique.
 *
 * \param e The expression to factor.
 *
 * \tparam N The number of times to factor each type out.
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <size_t N, typename C0, typename C1, typename... Cs, typename E>
auto factor(OpOperator<E> const& e) {
  return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
}

//! Factor an expression by the given terms the given number of times.
/*!
 * Factor an expression by the given terms which are represented by
 * the given types, not passed explicitly. Each term is factored out
 * the given number of times.
 * Returns the result of factoring the given term in the first entry of the
 * pair, and the product of all the factored terms in the second entry of the
 * pair.
 *
 * The types to factor by must all be unique.
 *
 * \param e The expression to factor.
 *
 * \tparam N The number of times to factor each type out.
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <size_t N, size_t Z0, size_t... Zs, typename E,
          typename std::enable_if_t<
              (!factor_pred<Variable<Z0>, Variable<Zs>...>::value), int> = 0>
auto factor(OpExpression<E> const& e) {
  return factor<N, Variable<Z0>, Variable<Zs>...>(*static_cast<E const*>(&e));
}

//! Factor an expression by the given terms the given number of times.
/*!
 * Factor an expression by the given terms which are represented by
 * the given types, not passed explicitly. Each term is factored out
 * the given number of times.
 * Returns the result of factoring the given term in the first entry of the
 * pair, and the product of all the factored terms in the second entry of the
 * pair.
 *
 * The types to factor by must all be unique.
 *
 * \param e The expression to factor.
 *
 * \tparam N The number of times to factor each type out.
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <size_t N, size_t Z0, size_t... Zs, typename E>
auto factor(OpOperator<E> const& e) {
  return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
}

//! Make a list of the result of each factorization.
/*!
 * For each given term, the original given expression is factored. The
 * results are concatenated into a list.
 * Returns the result of factoring the given term in the first entry of the
 * pair, and the product of all the factored terms in the second entry of the
 * pair.
 *
 * \param e The expression to factor.
 *
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <typename C0, typename... Cs, typename E,
          typename std::enable_if_t<
              (sizeof...(Cs) == 0 && !factor_pred<C0>::value), int> = 0>
auto factor_tuple(OpExpression<E> const& e) {
  return std::make_tuple(factor<C0>(*static_cast<const E*>(&e)));
}

//! Make a list of the result of each factorization.
/*!
 * For each given term, the original given expression is factored. The
 * results are concatenated into a list.
 * Returns the result of factoring the given term in the first entry of the
 * pair, and the product of all the factored terms in the second entry of the
 * pair.
 *
 * \param e The expression to factor.
 *
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <typename C0, typename... Cs, typename E,
          typename std::enable_if_t<
              (sizeof...(Cs) > 0 && !factor_pred<C0, Cs...>::value), int> = 0>
auto factor_tuple(OpExpression<E> const& e) {
  auto a = factor<C0>(*static_cast<const E*>(&e));
  auto b = factor_tuple<Cs...>(*static_cast<const E*>(&e));
  return std::tuple_cat(a, b);
}

//! Make a list of the result of each factorization.
/*!
 * For each given term, the original given expression is factored. The
 * results are concatenated into a list.
 * Returns the result of factoring the given term in the first entry of the
 * pair, and the product of all the factored terms in the second entry of the
 * pair.
 *
 * \param e The expression to factor.
 *
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <typename C0, typename... Cs, typename E>
auto factor_tuple(OpOperator<E> const& e) {
  auto a = factor<C0>(*static_cast<const E*>(&e));
  auto b = factor_tuple<Cs...>(*static_cast<const E*>(&e));
  return std::tuple_cat(a, b);
}

//! Keep factoring the given expression once by all terms.
/*!
 * The expression is factored only once by each given term.
 *
 * All the factors that are removed from the expression are in the first entry
 * of the returned pair, and the factored expression is given in the
 * second entry.
 *
 * That is: Given a list of factors and an expression *E*, this function returns
 * *A*, *B*:
 * - *A* is a product of the factors that have been taken from *E*.
 * - *B* is the factored expression, such that *E* = *A* * *B*.
 *
 * \param e The expression to factor.
 *
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <typename C0, typename... Cs, typename E,
          typename std::enable_if_t<(sizeof...(Cs) == 0), int> = 0>
auto factor_list(OpExpression<E> const& e) {
  constexpr size_t min_order = fixed_min<expr::factor_count<C0, E>::value, 1>;
  auto a = _factor<min_order, C0>(*static_cast<const E*>(&e));
  return std::make_pair(a.first, a.second);
}

//! Keep factoring the given expression once by all terms.
/*!
 * The expression is factored only once by each given term.
 *
 * All the factors that are removed from the expression are in the first entry
 * of the returned pair, and the factored expression is given in the
 * second entry.
 *
 * That is: Given a list of factors and an expression *E*, this function returns
 * *A*, *B*:
 * - *A* is a product of the factors that have been taken from *E*.
 * - *B* is the factored expression, such that *E* = *A* * *B*.
 *
 * \param e The expression to factor.
 *
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <typename C0, typename... Cs, typename E,
          typename std::enable_if_t<(sizeof...(Cs) > 0), int> = 0>
auto factor_list(OpExpression<E> const& e) {
  constexpr size_t min_order = fixed_min<expr::factor_count<C0, E>::value, 1>;

  auto a = _factor<min_order, C0>(*static_cast<const E*>(&e));
  auto b = factor_list<Cs...>(a.second);
  return std::make_pair(a.first * b.first, b.second);
}

//! Keep factoring the given expression once by all terms.
/*!
 * The expression is factored only once by each given term.
 *
 * All the factors that are removed from the expression are in the first entry
 * of the returned pair, and the factored expression is given in the
 * second entry.
 *
 * That is: Given a list of factors and an expression *E*, this function returns
 * *A*, *B*:
 * - *A* is a product of the factors that have been taken from *E*.
 * - *B* is the factored expression, such that *E* = *A* * *B*.
 *
 * \param e The expression to factor.
 *
 * \tparam Z0 The first variable index to factor.
 * \tparam Cs... The remaining variable indices to factor.
 */
template <size_t Z0, size_t... Zs, typename E,
          typename std::enable_if_t<(sizeof...(Zs) == 0), int> = 0>
auto factor_list(OpExpression<E> const& e) {
  return factor_list<Variable<Z0>, Variable<Zs>...>(*static_cast<const E*>(&e));
}

//! Keep factoring the given expression once by all terms.
/*!
 * The expression is factored only once by each given term.
 *
 * All the factors that are removed from the expression are in the first entry
 * of the returned pair, and the factored expression is given in the
 * second entry.
 *
 * That is: Given a list of factors and an expression *E*, this function returns
 * *A*, *B*:
 * - *A* is a product of the factors that have been taken from *E*.
 * - *B* is the factored expression, such that *E* = *A* * *B*.
 *
 * \param e The expression to factor.
 *
 * \tparam C0 The first term to factor.
 * \tparam Cs... The remaining terms to factor.
 */
template <typename C0, typename... Cs, typename E>
auto factor_list(OpOperator<E> const& e) {
  return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
}

//! Keep factoring the given expression once by all terms.
/*!
 * The expression is factored only once by each given term.
 *
 * All the factors that are removed from the expression are in the first entry
 * of the returned pair, and the factored expression is given in the
 * second entry.
 *
 * That is: Given a list of factors and an expression *E*, this function returns
 * *A*, *B*:
 * - *A* is a product of the factors that have been taken from *E*.
 * - *B* is the factored expression, such that *E* = *A* * *B*.
 *
 * \param e The expression to factor.
 *
 * \tparam Z0 The first variable index to factor.
 * \tparam Cs... The remaining variable indices to factor.
 */
template <size_t Z0, size_t... Zs, typename E>
auto factor_list(OpOperator<E> const& e) {
  return factor_list<Variable<Z0>, Variable<Zs>...>(*static_cast<const E*>(&e));
}

template <typename... condition_ts, typename E>
auto separate_by(OpExpression<E> const& e) {
  using namespace symphas::internal;
  if constexpr ((expression_satisfies_condition<E, condition_ts> || ... ||
                 false)) {
    return pack_left(*static_cast<E const*>(&e));
  } else {
    return pack_right(*static_cast<E const*>(&e));
  }
}

template <typename... condition_ts, typename... Es, size_t... Is>
auto separate_by(OpAdd<Es...> const& e, std::index_sequence<Is...>) {
  return adds_expand_pair(separate_by<condition_ts...>(expr::get<Is>(e))...);
}

template <typename... condition_ts, typename... Es>
auto separate_by(OpAdd<Es...> const& e) {
  return separate_by<condition_ts...>(
      e, std::make_index_sequence<sizeof...(Es)>{});
}

template <typename... condition_ts, typename E>
auto filter(OpExpression<E> const& e) {
  using namespace symphas::internal;
  if constexpr ((expression_satisfies_condition<E, condition_ts> || ... ||
                 false)) {
    return *static_cast<E const*>(&e);
  } else {
    return OpVoid{};
  }
}

template <typename... condition_ts, typename... Es, size_t... Is>
auto filter(OpAdd<Es...> const& e, std::index_sequence<Is...>) {
  return (expr::get<Is>(e) + ... + OpVoid{});
}

template <typename... condition_ts, typename... Es>
auto filter(OpAdd<Es...> const& e) {
  using namespace symphas::internal;
  using mask = std::integer_sequence<
      bool, expression_satisfies_condition<Es, expr::or_<condition_ts...>>...>;
  using seq = std::make_index_sequence<sizeof...(Es)>;
  return filter<condition_ts...>(e, symphas::lib::filter_seq_t<seq, mask>{});
}
}  // namespace expr::split

namespace expr {
template <size_t O, typename factor_t>
template <typename E>
auto compile_trait::factor_trait<O, factor_t>::operator()(E&& e) {
  return expr::split::factor<O, factor_t>(std::forward<E>(e));
}

template <typename E>
auto compile_trait::separate_operator_trait::operator()(E&& e) {
  return expr::split::separate_operator(std::forward<E>(e));
}

template <typename condition_t>
template <typename E>
auto compile_trait::separate_by_trait<condition_t>::operator()(E&& e) {
  return expr::split::separate_by<condition_t>(std::forward<E>(e));
}

}  // namespace expr
