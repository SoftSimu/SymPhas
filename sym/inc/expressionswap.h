
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

#include "expressionsplit.h"

namespace expr::transform {

namespace {

//! Used for the grid swapping routines; swaps all indicies with the given
//! replacement.
/*!
 * Defined using the expr::split::factor routine, implemented below. All the
 * matching indices are factored from the expression, and then the replacements
 * are combined into a single term and substituted with the replacement.
 */
template <typename E, int N0, int... P0s, typename G_F>
auto swap_matching_i(OpExpression<E> const& e,
                     symphas::lib::types_list<expr::symbols::i_<N0, P0s>...>,
                     G_F&& g);

template <int N0, int... P0s, expr::exp_key_t... Xs, typename G_F>
auto swap_matching_i(
    OpTerms<OpIdentity, Term<expr::symbols::i_<N0, P0s>, Xs>...> const& e,
    G_F&& g);

template <typename G_F>
auto swap_matching_i(OpIdentity, G_F&& g) {
  return OpIdentity{};
}

}  // namespace

// **************************************************************************************

/* algorithm that takes the id of a variable and swaps all instances of the
 * variable a different grid accounts for both being given a grid and being
 * given a usual variable; doesn't handle swapping an expression in place of a
 * variable
 */

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename E, typename G_F>
auto swap_grid(OpExpression<E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename E, typename G_F>
auto swap_grid(OpOperator<E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Implementation of a successful search, where the given variable term
 * associated with the prescribed index will be switched with the given
 * expression.
 *
 * \param v The term which is swapped.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename V, typename... Gs, exp_key_t... Xs, typename G_F>
auto swap_grid(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename... Es, typename G_F>
auto swap_grid(OpAdd<Es...> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename E1, typename E2, typename G_F>
auto swap_grid(OpBinaryMul<E1, E2> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename E1, typename E2, typename G_F>
auto swap_grid(OpBinaryDiv<E1, E2> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, size_t O, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorDerivative<O, V, Sp> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, Axis ax, size_t O, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, size_t... Os, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorMixedDerivative<V, Sp, Os...> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename A1, typename A2, typename G_F>
auto swap_grid(OpOperatorChain<A1, A2> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename A1, typename A2, typename E, typename G_F>
auto swap_grid(OpChain<A1, A2, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename A1, typename E, typename G_F>
auto swap_grid(OpChain<A1, OpIdentity, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename A1, typename A2, typename G_F>
auto swap_grid(OpOperatorCombination<A1, A2> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename A1, typename A2, typename E, typename G_F>
auto swap_grid(OpCombination<A1, A2, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename V, typename E1, typename E2, typename G_F>
auto swap_grid(OpConvolution<V, E1, E2> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename V, size_t D, typename E, typename G_F>
auto swap_grid(OpConvolution<V, GaussianSmoothing<D>, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename Dd, typename V, typename E, typename Sp,
          typename G_F>
auto swap_grid(OpDerivative<Dd, V, E, Sp> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, size_t O, typename V, typename E, typename GG, typename G_F>
auto swap_grid(
    OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& e,
    G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename V, typename E, typename T, typename G_F>
auto swap_grid(OpIntegral<V, E, T> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename G, typename V, typename E, typename G_F>
auto swap_grid(OpMap<G, V, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename G, typename E, typename G_F>
auto swap_grid(OpMap<G, OpIdentity, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, expr::exp_key_t X, typename V, typename E, typename G_F>
auto swap_grid(OpPow<X, V, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, typename V, typename E, typename F, typename Arg0,
          typename... Args, typename G_F>
auto swap_grid(OpFunction<V, E, F, Arg0, Args...> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \tparam Z The index of the variable to change.
 */
template <size_t Z, auto f, typename V, typename E, typename G_F>
auto swap_grid(OpFunctionApply<f, V, E> const& e, G_F&& g);

template <size_t Z, typename V, typename sub_t, typename E, typename... Ts,
          typename G_F>
auto swap_grid(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e,
               G_F&& g);

template <size_t Z, int N, int P, typename G_F>
auto swap_grid(expr::symbols::i_<N, P> const& e, G_F&& g) {
  return expr::symbols::i_<N, P>{};
}

template <size_t Z, typename T, typename G_F>
auto swap_grid(SymbolicData<T> const& e, G_F&& g) {
  return e;
}

namespace {

template <typename V, typename... Gs, exp_key_t... Xs, typename... Ts,
          size_t I0, size_t... I0s>
auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e,
                     std::index_sequence<I0, I0s...>,
                     std::tuple<Ts...> const& subbed, std::index_sequence<>);
template <typename V, typename... Gs, exp_key_t... Xs, typename... Ts,
          size_t I1, size_t... I1s>
auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e,
                     std::index_sequence<>, std::tuple<Ts...> const& subbed,
                     std::index_sequence<I1, I1s...>);
template <typename V, typename... Gs, exp_key_t... Xs, typename... Ts>
auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e,
                     std::index_sequence<>, std::tuple<Ts...> const& subbed,
                     std::index_sequence<>);
template <typename V, typename... Gs, exp_key_t... Xs, typename... Ts,
          size_t I0, size_t... I0s, size_t I1, size_t... I1s>
auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e,
                     std::index_sequence<I0, I0s...>,
                     std::tuple<Ts...> const& subbed,
                     std::index_sequence<I1, I1s...>);

template <typename V, typename... Gs, exp_key_t... Xs, typename... Ts,
          size_t I0, size_t... I0s>
auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e,
                     std::index_sequence<I0, I0s...>,
                     std::tuple<Ts...> const& subbed, std::index_sequence<>) {
  return expr::make_term(expr::get<I0 + 1>(e)) *
         recombine_terms(e, std::index_sequence<I0s...>{}, subbed,
                         std::index_sequence<>{});
}

template <typename V, typename... Gs, exp_key_t... Xs, typename... Ts,
          size_t I1, size_t... I1s>
auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e,
                     std::index_sequence<>, std::tuple<Ts...> const& subbed,
                     std::index_sequence<I1, I1s...>) {
  return std::get<I1>(subbed) * recombine_terms(e, std::index_sequence<>{},
                                                subbed,
                                                std::index_sequence<I1s...>{});
}

template <typename V, typename... Gs, exp_key_t... Xs, typename... Ts>
auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e,
                     std::index_sequence<>, std::tuple<Ts...> const& subbed,
                     std::index_sequence<>) {
  return OpIdentity{};
}

template <typename V, typename... Gs, exp_key_t... Xs, typename... Ts,
          size_t I0, size_t... I0s, size_t I1, size_t... I1s>
auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e,
                     std::index_sequence<I0, I0s...>,
                     std::tuple<Ts...> const& subbed,
                     std::index_sequence<I1, I1s...>) {
  if constexpr (I0 < I1) {
    return expr::make_term(expr::get<I0 + 1>(e)) *
           recombine_terms(e, std::index_sequence<I0s...>{}, subbed,
                           std::index_sequence<I1, I1s...>{});
  } else {
    return std::get<I1>(subbed) *
           recombine_terms(e, std::index_sequence<I0, I0s...>{}, subbed,
                           std::index_sequence<I1s...>{});
  }
}

template <typename V, typename... G1s, exp_key_t... X1s, size_t... Ns,
          size_t... Ms, exp_key_t X, typename G_F>
auto _pick_terms(OpTerms<V, Term<G1s, X1s>...> const& a,
                 std::index_sequence<Ns...>, std::index_sequence<Ms...>,
                 std::integer_sequence<exp_key_t, X>, G_F&& g) {
  if constexpr (!expr::is_expression<G_F>) {
    if constexpr (symphas::is_simple_data<G_F>) {
      return sift_term(
          sift_term(g), expr::get<Ns>(a)...,
          Term(expr::make_literal(std::forward<G_F>(g))).template pow<X>(),
          expr::get<Ms>(a)...);
    } else {
      return sift_term(OpIdentity{}, expr::get<Ns>(a)...,
                       Term(std::forward<G_F>(g)).template pow<X>(),
                       expr::get<Ms>(a)...);
    }
  } else {
    if constexpr (_Xk_t<X>::D == 1) {
      if constexpr (_Xk_t<X>::N == 1) {
        return sift_term(OpIdentity{}, expr::get<Ns>(a)...) *
               std::forward<G_F>(g) *
               sift_term(OpIdentity{}, expr::get<Ms>(a)...);
      } else {
        return sift_term(OpIdentity{}, expr::get<Ns>(a)...) *
               expr::pow<_Xk_t<X>::N>(std::forward<G_F>(g)) *
               sift_term(OpIdentity{}, expr::get<Ms>(a)...);
      }
    } else {
      return sift_term(OpIdentity{}, expr::get<Ns>(a)...) *
             expr::exp(((_Xk_t<X>::sign) ? OpNegIdentity{} : OpIdentity{}) *
                       expr::make_fraction<_Xk_t<X>::N, _Xk_t<X>::D>() *
                       expr::log(std::forward<G_F>(g))) *
             sift_term(OpIdentity{}, expr::get<Ms>(a)...);
    }
  }
}

template <typename V, typename... G1s, exp_key_t... X1s, typename G_F>
auto pick_terms(OpTerms<V, Term<G1s, X1s>...> const& e, std::index_sequence<>,
                G_F&& g) {
  return expr::make_term(expr::terms_after_first(e));
}

template <typename V, typename... G1s, exp_key_t... X1s, size_t N, typename G_F>
auto pick_terms(OpTerms<V, Term<G1s, X1s>...> const& e, std::index_sequence<N>,
                G_F&& g) {
  using namespace symphas::lib;

  using pick_terms1_t = seq_offset_t<1, std::make_index_sequence<N>>;
  using pick_terms2_t =
      seq_offset_t<N + 2, std::make_index_sequence<sizeof...(G1s) - N - 1>>;
  using pick_power_t = symphas::lib::direct_type_at_index<
      N, std::integer_sequence<exp_key_t, X1s>...>;

  return _pick_terms(e, pick_terms1_t{}, pick_terms2_t{}, pick_power_t{},
                     std::forward<G_F>(g));
}

template <typename V, typename... G1s, exp_key_t... X1s, size_t N0, size_t N1,
          size_t... Ns, typename G_F>
auto pick_terms(OpTerms<V, Term<G1s, X1s>...> const& e,
                std::index_sequence<N0, N1, Ns...>, G_F&& g) {
  auto list =
      std::make_tuple(std::forward<G_F>(g), std::forward<G_F>(g),
                      ([&](auto) { return std::forward<G_F>(g); })(Ns)...);
  using seq_t =
      symphas::lib::filter_seq_t<std::make_index_sequence<sizeof...(G1s)>,
                                 std::index_sequence<Ns...>>;
  return recombine_terms(e, seq_t{}, list,
                         std::make_index_sequence<sizeof...(Ns) + 2>{});
}

template <typename V, typename... Gs, exp_key_t... Xs, size_t... Ns>
auto select_index(OpTerms<V, Term<Gs, Xs>...> const& e,
                  std::index_sequence<Ns...>, DynamicIndexSet g) {
  using namespace symphas::lib;

  ((expr::get<Ns + 1>(const_cast<OpTerms<V, Term<Gs, Xs>...>&>(e))
        .data()
        .fix(g)),
   ...);
  return e;
}

template <typename V, typename... Gs, exp_key_t... Xs, size_t... Is, bool... fs,
          typename G_F>
decltype(auto) swap_terms(OpTerms<V, Term<Gs, Xs>...> const& e,
                          std::index_sequence<Is...>,
                          std::integer_sequence<bool, fs...>, G_F&& g) {
  using ::symphas::lib::seq_join_t;

  using swap_seq_t = seq_join_t<std::index_sequence<>,
                                std::conditional_t<fs, std::index_sequence<Is>,
                                                   std::index_sequence<>>...>;

  return pick_terms(e, swap_seq_t{}, std::forward<G_F>(g));
}

template <size_t Z, typename... Es, typename G_F, size_t... Is>
auto swap_grid_adds(OpAdd<Es...> const& e, G_F&& g,
                    std::index_sequence<Is...>) {
  return (swap_grid<Z>(expr::get<Is>(e), std::forward<G_F>(g)) + ... +
          OpVoid{});
}
}  // namespace

template <size_t Z, typename E, typename G_F>
auto swap_grid(OpExpression<E> const& e, G_F&&) {
  return *static_cast<E const*>(&e);
}

template <size_t Z, typename E, typename G_F>
auto swap_grid(OpOperator<E> const& e, G_F&&) {
  return *static_cast<E const*>(&e);
}

template <size_t Z, typename V, typename... Gs, exp_key_t... Xs, typename G_F>
auto swap_grid(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g) {
  using mask_t =
      std::integer_sequence<bool, (expr::factor_count<Variable<Z>, Gs>::value >
                                   0)...>;
  return expr::coeff(e) * swap_terms(e,
                                     std::make_index_sequence<sizeof...(Gs)>{},
                                     mask_t{}, std::forward<G_F>(g));
}

template <size_t Z, typename... Es, typename G_F>
auto swap_grid(OpAdd<Es...> const& e, G_F&& g) {
  return swap_grid_adds<Z>(e, std::forward<G_F>(g),
                           std::make_index_sequence<sizeof...(Es)>{});
}

template <size_t Z, typename E1, typename E2, typename G_F>
auto swap_grid(OpBinaryMul<E1, E2> const& e, G_F&& g) {
  return swap_grid<Z>(e.a, std::forward<G_F>(g)) *
         swap_grid<Z>(e.b, std::forward<G_F>(g));
}

template <size_t Z, typename E1, typename E2, typename G_F>
auto swap_grid(OpBinaryDiv<E1, E2> const& e, G_F&& g) {
  return swap_grid<Z>(e.a, std::forward<G_F>(g)) /
         swap_grid<Z>(e.b, std::forward<G_F>(g));
}

template <size_t Z, typename V, typename E1, typename E2, typename G_F>
auto swap_grid(OpConvolution<V, E1, E2> const& e, G_F&& g) {
  return expr::make_convolution(e.value,
                                swap_grid<Z>(e.a, std::forward<G_F>(g)),
                                swap_grid<Z>(e.b, std::forward<G_F>(g)));
}

template <size_t Z, typename V, size_t D, typename E, typename G_F>
auto swap_grid(OpConvolution<V, GaussianSmoothing<D>, E> const& e, G_F&& g) {
  return expr::make_convolution(
      e.value,
      swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)),
      e.smoother);
}

template <size_t Z, size_t O, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorDerivative<O, V, Sp> const& e, G_F&& g) {
  return e;
}

template <size_t Z, Axis ax, size_t O, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e,
               G_F&& g) {
  return e;
}

template <size_t Z, size_t... Os, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorMixedDerivative<V, Sp, Os...> const& e, G_F&& g) {
  return e;
}

template <size_t Z, typename A1, typename A2, typename G_F>
auto swap_grid(OpOperatorChain<A1, A2> const& e, G_F&& g) {
  return OpOperatorChain(swap_grid<Z>(e.f, std::forward<G_F>(g)),
                         swap_grid<Z>(e.g, std::forward<G_F>(g)));
}

template <size_t Z, typename A1, typename A2, typename E, typename G_F>
auto swap_grid(OpChain<A1, A2, E> const& e, G_F&& g) {
  return swap_grid<Z>(e.combination, std::forward<G_F>(g))(
      swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
}

template <size_t Z, typename A1, typename E, typename G_F>
auto swap_grid(OpChain<A1, OpIdentity, E> const& e, G_F&& g) {
  return OpChain(
      OpOperatorChain(swap_grid<Z>(e.combination.f, std::forward<G_F>(g)),
                      OpIdentity{}),
      swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
}

template <size_t Z, typename A1, typename A2, typename G_F>
auto swap_grid(OpOperatorCombination<A1, A2> const& e, G_F&& g) {
  return OpOperatorCombination(swap_grid<Z>(e.f, std::forward<G_F>(g)),
                               swap_grid<Z>(e.g, std::forward<G_F>(g)));
}

template <size_t Z, typename A1, typename A2, typename E, typename G_F>
auto swap_grid(OpCombination<A1, A2, E> const& e, G_F&& g) {
  return swap_grid<Z>(e.combination, std::forward<G_F>(g))(
      swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
}

template <size_t Z, typename Dd, typename V, typename E, typename Sp,
          typename G_F>
auto swap_grid(OpDerivative<Dd, V, E, Sp> const& e, G_F&& g) {
  constexpr size_t order = OpDerivative<Dd, V, E, Sp>::order;
  constexpr Axis axis = OpDerivative<Dd, V, E, Sp>::axis;
  return symphas::internal::nth_derivative_apply<axis, order, Sp>::get(
      e.value,
      swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)),
      e.solver);
}

template <size_t Z, size_t O, typename V, typename E, typename GG, typename G_F>
auto swap_grid(
    OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& e,
    G_F&& g) {
  return expr::make_derivative<O, GG>(
      e.value,
      swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)),
      e.solver);
}

template <size_t Z, typename V, typename E, typename T, typename G_F>
auto swap_grid(OpIntegral<V, E, T> const& e, G_F&& g) {
  return expr::make_integral(
      expr::coeff(e),
      swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)),
      e.domain);
}

template <size_t Z, typename G, typename V, typename E, typename G_F>
auto swap_grid(OpMap<G, V, E> const& e, G_F&& g) {
  auto eg = swap_grid<Z>(e.e, std::forward<G_F>(g));
  return OpMap<G, V, decltype(eg)>(expr::coeff(e), eg);
}

template <size_t Z, typename G, typename E, typename G_F>
auto swap_grid(OpMap<G, OpIdentity, E> const& e, G_F&& g) {
  auto eg = swap_grid<Z>(e.e, std::forward<G_F>(g));
  return OpMap<G, OpIdentity, decltype(eg)>(OpIdentity{}, eg);
}

template <size_t Z, expr::exp_key_t X, typename V, typename E, typename G_F>
auto swap_grid(OpPow<X, V, E> const& e, G_F&& g) {
  return expr::make_pow<X>(
      expr::coeff(e),
      swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
}

template <size_t Z, typename V, typename E, typename F, typename Arg0,
          typename... Args, typename G_F>
auto swap_grid(OpFunction<V, E, F, Arg0, Args...> const& e, G_F&& g) {
  return OpFunction(
      e.name, e.value,
      swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), e.f,
      e.tt);
}

template <size_t Z, auto f, typename V, typename E, typename G_F>
auto swap_grid(OpFunctionApply<f, V, E> const& e, G_F&& g) {
  auto eg = swap_grid<Z>(e.e, std::forward<G_F>(g));
  return expr::coeff(e) * expr::make_function<f>(eg);
}

namespace {

template <size_t Z, typename E, size_t... ArgNs, typename... Ts,
          size_t... ArgMs, size_t... Is, typename G_F>
auto swap_grid_symbolic(SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f,
                        SymbolicTemplate<E, ArgNs...> const& tmpl,
                        std::index_sequence<ArgMs...>,
                        std::index_sequence<Is...>, G_F&& g) {
  auto swapped = swap_grid<Z>(
      tmpl(swap_grid<Z>(std::get<Is>(f.data), std::forward<G_F>(g))...),
      std::forward<G_F>(g));
  return (expr::function_of(as_variable<ArgMs>(swap_grid<Z>(
              std::get<Is>(f.data), std::forward<G_F>(g)))...) = swapped);
}

template <size_t Z, typename E, size_t... ArgNs, typename... Ts,
          size_t... ArgMs, size_t... Is, typename G_F>
auto swap_grid_symbolic(SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f,
                        OpExpression<E> const& e, std::index_sequence<ArgMs...>,
                        std::index_sequence<Is...>, G_F&& g) {
  auto swapped = swap_grid<Z>(*static_cast<E const*>(&e), std::forward<G_F>(g));
  return (expr::function_of(as_variable<ArgMs>(swap_grid<Z>(
              std::get<Is>(f.data), std::forward<G_F>(g)))...) = swapped);
}

template <size_t Z, typename E, size_t... ArgNs, typename... Ts, typename G_F>
auto swap_grid_symbolic(SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f,
                        G_F&& g) {
  using var_g = decltype(expr::get_independent_variables(std::forward<G_F>(g)));
  using seq_t =
      std::make_index_sequence<fixed_max<sizeof...(Ts), var_g::size()> + 1>;
  using seq_filt_t = symphas::lib::filter_seq_t<seq_t, var_g>;
  using seq_cut_t = symphas::lib::seq_lt_t<sizeof...(Ts), seq_filt_t>;
  return swap_grid_symbolic<Z>(
      f, expr::template_of(Variable<ArgNs>{}...) = f.e, seq_cut_t{},
      std::make_index_sequence<sizeof...(Ts)>{}, std::forward<G_F>(g));
}

template <size_t Z, typename T, typename G_F>
auto swap_grid_symbolic(SymbolicData<T> const& e, G_F&& g) {
  return e;
}

template <size_t Z, typename T, typename G_F>
auto swap_grid_symbolic(SymbolicDataArray<T> const& e, G_F&& g) {
  return e;
}

template <bool flag, typename A, typename B>
decltype(auto) switch_tuple_element(A&& a, B&& b) {
  if constexpr (flag) {
    return std::forward<A>(a);
  } else {
    return std::forward<B>(b);
  }
}

template <typename... Ts, bool... Bs, size_t... Is, typename G_F>
auto swap_symbolic_data_array(
    SymbolicDataArray<std::tuple<Term<Ts>...>> const& data,
    std::integer_sequence<bool, Bs...>, std::index_sequence<Is...>, G_F&& g) {
  return std::make_tuple(switch_tuple_element<Bs>(
      std::forward<G_F>(g), std::get<Is>(data.get_data_tuple()))...);
}

template <size_t Z, size_t... Zs, typename... Ts, typename G_F>
auto swap_grid_symbolic(
    SymbolicDataArray<std::tuple<Term<Variable<Zs, Ts>>...>> const& data,
    G_F&& g) {
  using mask_seq = std::integer_sequence<bool, (Z == Zs)...>;
  return SymbolicDataArray(swap_symbolic_data_array(
      data, mask_seq{}, std::make_index_sequence<sizeof...(Ts)>{},
      std::forward<G_F>(g)));
}

template <size_t Z, typename... Ts, size_t... Is, typename G_F>
auto swap_grid_symbolic(Substitution<Ts...> const& data,
                        std::index_sequence<Is...>, G_F&& g) {
  return Substitution(
      swap_grid_symbolic<Z>(std::get<Is>(data), std::forward<G_F>(g))...);
}

template <size_t Z, typename... Ts, typename G_F>
auto swap_grid_symbolic(Substitution<Ts...> const& data, G_F&& g) {
  return swap_grid_symbolic<Z>(data, std::make_index_sequence<sizeof...(Ts)>{},
                               std::forward<G_F>(g));
}

template <size_t Z, expr::NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type, typename G_F>
auto swap_grid_symbolic(NoiseData<nt, T, D, grid_type> const& data, G_F&& g) {
  return data;
}

template <size_t Z, typename G_F>
auto swap_grid_symbolic(DynamicIndex const& data, G_F&& g) {
  return data;
}

template <size_t Z, int N, int P, typename G_F>
auto swap_grid_symbolic(expr::symbols::i_<N, P> const& data, G_F&& g) {
  return data;
}

template <size_t Z, typename T, typename E0, typename... T0s, typename G_F>
auto swap_grid_symbolic(T const& data, SymbolicFunction<E0, T0s...> const& f,
                        G_F&& g) {
  auto e = swap_grid_symbolic<Z>(f, std::forward<G_F>(g));
  auto substitution = swap_grid_symbolic<Z>(data, std::forward<G_F>(g));

  return symphas::internal::make_symbolic_eval(OpIdentity{}, substitution, f);
}

template <size_t Z, typename Op, typename... Ts, typename E, typename E0,
          typename... T0s, typename G_F>
auto swap_grid_symbolic(
    SymbolicSeries<Op, Substitution<Ts...>, E> const& series,
    SymbolicFunction<E0, T0s...> const& f, G_F&& g) {
  auto e = swap_grid_symbolic<Z>(series.e, std::forward<G_F>(g));
  auto substitution =
      swap_grid_symbolic<Z>(series.substitution, std::forward<G_F>(g));

  return expr::series<Op>(e)(series.limits, substitution);
}
}  // namespace

template <size_t Z, typename V, typename sub_t, typename E, typename... Ts,
          typename G_F>
auto swap_grid(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e,
               G_F&& g) {
  return expr::coeff(e) *
         swap_grid_symbolic<Z>(e.data, e.f, std::forward<G_F>(g));
}

template <size_t Z0, size_t Z1, size_t... Zs, typename E, typename G_F>
auto swap_grid(OpExpression<E> const& e, G_F&& g) {
  return swap_grid<Z1, Zs...>(
      swap_grid<Z0>(*static_cast<E const*>(&e), std::forward<G_F>(g)),
      std::forward<G_F>(g));
}

template <size_t Z0, size_t Z1, size_t... Zs, typename E, typename G_F>
auto swap_grid(OpOperator<E> const& e, G_F&& g) {
  return swap_grid<Z1, Zs...>(
      swap_grid<Z0>(*static_cast<E const*>(&e), std::forward<G_F>(g)),
      std::forward<G_F>(g));
}

template <size_t Z0, size_t Z1, size_t... Zs, typename E, typename G_F0,
          typename G_F1, typename... G_Fs>
auto swap_grid(OpExpression<E> const& e, G_F0&& g0, G_F1&& g1, G_Fs&&... gs) {
  return swap_grid<Z1, Zs...>(
      swap_grid<Z0>(*static_cast<E const*>(&e), std::forward<G_F0>(g0)),
      std::forward<G_F1>(g1), std::forward<G_Fs>(gs)...);
}

template <size_t Z0, size_t Z1, size_t... Zs, typename E, typename G_F0,
          typename G_F1, typename... G_Fs>
auto swap_grid(OpOperator<E> const& e, G_F0&& g0, G_F1&& g1, G_Fs&&... gs) {
  return swap_grid<Z1, Zs...>(
      swap_grid<Z0>(*static_cast<E const*>(&e), std::forward<G_F0>(g0)),
      std::forward<G_F1>(g1), std::forward<G_Fs>(gs)...);
}

namespace impl {
/* A development of the previous algorithm which will swap the opvariable
 * based on the grid type. An important distinction is that the grid
 * will not be swapped if it is nested within a variable.
 */

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename E, typename G_F>
decltype(auto) swap_grid(OpEvaluable<E> const& e, G_F&& g) {
  return *static_cast<E const*>(&e);
}

//! Swap a data term in the expression.
/*!
 * Implementation of a successful search, where the given variable term
 * associated with the prescribed type will be switched with the given
 * expression.
 *
 * \param v The term which is swapped.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
// template<typename Sg, typename T, typename G, typename G_F,
//	typename std::enable_if_t<(expr::is_same_base<Sg, G> &&
//! expr::is_expression<G_F>), int> = 0> decltype(auto) swap_grid(OpTerm<T, G>
// const& v, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Implementation of a successful search, where the given variable term
 * associated with the prescribed type will be switched with the given
 * expression.
 *
 * \param v The term which is swapped.
 * \param g The element which will replace the variable.
 *
 *\param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename V, typename... Gs, exp_key_t... Xs,
          typename G_F>
auto swap_grid(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename... Es, typename G_F>
auto swap_grid(OpAdd<Es...> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename E1, typename E2, typename G_F>
auto swap_grid(OpBinaryMul<E1, E2> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename E1, typename E2, typename G_F>
auto swap_grid(OpBinaryDiv<E1, E2> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename A1, typename A2, typename G_F>
auto swap_grid(OpOperatorCombination<A1, A2> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename A1, typename A2, typename E, typename G_F>
auto swap_grid(OpCombination<A1, A2, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, size_t O, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorDerivative<O, V, Sp> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, Axis ax, size_t O, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, size_t... Os, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorMixedDerivative<V, Sp, Os...> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename A1, typename A2, typename G_F>
auto swap_grid(OpOperatorChain<A1, A2> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename A1, typename A2, typename E, typename G_F>
auto swap_grid(OpChain<A1, A2, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename A1, typename E, typename G_F>
auto swap_grid(OpChain<A1, OpIdentity, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename V, typename E1, typename E2, typename G_F>
auto swap_grid(OpConvolution<V, E1, E2> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename V, size_t D, typename E, typename G_F>
auto swap_grid(OpConvolution<V, GaussianSmoothing<D>, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename Dd, typename V, typename E, typename Sp,
          typename G_F>
auto swap_grid(OpDerivative<Dd, V, E, Sp> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, size_t O, typename V, typename E, typename GG,
          typename G_F>
auto swap_grid(
    OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& e,
    G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename V, typename E, typename T, typename G_F>
auto swap_grid(OpIntegral<V, E, T> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename G, typename V, typename E, typename G_F>
auto swap_grid(OpMap<G, V, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename G, typename E, typename G_F>
auto swap_grid(OpMap<G, OpIdentity, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, expr::exp_key_t X, typename V, typename E, typename G_F>
auto swap_grid(OpPow<X, V, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename V, typename E, typename F, typename Arg0,
          typename... Args, typename G_F>
auto swap_grid(OpFunction<V, E, F, Arg0, Args...> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename V, typename sub_t, typename E, typename... Ts,
          typename G_F>
auto swap_grid(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e,
               G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename V, int N, int P, typename E, typename... Ts,
          typename G_F>
auto swap_grid(OpSymbolicEval<V, expr::symbols::i_<N, P>,
                              SymbolicFunction<E, Ts...>> const& e,
               G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename V, typename E0, typename K, typename E,
          typename... Ts, typename G_F>
auto swap_grid(OpSymbolicEval<V, SymbolicListIndex<E0, K>,
                              SymbolicFunction<E, Ts...>> const& e,
               G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, auto f, typename V, typename E, typename G_F>
auto swap_grid(OpFunctionApply<f, V, E> const& e, G_F&& g);

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the data type
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * \param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename T, typename I, typename G_F>
auto swap_grid(OpCoeff<T, I> const&, G_F&& g);

template <typename Sg, typename T, typename I, size_t... Ns, typename G_F>
auto swap_grid(OpTensor<OpCoeff<T, I>, Ns...> const& coeff, G_F&& g);

template <typename Sg, typename G_F>
auto swap_grid(DynamicIndex const& data, G_F&& g);

template <typename T>
using strip_qualifiers_t = std::remove_const_t<std::remove_reference_t<T>>;

namespace {

template <typename... As, typename... Bs, exp_key_t X, typename C,
          typename... Sgs, size_t... Is, typename... G_Fs>
auto handle_case(Term<SymbolicCase<expr::case_entry<As, Bs>...>, X> const& term,
                 SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>,
                 std::tuple<G_Fs...> const& gs, std::index_sequence<Is...>) {
  constexpr int N = symphas::lib::index_of_type<C, As...>;
  if constexpr (N >= 0) {
    if constexpr (sizeof...(G_Fs) > 0) {
      auto subbed_case = swap_grid<Sgs...>(
          std::get<size_t(N)>(term.data().cases), std::get<Is>(gs)...);
      return expr::pow_x<X>(subbed_case);
    } else {
      return expr::pow_x<X>(std::get<size_t(N)>(term.data().cases));
    }
  } else {
    return expr::make_term(term);
  }
}

template <typename V, typename... Gs, exp_key_t... Xs, typename C,
          typename... Sgs, size_t... Is, typename... G_Fs>
auto handle_cases(OpTerms<V, Term<Gs, Xs>...> const& e,
                  SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>,
                  std::index_sequence<Is...>, std::tuple<G_Fs...> const& gs) {
  return std::make_tuple(
      handle_case(expr::get<Is + 1>(e),
                  SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>{}, gs,
                  std::make_index_sequence<sizeof...(G_Fs)>{})...);
}

template <typename V, typename... Gs, exp_key_t... Xs, typename C,
          typename... Sgs, size_t... Is, bool... fs, typename... G_Fs>
decltype(auto) swap_terms_case(
    OpTerms<V, Term<Gs, Xs>...> const& e,
    SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>,
    std::index_sequence<Is...>, std::integer_sequence<bool, fs...>,
    std::tuple<G_Fs...> const& gs) {
  using ::symphas::lib::seq_join_t;

  using swap_seq_t = seq_join_t<std::index_sequence<>,
                                std::conditional_t<fs, std::index_sequence<Is>,
                                                   std::index_sequence<>>...>;

  if constexpr (swap_seq_t::size() > 0) {
    auto subbed_cases =
        handle_cases(e, SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>{},
                     swap_seq_t{}, gs);
    return recombine_terms(
        e,
        symphas::lib::filter_seq_t<std::make_index_sequence<sizeof...(Gs)>,
                                   swap_seq_t>{},
        subbed_cases, std::make_index_sequence<swap_seq_t::size()>{});
  } else {
    return expr::make_term(expr::terms_after_first(e));
  }
}

template <typename V, typename... Gs, exp_key_t... Xs, size_t... Is, bool... fs,
          typename G_F>
decltype(auto) swap_terms_case(OpTerms<V, Term<Gs, Xs>...> const& e,
                               SymbolicCaseSwap<>, std::index_sequence<Is...>,
                               std::integer_sequence<bool, fs...>, G_F&& g) {
  using namespace symphas::lib;

  using swap_seq_t = seq_join_t<std::index_sequence<>,
                                std::conditional_t<fs, std::index_sequence<Is>,
                                                   std::index_sequence<>>...>;

  if constexpr (swap_seq_t::size() > 0) {
    return pick_terms(e, swap_seq_t{}, std::forward<G_F>(g));
  } else {
    return expr::make_term(expr::terms_after_first(e));
  }
}

template <typename V, typename... Gs, exp_key_t... Xs, size_t... Is, bool... fs>
decltype(auto) swap_terms_case(OpTerms<V, Term<Gs, Xs>...> const& e,
                               DynamicIndexSet, std::index_sequence<Is...>,
                               std::integer_sequence<bool, fs...>,
                               DynamicIndexSet g) {
  using namespace symphas::lib;

  using swap_seq_t = seq_join_t<std::index_sequence<>,
                                std::conditional_t<fs, std::index_sequence<Is>,
                                                   std::index_sequence<>>...>;

  if constexpr (swap_seq_t::size() > 0) {
    return select_index(e, swap_seq_t{}, g);
  } else {
    return expr::make_term(expr::terms_after_first(e));
  }
}

template <typename I>
struct swap_index {
  template <int N, int P, typename G_F>
  auto operator()(expr::symbols::i_<N, P> const& e, G_F&& g) const {
    return expr::symbols::i_<N, P>{};
  }
};

template <int N0, int P0>
struct swap_index<expr::symbols::i_<N0, P0>> {
  template <int P, int N00, int P00>
  auto operator()(expr::symbols::i_<N0, P> const& e,
                  expr::symbols::i_<N00, P00> const&) const {
    return expr::symbols::i_<N00, P>{};
  }

  template <int P, typename G_F>
  auto operator()(expr::symbols::i_<N0, P> const& e, G_F const& g) const {
    return g;
  }

  template <int N, int P, typename G_F>
  auto operator()(expr::symbols::i_<N, P> const& e, G_F const& g) const {
    return expr::symbols::i_<N, P>{};
  }
};

}  // namespace

template <typename Sg, int N, int P, typename G_F>
auto swap_grid(expr::symbols::i_<N, P> const& e, G_F&& g) {
  return swap_index<Sg>{}(e, std::forward<G_F>(g) + DynamicIndex(P));
}

template <typename Sg, typename G_F>
auto swap_grid(int n, G_F&& g) {
  return n;
}

template <typename Sg, typename T, typename G_F>
auto swap_grid(SymbolicData<T> const& e, G_F&& g) {
  return e;
}

// template<typename Sg>
// struct swap_grid_terms_redirect<false, false, true, Sg>
//{
//	template<typename V, typename... Gs, exp_key_t... Xs, typename G_F>
//	auto operator()(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g)
//	{
//		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
//		return c * symphas::internal::fix_dynamic_indices(e,
// std::make_index_sequence<sizeof...(Gs)>{}, std::forward<G_F>(g));
//	}
// };

template <bool symbolic_flag, typename Sg, typename V, typename... Gs,
          exp_key_t... Xs, typename... factors_t, typename G_F>
auto swap_grid_terms_index_false_2(OpTerms<V, Term<Gs, Xs>...> const& e,
                                   symphas::lib::types_list<factors_t...>,
                                   G_F&& g) {
  using mask_t = std::integer_sequence<bool, (factors_t::value > 0)...>;
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  if constexpr (!symbolic_flag) {
    return c * swap_terms(e, std::make_index_sequence<sizeof...(Gs)>{},
                          mask_t{}, std::forward<G_F>(g));
  } else {
    return c * swap_terms_case(e, Sg{},
                               std::make_index_sequence<sizeof...(Gs)>{},
                               mask_t{}, std::forward<G_F>(g));
  }
}

template <bool symbolic_flag, typename Sg, typename V, typename... Gs,
          exp_key_t... Xs, typename G_F>
auto swap_grid_terms_index_false(OpTerms<V, Term<Gs, Xs>...> const& e,
                                 G_F&& g) {
  using factors_t = symphas::lib::types_list<expr::factor_count<Sg, Gs>...>;
  return swap_grid_terms_index_false_2<symbolic_flag, Sg>(e, factors_t{},
                                                          std::forward<G_F>(g));
}

template <typename Sg, typename V, typename... Gs, exp_key_t... Xs,
          typename G_F>
auto swap_grid_terms_index_true(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g) {
  using matching_ids_t =
      symphas::internal::select_all_i_<Sg,
                                       op_types_t<OpTerms<V, Term<Gs, Xs>...>>>;
  return swap_matching_i(e, matching_ids_t{}, std::forward<G_F>(g));
}

template <bool symbolic_case_flag, bool index_flag, typename Sg, typename V,
          typename... Gs, exp_key_t... Xs, typename G_F>
auto swap_grid_terms_2(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g) {
  constexpr bool symbolic_flag =
      symbolic_case_flag || std::is_same<Sg, DynamicIndexSet>::value;
  if constexpr (!index_flag) {
    return swap_grid_terms_index_false<symbolic_flag, Sg>(e,
                                                          std::forward<G_F>(g));
  } else if constexpr (!symbolic_flag && index_flag) {
    return swap_grid_terms_index_true<Sg>(e, std::forward<G_F>(g));
  }
}

template <typename E>
struct is_symbolic_case {
  static const bool value = false;
};

//! Specialization based on expr::factor_count;
template <typename... Es>
struct is_symbolic_case<SymbolicCase<Es...>> {
  static const bool value = true;
};

//! Specialization based on expr::factor_count;
template <typename C, typename T>
struct is_symbolic_case<SymbolicCaseSwap<C, T>> {
  static const bool value = true;
};

template <bool index_flag, typename Sg, typename E, typename G_F>
auto swap_grid_terms_1(E const& e, G_F&& g) {
  constexpr bool symbolic_case_flag = is_symbolic_case<Sg>::value;
  return swap_grid_terms_2<symbolic_case_flag, index_flag, Sg>(
      e, std::forward<G_F>(g));
}

template <typename Sg, typename E, typename G_F>
auto swap_grid_terms(E const& e, G_F&& g) {
  constexpr bool index_flag = expr::has_selected_index<Sg, E>;
  return swap_grid_terms_1<index_flag, Sg>(e, std::forward<G_F>(g));
}

template <typename Sg, typename V, typename... Gs, exp_key_t... Xs,
          typename G_F>
auto swap_grid(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g) {
  return swap_grid_terms<Sg>(e, std::forward<G_F>(g));
}

namespace {
template <typename Sg, typename... Es, typename G_F, size_t... Is>
auto swap_grid_adds(OpAdd<Es...> const& e, G_F&& g,
                    std::index_sequence<Is...>) {
  return (swap_grid<Sg>(expr::get<Is>(e), std::forward<G_F>(g)) + ... +
          OpVoid{});
}

template <typename Sg, typename G, typename G_F>
auto swap_grid_solver(
    SymbolicFunctionalDerivative<DynamicVariable<G>> const& solver, G_F&& g) {
  return SymbolicFunctionalDerivative<DynamicVariable<G>>(
      swap_grid<Sg>(solver.index, std::forward<G_F>(g)));
}

template <typename Sg, typename G, typename G_F>
auto swap_grid_solver(SymbolicDerivative<DynamicVariable<G>> const& solver,
                      G_F&& g) {
  return SymbolicDerivative<DynamicVariable<G>>(
      swap_grid<Sg>(solver.index, std::forward<G>(g)));
}

template <typename Sg, typename S, typename G_F>
decltype(auto) swap_grid_solver(S const& solver, G_F&& g) {
  return solver;
}
}  // namespace

template <typename Sg, typename... Es, typename G_F>
auto swap_grid(OpAdd<Es...> const& e, G_F&& g) {
  return swap_grid_adds<Sg>(e, std::forward<G_F>(g),
                            std::make_index_sequence<sizeof...(Es)>{});
}

template <typename Sg, typename E1, typename E2, typename G_F>
auto swap_grid(OpBinaryMul<E1, E2> const& e, G_F&& g) {
  return swap_grid<Sg>(e.a, std::forward<G_F>(g)) *
         swap_grid<Sg>(e.b, std::forward<G_F>(g));
}

template <typename Sg, typename E1, typename E2, typename G_F>
auto swap_grid(OpBinaryDiv<E1, E2> const& e, G_F&& g) {
  return swap_grid<Sg>(e.a, std::forward<G_F>(g)) /
         swap_grid<Sg>(e.b, std::forward<G_F>(g));
}

template <typename Sg, typename V, typename E1, typename E2, typename G_F>
auto swap_grid(OpConvolution<V, E1, E2> const& e, G_F&& g) {
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c * expr::make_convolution(swap_grid<Sg>(e.a, std::forward<G_F>(g)),
                                    swap_grid<Sg>(e.b, std::forward<G_F>(g)));
}

template <typename Sg, typename V, size_t D, typename E, typename G_F>
auto swap_grid(OpConvolution<V, GaussianSmoothing<D>, E> const& e, G_F&& g) {
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c *
         expr::make_convolution(swap_grid<Sg>(expr::get_enclosed_expression(e),
                                              std::forward<G_F>(g)),
                                e.smoother);
}

template <typename Sg, size_t O, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorDerivative<O, V, Sp> const& e, G_F&& g) {
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c * expr::make_operator_derivative<O>(
                 swap_grid_solver<Sg>(e.solver, std::forward<G_F>(g)));
}

template <typename Sg, Axis ax, size_t O, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e,
               G_F&& g) {
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c * expr::make_operator_directional_derivative<ax, O>(
                 swap_grid_solver<Sg>(e.solver, std::forward<G_F>(g)));
}

template <typename Sg, size_t... Os, typename V, typename Sp, typename G_F>
auto swap_grid(OpOperatorMixedDerivative<V, Sp, Os...> const& e, G_F&& g) {
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c * expr::make_operator_mixed_derivative<Os...>(
                 swap_grid_solver<Sg>(e.solver, std::forward<G_F>(g)));
}

template <typename Sg, typename A1, typename A2, typename G_F>
auto swap_grid(OpOperatorChain<A1, A2> const& e, G_F&& g) {
  return OpOperatorChain(swap_grid<Sg>(e.f, std::forward<G_F>(g)),
                         swap_grid<Sg>(e.g, std::forward<G_F>(g)));
}

template <typename Sg, typename A1, typename A2, typename E, typename G_F>
auto swap_grid(OpChain<A1, A2, E> const& e, G_F&& g) {
  return swap_grid<Sg>(e.combination, std::forward<G_F>(g))(
      swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
}

template <typename Sg, typename A1, typename E, typename G_F>
auto swap_grid(OpChain<A1, OpIdentity, E> const& e, G_F&& g) {
  return OpChain(
      OpOperatorChain(swap_grid<Sg>(e.combination.f, std::forward<G_F>(g)),
                      OpIdentity{}),
      swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
}

template <typename Sg, typename A1, typename A2, typename G_F>
auto swap_grid(OpOperatorCombination<A1, A2> const& e, G_F&& g) {
  return OpOperatorCombination(swap_grid<Sg>(e.f, std::forward<G_F>(g)),
                               swap_grid<Sg>(e.g, std::forward<G_F>(g)));
}

template <typename Sg, typename A1, typename A2, typename E, typename G_F>
auto swap_grid(OpCombination<A1, A2, E> const& e, G_F&& g) {
  return swap_grid<Sg>(e.combination, std::forward<G_F>(g))(
      swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
}

template <typename Sg, typename Dd, typename V, typename E, typename Sp,
          typename G_F>
auto swap_grid(OpDerivative<Dd, V, E, Sp> const& e, G_F&& g) {
  // constexpr size_t order = OpDerivative<Dd, V, E, Sp>::order;
  // constexpr Axis axis = OpDerivative<Dd, V, E, Sp>::axis;
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c * expr::make_derivative<Dd>(
                 swap_grid<Sg>(expr::get_enclosed_expression(e),
                               std::forward<G_F>(g)),
                 swap_grid_solver<Sg>(e.solver, std::forward<G_F>(g)));
}

template <typename Sg, size_t O, typename V, typename E, typename GG,
          typename G_F>
auto swap_grid(
    OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& e,
    G_F&& g) {
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c * expr::make_derivative<O, GG>(
                 swap_grid<Sg>(expr::get_enclosed_expression(e),
                               std::forward<G_F>(g)),
                 swap_grid_solver<Sg>(e.solver, std::forward<G_F>(g)));
}

template <typename Sg, typename V, typename E, typename T, typename G_F>
auto swap_grid(OpIntegral<V, E, T> const& e, G_F&& g) {
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return expr::make_integral(
      c, swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)),
      e.domain);
}

template <typename Sg, typename G, typename V, typename E, typename G_F>
auto swap_grid(OpMap<G, V, E> const& e, G_F&& g) {
  auto eg = swap_grid<Sg>(e.e, std::forward<G_F>(g));
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c * OpMap<G, OpIdentity, decltype(eg)>(OpIdentity{}, eg);
}

template <typename Sg, typename G, typename E, typename G_F>
auto swap_grid(OpMap<G, OpIdentity, E> const& e, G_F&& g) {
  auto eg = swap_grid<Sg>(e.e, std::forward<G_F>(g));
  return OpMap<G, OpIdentity, decltype(eg)>(OpIdentity{}, eg);
}

template <typename Sg, expr::exp_key_t X, typename V, typename E, typename G_F>
auto swap_grid(OpPow<X, V, E> const& e, G_F&& g) {
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c * expr::make_pow<X>(swap_grid<Sg>(expr::get_enclosed_expression(e),
                                             std::forward<G_F>(g)));
}

template <typename Sg, typename V, typename E, typename F, typename Arg0,
          typename... Args, typename G_F>
auto swap_grid(OpFunction<V, E, F, Arg0, Args...> const& e, G_F&& g) {
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c * ::OpFunction(e.name, OpIdentity{},
                          swap_grid<Sg>(e.e, std::forward<G_F>(g)), e.f, e.tt);
}

template <typename Sg, auto f, typename V, typename E, typename G_F>
auto swap_grid(OpFunctionApply<f, V, E> const& e, G_F&& g) {
  auto&& eg = swap_grid<Sg>(e.e, std::forward<G_F>(g));
  auto&& c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c * expr::make_function<f>(OpIdentity{}, eg);
}

namespace {

template <typename Sg, typename E, size_t... ArgNs, typename... Ts,
          size_t... ArgMs, size_t... Is, typename G_F>
auto swap_grid_symbolic(SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f,
                        SymbolicTemplate<E, ArgNs...> const& tmpl,
                        std::index_sequence<ArgMs...>,
                        std::index_sequence<Is...>, G_F&& g) {
  auto swapped = swap_grid<Sg>(
      tmpl(swap_grid<Sg>(std::get<Is>(f.data), std::forward<G_F>(g))...),
      std::forward<G_F>(g));
  return (expr::function_of(as_variable<ArgMs>(swap_grid<Sg>(
              std::get<Is>(f.data), std::forward<G_F>(g)))...) = swapped);
}

template <typename Sg, typename E, size_t... ArgNs, typename... Ts,
          size_t... ArgMs, size_t... Is, typename G_F>
auto swap_grid_symbolic(SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f,
                        OpExpression<E> const& e, std::index_sequence<ArgMs...>,
                        std::index_sequence<Is...>, G_F&& g) {
  auto swapped =
      swap_grid<Sg>(*static_cast<E const*>(&e), std::forward<G_F>(g));
  return (expr::function_of(as_variable<ArgMs>(swap_grid<Sg>(
              std::get<Is>(f.data), std::forward<G_F>(g)))...) = swapped);
}

template <typename Sg, typename E, size_t... ArgNs, typename... Ts,
          typename G_F>
auto swap_grid_symbolic(SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f,
                        G_F&& g) {
  using var_g = decltype(expr::get_independent_variables(std::forward<G_F>(g)));
  using seq_t = std::make_index_sequence<sizeof...(Ts) + var_g::size()>;
  using seq_filt_t = symphas::lib::filter_seq_t<seq_t, var_g>;
  using seq_cut_t = symphas::lib::seq_lt_t<sizeof...(Ts), seq_filt_t>;
  return swap_grid_symbolic<Sg>(
      f, expr::template_of(Variable<ArgNs>{}...) = f.e, seq_cut_t{},
      std::make_index_sequence<sizeof...(Ts)>{}, std::forward<G_F>(g));
}

template <typename Sg, typename T, typename G_F>
auto swap_grid_symbolic(SymbolicData<T> const& e, G_F&& g) {
  return e;
}

template <typename Sg, typename T, typename G_F>
auto swap_grid_symbolic(SymbolicDataArray<T> const& e, G_F&& g) {
  return e;
}

template <typename Sg, typename... Ts, typename G_F>
auto swap_grid_symbolic(SymbolicDataArray<std::tuple<Term<Ts>...>> const& e,
                        G_F&& g) {
  using mask_seq = std::integer_sequence<
      bool, (expr::factor_count<OpTerm<OpIdentity, Ts>, Sg>::value > 0)...>;
  return SymbolicDataArray(swap_symbolic_data_array(
      e, mask_seq{}, std::make_index_sequence<sizeof...(Ts)>{},
      std::forward<G_F>(g)));
}

template <typename Sg, typename... Ts, size_t... Is, typename G_F>
auto swap_grid_symbolic(Substitution<Ts...> const& e,
                        std::index_sequence<Is...>, G_F&& g) {
  return Substitution(
      swap_grid_symbolic<Sg>(std::get<Is>(e), std::forward<G_F>(g))...);
}

template <typename Sg, typename... Ts, typename G_F>
auto swap_grid_symbolic(Substitution<Ts...> const& e, G_F&& g) {
  return swap_grid_symbolic<Sg>(e, std::make_index_sequence<sizeof...(Ts)>{},
                                std::forward<G_F>(g));
}

template <typename Sg, expr::NoiseType nt, typename T, size_t D,
          template <typename, size_t> typename grid_type, typename G_F>
auto swap_grid_symbolic(NoiseData<nt, T, D, grid_type> const& data, G_F&& g) {
  return data;
}

template <typename Sg, typename G_F>
auto swap_grid_symbolic(DynamicIndex const& data, G_F&& g) {
  return swap_grid<Sg>(data, std::forward<G_F>(g));
}

template <typename Sg, typename E0, typename K, typename G_F>
auto swap_grid_symbolic(SymbolicListIndex<E0, K> const& data, G_F&& g) {
  return SymbolicListIndex{swap_grid<Sg>(data.e, std::forward<G_F>(g)), K{}};
}

template <typename Sg, int N, int P, typename G_F>
auto swap_grid_symbolic(expr::symbols::i_<N, P> const& data, G_F&& g) {
  return SymbolicListIndex{swap_grid<Sg>(data, std::forward<G_F>(g))};
}

template <typename Sg, typename A1, typename B1, typename A2, typename B2,
          typename G_F>
auto swap_grid_symbolic(
    expr::series_limits<std::pair<A1, B1>, std::pair<A2, B2>> const& limit,
    G_F&& g) {
  return expr::series_limits(
      std::make_pair(
          swap_grid<Sg>(expr::limit_0(limit).first, std::forward<G_F>(g)),
          swap_grid<Sg>(expr::limit_0(limit).second, std::forward<G_F>(g))),
      std::make_pair(
          swap_grid<Sg>(expr::limit_1(limit).first, std::forward<G_F>(g)),
          swap_grid<Sg>(expr::limit_1(limit).second, std::forward<G_F>(g))));
}

template <typename Sg, typename A, typename B, typename T1, typename G_F>
auto swap_grid_symbolic(expr::series_limits<T1, std::pair<A, B>> const& limit,
                        G_F&& g) {
  return expr::series_limits(
      swap_grid<Sg>(expr::limit_0(limit), std::forward<G_F>(g)),
      std::make_pair(
          swap_grid<Sg>(expr::limit_1(limit).first, std::forward<G_F>(g)),
          swap_grid<Sg>(expr::limit_1(limit).second, std::forward<G_F>(g))));
}

template <typename Sg, typename A, typename B, typename T2, typename G_F>
auto swap_grid_symbolic(expr::series_limits<std::pair<A, B>, T2> const& limit,
                        G_F&& g) {
  return expr::series_limits(
      std::make_pair(
          swap_grid<Sg>(expr::limit_0(limit).first, std::forward<G_F>(g)),
          swap_grid<Sg>(expr::limit_0(limit).second, std::forward<G_F>(g))),
      swap_grid<Sg>(expr::limit_1(limit), std::forward<G_F>(g)));
}

template <typename Sg, typename T1, typename T2, typename G_F>
auto swap_grid_symbolic(expr::series_limits<T1, T2> const& limit, G_F&& g) {
  return expr::series_limits(
      swap_grid<Sg>(expr::limit_0(limit), std::forward<G_F>(g)),
      swap_grid<Sg>(expr::limit_1(limit), std::forward<G_F>(g)));
}

template <typename Sg, typename... T1s, typename... T2s, size_t... Is,
          typename G_F>
auto swap_grid_symbolic(
    std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
    std::index_sequence<Is...>, G_F&& g) {
  return std::make_tuple(
      swap_grid_symbolic<Sg>(std::get<Is>(limits), std::forward<G_F>(g))...);
}

template <typename Sg, typename... T1s, typename... T2s, typename G_F>
auto swap_grid_symbolic(
    std::tuple<expr::series_limits<T1s, T2s>...> const& limits, G_F&& g) {
  return swap_grid_symbolic<Sg>(
      limits, std::make_index_sequence<sizeof...(T1s)>{}, std::forward<G_F>(g));
}

template <typename Sg, typename T, typename E0, typename... T0s, typename G_F>
auto swap_grid_symbolic(T const& data, SymbolicFunction<E0, T0s...> const& f,
                        G_F&& g) {
  auto e = swap_grid_symbolic<Sg>(f, std::forward<G_F>(g));
  auto substitution = swap_grid_symbolic<Sg>(data, std::forward<G_F>(g));

  return symphas::internal::make_symbolic_eval(OpIdentity{}, substitution, e);
}

template <typename Sg, typename Op, typename... Ts, typename E, typename E0,
          typename... T0s, typename G_F>
auto swap_grid_symbolic(
    SymbolicSeries<Op, Substitution<Ts...>, E> const& series,
    SymbolicFunction<E0, T0s...> const& f, G_F&& g) {
  auto limits = swap_grid_symbolic<Sg>(series.limits, std::forward<G_F>(g));
  auto e = swap_grid<Sg>(series.e, std::forward<G_F>(g));
  auto substitution =
      swap_grid_symbolic<Sg>(series.substitution, std::forward<G_F>(g));

  return expr::recreate_series(e, limits, series, substitution);
}
}  // namespace

template <typename Sg, typename V, int N, int P, typename E, typename... Ts,
          typename G_F>
auto swap_grid(OpSymbolicEval<V, expr::symbols::i_<N, P>,
                              SymbolicFunction<E, Ts...>> const& e,
               G_F&& g) {
  if constexpr (expr::factor_count<Sg, expr::symbols::i_<N, P>>::value > 0) {
    auto substitution = swap_grid_symbolic<Sg>(e.data, std::forward<G_F>(g));
    auto ev = symphas::internal::make_symbolic_eval(
        OpIdentity{},
        SymbolicListIndex{substitution, expr::symbols::i_<N, P>{}}, e.f);

    auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
    return c * ev;
  } else {
    auto swapped = swap_grid_symbolic<Sg>(e.data, e.f, std::forward<G_F>(g));
    auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
    return c * swapped;
  }
}

template <typename Sg, typename V, typename E0, typename K, typename E,
          typename... Ts, typename G_F>
auto swap_grid(OpSymbolicEval<V, SymbolicListIndex<E0, K>,
                              SymbolicFunction<E, Ts...>> const& e,
               G_F&& g) {
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  if constexpr (expr::factor_count<Sg, K>::value > 0 ||
                expr::factor_count<Sg, SymbolicListIndex<E0, K>>::value > 0) {
    auto substitution = swap_grid_symbolic<Sg>(e.data, std::forward<G_F>(g));
    auto ev =
        symphas::internal::make_symbolic_eval(OpIdentity{}, substitution, e.f);
    return c * ev;
  } else {
    auto swapped = swap_grid_symbolic<Sg>(e.data, e.f, std::forward<G_F>(g));
    return c * swapped;
  }
}

template <typename Sg, typename V, typename sub_t, typename E, typename... Ts,
          typename G_F>
auto swap_grid(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e,
               G_F&& g) {
  auto swapped = swap_grid_symbolic<Sg>(e.data, e.f, std::forward<G_F>(g));
  auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
  return c * swapped;
}

// template<typename... Sgs, typename E, typename... G_Fs>
// auto swap_grid(std::tuple<>, OpExpression<E> const& e, G_Fs&&... gs)
//{
//	return *static_cast<E const*>(&e);
// }

// template<typename Sg, typename... Sgs, typename C0, typename... Cs, typename
// E, typename G_F0, typename... G_Fs> auto
// swap_grid(std::tuple<SymbolicCaseSwap<C0>, SymbolicCaseSwap<Cs>...>,
// OpExpression<E> const& e, G_F0&& g0, G_Fs&&... gs)
//{
//	auto c0 = swap_grid<Sg>(SymbolicCaseSwap<C0>{}, * static_cast<E
// const*>(&e), std::forward<G_F0>(g0)); 	return
// swap_grid<Sgs...>(std::tuple<SymbolicCaseSwap<Cs>...>{}, * static_cast<E
// const*>(&e), std::forward<G_Fs>(gs)...);
// }

namespace {
/*template<typename T, typename I, typename E>
auto handle_coeff_swap(OpCoeff<T, I> const& coeff, E const& swap)
{
        return coeff[expr::eval(swap)];
}*/

template <typename T, size_t N0>
auto handle_coeff_swap(
    OpCoeff<T, expr::symbols::placeholder_N_symbol_<N0>> const& coeff,
    DynamicIndex const& index) {
  return coeff(index);
}

template <typename T, int N, int P>
auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, P>> const& coeff,
                       int index) {
  return OpCoeff(coeff.data.data + P, coeff.data.len - P)[index];
}

template <typename T, int N>
auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, 0>> const& coeff,
                       int index) {
  return coeff[index];
}

template <typename T, int N, int P, size_t N0>
auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, P>> const& coeff,
                       expr::symbols::placeholder_N_symbol_<N0>) {
  return OpCoeff(coeff.data.data + P, coeff.data.len - P)(
      expr::symbols::placeholder_N_symbol_<N0>{});
}

template <typename T, int N, size_t N0>
auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, 0>> const& coeff,
                       expr::symbols::placeholder_N_symbol_<N0>) {
  return coeff(expr::symbols::placeholder_N_symbol_<N0>{});
}

template <typename T, int N, int P, size_t N0>
auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, P>> const& coeff,
                       expr::symbols::placeholder_N_<N0>) {
  return handle_coeff_swap(coeff, expr::symbols::placeholder_N_symbol_<N0>{});
}

template <typename T, int N, int P>
auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, P>> const& coeff,
                       DynamicIndex const& index) {
  return OpCoeff(coeff.data.data + P, coeff.data.len - P)(index);
}

template <typename T, int N>
auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, 0>> const& coeff,
                       DynamicIndex const& index) {
  return coeff(index);
}

template <typename T>
auto handle_coeff_swap(OpCoeff<T, DynamicIndex> const& coeff,
                       OpCoeffSwap<DynamicIndex> const&) {
  return const_cast<OpCoeff<T, DynamicIndex>&>(coeff).fix();
}

template <typename T>
auto handle_coeff_swap(OpCoeff<T, DynamicIndex> const& coeff,
                       DynamicIndexSet const& set) {
  return const_cast<OpCoeff<T, DynamicIndex>&>(coeff).fix(set);
}

template <typename T, typename I>
auto handle_coeff_swap(OpCoeff<T, I> const& coeff, DynamicIndexSet const& set) {
  return coeff;
}

// template<typename T, int N, size_t N0>
// auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, 0>> const& coeff,
// expr::symbols::placeholder_N_symbol_<N0> const&)
//{
//	return coeff(expr::symbols::placeholder_N_symbol_<N0>{});
// }

// template<typename T, int N, int P, size_t N0>
// auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, P>> const& coeff,
// expr::symbols::placeholder_N_symbol_<N0> const&)
//{
//	return OpCoeff(coeff.data.data + P, coeff.data.len -
// P)(expr::symbols::placeholder_N_symbol_<N0>{});
// }
}  // namespace

template <typename Sg, typename G_F>
auto swap_grid(DynamicIndex const& data, G_F&& g) {
  if constexpr (std::is_same<Sg, DynamicIndexSet>::value) {
    const_cast<DynamicIndex&>(data).fix(std::forward<G_F>(g));
  }
  return data;
}

template <typename Sg, typename T, typename I, typename G_F>
auto swap_grid(OpCoeff<T, I> const& coeff, G_F&& g) {
  if constexpr (expr::factor_count<Sg, OpCoeff<T, I>>::value > 0 ||
                std::is_same<Sg, DynamicIndexSet>::value) {
    return handle_coeff_swap(coeff, std::forward<G_F>(g));
  } else {
    return coeff;
  }
}

template <typename Sg, typename T, typename I, size_t... Ns, typename G_F>
auto swap_grid(OpTensor<OpCoeff<T, I>, Ns...> const& coeff, G_F&& g) {
  if constexpr (expr::factor_count<Sg, OpCoeff<T, I>>::value > 0 ||
                std::is_same<Sg, DynamicIndexSet>::value) {
    return expr::make_tensor<Ns...>(handle_coeff_swap(
        static_cast<OpCoeff<T, I> const&>(coeff), std::forward<G_F>(g)));
  } else {
    return coeff;
  }
}

}  // namespace impl

/* A development of the previous algorithm which will swap the opvariable
 * based on the grid type. An important distinction is that the grid
 * will not be swapped if it is nested within a variable.
 */

//! Swap a data term in the expression.
/*!
 * Swaps the instance of the variable term which matches the given index
 * for a different term or expression. If the term is not found, the
 * expression is returned unchanged.
 *
 * \param e The expression to search for the term to swap.
 * \param g The element which will replace the variable.
 *
 * param Sg The type of the grid to match for the swap.
 */
template <typename Sg, typename E, typename G_F>
decltype(auto) swap_grid(E const& e, G_F&& g) {
  if constexpr (expr::satisfies<E, expr::matches_with<Sg>>) {
    return expr::coeff(e) * std::forward<G_F>(g);
  } else {
    return impl::swap_grid<Sg>(e, std::forward<G_F>(g));
  }
}

template <typename... Sgs, typename C, typename E, typename... G_Fs>
auto swap_grid(SymbolicCaseSwap<C>, OpExpression<E> const& e, G_Fs&&... gs) {
  return impl::swap_grid<SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>>(
      *static_cast<E const*>(&e), std::make_tuple(std::forward<G_Fs>(gs)...));
}

template <typename... Sgs, typename C, typename E, typename... G_Fs>
auto swap_grid(SymbolicCaseSwap<C>, OpOperator<E> const& e, G_Fs&&... gs) {
  return impl::swap_grid<SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>>(
      *static_cast<E const*>(&e), std::make_tuple(std::forward<G_Fs>(gs)...));
}

template <typename Sg0, typename Sg1, typename... Sgs, typename E, typename G_F>
auto swap_grid(OpEvaluable<E> const& e, G_F&& g) {
  return swap_grid<Sg1, Sgs...>(
      swap_grid<Sg0>(*static_cast<E const*>(&e), std::forward<G_F>(g)),
      std::forward<G_F>(g));
}

template <typename Sg0, typename Sg1, typename... Sgs, typename E,
          typename G_F0, typename G_F1, typename... G_Fs>
auto swap_grid(OpEvaluable<E> const& e, G_F0&& g0, G_F1&& g1, G_Fs&&... gs) {
  return swap_grid<Sg1, Sgs...>(
      swap_grid<Sg0>(*static_cast<E const*>(&e), std::forward<G_F0>(g0)),
      std::forward<G_F1>(g1), std::forward<G_Fs>(gs)...);
}

template <Axis ax, size_t Z, typename E, typename G_F>
auto swap_grid(OpEvaluable<E> const& e, G_F&& g) {
  return swap_grid<Variable<Z, VectorComponent<ax>>>(*static_cast<E const*>(&e),
                                                     std::forward<G_F>(g));
}

template <Axis ax, typename Sg, typename E, typename G_F>
auto swap_grid(OpEvaluable<E> const& e, G_F&& g) {
  return swap_grid<VectorComponent<ax, Sg>>(*static_cast<E const*>(&e),
                                            std::forward<G_F>(g));
}

template <typename... Sgs, typename E, typename G_F>
auto swap_grid(symphas::lib::types_list<Sgs...>, E&& e, G_F&& g) {
  return swap_grid<Sgs...>(std::forward<E>(e), std::forward<G_F>(g));
}

template <typename Eg, typename E, typename E_F>
auto swap_expression(OpEvaluable<E> const& e, E_F&& g) {
  return swap_grid<expr::matches_with<Eg>>(*static_cast<E const*>(&e),
                                           std::forward<E_F>(g));
}

namespace {

template <int N0, int... P0s, expr::exp_key_t... Xs, typename G_F>
auto swap_matching_i(
    OpTerms<OpIdentity, Term<expr::symbols::i_<N0, P0s>, Xs>...> const& e,
    G_F&& g) {
  return (
      expr::pow_x<Xs>(expr::transform::swap_grid<expr::symbols::i_<N0, P0s>>(
          expr::symbols::i_<N0, P0s>{}, std::forward<G_F>(g))) *
      ...);
  // return swap_matching_i(*static_cast<OpTermsList<Term<expr::symbols::i_<N0,
  // P0s>, Xs>...> const*>(&e), std::forward<G_F>(g));
}

template <typename E, int N0, int... P0s, typename G_F>
auto swap_matching_i(OpExpression<E> const& e,
                     symphas::lib::types_list<expr::symbols::i_<N0, P0s>...>,
                     G_F&& g) {
  auto [ind, r] = expr::split::factor<expr::symbols::i_<N0, P0s>...>(
      *static_cast<E const*>(&e));
  return expr::coeff(r) * swap_matching_i(ind, std::forward<G_F>(g)) * r;
}
}  // namespace

}  // namespace expr::transform
