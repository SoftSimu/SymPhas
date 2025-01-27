
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

#include "expressions.h"

namespace symphas::internal {

template <typename expression_condition_t>
struct expression_condition_impl {
  using this_type = expression_condition_impl<expression_condition_t>;

  template <typename... Es>
  auto operator()(OpAddList<Es...> const& e) const {
    return this_type{};
  }

  template <typename... Gs>
  auto operator()(OpTermsList<Gs...> const& e) const {
    return this_type{};
  }

  template <typename... Es>
  auto operator()(OpAdd<Es...> const& e) const {
    return this_type{};
  }

  template <typename... Gs>
  auto operator()(OpTerms<Gs...> const& e) const {
    return this_type{};
  }

  template <typename E>
  auto operator()(OpEvaluable<E> const& e) const {
    return this_type{};
  }

  auto operator()(...) const { return this_type{}; }

  template <typename E>
  auto get_value() const {
    return cast().operator()(E{});
  }

  expression_condition_t& cast() {
    return *static_cast<expression_condition_t*>(this);
  }

  const expression_condition_t& cast() const {
    return *static_cast<expression_condition_t const*>(this);
  }
};

template <typename... expression_condition_ts>
struct not_expression_condition {};

template <typename... expression_condition_ts>
struct and_expression_condition {};

template <typename... expression_condition_ts>
struct or_expression_condition {};

template <typename E, typename expression_condition_t>
struct expression_satisfies_condition_impl {
  using invoke_type = expression_condition_impl<expression_condition_t>;
  using type = decltype(expression_condition_t{}(E{}));
  static const bool value = !std::is_same<type, invoke_type>::value;
};

template <typename E>
struct expression_satisfies_condition_impl<E, void> {
  static const bool value = false;
};

template <typename E, typename... expression_condition_ts>
struct expression_satisfies_condition_impl<
    E, not_expression_condition<expression_condition_ts...>> {
  static const bool value =
      (!expression_satisfies_condition_impl<E,
                                            expression_condition_ts>::value &&
       ... && true);
};

template <typename E, typename... expression_condition_ts>
struct expression_satisfies_condition_impl<
    E, and_expression_condition<expression_condition_ts...>> {
  static const bool value =
      (expression_satisfies_condition_impl<E, expression_condition_ts>::value &&
       ... && true);
};

template <typename E, typename... expression_condition_ts>
struct expression_satisfies_condition_impl<
    E, or_expression_condition<expression_condition_ts...>> {
  static const bool value =
      (expression_satisfies_condition_impl<E, expression_condition_ts>::value ||
       ... || (sizeof...(expression_condition_ts) == 0));
};

template <typename E, typename expression_condition_t>
constexpr bool expression_satisfies_condition =
    expression_satisfies_condition_impl<E, expression_condition_t>::value;

}  // namespace symphas::internal

// **************************************************************************************
// Split an expression into groups for evaluating by individual regions
// **************************************************************************************

namespace expr {
template <typename condition_t>
using expr_cond = symphas::internal::expression_condition_impl<condition_t>;
template <typename... condition_ts>
using not_ = symphas::internal::not_expression_condition<condition_ts...>;
template <typename... condition_ts>
using and_ = symphas::internal::and_expression_condition<condition_ts...>;
template <typename... condition_ts>
using or_ = symphas::internal::or_expression_condition<condition_ts...>;

template <typename E, typename expression_condition_t>
constexpr bool satisfies =
    symphas::internal::expression_satisfies_condition_impl<
        E, expression_condition_t>::value;

struct matches_any : expr_cond<matches_any> {
  using expr_cond<matches_any>::operator();
  template <typename E>
  auto operator()(OpEvaluable<E> const& e) const {
    return matches_any{};
  }
};

struct matches_series : expr_cond<matches_series> {
  using expr_cond<matches_series>::operator();
  template <typename V0, typename E, typename... Ts, int... I0s, int... P0s,
            typename A, typename B, typename... Vs>
  auto operator()(
      OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
            symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
            symphas::lib::types_list<Vs...>> const& sum) const {
    return matches_series{};
  }
};

struct matches_mul : expr_cond<matches_mul> {
  using expr_cond<matches_mul>::operator();
  template <typename A, typename B>
  auto operator()(OpBinaryMul<A, B> const& e) const {
    return matches_mul{};
  }
};

struct matches_div : expr_cond<matches_div> {
  using expr_cond<matches_div>::operator();
  template <typename A, typename B>
  auto operator()(OpBinaryDiv<A, B> const& e) const {
    return matches_div{};
  }
};

struct matches_term : expr_cond<matches_term> {
  using expr_cond<matches_term>::operator();
  template <typename V, typename... Gs, expr::exp_key_t... Xs>
  auto operator()(OpTerms<V, Term<Gs, Xs>...> const& e) const {
    return matches_term{};
  }
};

struct matches_integral : expr_cond<matches_integral> {
  using expr_cond<matches_integral>::operator();
  template <typename V, typename E, typename T>
  auto operator()(OpIntegral<V, E, T> const& e) const {
    return matches_integral{};
  }
};

struct matches_derivative : expr_cond<matches_derivative> {
  using expr_cond<matches_derivative>::operator();
  template <typename Dd, typename V, typename E, typename Sp>
  auto operator()(OpDerivative<Dd, V, E, Sp> const& e) const {
    return matches_derivative{};
  }
};

struct matches_operator : expr_cond<matches_operator> {
  using expr_cond<matches_operator>::operator();
  template <typename E>
  auto operator()(OpOperator<E> const& e) const {
    return matches_operator{};
  }
};

template <typename matches_t, typename A, typename B>
struct matching_in_mul_apply {
  static const bool value =
      symphas::internal::expression_satisfies_condition<A, matches_t> ||
      symphas::internal::expression_satisfies_condition<B, matches_t>;
};

template <typename matches_t, typename A, typename B, typename C>
struct matching_in_mul_apply<matches_t, A, OpBinaryMul<B, C>> {
  static const bool value =
      symphas::internal::expression_satisfies_condition<A, matches_t> ||
      matching_in_mul_apply<matches_t, B, C>::value;
};

template <typename matches_t, typename A, typename B, typename C>
struct matching_in_mul_apply<matches_t, OpBinaryMul<A, B>, C> {
  static const bool value =
      matching_in_mul_apply<matches_t, A, B>::value ||
      symphas::internal::expression_satisfies_condition<C, matches_t>;
};

template <typename matches_t, typename A, typename B, typename C, typename D>
struct matching_in_mul_apply<matches_t, OpBinaryMul<A, B>, OpBinaryMul<C, D>> {
  static const bool value = matching_in_mul_apply<matches_t, A, B>::value ||
                            matching_in_mul_apply<matches_t, C, D>::value;
};

template <typename matches_t>
struct matching_in_mul : expr_cond<matching_in_mul<matches_t>> {
  using expr_cond<matching_in_mul<matches_t>>::operator();
  template <typename A, typename B>
  auto operator()(OpBinaryMul<A, B> const& e) const {
    if constexpr (matching_in_mul_apply<matches_t, A, B>::value) {
      return matches_operator{};
    } else {
      return expr_cond<matching_in_mul<matches_t>>::operator()(OpVoid{});
    }
  }
};

template <typename E>
struct matches_with : expr_cond<matches_with<E>> {
  using expr_cond<matches_with<E>>::operator();
  auto operator()(E const& e) const { return matches_with<E>{}; }
};

template <typename V, typename... Gs, expr::exp_key_t... Xs>
struct matches_with<OpTerms<V, Term<Gs, Xs>...>>
    : expr_cond<matches_with<OpTerms<V, Term<Gs, Xs>...>>> {
  using expr_cond<matches_with<OpTerms<V, Term<Gs, Xs>...>>>::operator();
  template <typename V0, typename... G0s, expr::exp_key_t... X0s>
  auto operator()(OpTerms<V0, Term<G0s, X0s>...> const& e) const {
    using namespace symphas::lib;

    using G_list = types_list<Gs...>;
    using G0_list = types_list<G0s...>;
    using X_seq = std::integer_sequence<expr::exp_key_t, Xs...>;
    using X0_seq = std::integer_sequence<expr::exp_key_t, X0s...>;
    if constexpr (types_list_size<expand_types_list<
                          filter_types<G_list, G0_list>,
                          filter_types<G0_list, G_list>>>::value == 0 &&
                  seq_join_t<filter_seq_t<X_seq, X0_seq>,
                             filter_seq_t<X0_seq, X_seq>>::size() == 0 &&
                  sizeof...(Gs) == sizeof...(G0s)) {
      return matches_with<OpTerms<V, Term<Gs, Xs>...>>{};
    } else {
      return expr_cond<matches_with<OpTerms<V, Term<Gs, Xs>...>>>::operator()(
          OpVoid{});
    }
  }
};

template <typename M>
struct contains_satisfying : expr_cond<contains_satisfying<M>> {
  using self_match_t = contains_satisfying<M>;

  template <typename E, std::enable_if_t<satisfies<E, M>, int> = 0>
  auto operator()(OpEvaluable<E> const& e) const {
    return contains_satisfying<M>{};
  }

  template <typename E, std::enable_if_t<!satisfies<E, M>, int> = 0>
  auto operator()(OpEvaluable<E> const& e) const {
    return expr_cond<contains_satisfying<M>>::operator()(OpVoid{});
  }

  template <typename... Es,
            std::enable_if_t<symphas::disjunction_values<satisfies<Es, M>...>,
                             int> = 0>
  auto operator()(OpAdd<Es...> const& e) const {
    return contains_satisfying<M>{};
  }
};

template <typename M>
struct contains_satisfying_anywhere
    : expr_cond<contains_satisfying_anywhere<M>> {
  using self_match_t = contains_satisfying_anywhere<M>;

  template <typename E, std::enable_if_t<satisfies<E, M>, int> = 0>
  auto operator()(OpEvaluable<E> const& e) const {
    return contains_satisfying_anywhere<M>{};
  }

  template <typename E, std::enable_if_t<!satisfies<E, M>, int> = 0>
  auto operator()(OpEvaluable<E> const& e) const {
    return expr_cond<contains_satisfying_anywhere<M>>::operator()(OpVoid{});
  }

  template <
      typename A, typename B, typename E,
      std::enable_if_t<(satisfies<OpOperatorCombination<A, B>, self_match_t> ||
                        satisfies<E, self_match_t>),
                       int> = 0>
  auto operator()(OpCombination<A, B, E> const& e) const {
    return contains_satisfying_anywhere<M>{};
  }

  template <typename A, typename B, typename E,
            std::enable_if_t<(satisfies<OpOperatorChain<A, B>, self_match_t> ||
                              satisfies<E, self_match_t>),
                             int> = 0>
  auto operator()(OpChain<A, B, E> const& e) const {
    return contains_satisfying_anywhere<M>{};
  }

  template <
      typename A, typename B,
      std::enable_if_t<
          (satisfies<A, self_match_t> || satisfies<B, self_match_t>), int> = 0>
  auto operator()(OpOperatorCombination<A, B> const& e) const {
    return contains_satisfying_anywhere<M>{};
  }

  template <
      typename A, typename B,
      std::enable_if_t<
          (satisfies<A, self_match_t> || satisfies<B, self_match_t>), int> = 0>
  auto operator()(OpOperatorChain<A, B> const& e) const {
    return contains_satisfying_anywhere<M>{};
  }

  template <
      typename A, typename B,
      std::enable_if_t<
          (satisfies<A, self_match_t> || satisfies<B, self_match_t>), int> = 0>
  auto operator()(OpBinaryDiv<A, B> const& e) const {
    return contains_satisfying_anywhere<M>{};
  }

  template <
      typename A, typename B,
      std::enable_if_t<
          (satisfies<A, self_match_t> || satisfies<B, self_match_t>), int> = 0>
  auto operator()(OpBinaryMul<A, B> const& e) const {
    return contains_satisfying_anywhere<M>{};
  }

  template <
      typename... Es,
      std::enable_if_t<
          symphas::disjunction_values<satisfies<Es, self_match_t>...>, int> = 0>
  auto operator()(OpAdd<Es...> const& e) const {
    return contains_satisfying_anywhere<M>{};
  }
};

template <typename E>
using contains_matching = contains_satisfying<matches_with<E>>;
template <typename E>
using contains_matching_anywhere =
    contains_satisfying_anywhere<matches_with<E>>;

// I can introduce other criteria for terms like, is linear, or is operator
}  // namespace expr

namespace symphas::internal {

template <typename E, typename M>
struct expression_satisfies_condition_impl<
    E, expr::matches_with<expr::matches_with<M>>> {
  static const bool value =
      (expression_satisfies_condition_impl<E, expr::matches_with<M>>::value);
};
}  // namespace symphas::internal
