
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

#include "expressionconditions.h"
#include "expressionconvolution.h"
#include "expressionderivatives.h"
#include "expressionexponentials.h"
#include "expressionfunctions.h"
#include "expressionintegrals.h"
#include "expressionrules.h"
#include "symbolicdata.h"

// ******************************************************************************************

namespace expr {

namespace {
/* returns whether a derivative is present in the expression by
 * checking that the derivative index is greater than the index associated
 * with a constant
 */

template <typename E>
struct expr_has_deriv {
 protected:
  static const size_t deriv_order_1 = derivative_index_raw<1>::value;
  static const size_t deriv_order_0 = derivative_index_raw<0>::value;

 public:
  static const bool value =
      expr::derivative_index<deriv_order_1, E>::value > deriv_order_0;
};

template <typename T>
constexpr bool is_term_type = false;

template <typename G>
constexpr bool is_term_type<OpTerm<OpIdentity, G>> = true;

}  // namespace

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename E>
auto apply_operators(OpExpression<E> const& e) {
  return *static_cast<E const*>(&e);
}

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The operator which is distributed.
 */
template <typename E>
auto apply_operators(OpOperator<E> const& e) {
  return *static_cast<E const*>(&e);
}

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename A1, typename A2>
auto apply_operators(OpOperatorChain<A1, A2> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename A2>
auto apply_operators(OpOperatorChain<OpIdentity, A2> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename A1, typename A2>
auto apply_operators(OpOperatorCombination<A1, A2> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename A1, typename A2, typename E>
auto apply_operators(OpChain<A1, A2, E> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename A1, typename A2, typename E>
auto apply_operators(OpCombination<A1, A2, E> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename Dd, typename V, typename E, typename Sp,
          typename = std::enable_if_t<expr_has_deriv<E>::value, int>>
auto apply_operators(OpDerivative<Dd, V, E, Sp> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename Dd1, typename Dd2, typename V1, typename V2, typename E,
          typename Sp1, typename Sp2>
auto apply_operators(
    OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp1>, Sp2> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename Dd1, typename Dd2, typename V1, typename V2, typename E,
          typename Sp>
auto apply_operators(
    OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp>, Sp> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <size_t O1, typename V, size_t O2, typename V1, typename E,
          typename G00>
auto apply_operators(OpDerivative<std::index_sequence<O1>, V,
                                  OpDerivative<std::index_sequence<O2>, V1, E,
                                               SymbolicDerivative<G00>>,
                                  SymbolicDerivative<G00>> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <size_t O1, typename V, size_t O2, typename V1, typename E,
          typename G00, typename G01,
          typename = std::enable_if_t<
              !std::is_same<G01, OpDerivative<std::index_sequence<O2>, V1, E,
                                              SymbolicDerivative<G00>>>::value,
              int>>
auto apply_operators(OpDerivative<std::index_sequence<O1>, V,
                                  OpDerivative<std::index_sequence<O2>, V1, E,
                                               SymbolicDerivative<G00>>,
                                  SymbolicDerivative<G01>> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename Dd, typename V, typename... Es, typename Sp>
auto apply_operators(OpDerivative<Dd, V, OpAdd<Es...>, Sp> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename... Es>
auto apply_operators(OpAdd<Es...> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename E1, typename E2>
auto apply_operators(OpBinaryMul<E1, E2> const& e);

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
template <typename E1, typename E2>
auto apply_operators(OpBinaryDiv<E1, E2> const& e);

//! Implementation of the division rule.
template <
    typename Dd, typename Sp, typename V, typename E1, typename E2,
    typename std::enable_if_t<
        (OpDerivative<Dd, V, OpBinaryMul<E1, E2>, Sp>::order > 0), int> = 0>
auto apply_operators(OpDerivative<Dd, V, OpBinaryMul<E1, E2>, Sp> const& e);

//! Implementation of the product rule.
template <
    typename Dd, typename Sp, typename V, typename E1, typename E2,
    typename std::enable_if_t<
        (OpDerivative<Dd, V, OpBinaryDiv<E1, E2>, Sp>::order > 0), int> = 0>
auto apply_operators(OpDerivative<Dd, V, OpBinaryDiv<E1, E2>, Sp> const& e);

////! Implementation of the product rule for symbolic derivatives.
// template<size_t O, typename V, typename A1, typename A2, typename G>
// auto apply_operators(OpDerivative<std::index_sequence<O>, V,
// OpOperatorChain<A1, A2>, SymbolicDerivative<G>> const& e);

////! Implementation of the product rule for symbolic derivatives.
// template<size_t O, typename V, typename A1, typename A2, typename E, typename
// G> auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpChain<A1,
// A2, E>, SymbolicDerivative<G>> const& e);

////! Implementation of the product rule for symbolic derivatives.
// template<size_t O, typename V, typename A1, typename A2, typename G>
// auto apply_operators(OpDerivative<std::index_sequence<O>, V,
// OpOperatorCombination<A1, A2>, SymbolicDerivative<G>> const& e);

////! Implementation of the product rule for symbolic derivatives.
// template<size_t O, typename V, typename A1, typename A2, typename E, typename
// G> auto apply_operators(OpDerivative<std::index_sequence<O>, V,
// OpCombination<A1, A2, E>, SymbolicDerivative<G>> const& e);

//! Specialization based on expr::grid_dim.
template <typename V, typename E, typename F>
auto apply_operators(OpFunction<V, E, F, void> const& e);

//! Specialization based on expr::grid_dim.
template <typename V, typename E, typename F, typename Arg0, typename... Args>
auto apply_operators(OpFunction<V, E, F, Arg0, Args...> const& e);

//! Specialization based on expr::grid_dim.
template <auto f, typename V, typename E>
auto apply_operators(OpFunctionApply<f, V, E> const& e);

template <typename V, typename V0, typename E, typename T, typename G>
auto apply_operators(
    OpFunctionalDerivative<V, OpDomainIntegral<V0, E, T>, G> const& e);

template <size_t O, typename V, size_t Z, int I0, int P0, typename V0,
          typename E, typename... Ts, int... I0s, int... P0s, typename A,
          typename B, typename... Vs>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
                       symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
                       A, B, symphas::lib::types_list<Vs...>>,
                 SymbolicDerivative<expr::symbols::v_id_type<
                     expr::symbols::i_<I0, P0>>>> const& e);

template <size_t O, typename V, typename V0, typename E, typename... Ts,
          int... I0s, int... P0s, typename A, typename B, typename... Vs,
          typename GG>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
                       symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
                       A, B, symphas::lib::types_list<Vs...>>,
                 SymbolicDerivative<GG>> const& e);

template <typename V, typename E, typename... Ts, int... I0s, int... P0s,
          typename A, typename B, typename... Vs>
auto apply_operators(
    OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>,
          symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
          symphas::lib::types_list<Vs...>> const& e);

//! Implementation of the product rule for terms.
// template<typename Dd, typename Sp, typename V, typename V0, typename G0,
// expr::exp_key_t X0, typename G1, expr::exp_key_t X1, typename... Gs,
// expr::exp_key_t... Xs, 	typename std::enable_if_t<(Dd::order > 0), int>
// = 0> auto apply_operators(OpDerivative<Dd, V, OpTerms<V0, Term<G0, X0>,
// Term<G1, X1>, Term<Gs, Xs>...>, Sp> const& e);

template <
    typename Dd, typename Sp, typename V, typename V0, typename E,
    typename std::enable_if_t<
        (OpDerivative<Dd, V, OpExponential<V0, E>, Sp>::order > 0), int> = 0>
auto apply_operators(OpDerivative<Dd, V, OpExponential<V0, E>, Sp> const& e);

template <size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0,
          typename GG>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<expr::symbols::i_<I0, P0>, X0>>,
                 SymbolicDerivative<GG>> const& e);

template <size_t O, typename V, typename V1, int I0, int P0, typename G1,
          typename... Gs, expr::exp_key_t X0, expr::exp_key_t X1,
          expr::exp_key_t... Xs, typename GG>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V1, Term<expr::symbols::i_<I0, P0>, X0>, Term<G1, X1>,
                         Term<Gs, Xs>...>,
                 SymbolicDerivative<GG>> const& e);

template <size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0,
          size_t Z, typename GG, size_t N = expr::_Xk_t<X0>::N,
          typename = std::enable_if_t<(O > 0 && N >= O), int>>
auto apply_operators(
    OpDerivative<
        std::index_sequence<O>, V,
        OpTerms<V0,
                Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>,
        SymbolicDerivative<Variable<Z, GG>>> const& e);

template <size_t O, typename V, typename V0, int I0, int P0, size_t D,
          expr::exp_key_t X0, size_t Z, typename GG,
          size_t N = expr::_Xk_t<X0>::N,
          typename = std::enable_if_t<(O > 0 && N >= O), int>>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<
                                                 expr::symbols::i_<I0, P0>>,
                                             D>,
                                  X0>>,
                 SymbolicDerivative<Variable<Z, GG>>> const& e);

template <size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0,
          typename GG, size_t N = expr::_Xk_t<X0>::N,
          typename = std::enable_if_t<(O > 0 && N >= O), int>>
auto apply_operators(
    OpDerivative<
        std::index_sequence<O>, V,
        OpTerms<V0,
                Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>,
        SymbolicDerivative<DynamicVariable<GG>>> const& e);

template <size_t O, typename V, typename V0, int I0, int P0, size_t D,
          expr::exp_key_t X0, typename GG, size_t N = expr::_Xk_t<X0>::N,
          typename = std::enable_if_t<(O > 0 && N >= O), int>>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<
                                                 expr::symbols::i_<I0, P0>>,
                                             D>,
                                  X0>>,
                 SymbolicDerivative<DynamicVariable<GG>>> const& e);

template <size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0,
          typename GG, size_t N = expr::_Xk_t<X0>::N,
          typename = std::enable_if_t<
              (O > 0 && N >= O && !expr::is_expression<GG>), int>>
auto apply_operators(
    OpDerivative<
        std::index_sequence<O>, V,
        OpTerms<V0,
                Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>,
        SymbolicDerivative<GG>> const& e);

template <size_t O, typename V, typename V0, int I0, int P0, size_t D,
          expr::exp_key_t X0, typename GG, size_t N = expr::_Xk_t<X0>::N,
          typename = std::enable_if_t<
              (O > 0 && N >= O && !expr::is_expression<GG>), int>>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<
                                                 expr::symbols::i_<I0, P0>>,
                                             D>,
                                  X0>>,
                 SymbolicDerivative<GG>> const& e);

template <
    size_t O, typename V, typename V0, typename G0, expr::exp_key_t X0,
    typename GG, size_t N = expr::factor_count_list<GG, Term<G0, X0>>::value,
    typename =
        std::enable_if_t<(O > 0 && N >= O && !expr::is_expression<GG>), int>>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>>,
                 SymbolicDerivative<GG>> const& e);

template <size_t O, typename V, typename V0, typename G0, typename G1,
          typename... Gs, expr::exp_key_t X0, expr::exp_key_t X1,
          expr::exp_key_t... Xs, typename GG,
          typename = std::enable_if_t<!is_term_type<GG>, int>>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>,
                 SymbolicDerivative<GG>> const& e);

template <size_t O, typename V, typename V0, typename G0, expr::exp_key_t X0,
          typename GG>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>>,
                 SymbolicDerivative<OpTerm<OpIdentity, GG>>> const& e);

template <size_t O, typename V, typename V0, typename G0, typename G1,
          typename... Gs, expr::exp_key_t X0, expr::exp_key_t X1,
          expr::exp_key_t... Xs, typename GG>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>,
                 SymbolicDerivative<OpTerm<OpIdentity, GG>>> const& e);

template <size_t O, typename V0, auto f, typename V1, typename E, typename G0>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V0, OpFunctionApply<f, V1, E>,
                 SymbolicDerivative<G0>> const& e);

template <size_t O, typename V, typename... Es, typename G0>
auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpAdd<Es...>,
                                  SymbolicDerivative<G0>> const& e);

template <size_t O, typename G, typename V, typename E,
          typename std::enable_if_t<!expr_has_deriv<E>::value, int> = 0>
auto apply_operators(OpDerivative<std::index_sequence<O>, V, E,
                                  SymbolicDerivative<G>> const& e) {
  return OpVoid{};
}

template <typename V, expr::exp_key_t X0, typename V0, typename E0, typename GG>
auto apply_operators(OpDerivative<std::index_sequence<1>, V, OpPow<X0, V0, E0>,
                                  SymbolicDerivative<GG>> const& e);

template <size_t O, typename V, expr::exp_key_t X0, typename V0, typename E0,
          typename GG>
auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpPow<X0, V0, E0>,
                                  SymbolicDerivative<GG>> const& e);

template <expr::exp_key_t X, typename V, typename E>
auto apply_operators(OpPow<X, V, E> const& e);

template <size_t O, typename V, typename E1, typename E2, typename G>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V, OpBinaryDiv<E1, E2>,
                 SymbolicDerivative<G>> const& e);

//! Implementation of the product rule for symbolic derivatives.
template <size_t O, typename V, typename E1, typename E2, typename G>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V, OpBinaryMul<E1, E2>,
                 SymbolicDerivative<G>> const& e);

//! Implementation of the product rule for symbolic derivatives.
template <size_t O, typename V, typename A1, typename E, typename G,
          std::enable_if_t<expr::is_coeff<A1>, int> = 0>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V, OpChain<A1, OpIdentity, E>,
                 SymbolicDerivative<G>> const& e);

template <typename V, typename V1, typename E, typename Dd, typename Sp>
auto apply_operators(
    OpDerivative<std::index_sequence<1>, V, OpDerivative<Dd, V1, E, Sp>,
                 SymbolicDerivative<OpDerivative<Dd, OpIdentity, E, Sp>>> const&
        e) {
  auto expr = expr::get_enclosed_expression(e);
  return expr::coeff(e) * expr::coeff(expr);
}

template <typename V, typename E, typename Dd, typename Sp>
auto apply_operators(
    OpDerivative<std::index_sequence<1>, V, OpDerivative<Dd, OpIdentity, E, Sp>,
                 SymbolicDerivative<OpDerivative<Dd, OpIdentity, E, Sp>>> const&
        e) {
  // auto expr = expr::get_enclosed_expression(e);
  return expr::coeff(e);
}

template <size_t O, typename V, typename V1, typename E, typename Dd,
          typename Sp, typename G0,
          typename std::enable_if_t<
              !std::is_same<G0, OpDerivative<Dd, V1, E, Sp>>::value, int> = 0>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V, OpDerivative<Dd, V1, E, Sp>,
                 SymbolicDerivative<G0>> const& e) {
  return OpVoid{};
}

template <size_t O1, typename V, size_t O2, typename V1, typename E,
          typename G00>
auto apply_operators(OpDerivative<std::index_sequence<O1>, V,
                                  OpDerivative<std::index_sequence<O2>, V1, E,
                                               SymbolicDerivative<G00>>,
                                  SymbolicDerivative<G00>> const& e) {
  auto expr = expr::get_enclosed_expression(e);
  return expr::coeff(e) * expr::coeff(expr) *
         apply_operators(expr::make_derivative<O1 + O2, G00>(
             expr::get_enclosed_expression(expr), e.solver));
}

template <size_t O1, typename V, size_t O2, typename V1, typename E,
          typename G00, typename G01, typename>
auto apply_operators(OpDerivative<std::index_sequence<O1>, V,
                                  OpDerivative<std::index_sequence<O2>, V1, E,
                                               SymbolicDerivative<G00>>,
                                  SymbolicDerivative<G01>> const& e) {
  auto expr = expr::get_enclosed_expression(e);
  return expr::coeff(e) * apply_operators(expr::make_derivative<O1, G01>(
                              apply_operators(expr), e.solver));
}

template <typename G, typename V, typename E>
auto apply_operators(OpDerivative<std::index_sequence<0>, V, E,
                                  SymbolicDerivative<G>> const& e) {
  return expr::coeff(e) * expr::get_enclosed_expression(e);
}

template <typename V, typename E>
auto apply_operators(OpDerivative<std::index_sequence<1>, V, E,
                                  SymbolicDerivative<E>> const& e) {
  return expr::coeff(e);
}

namespace {

template <typename Dd, typename Sp, typename E1, typename E2>
auto handle_apply_mul(Sp const& solver, OpExpression<E1> const& lhs0,
                      OpExpression<E2> const& rhs0) {
  auto lhs = apply_operators(
      expr::make_derivative<Dd>(*static_cast<E1 const*>(&lhs0), solver) *
      (*static_cast<E2 const*>(&rhs0)));
  auto rhs = apply_operators(
      (*static_cast<E1 const*>(&lhs0)) *
      expr::make_derivative<Dd>(*static_cast<E2 const*>(&rhs0), solver));
  return lhs + rhs;
}

template <typename Dd, typename Sp, typename E1, typename E2>
auto handle_apply_mul(Sp const& solver, OpOperator<E1> const& lhs0,
                      OpExpression<E2> const& rhs0) {
  return apply_operators(expr::make_derivative<Dd>(
      apply_operators((*static_cast<E1 const*>(&lhs0)) *
                      (*static_cast<E2 const*>(&rhs0))),
      solver));
}

template <typename Dd, typename Sp, typename E1, typename E2>
auto handle_apply_mul(Sp const& solver, OpExpression<E1> const& lhs0,
                      OpOperator<E2> const& rhs0) {
  auto lhs = apply_operators(
      expr::make_derivative<Dd>(*static_cast<E1 const*>(&lhs0), solver) *
      (*static_cast<E2 const*>(&rhs0)));
  auto rhs = apply_operators(
      (*static_cast<E1 const*>(&lhs0)) *
      expr::make_derivative<Dd>(*static_cast<E2 const*>(&rhs0), solver));
  return lhs + rhs;
}

template <typename Dd, typename Sp, typename E1, typename E2>
auto handle_apply_mul(Sp const& solver, OpOperator<E1> const& lhs0,
                      OpOperator<E2> const& rhs0) {
  return apply_operators(expr::make_derivative<Dd>(
      apply_operators((*static_cast<E1 const*>(&lhs0)) *
                      (*static_cast<E2 const*>(&rhs0))),
      solver));
}

template <typename Dd, typename Sp, typename E10, typename E11, typename E2>
auto handle_apply_mul(Sp const& solver,
                      OpOperatorCombination<E10, E11> const& lhs0,
                      OpExpression<E2> const& rhs0) {
  return apply_operators(
      expr::make_derivative<Dd>(lhs0.f * (*static_cast<E2 const*>(&rhs0)) +
                                    lhs0.g * (*static_cast<E2 const*>(&rhs0)),
                                solver));
}

template <typename Dd, typename Sp, typename E10, typename E11, typename E2>
auto handle_apply_mul(Sp const& solver, OpOperatorChain<E10, E11> const& lhs0,
                      OpExpression<E2> const& rhs0) {
  return apply_operators(expr::make_derivative<Dd>(
      lhs0.f(lhs0.g * (*static_cast<E2 const*>(&rhs0))), solver));
}

template <typename Dd, typename Sp, typename E10, typename E11, typename E2>
auto handle_apply_mul(Sp const& solver,
                      OpOperatorCombination<E10, E11> const& lhs0,
                      OpOperator<E2> const& rhs0) {
  return apply_operators(
      expr::make_derivative<Dd>(lhs0.f * (*static_cast<E2 const*>(&rhs0)) +
                                    lhs0.g * (*static_cast<E2 const*>(&rhs0)),
                                solver));
}

template <typename Dd, typename Sp, typename E10, typename E11, typename E2>
auto handle_apply_mul(Sp const& solver, OpOperatorChain<E10, E11> const& lhs0,
                      OpOperator<E2> const& rhs0) {
  return apply_operators(expr::make_derivative<Dd>(
      lhs0.f(lhs0.g * (*static_cast<E2 const*>(&rhs0))), solver));
}

template <typename E>
constexpr bool is_op_operator_v = std::is_base_of<OpOperator<E>, E>::value;

//! Sum over axes of (d_i a)(d_i b), the cross term of the laplacian rule.
template <typename Sp, typename E1, typename E2, Axis... axs>
auto product_rule_cross(Sp const& solver, E1 const& a, E2 const& b,
                        symphas::lib::axis_list<axs...>) {
  return (OpVoid{} + ... +
          (apply_operators(expr::make_derivative<
                           typename Solver<Sp>::template derivative<axs, 1>>(
               a, solver)) *
           apply_operators(expr::make_derivative<
                           typename Solver<Sp>::template derivative<axs, 1>>(
               b, solver))));
}

template <size_t... Os>
constexpr size_t first_nonzero_order() {
  constexpr size_t orders[]{Os...};
  size_t i = 0;
  while (i < sizeof...(Os) && orders[i] == 0) ++i;
  return i;
}

//! Mixed derivative of a product: takes one first-order derivative along the
//! first axis with a nonzero order, then applies the rest. The rest is built
//! with one order per grid axis (axs), the form the derivative combination
//! rules expect.
template <typename Sp, typename E1, typename E2, size_t... Os, Axis... axs>
auto apply_mixed_product_rule(Sp const& solver, E1 const& a, E2 const& b,
                              std::index_sequence<Os...>,
                              symphas::lib::axis_list<axs...>) {
  constexpr size_t orders[]{Os...};
  constexpr Axis axes[]{axs...};
  constexpr size_t i0 = first_nonzero_order<Os...>();
  static_assert(i0 < sizeof...(axs), "mixed derivative along a missing axis");
  using D1 = typename Solver<Sp>::template directional_derivative<axes[i0], 1>;
  using Dr = typename Solver<Sp>::template mixed_derivative<(
      orders[static_cast<size_t>(axs)] - ((axs == axes[i0]) ? 1 : 0))...>;
  return apply_operators(expr::make_derivative<Dr>(
      handle_apply_mul<D1>(solver, a, b), solver));
}

//! Recovers the per-axis orders of a mixed derivative type, which cannot be
//! deduced from the nested template. Returns std::false_type if no match.
template <typename Sp, typename Dd, size_t O, size_t Rem, size_t... Os>
auto mixed_derivative_orders();

template <typename Sp, typename Dd, size_t O, size_t Rem, size_t K,
          size_t... Os>
auto mixed_derivative_orders_from() {
  using r = decltype(mixed_derivative_orders<Sp, Dd, O, Rem - 1, Os..., K>());
  if constexpr (!std::is_same<r, std::false_type>::value) {
    return r{};
  } else if constexpr (K < O) {
    return mixed_derivative_orders_from<Sp, Dd, O, Rem, K + 1, Os...>();
  } else {
    return std::false_type{};
  }
}

template <typename Sp, typename Dd, size_t O, size_t Rem, size_t... Os>
auto mixed_derivative_orders() {
  if constexpr (Rem > 0) {
    return mixed_derivative_orders_from<Sp, Dd, O, Rem, 0, Os...>();
  } else if constexpr (std::is_same<Dd, typename Solver<Sp>::template
                                            mixed_derivative<Os...>>::value) {
    return std::index_sequence<Os...>{};
  } else {
    return std::false_type{};
  }
}

//! Product rule for a spatial derivative Dd of any order.
/*!
 * handle_apply_mul applies D(ab) = (Da) b + a (Db), which holds only for
 * first-order D. A higher-order derivative is split into first-order steps
 * and the laplacian:
 * - laplacian: lap(ab) = (lap a) b + 2 grad a . grad b + a lap b;
 * - derivative<ax, O>, O even: lap^(O/2), one laplacian at a time;
 * - derivative<ax, O>, O odd: d_ax lap^((O-1)/2), d_ax applied first;
 * - directional and mixed derivatives: one first-order derivative at a time.
 */
template <typename Dd, typename Sp, typename E1, typename E2>
auto apply_product_rule(Sp const& solver, E1 const& a, E2 const& b) {
  using solver_type = Solver<Sp>;
  constexpr size_t D = expr::grid_dim<OpBinaryMul<E1, E2>>::value;
  constexpr size_t O = Dd::order;
  constexpr Axis ax = Dd::axis;

  // An operator on the left is composed with the right side, not multiplied
  // (see handle_apply_mul). An operator on the right, such as the row vector
  // in a dot product, is multiplied like any expression.
  if constexpr (O <= 1 || D == 0 || is_op_operator_v<E1>) {
    return handle_apply_mul<Dd>(solver, a, b);
  } else if constexpr (Dd::is_mixed) {
    // The mixed derivative factories list D orders or, for the axis-named
    // forms, always 3.
    using orders_d = decltype(mixed_derivative_orders<Sp, Dd, O, D>());
    using orders_t =
        std::conditional_t<std::is_same<orders_d, std::false_type>::value,
                           decltype(mixed_derivative_orders<Sp, Dd, O, 3>()),
                           orders_d>;
    static_assert(!std::is_same<orders_t, std::false_type>::value,
                  "unrecognized mixed derivative");
    return apply_mixed_product_rule(
        solver, a, b, orders_t{},
        symphas::lib::make_axis_list<D>());
  } else if constexpr (Dd::is_directional) {
    using D1 = typename solver_type::template directional_derivative<ax, 1>;
    using Dr = typename solver_type::template directional_derivative<ax, O - 1>;
    return apply_operators(expr::make_derivative<Dr>(
        handle_apply_mul<D1>(solver, a, b), solver));
  } else if constexpr (O % 2 == 1) {
    using D1 = typename solver_type::template derivative<ax, 1>;
    using Dr = typename solver_type::template derivative<ax, O - 1>;
    return apply_operators(expr::make_derivative<Dr>(
        handle_apply_mul<D1>(solver, a, b), solver));
  } else {
    using L = typename solver_type::template derivative<ax, 2>;
    auto lap = apply_operators(expr::make_derivative<L>(a, solver)) * b +
               expr::make_integer<2>() *
                   product_rule_cross(solver, a, b,
                                      symphas::lib::make_axis_list<D>()) +
               a * apply_operators(expr::make_derivative<L>(b, solver));
    if constexpr (O == 2) {
      return lap;
    } else {
      using Dr = typename solver_type::template derivative<ax, O - 2>;
      return apply_operators(expr::make_derivative<Dr>(lap, solver));
    }
  }
}

}  // namespace

//! Implementation of the quotient rule for symbolic derivatives.
template <size_t O, typename V, typename E1, typename E2, typename G>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V, OpBinaryDiv<E1, E2>,
                 SymbolicDerivative<G>> const& e) {
  // auto expr = expr::get_enclosed_expression(e);
  // return expr::coeff(e) * handle_apply_mul(expr::make_operator_derivative<O,
  // G>(e.solver), apply_operators(expr.a), apply_operators(expr.b));

  auto expr = expr::get_enclosed_expression(e);
  auto a = apply_operators(expr.a);
  auto b = apply_operators(expr.b);

  auto lhs = apply_operators(expr::make_derivative<O, G>(a, e.solver) * b);
  auto rhs = apply_operators(a * expr::make_derivative<O, G>(b, e.solver));
  return expr::coeff(e) * (lhs - rhs) / (b * b);
}

//! Implementation of the quotient rule for symbolic derivatives.
template <size_t O, typename V, typename E1, typename E2, typename G>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V, OpBinaryMul<E1, E2>,
                 SymbolicDerivative<G>> const& e) {
  auto expr = expr::get_enclosed_expression(e);
  return expr::coeff(e) *
         handle_apply_mul<std::index_sequence<O>>(
             e.solver, apply_operators(expr.a), apply_operators(expr.b));
}

//! Implementation of the product rule for symbolic derivatives.
template <size_t O, typename V, typename A1, typename E, typename G,
          std::enable_if_t<expr::is_coeff<A1>, int>>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V, OpChain<A1, OpIdentity, E>,
                 SymbolicDerivative<G>> const& e) {
  auto expr = expr::get_enclosed_expression(e);
  auto applied = apply_operators(expr::make_derivative<O, G>(
      expr::get_enclosed_expression(expr), e.solver));
  return expr::make_operator_chain(expr::coeff(e) * expr.combination.f,
                                   OpIdentity{})(applied);
}

//! Implementation of the product rule.
template <typename Dd, typename Sp, typename V, typename E1, typename E2,
          typename std::enable_if_t<
              (OpDerivative<Dd, V, OpBinaryMul<E1, E2>, Sp>::order > 0), int>>
auto apply_operators(OpDerivative<Dd, V, OpBinaryMul<E1, E2>, Sp> const& e) {
  auto expr = expr::get_enclosed_expression(e);
  return expr::coeff(e) * apply_product_rule<Dd>(e.solver,
                                                 apply_operators(expr.a),
                                                 apply_operators(expr.b));

  // auto expr = expr::get_enclosed_expression(e);
  // auto a = apply_operators(expr.a);
  // auto b = apply_operators(expr.b);
  //
  // auto lhs = apply_operators(expr::make_derivative<Dd>(a, e.solver) * b);
  // auto rhs = apply_operators(a * expr::make_derivative<Dd>(b, e.solver));
  // return expr::coeff(e) * (lhs + rhs);
}

//! Implementation of the quotient rule.
template <typename Dd, typename Sp, typename V, typename E1, typename E2,
          typename std::enable_if_t<
              (OpDerivative<Dd, V, OpBinaryDiv<E1, E2>, Sp>::order > 0), int>>
auto apply_operators(OpDerivative<Dd, V, OpBinaryDiv<E1, E2>, Sp> const& e) {
  auto expr = expr::get_enclosed_expression(e);
  auto a = apply_operators(expr.a);
  auto b = apply_operators(expr.b);

  auto lhs = apply_operators(expr::make_derivative<Dd>(a, e.solver) * b);
  auto rhs = apply_operators(a * expr::make_derivative<Dd>(b, e.solver));
  return expr::coeff(e) * (lhs - rhs) / (b * b);
}

////! Implementation of the product rule for symbolic derivatives.
// template<size_t O, typename V, typename A1, typename A2, typename G>
// auto apply_operators(OpDerivative<std::index_sequence<O>, V,
// OpOperatorChain<A1, A2>, SymbolicDerivative<G>> const& e)
//{
//	auto&& expr = expr::get_enclosed_expression(e);
//	auto lhs = apply_operators(expr::make_derivative<O, G>(expr.f,
// e.solver))(expr.g); 	auto rhs =
// expr.f(apply_operators(expr::make_derivative<O, G>(expr.g, e.solver)));
// return expr::coeff(e) * (lhs + rhs);
// }

////! Implementation of the product rule for symbolic derivatives.
// template<size_t O, typename V, typename A1, typename A2, typename E, typename
// G> auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpChain<A1,
// A2, E>, SymbolicDerivative<G>> const& e)
//{
//	auto&& expr = apply_operators(expr::get_enclosed_expression(e));
//	auto lhs = apply_operators(expr::make_derivative<O, G>(expr.combination,
// e.solver)) * expr::get_enclosed_expression(expr); 	auto rhs =
// expr.combination
//* apply_operators(expr::make_derivative<O,
// G>(expr::get_enclosed_expression(expr), e.solver)); 	return expr::coeff(e) *
//(lhs + rhs);
// }

////! Implementation of the product rule for symbolic derivatives.
// template<size_t O, typename V, typename A1, typename A2, typename G>
// auto apply_operators(OpDerivative<std::index_sequence<O>, V,
// OpOperatorCombination<A1, A2>, SymbolicDerivative<G>> const& e)
//{
//	auto&& expr = expr::get_enclosed_expression(e);
//	return expr::coeff(e) * (apply_operators(expr::make_derivative<O,
// G>(expr.f, e.solver)) + apply_operators(expr::make_derivative<O, G>(expr.g,
// e.solver)));
// }

////! Implementation of the product rule for symbolic derivatives.
// template<size_t O, typename V, typename A1, typename A2, typename E, typename
// G> auto apply_operators(OpDerivative<std::index_sequence<O>, V,
// OpCombination<A1, A2, E>, SymbolicDerivative<G>> const& e)
//{
//	auto&& expr = apply_operators(expr::get_enclosed_expression(e));
//	auto lhs = apply_operators(expr::make_derivative<O, G>(expr.combination,
// e.solver)) * expr::get_enclosed_expression(expr); 	auto rhs =
// expr.combination
//* apply_operators(expr::make_derivative<O,
// G>(expr::get_enclosed_expression(expr), e.solver)); 	return expr::coeff(e) *
//(lhs + rhs);
// }

//! Specialization based on expr::grid_dim.
template <typename V, typename E, typename F>
auto apply_operators(OpFunction<V, E, F, void> const& e) {
  return expr::coeff(e) *
         expr::make_function(apply_operators(expr::get_enclosed_expression(e)),
                             e.name, e.f);
};

//! Specialization based on expr::grid_dim.
template <typename V, typename E, typename F, typename Arg0, typename... Args>
auto apply_operators(OpFunction<V, E, F, Arg0, Args...> const& e) {
  return expr::coeff(e) *
         expr::make_function(apply_operators(expr::get_enclosed_expression(e)),
                             e.name, e.f, e.args);
};

//! Specialization based on expr::grid_dim.
template <auto f, typename V, typename E>
auto apply_operators(OpFunctionApply<f, V, E> const& e) {
  return expr::coeff(e) * expr::make_function<f>(apply_operators(
                              expr::get_enclosed_expression(e)));
};

//! Specialization based on expr::grid_dim.
template <typename V, typename E, typename T>
auto apply_operators(OpIntegral<V, E, T> const& e) {
  return expr::coeff(e) *
         expr::make_integral(apply_operators(expr::get_enclosed_expression(e)),
                             e.domain);
};

//! Implementation of the product rule for terms.
// template<typename Dd, typename Sp, typename V, typename V0, typename G0,
// expr::exp_key_t X0, typename G1, expr::exp_key_t X1, typename... Gs,
// expr::exp_key_t... Xs, 	typename std::enable_if_t<(OpDerivative<Dd, V,
// OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>, Sp>::order > 0),
// int>> auto apply_operators(OpDerivative<Dd, V, OpTerms<V0, Term<G0, X0>,
// Term<G1, X1>, Term<Gs, Xs>...>, Sp> const& e)
//{
//	auto&& terms = expr::get_enclosed_expression(e);
//	auto coeff = expr::coeff(terms);
//	auto a = OpTerms(OpIdentity{}, expr::get<1>(terms));
//	auto b = OpTerms(OpIdentity{}, *static_cast<OpTerms<Term<G1, X1>,
// Term<Gs, Xs>...> const*>(&terms));
//
//	auto lhs = apply_operators(expr::make_derivative<Dd>(a, e.solver)) * b;
//	auto rhs = a * apply_operators(expr::make_derivative<Dd>(b, e.solver));
//	return (coeff * e.value) * (lhs + rhs);
// }

/*
//! Implementation of the product rule for terms.
template<typename Dd, typename Sp, typename V, typename V0, typename G0,
expr::exp_key_t X0, typename std::enable_if_t<(OpDerivative<Dd, V, OpTerms<V0,
Term<G0, X0>>, Sp>::order > 0
                && (expr::Xk_t<X0>::N > expr::Xk_t<X0>::D &&
!expr::Xk_t<X0>::sign)), int>> auto apply_operators(OpDerivative<Dd, V,
OpTerms<V0, Term<G0, X0>>, Sp> const& e)
{
        auto&& term = expr::get_enclosed_expression(e);

        constexpr size_t N = expr::_Xk_t<X0>::N;
        constexpr size_t D = expr::_Xk_t<X0>::D;
        constexpr bool sign = expr::_Xk_t<X0>::sign;

        auto dvar = symphas::internal::to_term_element<expr::Xk<N - D, D,
sign>>(expr::get<1>(term).data()); auto var1 =
symphas::internal::to_term_element(expr::get<1>(term).data());

        auto coeff = expr::coeff(term);
        auto dterm = OpTerms(coeff, dvar);

        auto lhs = apply_operators(expr::make_derivative<Dd>(dvar, e.solver)) *
var1; auto rhs = dvar * apply_operators(expr::make_derivative<Dd>(var1,
e.solver)); return (lhs + rhs);
}*/

template <typename Dd, typename Sp, typename V, typename V0, typename E,
          typename std::enable_if_t<
              (OpDerivative<Dd, V, OpExponential<V0, E>, Sp>::order > 0), int>>
auto apply_operators(OpDerivative<Dd, V, OpExponential<V0, E>, Sp> const& e) {
  auto&& expr = expr::get_enclosed_expression(e);
  auto&& pow = expr::get_enclosed_expression(expr);
  return expr::coeff(e) * expr *
         apply_operators(expr::make_derivative<Dd>(pow, e.solver));
}

namespace {
template <typename A1, typename A2, typename E,
          typename std::enable_if_t<expr_has_deriv<E>::value, int> = 0>
auto apply_operators_chain(OpChain<A1, A2, E> const& e) {
  return apply_operators(expr::distribute_operator(
      e.combination.f,
      apply_operators(expr::distribute_operator(
          e.combination.g,
          apply_operators(expr::get_enclosed_expression(e))))));
}

template <typename A1, typename A2, typename E,
          typename std::enable_if_t<!expr_has_deriv<E>::value, int> = 0>
auto apply_operators_chain(OpChain<A1, A2, E> const& e) {
  return apply_operators(
      expr::expand_operator(e.combination, expr::get_enclosed_expression(e)));
}

template <typename V, typename E>
auto apply_operators_chain(OpChain<V, OpIdentity, E> const& e) {
  return expr::make_operator_chain(
      apply_operators(e.combination.f),
      OpIdentity{})(apply_operators(expr::get_enclosed_expression(e)));
}

template <typename A1, typename A2, typename E,
          typename std::enable_if_t<expr_has_deriv<E>::value, int> = 0>
auto apply_operators_combination(OpCombination<A1, A2, E> const& e) {
  return apply_operators(expr::distribute_operator(
             e.combination.f,
             apply_operators(expr::get_enclosed_expression(e)))) +
         apply_operators(expr::distribute_operator(
             e.combination.g,
             apply_operators(expr::get_enclosed_expression(e))));
}

template <typename A1, typename A2, typename E,
          typename std::enable_if_t<!expr_has_deriv<E>::value, int> = 0>
auto apply_operators_combination(OpCombination<A1, A2, E> const& e) {
  return apply_operators(
      expr::expand_operator(e.combination, expr::get_enclosed_expression(e)));
}

template <
    Axis ax1, Axis ax2, size_t O, typename V, typename Sp, typename E,
    typename std::enable_if_t<(!expr_has_deriv<E>::value && expr::is_const<V> &&
                               expr::grid_dim<E>::value == 3),
                              int> = 0>
auto apply_operators_combination(
    OpCombination<OpOperatorDirectionalDerivative<ax1, O, V, Sp>,
                  OpOperatorDirectionalDerivative<ax2, O, V, Sp>, E> const& e) {
  return e;
}

template <typename... Es, size_t... Is>
auto apply_operators_adds(OpAdd<Es...> const& e, std::index_sequence<Is...>) {
  return (apply_operators(expr::get<Is>(e)) + ...);
}

template <typename Dd, typename Sp, typename... Es, size_t... Is>
auto apply_operators_adds(Sp const& solver, OpAdd<Es...> const& e,
                          std::index_sequence<Is...>) {
  return (apply_operators(expr::make_derivative<Dd>(expr::get<Is>(e), solver)) +
          ...);
}

template <size_t O, typename G0, typename... Es, size_t... Is>
auto apply_operators_adds(SymbolicDerivative<G0> const& solver,
                          OpAdd<Es...> const& e, std::index_sequence<Is...>) {
  return (
      apply_operators(expr::make_derivative<O, G0>(expr::get<Is>(e), solver)) +
      ...);
}

template <typename A1, typename A2>
auto apply_operators_chain(OpOperatorChain<A1, A2> const& e) {
  return e;
}

template <typename A1, typename B1, typename B2>
auto apply_operators_chain(
    OpOperatorChain<A1, OpOperatorCombination<B1, B2>> const& e) {
  return apply_operators(e.f(e.g.f)) + apply_operators(e.f(e.g.g));
}

template <typename A1, typename B1, typename B2, typename E>
auto apply_operators_chain(
    OpOperatorChain<A1, OpCombination<B1, B2, E>> const& e) {
  return apply_operators(distribute_operator(e.f, e.g));
}

template <typename A1, typename A2>
auto apply_operators_chain(
    OpOperatorChain<A1, OpOperatorChain<OpIdentity, A2>> const& e) {
  return expr::make_operator_chain(OpIdentity{}, apply_operators(e.f(e.g.g)));
}

template <typename A1, typename A2>
auto apply_operators_chain(
    OpOperatorChain<OpOperatorChain<OpIdentity, A2>,
                    OpOperatorChain<OpIdentity, A2>> const& e) {
  return expr::make_operator_chain(OpIdentity{},
                                   apply_operators((e.f.g)(e.g.g)));
}

template <typename A1, typename A2>
auto apply_operators_chain(
    OpOperatorChain<OpOperatorChain<OpIdentity, A2>, A2> const& e) {
  return expr::make_operator_chain(OpIdentity{}, apply_operators((e.f.g)(e.g)));
}

template <typename E>
auto apply_operators_chain(OpExpression<E> const& e) {
  return apply_operators(*static_cast<E const*>(&e));
}

template <typename E>
auto apply_operators_chain(OpOperator<E> const& e) {
  return apply_operators(*static_cast<E const*>(&e));
}
}  // namespace

template <typename A1, typename A2>
auto apply_operators(OpOperatorChain<A1, A2> const& e) {
  return apply_operators_chain(apply_operators(e.f)(apply_operators(e.g)));
}

template <typename A2>
auto apply_operators(OpOperatorChain<OpIdentity, A2> const& e) {
  return expr::make_operator_chain(OpIdentity{}, apply_operators(e.g));
}

// template<typename A1, typename B1, typename B2>
// auto apply_operators(OpOperatorChain<A1, OpOperatorChain<B1, B2>> const& e)
//{
//	return apply_operators(((apply_operators(e.f)(apply_operators(e.g.f))) *
// apply_operators(e.g.g)))
//		+ apply_operators((apply_operators(e.g.f) *
//(apply_operators(e.f)(apply_operators(e.g.g)))));
// }

// template<typename A1, typename B1, typename B2>
// auto apply_operators(OpOperatorChain<A1, OpBinaryMul<B1, B2>> const& e)
//{
//	return apply_operators(((apply_operators(e.f)(apply_operators(e.g.a))) *
// apply_operators(e.g.b)))
//		+ apply_operators((apply_operators(e.g.a) *
//(apply_operators(e.f)(apply_operators(e.g.b)))));
// }

template <typename A1, typename A2>
auto apply_operators(OpOperatorCombination<A1, A2> const& e) {
  return apply_operators(e.f) + apply_operators(e.g);
}

template <typename A1, typename A2, typename E>
auto apply_operators(OpChain<A1, A2, E> const& e) {
  return apply_operators_chain(e);
}

template <typename A1, typename A2, typename E>
auto apply_operators(OpCombination<A1, A2, E> const& e) {
  return apply_operators_combination(e);
}

namespace {

template <typename Sp, typename Seq, typename Axes>
struct mixed_derivative_type {};

template <typename Sp, size_t... Os, Axis... axs>
struct mixed_derivative_type<Sp, std::index_sequence<Os...>,
                             symphas::lib::axis_list<axs...>> {
  template <Axis ax0, Axis ax1, size_t O>
  static const size_t p1 = (ax0 == ax1) ? (O % 2) : 0;

  template <Axis ax0, Axis ax1, size_t O>
  static const size_t pO = (ax0 == ax1) ? (O - (O % 2)) : 0;

  template <size_t O, Axis axd>
  // using type = typename Solver<Sp>::template mixed_derivative<(Os + pO<axm,
  // axs, O> + p1<axs, O>)...>;
  using type = typename Solver<Sp>::template mixed_derivative<(
      Os + pO<axd, axs, O> + p1<axd, axs, O>)...>;
};

template <typename Sp>
struct combine_mixed_derivatives {
  using solver_type = Solver<Sp>;

  template <size_t... Os>
  using mixed_derivative =
      typename solver_type::template mixed_derivative<Os...>;
  template <Axis ax, size_t O>
  using directional_derivative =
      typename solver_type::template directional_derivative<ax, O>;
  template <Axis ax, size_t O>
  using derivative = typename solver_type::template derivative<ax, O>;

  template <typename E, Axis ax1, size_t O1, Axis ax2, size_t O2, Axis... axs>
  auto operator()(OpExpression<E> const& enclosed, Sp const& solver,
                  directional_derivative<ax1, O1>,
                  directional_derivative<ax2, O2>,
                  symphas::lib::axis_list<axs...>) const {
    if constexpr (ax1 == ax2) {
      using Dd =
          typename Solver<Sp>::template directional_derivative<ax1, O1 + O2>;
      return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed),
                                       solver);
    } else {
      using Dd = typename Solver<Sp>::template mixed_derivative<(
          (ax1 == axs)   ? O1
          : (ax2 == axs) ? O2
                         : 0)...>;
      return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed),
                                       solver);
    }
  }

  template <typename E, Axis ax1, Axis ax2, Axis... axs>
  auto operator()(OpExpression<E> const& enclosed, Sp const& solver,
                  derivative<ax1, 1>, derivative<ax2, 1>,
                  symphas::lib::axis_list<axs...>) const {
    if constexpr (ax1 == ax2) {
      using Dd = typename Solver<Sp>::template directional_derivative<ax1, 2>;
      return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed),
                                       solver);
    } else {
      using Dd = typename Solver<Sp>::template mixed_derivative<(
          (ax1 == axs || ax2 == axs) ? 1 : 0)...>;
      return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed),
                                       solver);
    }
  }

  //! derivative<ax, O> composed with another derivative D2.
  /*!
   * derivative<ax, O> is lap^(O/2) for even O and d/dax lap^((O-1)/2) for
   * odd O. One laplacian, sum_i d2/dx_i2, is composed with D2 at a time and
   * the rest, derivative<ax, O - 2>, is applied to the result; d/dax is
   * composed last. Composing every axis with the full order O at once, as
   * was done before, is only correct for O = 2, and folding an order-0
   * directional derivative over the axes counted the term once per axis.
   */
  template <Axis ax, size_t O, typename E, typename D2, Axis... axs>
  auto compose_derivative(OpExpression<E> const& enclosed, Sp const& solver,
                          D2, symphas::lib::axis_list<axs...>) const {
    if constexpr (O == 1) {
      return this->operator()(*static_cast<E const*>(&enclosed), solver,
                              directional_derivative<ax, 1>{}, D2{},
                              symphas::lib::axis_list<axs...>{});
    } else {
      auto lap = (this->operator()(*static_cast<E const*>(&enclosed), solver,
                                   directional_derivative<axs, 2>{}, D2{},
                                   symphas::lib::axis_list<axs...>{}) +
                  ...);
      if constexpr (O == 2) {
        return lap;
      } else {
        return apply_operators(
            expr::make_derivative<derivative<ax, O - 2>>(lap, solver));
      }
    }
  }

  template <typename E, Axis ax1, size_t O1, Axis ax2, size_t O2, Axis... axs>
  auto operator()(OpExpression<E> const& enclosed, Sp const& solver,
                  derivative<ax1, O1>, derivative<ax2, O2>,
                  symphas::lib::axis_list<axs...>) const {
    return compose_derivative<ax1, O1>(*static_cast<E const*>(&enclosed),
                                       solver, derivative<ax2, O2>{},
                                       symphas::lib::axis_list<axs...>{});
  }

  template <size_t... O1s>
  struct deduction_redirection {
    template <typename E, size_t... O2s, Axis... axs>
    auto operator()(OpExpression<E> const& enclosed, Sp const& solver,
                    mixed_derivative<O2s...>,
                    symphas::lib::axis_list<axs...>) const {
      using Dd = typename Solver<Sp>::template mixed_derivative<(O1s + O2s)...>;
      return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed),
                                       solver);
    }
  };

  template <typename E, size_t... O1s, typename mixed_other, Axis... axs>
  auto operator()(OpExpression<E> const& enclosed, Sp const& solver,
                  mixed_derivative<O1s...>, mixed_other,
                  symphas::lib::axis_list<axs...>) const {
    return deduction_redirection<O1s...>{}(*static_cast<E const*>(&enclosed),
                                           solver, mixed_other{},
                                           symphas::lib::axis_list<axs...>{});
  }

  template <typename E, Axis ax, size_t O, size_t... O2s, Axis... axs>
  auto operator()(OpExpression<E> const& enclosed, Sp const& solver,
                  directional_derivative<ax, O>, mixed_derivative<O2s...>,
                  symphas::lib::axis_list<axs...>) const {
    using Dd = typename Solver<Sp>::template mixed_derivative<(
        (ax == axs) ? O + O2s : O2s)...>;
    return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver);
  }

  template <typename E, Axis ax, size_t O, size_t... O1s, Axis... axs>
  auto operator()(OpExpression<E> const& enclosed, Sp const& solver,
                  mixed_derivative<O1s...>, directional_derivative<ax, O>,
                  symphas::lib::axis_list<axs...>) const {
    using Dd = typename Solver<Sp>::template mixed_derivative<(
        (ax == axs) ? O + O1s : O1s)...>;
    return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver);
  }

  template <typename E, Axis ax, size_t O, size_t... O2s, Axis... axs>
  auto operator()(OpExpression<E> const& enclosed, Sp const& solver,
                  derivative<ax, O>, mixed_derivative<O2s...>,
                  symphas::lib::axis_list<axs...>) const {
    return compose_derivative<ax, O>(*static_cast<E const*>(&enclosed), solver,
                                     mixed_derivative<O2s...>{},
                                     symphas::lib::axis_list<axs...>{});
  }

  template <typename E, Axis ax, size_t O, size_t... O1s, Axis... axs>
  auto operator()(OpExpression<E> const& enclosed, Sp const& solver,
                  mixed_derivative<O1s...>, derivative<ax, O>,
                  symphas::lib::axis_list<axs...>) const {
    return this->operator()(*static_cast<E const*>(&enclosed), solver,
                            derivative<ax, O>{}, mixed_derivative<O1s...>{},
                            symphas::lib::axis_list<axs...>{});
  }

  template <typename E, Axis ax1, size_t O1, Axis ax2, size_t O2, Axis... axs>
  auto operator()(OpExpression<E> const& enclosed, Sp const& solver,
                  derivative<ax1, O1>, directional_derivative<ax2, O2>,
                  symphas::lib::axis_list<axs...>) const {
    return compose_derivative<ax1, O1>(*static_cast<E const*>(&enclosed),
                                       solver, directional_derivative<ax2, O2>{},
                                       symphas::lib::axis_list<axs...>{});
  }

  template <typename E, Axis ax1, size_t O1, Axis ax2, size_t O2, Axis... axs>
  auto operator()(OpExpression<E> const& enclosed, Sp const& solver,
                  directional_derivative<ax1, O1>, derivative<ax2, O2>,
                  symphas::lib::axis_list<axs...>) const {
    // using mixed_type = typename Solver<Sp>::template mixed_derivative<((axs
    // == ax1) ? O1 : (axs == ax2) ? O2 : 0)...>; return
    // expr::make_derivative<mixed_type>(*static_cast<E const*>(&enclosed),
    // solver);
    return this->operator()(
        *static_cast<E const*>(&enclosed), solver, derivative<ax2, O2>{},
        directional_derivative<ax1, O1>{}, symphas::lib::axis_list<axs...>{});
  }
};

template <typename Dd1, typename Dd2, typename V1, typename V2, typename E,
          typename Sp1, typename Sp2>
auto apply_operator_derivative_nested(
    OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp1>, Sp2> const& e) {
  return e;
}

template <typename Dd1, typename V1, typename E, typename Sp2>
auto apply_operator_derivative_nested(OpDerivative<Dd1, V1, E, Sp2> const& e) {
  return expr::apply_operators(e);
}

// FE models with reversible vector coupling (e.g. MF_FE's `grad * DF(2)`)
// produce an OpAdd of per-axis OpDerivatives (one for each vector
// component) before the outer apply_operators reaches here. The two
// overloads above are written for a single OpDerivative; without an
// OpAdd specialization the dispatch fails. Distribute over the sum.
template <typename... Es>
auto apply_operator_derivative_nested(OpAdd<Es...> const& e) {
  return expr::apply_operators(e);
}
}  // namespace

//! Distribute operators so they are applied to individual expressions.
/*!
 * For expressions that are derivatives of derivatives, the outermost
 * derivatives might need to be distributed to the rest of the expression.
 * Additionally, all derivatives which are applied to expressions that
 * are linear combinations are distributed.
 *
 * \param e The expression which is distributed.
 */
// template<typename Dd1, size_t N, typename V1, typename V2, typename E,
// typename Sp1, typename Sp2> auto apply_operators(OpDerivative<Dd1, V1,
// OpDerivative<std::index_sequence<N>, V2, E, SymbolicDerivative<Sp1>>, Sp2>
// const& e);

// template<typename Dd1, size_t N, typename V1, typename V2, typename E,
// typename Sp1, typename Sp2> auto apply_operators(OpDerivative<Dd1, V1,
// OpDerivative<std::index_sequence<N>, V2, E, SymbolicDerivative<Sp1>>, Sp2>
// const& e)
//{
//	return e;
// }

template <typename Dd1, typename Dd2, typename V1, typename V2, typename E,
          typename Sp1, typename Sp2>
auto apply_operators(
    OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp1>, Sp2> const& e) {
  return expr::coeff(e) *
         apply_operator_derivative_nested(expr::make_derivative<Dd1>(
             expr::apply_operators(expr::get_enclosed_expression(e)),
             e.solver));
}

template <bool parity1, bool parity2, size_t order1, size_t order2, Axis ax1,
          Axis ax2>
struct combine_derivatives {
  template <typename Sp, typename E, Axis... axs>
  auto operator()(Sp const& solver, OpExpression<E> const& enclosed,
                  symphas::lib::axis_list<axs...>) {
    auto dd1 =
        expr::make_operator_directional_derivative<ax1, parity1>(solver) *
        break_up_derivative<order1 - parity1>(
            solver, symphas::lib::axis_list<axs...>{});
    auto dd2 =
        expr::make_operator_directional_derivative<ax2, parity2>(solver) *
        break_up_derivative<order2 - parity2>(
            solver, symphas::lib::axis_list<axs...>{});

    return apply_operators((dd1 * dd2) * *static_cast<E const*>(&enclosed));
  }
};

template <size_t order1, size_t order2, Axis ax1, Axis ax2>
struct combine_derivatives<0, 0, order1, order2, ax1, ax2> {
  template <typename Sp, typename E, Axis... axs>
  auto operator()(Sp const& solver, OpExpression<E> const& enclosed,
                  symphas::lib::axis_list<axs...>) {
    using Dd =
        typename Solver<Sp>::template derivative<Axis::X, order1 + order2>;
    return apply_operators(
        expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver));
    // return apply_operators(expr::make_operator_derivative<order1 +
    // order2>(solver)(OpVoid{}));
  }
};

template <size_t order1, size_t order2, Axis ax1, Axis ax2>
struct combine_derivatives<1, 0, order1, order2, ax1, ax2> {
  template <typename Sp, typename E, Axis... axs>
  auto operator()(Sp const& solver, OpExpression<E> const& enclosed,
                  symphas::lib::axis_list<axs...>) {
    using Dd = typename Solver<Sp>::template derivative<ax1, order1 + order2>;
    return apply_operators(
        expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver));
  }
};

template <size_t order1, size_t order2, Axis ax1, Axis ax2>
struct combine_derivatives<0, 1, order1, order2, ax1, ax2> {
  template <typename Sp, typename E, Axis... axs>
  auto operator()(Sp const& solver, OpExpression<E> const& enclosed,
                  symphas::lib::axis_list<axs...>) {
    using Dd = typename Solver<Sp>::template derivative<ax2, order1 + order2>;
    return apply_operators(
        expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver));
  }
};

template <typename Dd1, typename Dd2, typename V1, typename V2, typename E,
          typename Sp>
auto apply_operators(
    OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp>, Sp> const& e) {
  using d1t = OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp>, Sp>;
  using d2t = OpDerivative<Dd2, V2, E, Sp>;

  constexpr size_t D = expr::grid_dim<d1t>::value;

  auto enclosed1 = expr::get_enclosed_expression(e);
  auto enclosed2 = expr::get_enclosed_expression(enclosed1);
  auto enclosed = apply_operators(enclosed2);
  auto coeff = expr::coeff(e) * expr::coeff(enclosed1);

  auto axis_list = symphas::lib::make_axis_list<D>();

  if constexpr (!Dd1::is_mixed && !Dd2::is_mixed && !Dd1::is_directional &&
                !Dd2::is_directional) {
    return coeff *
           combine_derivatives<d1t::order % 2, d2t::order % 2, d1t::order,
                               d2t::order, d1t::axis, d2t::axis>{}(
               e.solver, enclosed, axis_list);
  } else {
    return coeff * apply_operators(combine_mixed_derivatives<Sp>{}(
                       enclosed, e.solver, Dd1{}, Dd2{}, axis_list));
  }
}

namespace {
template <size_t O, typename V, typename E, typename GG, typename E0>
auto apply_operators_deriv(
    OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& d,
    OpExpression<E0> const& e) {
  return apply_operators(
      expr::make_derivative<O, GG>(*static_cast<E0 const*>(&e), d.solver));
}

template <size_t O, typename V, typename E, typename GG, typename A1,
          typename A2, typename E0>
auto apply_operators_deriv(
    OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& d,
    OpCombination<A1, A2, E0> const& e) {
  return apply_operators(
             expr::make_derivative<O, GG>(e.combination.f(e.e), d.solver)) +
         apply_operators(
             expr::make_derivative<O, GG>(e.combination.g(e.e), d.solver));
}

template <size_t O, typename V, typename E, typename GG, typename E0>
auto apply_operators_deriv(
    OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& d,
    OpOperatorChain<OpIdentity, E0> const& e) {
  return apply_operators(expr::make_derivative<O, GG>(e.g));
}

template <size_t O, typename V1, typename V2, typename E, typename GG,
          typename E0>
auto apply_operators_deriv(OpDerivative<std::index_sequence<O>, V1, E,
                                        SymbolicDerivative<GG>> const& d,
                           OpChain<V2, OpIdentity, E0> const& e) {
  return apply_operators(
      expr::make_derivative<O, GG>(expr::get_enclosed_expression(e)));
}

template <size_t O, typename V, typename E, typename GG, typename E0>
auto apply_operators_deriv(
    OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& d,
    OpOperator<E0> const& e) {
  return OpVoid{};
}

// template<typename Dd, typename V, typename E, typename Sp, typename E0>
// auto apply_operators_deriv(OpDerivative<Dd, V, E, Sp> const& d,
// OpOperator<E0> const& e)
//{
//	return apply_operators(expr::make_derivative<Dd>(*static_cast<E0
// const*>(&e), d.solver));
// }

template <typename Dd, typename V, typename E, typename Sp, typename E0>
auto apply_operators_deriv(OpDerivative<Dd, V, E, Sp> const& d,
                           OpOperator<E0> const& e);

template <typename Dd, typename V, typename E, typename Sp, typename E0,
          typename std::enable_if_t<expr_has_deriv<E>::value, int> = 0>
auto apply_operators_deriv(OpDerivative<Dd, V, E, Sp> const& d,
                           OpExpression<E0> const& e) {
  return apply_operators(
      expr::make_derivative<Dd>(*static_cast<E0 const*>(&e), d.solver));
}
//
// template <typename Dd, typename V, typename E, typename Sp, typename E0,
//          typename std::enable_if_t<expr_has_deriv<E>::value, int> = 0>
// auto apply_operators_deriv(OpDerivative<Dd, V, E, Sp> const& d,
//                           OpCombination<A, B, E> const& e) {
//  return apply_operators(
//      expr::make_derivative<Dd>(*static_cast<E0 const*>(&e), d.solver));
//}

template <typename Dd, typename V, typename E, typename Sp, typename E0,
          typename std::enable_if_t<!expr_has_deriv<E>::value, int> = 0>
auto apply_operators_deriv(OpDerivative<Dd, V, E, Sp> const& d,
                           OpExpression<E0> const& e) {
  return expr::make_derivative<Dd>(*static_cast<E0 const*>(&e), d.solver);
}
}  // namespace

template <typename Dd, typename V, typename E, typename Sp, typename>
auto apply_operators(OpDerivative<Dd, V, E, Sp> const& e) {
  return expr::coeff(e) *
         apply_operators_deriv(
             e, apply_operators(expr::get_enclosed_expression(e)));
}

template <typename Dd, typename V, typename... Es, typename Sp>
auto apply_operators(OpDerivative<Dd, V, OpAdd<Es...>, Sp> const& e) {
  return expr::coeff(e) *
         apply_operators_adds<Dd>(e.solver, expr::get_enclosed_expression(e),
                                  std::make_index_sequence<sizeof...(Es)>{});
}

template <size_t O, typename V, typename... Es, typename G0>
auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpAdd<Es...>,
                                  SymbolicDerivative<G0>> const& e) {
  return expr::coeff(e) *
         apply_operators_adds<O>(e.solver, expr::get_enclosed_expression(e),
                                 std::make_index_sequence<sizeof...(Es)>{});
}

template <typename... Es>
auto apply_operators(OpAdd<Es...> const& e) {
  return apply_operators_adds(e, std::make_index_sequence<sizeof...(Es)>{});
}

namespace {

template <typename A, typename B, typename E2>
auto apply_operators_mul(OpOperatorCombination<A, B> const& combination,
                         OpExpression<E2> const& b) {
  return apply_operators(combination.f * (*static_cast<E2 const*>(&b))) +
         apply_operators(combination.g * (*static_cast<E2 const*>(&b)));
}

template <typename A, typename B, typename E2>
auto apply_operators_mul(OpOperatorCombination<A, B> const& combination,
                         OpOperator<E2> const& b) {
  return apply_operators(combination.f * (*static_cast<E2 const*>(&b))) +
         apply_operators(combination.g * (*static_cast<E2 const*>(&b)));
}

//! Apply the chain operation to an expression.
template <typename A1, typename A2, typename E>
auto apply_operators_mul(OpOperatorChain<A1, A2> const& combination,
                         OpExpression<E> const& b) {
  return apply_operators(combination.f(
      apply_operators(combination.g * *static_cast<E const*>(&b))));
}

//! Apply the chain operation to an expression.
template <typename A1, typename A2, typename E>
auto apply_operators_mul(OpOperatorChain<A1, A2> const& combination,
                         OpOperator<E> const& b) {
  return apply_operators(combination.f(
      apply_operators(combination.g * *static_cast<E const*>(&b))));
}

template <typename E1, typename E2>
auto apply_operators_mul(OpOperator<E1> const& a, OpExpression<E2> const& b) {
  return apply_operators((*static_cast<E1 const*>(&a)) *
                         (*static_cast<E2 const*>(&b)));
}

template <typename E1, typename E2>
auto apply_operators_mul(OpExpression<E1> const& a, OpOperator<E2> const& b) {
  return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
}

template <typename E1, typename E2>
auto apply_operators_mul(OpOperator<E1> const& a, OpOperator<E2> const& b) {
  return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
}

template <typename E1, typename E2>
auto apply_operators_mul(OpExpression<E1> const& a, OpExpression<E2> const& b) {
  return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
}

template <typename... Es, typename E2>
auto apply_operators_mul(OpAdd<Es...> const& a, OpExpression<E2> const& b) {
  return apply_operators(a * *static_cast<E2 const*>(&b));
}

template <typename... Es, typename E2>
auto apply_operators_mul(OpAdd<Es...> const& a, OpOperator<E2> const& b) {
  return apply_operators(a * *static_cast<E2 const*>(&b));
}

template <typename... Es, typename E2>
auto apply_operators_mul(OpAdd<Es...> const& a,
                         OpOperatorChain<OpIdentity, E2> const& b) {
  return (a * b);
}

template <size_t O, typename V, typename Sp, typename E2>
auto apply_operators_mul(OpOperatorDerivative<O, V, Sp> const& a,
                         OpOperatorChain<OpIdentity, E2> const& b) {
  return expr::make_operator_chain(OpIdentity{}, apply_operators(a * b.g));
}

//! Apply the chain operation to an expression.
template <typename A2, typename E>
auto apply_operators_mul(OpOperatorChain<OpIdentity, A2> const& a,
                         OpExpression<E> const& b) {
  return (a * *static_cast<E const*>(&b));
}

template <typename E1, typename E2>
auto apply_operators_mul(OpOperatorChain<OpIdentity, E1> const& a,
                         OpOperatorChain<OpIdentity, E2> const& b) {
  return expr::make_operator_chain(OpIdentity{}, apply_operators(a.g * b.g));
}

template <typename V, typename E1, typename E2>
auto apply_operators_mul(OpOperator<E1> const& a,
                         OpChain<V, OpIdentity, E2> const& b) {
  return b.combination.f * apply_operators((*static_cast<E1 const*>(&a)) *
                                           expr::get_enclosed_expression(b));
}
}  // namespace

template <typename E1, typename E2>
auto apply_operators(OpBinaryMul<E1, E2> const& e) {
  return apply_operators_mul(apply_operators(e.a), apply_operators(e.b));
}

template <typename E1, typename E2>
auto apply_operators(OpBinaryDiv<E1, E2> const& e) {
  return apply_operators(e.a) / apply_operators(e.b);
}

// precondition dot product
// the dot product is constructed as an OpChain(OpOperatorChain(V, OpIdentity),
// E). when a dot product appears in the addition below, we need to combine them
// when adding.

template <typename V1, typename V2, typename Dd1, typename Dd2, typename E,
          typename Sp>
auto operator+(OpChain<V1, OpIdentity,
                       OpDerivative<Dd1, OpTensor<V2, 0, 2>, E, Sp>> const& a,
               OpChain<V1, OpIdentity,
                       OpDerivative<Dd2, OpTensor<V2, 1, 2>, E, Sp>> const& b) {
  return OpChain(a.combination, expr::get_enclosed_expression(a) +
                                    expr::get_enclosed_expression(b));
}

template <typename V1, typename V2, typename Dd1, typename Dd2, typename E,
          typename Sp>
auto operator+(OpChain<V1, OpIdentity,
                       OpDerivative<Dd1, OpTensor<V2, 1, 3>, E, Sp>> const& a,
               OpChain<V1, OpIdentity,
                       OpDerivative<Dd2, OpTensor<V2, 2, 3>, E, Sp>> const& b) {
  return OpChain(a.combination, expr::get_enclosed_expression(a) +
                                    expr::get_enclosed_expression(b));
}

template <typename V1, typename V2, typename Dd1, typename Dd2, typename Dd3,
          typename E, typename Sp>
auto operator+(
    OpChain<V1, OpIdentity, OpDerivative<Dd1, OpTensor<V2, 0, 3>, E, Sp>> const&
        a,
    OpChain<V1, OpIdentity,
            OpAdd<OpDerivative<Dd2, OpTensor<V2, 1, 3>, E, Sp>,
                  OpDerivative<Dd3, OpTensor<V2, 2, 3>, E, Sp>>> const& b) {
  return OpChain(a.combination, expr::get_enclosed_expression(a) +
                                    expr::get_enclosed_expression(b));
}

namespace {
template <Axis ax, typename G, typename Sp>
using grad_term_t =
    OpDerivative<typename Solver<Sp>::template derivative<ax, 1>, OpIdentity,
                 OpTerm<OpIdentity, G>, Sp>;

template <typename E, typename Sp, typename G, Axis... axs, size_t... Is>
auto euler_lagrange_deriv(
    E const& e, Sp const& solver,
    symphas::lib::types_list<G, symphas::lib::axis_list<axs...>>,
    std::index_sequence<Is...>) {
  return expr::apply_operators(
      expr::make_operator_derivative<1>(solver) *
      (expr::apply_operators(
           expr::make_column_vector<Is, sizeof...(Is)>() *
           expr::make_derivative<1, grad_term_t<axs, G, Sp>>(e)) +
       ...));
}

template <typename G, size_t D, typename E, typename Sp>
auto euler_lagrange_deriv(E const& e, Sp const& solver) {
  return euler_lagrange_deriv(e, solver, symphas::lib::make_axis_list<D, G>(),
                              std::make_index_sequence<D>{});
}

template <typename G, size_t D, typename E>
auto euler_lagrange_deriv(E const& e, int) {
  return OpVoid{};
}
}  // namespace

template <typename G, typename E, typename Sp>
auto euler_lagrange_apply(SymbolicFunctionalDerivative<G> const& symbol,
                          OpExpression<E> const& e, Sp const& solver) {
  constexpr size_t D = expr::grid_dim<G>::value;

  auto interface_term =
      euler_lagrange_deriv<G, D>(*static_cast<E const*>(&e), solver);
  auto bulk_term = expr::apply_operators(
      expr::make_derivative<1, G>(*static_cast<E const*>(&e)));
  return bulk_term - interface_term;
}

template <typename G, typename E, typename Sp>
auto euler_lagrange_apply(
    SymbolicFunctionalDerivative<DynamicVariable<G>> const& symbol,
    OpExpression<E> const& e, Sp const& solver) {
  constexpr size_t D = expr::grid_dim<G>::value;

  auto interface_term = euler_lagrange_deriv<DynamicVariable<G>, D>(
      *static_cast<E const*>(&e), solver);
  auto bulk_term =
      expr::apply_operators(expr::make_derivative<1, DynamicVariable<G>>(
          *static_cast<E const*>(&e), symbol.index));
  return bulk_term - interface_term;
}

template <typename G, size_t D, typename E, typename Sp>
auto euler_lagrange_apply(SymbolicFunctionalDerivative<GridSymbol<G, D>>,
                          OpExpression<E> const& e, Sp const& solver) {
  auto interface_term = euler_lagrange_deriv<GridSymbol<G, D>, D>(
      *static_cast<E const*>(&e), solver);
  auto bulk_term = expr::apply_operators(
      expr::make_derivative<1, GridSymbol<G, D>>(*static_cast<E const*>(&e)));
  return bulk_term - interface_term;
}

template <typename V, typename V0, typename E, typename T, typename G>
auto apply_operators(
    OpFunctionalDerivative<V, OpDomainIntegral<V0, E, T>, G> const& e) {
  auto integral = expr::get_enclosed_expression(e);
  auto expr = expr::get_enclosed_expression(integral);
  return expr::coeff(e) * expr::coeff(integral) *
         euler_lagrange_apply(e.solver, expr, e.implicit_solver);
}

namespace compile_trait {
template <size_t O, typename factor_t>
struct factor_trait {
  template <typename E>
  auto operator()(E&& e);
};
struct separate_operator_trait {
  template <typename E>
  auto operator()(E&& e);
};

}  // namespace compile_trait

namespace {

template <typename Dd, typename V, typename E, typename Sp, typename E0>
auto apply_operators_deriv(OpDerivative<Dd, V, E, Sp> const& d,
                           OpOperator<E0> const& e) {
  auto [dd, ee] = compile_trait::separate_operator_trait{}(d);
  return OpOperatorChain(dd, ee);
}

template <typename G0, typename E>
struct symbolic_derivative_function {
  template <typename Sp>
  auto operator()(symphas::internal::wrap_base, E const& e) {
    return OpVoid{};
  }

  template <typename Sp>
  auto operator()(symphas::internal::wrap_f<func_cos<E>>, E const& e,
                  Sp const& solver) {
    return -OpFunctionApply<func_sin<E>, OpIdentity, E>(e) *
           apply_operators(expr::make_derivative<1, G0>(e, solver));
  }

  template <typename Sp>
  auto operator()(symphas::internal::wrap_f<func_sin<E>>, E const& e,
                  Sp const& solver) {
    return OpFunctionApply<func_cos<E>, OpIdentity, E>(e) *
           apply_operators(expr::make_derivative<1, G0>(e, solver));
  }

  template <typename Sp>
  auto operator()(symphas::internal::wrap_f<func_tan<E>>, E const& e,
                  Sp const& solver) {
    return OpFunctionApply<func_sec<E>, OpIdentity, E>(e) *
           OpFunctionApply<func_sec<E>, OpIdentity, E>(e) *
           apply_operators(expr::make_derivative<1, G0>(e, solver));
  }

  template <typename Sp>
  auto operator()(symphas::internal::wrap_f<func_csc<E>>, E const& e,
                  Sp const& solver) {
    return -OpFunctionApply<func_cot<E>, OpIdentity, E>(e) *
           OpFunctionApply<func_csc<E>, OpIdentity, E>(e) *
           apply_operators(expr::make_derivative<1, G0>(e, solver));
  }

  template <typename Sp>
  auto operator()(symphas::internal::wrap_f<func_sec<E>>, E const& e,
                  Sp const& solver) {
    return OpFunctionApply<func_tan<E>, OpIdentity, E>(e) *
           OpFunctionApply<func_sec<E>, OpIdentity, E>(e) *
           apply_operators(expr::make_derivative<1, G0>(e, solver));
  }

  template <typename Sp>
  auto operator()(symphas::internal::wrap_f<func_cot<E>>, E const& e,
                  Sp const& solver) {
    return -OpFunctionApply<func_csc<E>, OpIdentity, E>(e) *
           OpFunctionApply<func_csc<E>, OpIdentity, E>(e) *
           apply_operators(expr::make_derivative<1, G0>(e, solver));
  }

  template <typename Sp>
  auto operator()(symphas::internal::wrap_f<func_sqrt<E>>, E const& e,
                  Sp const& solver) {
    return expr::make_fraction<1, 2>() *
           expr::inverse(OpFunctionApply<func_sqrt<E>, OpIdentity, E>(e)) *
           apply_operators(expr::make_derivative<1, G0>(e, solver));
  }

  template <typename Sp>
  auto operator()(symphas::internal::wrap_f<func_log<E>>, E const& e,
                  Sp const& solver) {
    return expr::inverse(e);
  }
};
}  // namespace

template <size_t O, typename V0, auto f, typename V1, typename E, typename G0>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V0, OpFunctionApply<f, V1, E>,
                 SymbolicDerivative<G0>> const& e) {
  auto&& function = expr::get_enclosed_expression(e);
  auto&& expr = expr::get_enclosed_expression(function);

  return apply_operators(expr::make_derivative<O - 1, G0>(
      expr::coeff(e) * expr::coeff(function),
      symbolic_derivative_function<G0, E>{}(symphas::internal::wrap_f<f>{},
                                            expr, e.solver),
      e.solver));
}

template <size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0,
          typename GG>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<expr::symbols::i_<I0, P0>, X0>>,
                 SymbolicDerivative<GG>> const& e) {
  return OpVoid{};
}

template <size_t O, typename V, typename V1, int I0, int P0, typename G1,
          typename... Gs, expr::exp_key_t X0, expr::exp_key_t X1,
          expr::exp_key_t... Xs, typename GG>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V1, Term<expr::symbols::i_<I0, P0>, X0>, Term<G1, X1>,
                         Term<Gs, Xs>...>,
                 SymbolicDerivative<GG>> const& e) {
  auto&& terms = expr::get_enclosed_expression(e);
  auto coeff = expr::coeff(terms);
  auto a = OpTerms(OpIdentity{}, expr::get<1>(terms));
  auto b = expr::terms_after_n<1>(
      terms);  // OpTerms(OpIdentity{}, *static_cast<OpTerms<Term<G1, X1>,
               // Term<Gs, Xs>...> const*>(&terms));

  return (coeff * expr::coeff(e) * a) *
         expr::apply_operators(
             expr::make_derivative<O, GG>(expr::make_term(b), e.solver));
}

template <size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0,
          size_t Z, typename GG, size_t N, typename>
auto apply_operators(
    OpDerivative<
        std::index_sequence<O>, V,
        OpTerms<V0,
                Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>,
        SymbolicDerivative<Variable<Z, GG>>> const& e) {
  using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
  auto f = compile_trait::factor_trait<O, factor_t>{}(
      expr::get_enclosed_expression(e));

  SymbolicCase c(expr::symbols::i_<I0, P0>{} = Variable<Z, GG>{}, f.second,
                 OpVoid{});

  return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
}

template <size_t O, typename V, typename V0, int I0, int P0, size_t D,
          expr::exp_key_t X0, size_t Z, typename GG, size_t N, typename>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<
                                                 expr::symbols::i_<I0, P0>>,
                                             D>,
                                  X0>>,
                 SymbolicDerivative<Variable<Z, GG>>> const& e) {
  using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
  auto f = compile_trait::factor_trait<O, factor_t>{}(
      expr::get_enclosed_expression(e));

  SymbolicCase c(expr::symbols::i_<I0, P0>{} = Variable<Z, GG>{}, f.second,
                 OpVoid{});

  return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
}

template <size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0,
          typename GG, size_t N, typename>
auto apply_operators(
    OpDerivative<
        std::index_sequence<O>, V,
        OpTerms<V0,
                Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>,
        SymbolicDerivative<DynamicVariable<GG>>> const& e) {
  using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
  auto f = compile_trait::factor_trait<O, factor_t>{}(
      expr::get_enclosed_expression(e));

  SymbolicCase c(expr::symbols::i_<I0, P0>{} = DynamicVariable<GG>{}, f.second,
                 OpVoid{});

  return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
}

template <size_t O, typename V, typename V0, int I0, int P0, size_t D,
          expr::exp_key_t X0, typename GG, size_t N, typename>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<
                                                 expr::symbols::i_<I0, P0>>,
                                             D>,
                                  X0>>,
                 SymbolicDerivative<DynamicVariable<GG>>> const& e) {
  using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
  auto f = compile_trait::factor_trait<O, factor_t>{}(
      expr::get_enclosed_expression(e));

  SymbolicCase c(expr::symbols::i_<I0, P0>{} = DynamicVariable<GG>{}, f.second,
                 OpVoid{});

  return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
}

template <size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0,
          typename GG, size_t N, typename>
auto apply_operators(
    OpDerivative<
        std::index_sequence<O>, V,
        OpTerms<V0,
                Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>,
        SymbolicDerivative<GG>> const& e) {
  using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
  auto f = compile_trait::factor_trait<O, factor_t>{}(
      expr::get_enclosed_expression(e));

  SymbolicCase c(expr::symbols::i_<I0, P0>{} = GG{}, f.second, OpVoid{});

  return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
}

template <size_t O, typename V, typename V0, int I0, int P0, size_t D,
          expr::exp_key_t X0, typename GG, size_t N, typename>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<
                                                 expr::symbols::i_<I0, P0>>,
                                             D>,
                                  X0>>,
                 SymbolicDerivative<GG>> const& e) {
  using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
  auto f = compile_trait::factor_trait<O, factor_t>{}(
      expr::get_enclosed_expression(e));

  SymbolicCase c(expr::symbols::i_<I0, P0>{} = GG{}, f.second, OpVoid{});

  return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
}

template <size_t O, typename V, typename V0, typename G0, expr::exp_key_t X0,
          typename GG, size_t N, typename>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>>,
                 SymbolicDerivative<GG>> const& e) {
  auto f =
      compile_trait::factor_trait<O, GG>{}(expr::get_enclosed_expression(e));
  return expr::factorial<N, N - O>() * expr::coeff(e) * f.second;
}

template <size_t O, typename V, typename V0, typename G0, typename G1,
          typename... Gs, expr::exp_key_t X0, expr::exp_key_t X1,
          expr::exp_key_t... Xs, typename GG, typename>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>,
                 SymbolicDerivative<GG>> const& e) {
  auto&& terms = expr::get_enclosed_expression(e);
  auto coeff = expr::coeff(terms);
  auto a = OpTerms(OpIdentity{}, expr::get<1>(terms));
  auto b = expr::terms_after_n<1>(terms);

  auto lhs = apply_operators(expr::make_derivative<O, GG>(a, e.solver)) *
             expr::make_term(b);
  auto rhs = a * apply_operators(expr::make_derivative<O, GG>(
                     expr::make_term(b), e.solver));
  return (coeff * expr::coeff(e)) * (lhs + rhs);
}

template <size_t O, typename V, typename V0, typename G0, expr::exp_key_t X0,
          typename GG>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>>,
                 SymbolicDerivative<OpTerm<OpIdentity, GG>>> const& e) {
  return apply_operators(
      expr::coeff(e) *
      expr::make_derivative<O, GG>(expr::get_enclosed_expression(e)));
}

template <size_t O, typename V, typename V0, typename G0, typename G1,
          typename... Gs, expr::exp_key_t X0, expr::exp_key_t X1,
          expr::exp_key_t... Xs, typename GG>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>,
                 SymbolicDerivative<OpTerm<OpIdentity, GG>>> const& e) {
  return apply_operators(
      expr::coeff(e) *
      expr::make_derivative<O, GG>(expr::get_enclosed_expression(e)));
}

template <typename V, expr::exp_key_t X0, typename V0, typename E0, typename GG>
auto apply_operators(OpDerivative<std::index_sequence<1>, V, OpPow<X0, V0, E0>,
                                  SymbolicDerivative<GG>> const& e) {
  auto power = expr::get_enclosed_expression(e);
  auto p = apply_operators(expr::get_enclosed_expression(power));

  constexpr size_t N0 = expr::_Xk_t<X0>::N;
  constexpr size_t D0 = expr::_Xk_t<X0>::D;
  constexpr bool sign = expr::_Xk_t<X0>::sign;

  constexpr int N1 = (sign) ? N0 + D0 : (N0 < D0) ? D0 - N0 : N0 - D0;
  constexpr bool _sign = (sign) ? sign : (N0 < D0) ? true : false;

  constexpr expr::exp_key_t _X = expr::Xk<N1, D0, _sign>;

  auto result =
      expr::coeff(e) * expr::coeff(power) * expr::make_fraction<N0, D0>() *
      expr::dot(expr::make_pow<_X>(p),
                apply_operators(expr::make_derivative<1, GG>(p, e.solver)));

  if constexpr (sign) {
    return -result;
  } else {
    return result;
  }
}

template <size_t O, typename V, expr::exp_key_t X0, typename V0, typename E0,
          typename GG>
auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpPow<X0, V0, E0>,
                                  SymbolicDerivative<GG>> const& e) {
  return expr::coeff(e) *
         apply_operators(expr::make_derivative<O - 1, GG>(
             apply_operators(expr::make_derivative<1, GG>(
                 expr::get_enclosed_expression(e), e.solver))));
}

template <expr::exp_key_t X, typename V, typename E>
auto apply_operators(OpPow<X, V, E> const& e) {
  auto p = apply_operators(expr::get_enclosed_expression(e));
  return expr::coeff(e) * expr::make_pow<X>(p);
  // if constexpr (expr::_Xk<X> > 1)
  //{
  //	constexpr size_t N0 = expr::_Xk_t<X>::N;
  //	constexpr size_t D0 = expr::_Xk_t<X>::D;
  //	constexpr bool sign = expr::_Xk_t<X>::sign;
  //
  //	return expr::coeff(e) * apply_operators(expr::make_pow<expr::Xk<N0 - D0,
  // D0, sign>>(p)) * p;
  // }
  // else
  //{
  //	return expr::coeff(e) * expr::make_pow<X>(p);
  // }
}

}  // namespace expr

template <typename T, typename E>
auto operator<=(OpTerms<OpIdentity, T> const term, OpOperator<E> const& e) {
  return expr::make_lhs(term) =
             expr::apply_operators(*static_cast<E const*>(&e));
}

namespace expr {

//! Apply the div of a vector expression.
/*!
 * Calculate the divergence of a vector.
 */
template <typename Sp, typename E,
          typename = std::enable_if_t<(expr::eval_type<E>::rank > 0), int>>
auto divergence_of(OpExpression<E> const& e, solver_op_type<Sp> solver) {
  auto op = expr::make_operator_derivative<1>(solver);
  auto result = ::expr::dot(op, *static_cast<E const*>(&e));
  return result;
}
template <typename Sp, typename E,
          typename = std::enable_if_t<(expr::eval_type<E>::rank > 0), int>>
auto divergence_of(OpOperator<E> const& e, solver_op_type<Sp> solver) {
  auto op = expr::make_operator_derivative<1>(solver);
  auto result = ::expr::dot(op, *static_cast<E const*>(&e));
  return apply_operators(result);
}

//! Apply the curl of a vector.
/*!
 * Calculate the curl of a vector.
 */
template <typename Sp, typename E,
          typename = std::enable_if_t<(expr::eval_type<E>::rank > 0), int>>
auto curl_of(OpExpression<E> const& e, solver_op_type<Sp> solver) {
  if constexpr (expr::eval_type<E>::rank == 3) {
    auto op = expr::make_operator_derivative<1>(solver);
    auto result = cross(op, *static_cast<E const*>(&e));
    return result;
  }
  auto x = expr::make_row_vector<0, 2>() * (*static_cast<E const*>(&e));
  auto y = expr::make_row_vector<1, 2>() * (*static_cast<E const*>(&e));
  auto opx = expr::make_operator_directional_derivative<Axis::X, 1>(solver);
  auto opy = expr::make_operator_directional_derivative<Axis::Y, 1>(solver);
  // 2D curl (z-component): dm_y/dx - dm_x/dy
  return opx(y) - opy(x);
}

template <typename Sp, typename E,
          typename = std::enable_if_t<(expr::eval_type<E>::rank > 0), int>>
auto curl_of(OpOperator<E> const& e, solver_op_type<Sp> solver) {
  if constexpr (expr::eval_type<E>::rank == 3) {
    auto op = expr::make_operator_derivative<1>(solver);
    auto result = cross(op, *static_cast<E const*>(&e));
    return result;
  }
  auto x = expr::make_row_vector<0, 2>() * (*static_cast<E const*>(&e));
  auto y = expr::make_row_vector<1, 2>() * (*static_cast<E const*>(&e));
  auto opx = expr::make_operator_directional_derivative<Axis::X, 1>(solver);
  auto opy = expr::make_operator_directional_derivative<Axis::Y, 1>(solver);
  // 2D curl (z-component): dm_y/dx - dm_x/dy
  return opx(y) - opy(x);
}
}  // namespace expr

namespace expr::transform {

template <typename V>
auto sift_term(V value) {
  return value;
}

template <typename V, typename G0, exp_key_t X0, typename... G1s,
          exp_key_t... X1s>
auto sift_term(V value, Term<G0, X0> const& term0,
               Term<G1s, X1s> const&... rest) {
  return OpTerms(value, term0, rest...);
}

}  // namespace expr::transform
