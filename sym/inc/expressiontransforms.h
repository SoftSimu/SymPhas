
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

#include "expressionprune.h"

// ******************************************************************************************

namespace symphas::internal {
//! Constructs a grid of the prescribed dimension.
/*!
 * A grid is constructed using the given dimension, and the expression
 * evaluation type is used as the grid underlying type.
 *
 * The expression is evaluated into the grid.
 *
 * \tparam D The dimension of the grid.
 */
template <size_t D>
struct construct_grid_of_dimension {
  //! Create a grid and populate it with the evaluated expression.
  /*!
   * The grid is created using the dimensions of the underlying
   * expression data, and then populated by evaluating the expression
   * and storing the result in the grid.
   *
   * \param e The expression with which to populate the grid.
   */
  template <typename E, typename T = typename expr::eval_type<E>::type>
  auto operator()(OpExpression<E> const& e) {
    expr::prune::update(*const_cast<E*>(static_cast<E const*>(&e)));
    Grid<T, D> result(expr::data_dimensions(*static_cast<E const*>(&e)));
    expr::result(*static_cast<E const*>(&e), result.values, result.len);
    return result;
  }
};

//! Constructs an array for an expression.
/*!
 * A ::Block is constructed, and the expression
 * evaluation type is used as the grid underlying type. This specialization
 * is used when the dimension of the expression is not a positive nonzero
 * value.
 */
template <>
struct construct_grid_of_dimension<0> {
  //! Create a grid and populate it with the evaluated expression.
  /*!
   * The grid is created using the dimensions of the underlying
   * expression data, and then populated by evaluating the expression
   * and storing the result in the grid.
   *
   * \param e The expression with which to populate the grid.
   */
  template <typename E, typename T = typename expr::eval_type<E>::type>
  auto operator()(OpExpression<E> const& e) {
    expr::prune::update(*const_cast<E*>(static_cast<E const*>(&e)));
    Block<T> result(expr::data_length(*static_cast<E const*>(&e)));
    expr::result(*static_cast<E const*>(&e), result.values, result.len);
    return result;
  }
};

}  // namespace symphas::internal

namespace expr::transform {

//! Create a grid with values initialized to the result of the expression.
/*!
 * Each value is initialized to the final value of the expression. If the
 * dimension is not a positive non zero value, then ::Block is returned
 * instead.
 *
 * \param e The expression that is evaluated.
 *
 * \tparam E The expression type which is evaluated into the grid.
 * \tparam D The dimension of the grid to create.
 */
template <typename E, size_t D = expr::grid_dim<E>::value>
auto to_grid(OpExpression<E> const& e) {
  return symphas::internal::construct_grid_of_dimension<D>{}(
      *static_cast<E const*>(&e));
}

//! Evaluating a literal as a grid returns just the value of the literal.
/*!
 * Evaluating a literal as a grid returns just the value of the literal.
 * See expr::transform::to_grid(OpExpressioN<E>&)
 */
template <typename T>
auto to_grid(OpLiteral<T> const e) {
  return e.eval();
}

//! Evaluating a literal as a grid returns just the value of the literal.
/*!
 * Evaluating a literal as a grid returns just the value of the literal.
 * See expr::transform::to_grid(OpExpressioN<E>&)
 */
template <typename T, typename G>
auto to_grid(OpTerm<T, G> const e) {
  return expr::BaseData<G>::get(expr::get<1>(e).data());
}

//! Evaluating an identity or fraction as a grid returns just the value of the
//! literal.
/*!
 * Evaluating an identity or fraction as a grid returns just the value of the
 * literal. See expr::transform::to_grid(OpExpression<E>&)
 */
template <typename coeff_t,
          typename = std::enable_if_t<
              (expr::is_identity<coeff_t> || expr::is_fraction<coeff_t>), int>>
auto to_grid(coeff_t) {
  return coeff_t{}.eval();
}

namespace {

//! A list of expressions is evaluated into a list of grids.
template <typename... Es, size_t... Is>
decltype(auto) to_grid(std::tuple<Es...>& e, std::index_sequence<Is...>) {
  return std::make_tuple(expr::transform::to_grid(std::get<Is>(e))...);
}
}  // namespace

//! A list of expressions is evaluated into a list of grids.
/*!
 * Converts the expressions into grids, the order of the returned grid list
 * is consistent with the list of expressions. Constructs the list by
 * applying expr::to_grid(OpExpression<E>&) to each element of the list.
 *
 * \param es The list of expressions to evaluate into grids.
 */
template <typename... Es>
decltype(auto) to_grid(std::tuple<Es...>& es) {
  return expr::transform::to_grid(es,
                                  std::make_index_sequence<sizeof...(Es)>{});
}

//! Remove the terms for which the base_data has a match in the given list.
template <typename G0, expr::exp_key_t X0, typename... Gs,
          expr::exp_key_t... Xs, typename... Is>
auto filter_from_term(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& e,
                      symphas::lib::types_list<Is...>) {
  if constexpr ((std::is_same<expr::base_data_t<G0>, Is>::value || ... ||
                 false)) {
    return filter_from_term(expr::terms_after_first(e),
                            symphas::lib::types_list<Is...>{});
  } else {
    return expr::make_term(expr::get<1>(e)) *
           filter_from_term(expr::terms_after_first(e),
                            symphas::lib::types_list<Is...>{});
  }
}

//! Remove the terms for which the base_data has a match in the given list.
template <typename V, typename G0, expr::exp_key_t X0, typename... Gs,
          expr::exp_key_t... Xs, typename... Is>
auto filter_from_term(OpTerms<V, Term<G0, X0>, Term<Gs, Xs>...> const& e,
                      symphas::lib::types_list<Is...>) {
  return expr::coeff(e) * filter_from_term(expr::terms_after_first(e),
                                           symphas::lib::types_list<Is...>{});
}

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename A1, typename A2>
auto to_ft(OpOperatorCombination<A1, A2> const& e, double const* h,
           const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename A1, typename A2>
auto to_ft(OpOperatorChain<A1, A2> const& e, double const* h,
           const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename T0, typename T1, typename V, typename E>
auto to_ft(OpMap<MapGridFourier<T0, T1, D>, V, E> const& e, double const*,
           const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, size_t O, typename V, typename Sp>
auto to_ft(OpOperatorDerivative<O, V, Sp> const& e, double const* h,
           const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, Axis ax, size_t O, typename V, typename Sp>
auto to_ft(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e,
           double const* h, const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename V, typename Sp, size_t... Os>
auto to_ft(OpOperatorMixedDerivative<V, Sp, Os...> const& e, double const* h,
           const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename V, typename E1, typename E2>
auto to_ft(OpConvolution<V, E1, E2> const& e, double const* h,
           const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename V, typename E>
auto to_ft(OpConvolution<V, GaussianSmoothing<D>, E> const& e, double const* h,
           const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename... Es>
auto to_ft(OpAdd<Es...> const& e, double const* h, const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename E1, typename E2>
auto to_ft(OpBinaryMul<E1, E2> const& e, double const* h, const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename E>
auto to_ft(OpBinaryDiv<OpIdentity, E> const& e, double const* h,
           const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename E1, typename E2>
auto to_ft(OpBinaryDiv<E1, E2> const& e, double const* h, const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename A1, typename A2, typename E>
auto to_ft(OpChain<A1, A2, E> const& e, double const* h, const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename A1, typename A2, typename E>
auto to_ft(OpCombination<A1, A2, E> const& e, double const* h,
           const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename Dd, typename V, typename E, typename Sp>
auto to_ft(OpDerivative<Dd, V, E, Sp> const& e, double const* h,
           const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D>
auto to_ft(GaussianSmoothing<D> const& e, double const*, const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename V, typename sub_t, typename E, typename... Ts>
auto to_ft(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e,
           double const* h, const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename V, typename sub_t, typename E, typename... Ts>
auto to_ft(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e,
           double const* h, const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename V, typename E, typename... Ts, int... I0s,
          int... P0s, typename A, typename... As, typename B, typename C>
auto to_ft(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>,
                 symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
                 SymbolicFunction<A, As...>, B, C> const& series,
           double const* h, const len_type* dims);

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename V, size_t O>
auto to_ft(OpTerm<V, k_grid_type<O, D>> const& e, double const* h,
           const len_type*) {
  // return expr::make_operator_derivative<O>(Solver<void>{});
}

template <size_t D, typename V, Axis ax, size_t O>
auto to_ft(OpTerm<V, k_grid_axis_type<ax, O, D>> const& e, double const* h,
           const len_type*) {
  // return symphas::internal::nth_symbolic_derivative_function<ax, O,
  // Sp>::template get(Solver<void>{});
}

template <size_t D, typename V, Axis ax, size_t O>
auto to_ft(OpTerm<V, k_grid_component_type<ax, O, D>> const& e, double const* h,
           const len_type*) {
  // return
  // symphas::internal::make_directional_operator_derivative<O>(Solver<void>{});
}

//! Converts the expression to Fourier space.
/*!
 * Convert the given expression to the Fourier space equivalent.
 *
 * \param e The given expression.
 * \param h The spatial discretization of real space.
 * \param dims The dimensions of real space.
 *
 * \tparam D The real space dimension.
 */
template <size_t D, typename V, typename T,
          template <typename, size_t> typename grid_type, typename E,
          size_t... Ns, typename... Ts>
auto to_ft(OpSymbolicEval<V, NoiseData<NoiseType::WHITE, T, D, grid_type>,
                          SymbolicFunction<E, Variable<Ns, Ts>...>> const& e,
           double const* h, const len_type* dims);

template <size_t D, typename E>
auto to_ft(OpExpression<E> const& e, double const* h, const len_type*) {
  return expr::make_fourier_map(*static_cast<E const*>(&e));
}

template <size_t D, typename E>
auto to_ft(OpOperator<E> const& e, double const* h, const len_type*) {
  return expr::make_fourier_map(*static_cast<E const*>(&e));
}

template <size_t D>
auto to_ft(OpVoid const, double const* h, const len_type*) {
  return OpVoid{};
}

template <size_t D>
auto to_ft(OpIdentity const, double const* h, const len_type*) {
  return OpIdentity{};
}

template <size_t D>
auto to_ft(OpNegIdentity const, double const* h, const len_type* dims) {
  return -to_ft<D>(OpIdentity{}, h, dims);
}

template <size_t D, typename T>
auto to_ft(OpLiteral<T> const e, double const* h, const len_type* dims) {
  return e * to_ft<D>(OpIdentity{}, h, dims);
}

template <size_t D, size_t N, size_t D0>
auto to_ft(OpFractionLiteral<N, D>, double const* h, const len_type* dims) {
  return OpFractionLiteral<N, D>{} * to_ft<D>(OpIdentity{}, h, dims);
}

template <size_t D, size_t N, size_t D0>
auto to_ft(OpNegFractionLiteral<N, D>, double const* h, const len_type* dims) {
  return OpNegFractionLiteral<N, D>{} * to_ft<D>(OpIdentity{}, h, dims);
}

template <size_t D, typename T, size_t... Ns>
auto to_ft(OpTensor<T, Ns...> const& tensor, double const* h,
           const len_type* dims) {
  return tensor * to_ft<D>(OpIdentity{}, h, dims);
}

template <size_t D, typename A1, typename A2>
auto to_ft(OpOperatorCombination<A1, A2> const& e, double const* h,
           const len_type* dims) {
  return to_ft<D>(e.f, h, dims) + to_ft<D>(e.g, h, dims);
}

template <size_t D, typename A1, typename A2>
auto to_ft(OpOperatorChain<A1, A2> const& e, double const* h,
           const len_type* dims) {
  return to_ft<D>(e.f, h, dims) * to_ft<D>(e.g, h, dims);
}

template <size_t D, typename T0, typename T1, typename V, typename E>
auto to_ft(OpMap<MapGridInverseFourier<T0, T1, D>, V, E> const& e,
           double const*, const len_type* dims) {
  return expr::coeff(e) * expr::get_enclosed_expression(e);
}

template <size_t D, size_t O, typename V, typename Sp>
auto to_ft(OpOperatorDerivative<O, V, Sp> const& e, double const* h,
           const len_type* dims) {
  return expr::coeff(e) * expr::make_term(k_grid_type<O, D>(dims, h));
}

template <size_t D, Axis ax, size_t O, typename V, typename Sp>
auto to_ft(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e,
           double const* h, const len_type* dims) {
  return expr::coeff(e) * expr::make_term(k_grid_axis_type<ax, O, D>(dims, h));
}

namespace {

template <size_t D, typename V, typename Sp, size_t... Os, Axis... axs>
auto to_ft_mixed(OpOperatorMixedDerivative<V, Sp, Os...> const& e,
                 double const* h, const len_type* dims,
                 symphas::internal::axis_list<axs...>) {
  return expr::coeff(e) *
         (expr::make_term(k_grid_axis_type<axs, Os, D>(dims, h)) * ...);
}
}  // namespace

template <size_t D, typename V, typename Sp, size_t... Os>
auto to_ft(OpOperatorMixedDerivative<V, Sp, Os...> const& e, double const* h,
           const len_type* dims) {
  return to_ft_mixed(e, h, dims, symphas::lib::make_axis_list<sizeof...(Os)>());
}

template <size_t D, typename V, typename E1, typename E2>
auto to_ft(OpConvolution<V, E1, E2> const& e, double const* h,
           const len_type* dims) {
  return expr::coeff(e) * to_ft<D>(e.a, h, dims) * to_ft<D>(e.b, h, dims);
}

template <size_t D>
auto to_ft(GaussianSmoothing<D> const& e, double const* h,
           const len_type* dims) {
  return GaussianSmoothing<D>(e.data.dims, e.h, e.sigma, !e.fourier_space);
}

template <size_t D, typename V, typename sub_t, typename E, typename... Ts>
auto to_ft(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e,
           double const* h, const len_type* dims) {
  auto ft = to_ft<D>(e.f.e, h, dims);
  auto f = (expr::function_of(Ts{}...) = ft);
  f.set_data_tuple(e.data);
  return symphas::internal::make_symbolic_eval(expr::coeff(e), e.data, f);
}

template <size_t D, typename V, typename E, typename... Ts, int... I0s,
          int... P0s, typename A, typename... As, typename B, typename C>
auto to_ft(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>,
                 symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
                 SymbolicFunction<A, As...>, B, C> const& series,
           double const* h, const len_type* dims) {
  auto ft = to_ft<D>(series.data.e, h, dims);
  return expr::coeff(series) * expr::recreate_series(ft, series.data);
}

template <size_t D, typename V, typename E>
auto to_ft(OpConvolution<V, GaussianSmoothing<D>, E> const& e, double const* h,
           const len_type* dims) {
  return expr::coeff(e) * to_ft<D>(e.smoother, h, dims) *
         to_ft<D>(expr::get_enclosed_expression(e), h, dims);
}

namespace {
template <size_t D, typename... Es, size_t... Is>
auto to_ft_adds(OpAdd<Es...> const& e, double const* h, const len_type* dims,
                std::index_sequence<Is...>) {
  return (to_ft<D>(expr::get<Is>(e), h, dims) + ...);
}
}  // namespace

template <size_t D, typename... Es>
auto to_ft(OpAdd<Es...> const& e, double const* h, const len_type* dims) {
  return to_ft_adds<D>(e, h, dims, std::make_index_sequence<sizeof...(Es)>{});
}

template <size_t D, typename E1, typename E2>
auto to_ft(OpBinaryMul<E1, E2> const& e, double const* h,
           const len_type* dims) {
  return expr::make_convolution(to_ft<D>(e.a, h, dims), to_ft<D>(e.b, h, dims));
}

template <size_t D, typename E>
auto to_ft(OpBinaryDiv<OpIdentity, E> const& e, double const* h,
           const len_type* dims) {
  return expr::make_fourier_map(e);
}

template <size_t D, typename E1, typename E2>
auto to_ft(OpBinaryDiv<E1, E2> const& e, double const* h,
           const len_type* dims) {
  return expr::make_convolution(to_ft<D>(e.a, h, dims),
                                to_ft<D>(expr::inverse(e.b), h, dims));
}

/* conversion of operators into fourier space
 */

template <size_t D, typename A1, typename A2, typename E>
auto to_ft(OpChain<A1, A2, E> const& e, double const* h, const len_type* dims) {
  return to_ft<D>(e.combination, h, dims) *
         to_ft<D>(expr::get_enclosed_expression(e), h, dims);
}

template <size_t D, typename A1, typename A2, typename E>
auto to_ft(OpCombination<A1, A2, E> const& e, double const* h,
           const len_type* dims) {
  return to_ft<D>(e.combination, h, dims) *
         to_ft<D>(expr::get_enclosed_expression(e), h, dims);
}

namespace {
template <size_t D, typename Sp, size_t... Os, Axis... axs>
auto to_ft_mixed(Sp const& solver,
                 typename Sp::template mixed_derivative<Os...>,
                 symphas::internal::axis_list<axs...>, double const* h,
                 const len_type* dims) {
  return (OpIdentity{} * ... *
          to_ft<D>(expr::make_operator_directional_derivative<axs, Os>(solver),
                   h, dims));
}

template <size_t D, typename Sp, size_t... Os>
auto to_ft_mixed(Sp const& solver,
                 typename Sp::template mixed_derivative<Os...>, double const* h,
                 const len_type* dims) {
  return to_ft_mixed<D>(solver, typename Sp::template mixed_derivative<Os...>{},
                        symphas::lib::make_axis_list<sizeof...(Os)>(), h, dims);
}
}  // namespace

template <size_t D, typename Dd, typename V, typename E, typename Sp>
auto to_ft(OpDerivative<Dd, V, E, Sp> const& e, double const* h,
           const len_type* dims) {
  constexpr Axis axis = OpDerivative<Dd, V, E, Sp>::axis;
  constexpr size_t order = OpDerivative<Dd, V, E, Sp>::order;

  if constexpr (Dd::is_directional) {
    if constexpr (Dd::is_mixed) {
      return expr::coeff(e) * to_ft_mixed<D, Sp>(e.solver, Dd{}, h, dims) *
             to_ft<D>(expr::get_enclosed_expression(e), h, dims);
    } else {
      return expr::coeff(e) *
             expr::make_term(k_grid_axis_type<axis, order, D>(dims, h)) *
             to_ft<D>(expr::get_enclosed_expression(e), h, dims);
    }
  } else {
    if constexpr (order % 2 == 0) {
      return expr::coeff(e) * expr::make_term(k_grid_type<order, D>(dims, h)) *
             to_ft<D>(expr::get_enclosed_expression(e), h, dims);
    } else {
      return expr::coeff(e) *
             expr::make_term(k_grid_component_type<axis, order, D>(dims, h)) *
             to_ft<D>(expr::get_enclosed_expression(e), h, dims);
    }
  }
}

}  // namespace expr::transform

namespace expr {

template <typename E>
auto fix_index(OpExpression<E>& e, DynamicIndexSet set) {
  expr::transform::swap_grid<DynamicIndexSet>(*static_cast<E const*>(&e), set);
}

template <typename E>
auto fix_index(OpOperator<E>& e, DynamicIndexSet set) {
  expr::transform::swap_grid<DynamicIndexSet>(*static_cast<E const*>(&e), set);
}

template <typename E>
auto fix_coeff(OpExpression<E>& e) {
  expr::transform::swap_grid<OpCoeffSwap<DynamicIndex>>(
      *static_cast<E const*>(&e));
}

template <typename E>
auto fix_coeff(OpOperator<E>& e) {
  expr::transform::swap_grid<OpCoeffSwap<DynamicIndex>>(
      *static_cast<E const*>(&e));
}
}  // namespace expr

namespace symphas::internal {
template <typename E1, typename E2>
auto remove_factors(E1 const& e1, E2 const& e2, symphas::lib::types_list<>) {
  return std::make_pair(e1, e2);
}

template <size_t N01, typename G01, size_t... N1s, typename... G1s, typename E1,
          typename E2>
auto remove_factors(
    E1 const& e1, E2 const& e2,
    symphas::lib::types_list<std::pair<std::index_sequence<N01>, G01>,
                             std::pair<std::index_sequence<N1s>, G1s>...>) {
  auto f = expr::split::factor<N01, G01>(e1);
  auto g = expr::split::factor<N01, G01>(e2);
  return remove_factors(
      f.second, g.second,
      symphas::lib::types_list<std::pair<std::index_sequence<N1s>, G1s>...>{});
}

}  // namespace symphas::internal

namespace expr {

//! Applies division using the provided types as factors.
/*!
 * In the general case, no factoring is performed, but when some factors are
 * passed (and present), then the division will factor those provided elements.
 * The given types are provided in a tuple and have to be factors. If no
 * tuple is provided, a regular division will be returned.
 */
template <typename G>
struct divide_with_factors;

//! Specialization of expr::divide_with_factors with no factors given.
template <>
struct divide_with_factors<symphas::lib::types_list<>> {
  template <typename E1, typename E2>
  auto operator()(OpExpression<E1> const& a, OpExpression<E2> const& b) {
    return symphas::internal::terminate_div(
        *static_cast<E1 const*>(&a),
        expr::inverse(*static_cast<E2 const*>(&b)));
  }

  template <typename E1, typename E2>
  auto operator()(OpExpression<E1> const& a, OpOperator<E2> const& b) {
    return expr::make_div(*static_cast<E1 const*>(&a),
                          *static_cast<E2 const*>(&b));
  }

  template <typename E1, typename E2>
  auto operator()(OpOperator<E1> const& a, OpExpression<E2> const& b) {
    return expr::make_div(*static_cast<E1 const*>(&a),
                          *static_cast<E2 const*>(&b));
  }

  template <typename E1, typename E2>
  auto operator()(OpOperator<E1> const& a, OpOperator<E2> const& b) {
    return expr::make_div(*static_cast<E1 const*>(&a),
                          *static_cast<E2 const*>(&b));
  }
};

//! Specialization of expr::divide_with_factors a list of factors given.
template <size_t N0, typename G0, size_t... Ns, typename... Gs>
struct divide_with_factors<
    symphas::lib::types_list<std::pair<std::index_sequence<N0>, G0>,
                             std::pair<std::index_sequence<Ns>, Gs>...>> {
  template <typename E1, typename E2>
  auto operator()(OpExpression<E1> const& a, OpExpression<E2> const& b) {
    using factors_t =
        symphas::lib::types_list<std::pair<std::index_sequence<N0>, G0>,
                                 std::pair<std::index_sequence<Ns>, Gs>...>;
    auto [numerator, denominator] = symphas::internal::remove_factors(
        *static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b), factors_t{});

    return symphas::internal::terminate_div(numerator,
                                            expr::inverse(denominator));
  }

  template <typename E1, typename E2>
  auto operator()(OpExpression<E1> const& a, OpOperator<E2> const& b) {
    return expr::make_div(*static_cast<E1 const*>(&a),
                          *static_cast<E2 const*>(&b));
  }

  template <typename E1, typename E2>
  auto operator()(OpOperator<E1> const& a, OpExpression<E2> const& b) {
    return expr::make_div(*static_cast<E1 const*>(&a),
                          *static_cast<E2 const*>(&b));
  }

  template <typename E1, typename E2>
  auto operator()(OpOperator<E1> const& a, OpOperator<E2> const& b) {
    return expr::make_div(*static_cast<E1 const*>(&a),
                          *static_cast<E2 const*>(&b));
  }
};

}  // namespace expr

/*
 *
 *
 * Division rules
 *
 ******************************************************************************/

//! The division operator is overloaded to apply factoring in general.
template <typename E1, typename E2>
auto operator/(OpExpression<E1> const& a, OpExpression<E2> const& b) {
  return expr::divide_with_factors<
      typename expr::factor_list_all<E1, E2>::type>{}(
      *static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
}

/*
 *
 *
 * Additional overloads to avoid ambiguous overloads
 *
 ******************************************************************************/

namespace expr {
namespace {

template <typename... Ts>
struct filter_for_case;

template <typename... Seqs, typename... Ts>
struct filter_for_case<symphas::lib::types_list<Seqs...>,
                       symphas::lib::types_list<Ts...>> {
  using type = symphas::lib::types_list<Ts...>;
};

template <typename... Seqs, typename... Ts, typename G0, typename... Gs>
struct filter_for_case<symphas::lib::types_list<Seqs...>,
                       symphas::lib::types_list<Ts...>, G0, Gs...> {
  using type =
      typename filter_for_case<symphas::lib::types_list<Seqs...>,
                               symphas::lib::types_list<Ts...>, Gs...>::type;
};

template <typename... Seqs, typename... Ts, int I0, int P0, size_t Z0,
          typename G0, typename L, typename R, typename... Gs>
struct filter_for_case<
    std::integer_sequence<int, -1>, symphas::lib::types_list<Seqs...>,
    symphas::lib::types_list<Ts...>,
    SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>,
    Gs...> {
  using type = typename filter_for_case<
      symphas::lib::types_list<Seqs..., expr::symbols::i_<I0, P0>>,
      symphas::lib::types_list<
          Ts..., symphas::lib::types_list<SymbolicTernaryCase<
                     expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>>>,
      Gs...>::type;
};

template <int I, typename... Seqs, typename... Ts, int I0, int P0, size_t Z0,
          typename G0, typename L, typename R, typename... Gs>
struct filter_for_case<
    std::integer_sequence<int, I>, symphas::lib::types_list<Seqs...>,
    symphas::lib::types_list<Ts...>,
    SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>,
    Gs...> {
  static const size_t N = size_t(I);
  using type_N = symphas::lib::direct_type_at_index<N, Ts...>;

  using type_combined = symphas::lib::expand_types_list<
      type_N,
      SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>>;
  using seq_combined = symphas::lib::expand_types_list<
      symphas::lib::types_before_index<N, Seqs...>, expr::symbols::i_<I0, P0>,
      symphas::lib::types_after_at_index<N + 1, Seqs...>>;

  using type = typename filter_for_case<
      seq_combined,
      symphas::lib::expand_types_list<
          symphas::lib::types_before_index<N, Ts...>,
          symphas::lib::types_list<type_combined>,
          symphas::lib::types_after_at_index<N + 1, Ts...>>,
      Gs...>::type;
};

template <typename... Seqs, typename... Ts, int I0, int P0, size_t Z0,
          typename G0, typename L, typename R, typename... Gs>
struct filter_for_case<
    symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
    SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>,
    Gs...> {
  static const int ind =
      symphas::lib::index_of_type<expr::symbols::i_<I0, P0>, Seqs...>;

  using type = typename filter_for_case<
      std::integer_sequence<int, ind>, symphas::lib::types_list<Seqs...>,
      symphas::lib::types_list<Ts...>,
      SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>,
      Gs...>::type;
};

// template<typename... Seqs, typename... Ts, int I0, int P0, typename G0,
// typename L, typename R, typename... Gs> struct filter_for_case<
//	std::integer_sequence<int, -1>,
//	symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
//	SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, L,
// R>, Gs...>
//{
//	using type = typename filter_for_case<
//		symphas::lib::types_list<Seqs..., expr::symbols::i_<I0, P0>>,
//		symphas::lib::types_list<Ts...,
// symphas::lib::types_list<SymbolicTernaryCase<expr::symbols::i_<I0, P0>,
// DynamicVariable<G0>, L, R>>>, Gs...>::type;
// };

// template<int I, typename... Seqs, typename... Ts, int I0, int P0, typename
// G0, typename L, typename R, typename... Gs> struct filter_for_case<
//	std::integer_sequence<int, I>,
//	symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
//	SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, L,
// R>, Gs...>
//{
//	static const size_t N = size_t(I);
//	using type_N = symphas::lib::direct_type_at_index<N, Ts...>;

//	using type_combined = symphas::lib::expand_types_list<type_N,
// SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, L, R>>;
//	using seq_combined =
// symphas::lib::expand_types_list<symphas::lib::types_before_index<N, Seqs...>,
// expr::symbols::i_<I0, P0>, symphas::lib::types_after_at_index<N + 1,
// Seqs...>>;

//	using type = typename filter_for_case<
//		seq_combined,
//		symphas::lib::expand_types_list<symphas::lib::types_before_index<N,
// Ts...>, symphas::lib::types_list<type_combined>,
// symphas::lib::types_after_at_index<N + 1, Ts...>>, 		Gs...>::type;
//};

// template<typename... Seqs, typename... Ts, int I0, int P0, typename G0,
// typename L, typename R, typename... Gs> struct filter_for_case<
//	symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
//	SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, L,
// R>, Gs...>
//{
//	static const int ind = symphas::lib::index_of_type<expr::symbols::i_<I0,
// P0>, Seqs...>;

//	using type = typename filter_for_case<
//		std::integer_sequence<int, ind>,
//		symphas::lib::types_list<Seqs...>,
//		symphas::lib::types_list<Ts...>,
// SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, L, R>,
// Gs...>::type;
//};

template <typename... Seqs, typename... Ts, int I0, int P0, typename G0,
          typename L, typename R, typename... Gs>
struct filter_for_case<
    std::integer_sequence<int, -1>, symphas::lib::types_list<Seqs...>,
    symphas::lib::types_list<Ts...>,
    SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>, Gs...> {
  using type = typename filter_for_case<
      symphas::lib::types_list<Seqs..., expr::symbols::i_<I0, P0>>,
      symphas::lib::types_list<
          Ts..., symphas::lib::types_list<
                     SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>>>,
      Gs...>::type;
};

template <int I, typename... Seqs, typename... Ts, int I0, int P0, typename G0,
          typename L, typename R, typename... Gs>
struct filter_for_case<
    std::integer_sequence<int, I>, symphas::lib::types_list<Seqs...>,
    symphas::lib::types_list<Ts...>,
    SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>, Gs...> {
  static const size_t N = size_t(I);
  using type_N = symphas::lib::direct_type_at_index<N, Ts...>;

  using type_combined = symphas::lib::expand_types_list<
      type_N, SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>>;
  using seq_combined = symphas::lib::expand_types_list<
      symphas::lib::types_before_index<N, Seqs...>, expr::symbols::i_<I0, P0>,
      symphas::lib::types_after_at_index<N + 1, Seqs...>>;

  using type = typename filter_for_case<
      seq_combined,
      symphas::lib::expand_types_list<
          symphas::lib::types_before_index<N, Ts...>,
          symphas::lib::types_list<type_combined>,
          symphas::lib::types_after_at_index<N + 1, Ts...>>,
      Gs...>::type;
};

template <typename... Seqs, typename... Ts, int I0, int P0, typename G0,
          typename L, typename R, typename... Gs>
struct filter_for_case<
    symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
    SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>, Gs...> {
  static const int ind =
      symphas::lib::index_of_type<expr::symbols::i_<I0, P0>, Seqs...>;

  using type = typename filter_for_case<
      std::integer_sequence<int, ind>, symphas::lib::types_list<Seqs...>,
      symphas::lib::types_list<Ts...>,
      SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>, Gs...>::type;
};

template <typename... Seqs, typename... Ts, typename... Gs>
struct filter_for_case<symphas::lib::types_list<Seqs...>,
                       symphas::lib::types_list<Ts...>,
                       symphas::lib::types_list<Gs...>> {
  using type =
      typename filter_for_case<symphas::lib::types_list<Ts...>, Gs...>::type;
};

template <typename... Gs>
struct filter_for_case<symphas::lib::types_list<Gs...>> {
  using type =
      typename filter_for_case<symphas::lib::types_list<>,
                               symphas::lib::types_list<>, Gs...>::type;
};

template <size_t N0, typename T>
struct pos_in_substitution;

template <size_t N0, size_t... Ns, typename... Gs>
struct pos_in_substitution<N0, symphas::lib::types_list<Variable<Ns, Gs>...>> {
  static const int value = symphas::lib::index_of_value<size_t, N0, Ns...>;
};

template <size_t N0, size_t... Ns, typename... Gs>
struct pos_in_substitution<
    N0, symphas::lib::types_list<NamedData<Variable<Ns, Gs>>...>> {
  static const int value = symphas::lib::index_of_value<size_t, N0, Ns...>;
};

template <size_t N0, typename... Gs>
struct pos_in_substitution<N0, SymbolicDataArray<std::tuple<Gs...>>> {
  static const int value =
      pos_in_substitution<N0, symphas::lib::types_list<Gs...>>::value;
};

template <size_t N0, int N, int P>
struct pos_in_substitution<
    N0, SymbolicDataArray<expr::symbols::v_id_type<expr::symbols::i_<N, P>>>> {
  static const int value = -1;
};

template <size_t N0, typename T>
struct pos_in_substitution<N0, SymbolicDataArray<T>> {
  static const int value = N0;
};

template <size_t N0, typename T>
struct type_in_substitution;

template <size_t N0, size_t... Ns, typename... Gs>
struct type_in_substitution<
    N0, SymbolicDataArray<std::tuple<Variable<Ns, Gs>...>>> {
  using type =
      typename SymbolicDataArray<std::tuple<Variable<Ns, Gs>...>>::storage_type;
};

template <size_t N0, size_t... Ns, typename... Gs>
struct type_in_substitution<
    N0, SymbolicDataArray<std::tuple<NamedData<Variable<Ns, Gs>>...>>> {
  using type = typename SymbolicDataArray<
      std::tuple<NamedData<Variable<Ns, Gs>>...>>::storage_type;
};

template <size_t N0, int N, int P>
struct type_in_substitution<
    N0, SymbolicDataArray<expr::symbols::v_id_type<expr::symbols::i_<N, P>>>> {
  using type = expr::symbols::Symbol;
};

template <size_t N0, typename T>
struct type_in_substitution<N0, SymbolicDataArray<T>> {
  using type = T;
};

template <size_t Z0, size_t... Ns, typename... Gs>
auto get_arg(
    SymbolicDataArray<std::tuple<Variable<Ns, Gs>...>> const& substitution) {
  return expr::make_term<Z0>(
      *substitution.data[symphas::lib::index_of_value<size_t, Z0, Ns...>].data);
}

template <size_t Z0, size_t... Ns, typename... Gs>
auto get_arg(
    SymbolicDataArray<std::tuple<NamedData<Variable<Ns, Gs>>...>> const&
        substitution) {
  auto el = substitution.data[symphas::lib::index_of_value<size_t, Z0, Ns...>];
#ifdef PRINTABLE_EQUATIONS
  return expr::make_term<Z0>(NamedData(std::ref(*el.data), el.name));
#else
  return expr::make_term<Z0>(NamedData(std::ref(*el.data)));
#endif
}

template <size_t Z0, int N, int P>
auto get_arg(
    SymbolicDataArray<expr::symbols::v_id_type<expr::symbols::i_<N, P>>> const&
        substitution) {
  return expr::symbols::v_<expr::symbols::i_<N, P>>{};
}

template <size_t Z0, typename T>
auto get_arg(SymbolicDataArray<T> const& substitution) {
  // auto el = substitution.data[Z0];
  // return expr::make_term(NamedData(std::ref(*el.data), el.name));
  constexpr size_t D = expr::grid_dim<T>::value;
  return expr::make_term<Z0>(
      GridSymbol<expr::eval_type_t<OpTerm<OpIdentity, T>>, D>{});
}

template <size_t Z0, typename... Ts, size_t... Ns>
auto get_args(Substitution<Ts...> const& substitution,
              std::index_sequence<Ns...>) {
  return std::make_tuple(get_arg<Z0>(std::get<Ns>(substitution))...);
}

template <typename T>
auto get_arg(SymbolicDataArray<T> const& substitution) {
  constexpr size_t D = expr::grid_dim<T>::value;
  return expr::make_term(
      GridSymbol<expr::eval_type_t<OpTerm<OpIdentity, T>>, D>{});
}

template <typename T>
auto get_arg(SymbolicDataArray<T> const& substitution,
             DynamicIndex const& index) {
  return expr::make_term_dynamic(index, substitution.data);
}

template <typename T>
auto get_arg(SymbolicDataArray<NamedData<T*>> const& substitution,
             DynamicIndex const& index) {
  return expr::make_term_dynamic(
      index, NamedData(substitution.data, substitution.name));
}

template <typename... Ts, size_t... Ns>
auto get_args(Substitution<Ts...> const& substitution,
              std::index_sequence<Ns...>) {
  return std::make_tuple(get_arg(std::get<Ns>(substitution))...);
}

template <typename... Ts, size_t... Ns>
auto get_args(Substitution<Ts...> const& substitution,
              DynamicIndex const& index, std::index_sequence<Ns...>) {
  return std::make_tuple(get_arg(std::get<Ns>(substitution), index)...);
}

template <typename T, size_t D>
auto get_arg(GridSymbol<T, D> const& substitution) {
  return substitution;
}

template <typename T, size_t D>
auto get_arg(T const& substitution) {
  return substitution;
}

template <typename E, typename... Ts, int... I0s, int... P0s, typename B,
          typename C, typename E0, int I0, int P0, size_t Z0, typename G0,
          typename... Ls, typename... Rs>
auto apply_operators_sum_case(
    SymbolicSum<E, Substitution<Ts...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B,
                C> const& series,
    OpExpression<E0> const& applied,
    SymbolicDerivative<Variable<Z0, G0>> const& symbol,
    symphas::lib::types_list<SymbolicTernaryCase<
        expr::symbols::i_<I0, P0>, Variable<Z0, G0>, Ls, Rs>...>) {
  constexpr int ind_N = symphas::lib::index_of_value<int, I0, I0s...>;
  using pos_t =
      pos_in_substitution<Z0,
                          symphas::lib::type_at_index<size_t(ind_N), Ts...>>;
  constexpr int arg_N = pos_t::value;

  using v_type = expr::symbols::v_<expr::symbols::i_<I0, P0>>;

  auto arg = get_arg<Z0>(std::get<size_t(ind_N)>(series.substitution));

  auto left_case = expr::transform::swap_grid<v_type>(
      SymbolicCase(expr::symbols::i_<I0, P0>{} = Variable<Z0, G0>{}),
      *static_cast<E0 const*>(&applied), arg);
  auto right_case = expr::transform::swap_grid<v_type>(
      SymbolicCase(expr::symbols::i_<I0, P0>{} != Variable<Z0, G0>{}),
      *static_cast<E0 const*>(&applied), arg);

  auto left =
      expr::transform::swap_grid<SymbolicCaseSwap<>>(left_case, OpVoid{});
  auto right =
      expr::transform::swap_grid<SymbolicCaseSwap<>>(right_case, OpVoid{});

  auto series_left = expr::recreate_series(
      expr::symbols::i_<I0, P0>{} = expr::val<arg_N>, left, series);
  auto series_right =
      OpVoid{};  // expr::recreate_series(expr::symbols::i_<I0, 0>{} !=
                 // expr::val<arg_N>, expr_right, series);
  return series_left + series_right;
}

template <typename E, typename... Ts, int... I0s, int... P0s, typename B,
          typename C, typename E0, int I0, int P0, typename G0, typename... Ls,
          typename... Rs, std::enable_if_t<!expr::is_id_variable<G0>, int> = 0>
auto apply_operators_sum_case(
    SymbolicSum<E, Substitution<Ts...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B,
                C> const& series,
    OpExpression<E0> const& applied, SymbolicDerivative<G0> const& symbol,
    symphas::lib::types_list<
        SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, Ls, Rs>...>) {
  // constexpr int ind_N = symphas::lib::index_of_value<int, I0, I0s...>;

  using v_type = expr::symbols::v_<expr::symbols::i_<I0, P0>>;

  // auto arg = get_arg(std::get<size_t(ind_N)>(series.substitution));
  auto arg = get_arg(G0{});

  auto left_case = expr::transform::swap_grid<v_type>(
      SymbolicCase(expr::symbols::i_<I0, P0>{} = G0{}),
      *static_cast<E0 const*>(&applied), arg);
  auto right_case = expr::transform::swap_grid<v_type>(
      SymbolicCase(expr::symbols::i_<I0, P0>{} != G0{}),
      *static_cast<E0 const*>(&applied), arg);

  auto left =
      expr::transform::swap_grid<SymbolicCaseSwap<>>(left_case, OpVoid{});
  auto right =
      expr::transform::swap_grid<SymbolicCaseSwap<>>(right_case, OpVoid{});

  auto series_left = expr::recreate_series(
      expr::symbols::i_<I0, P0>{} = expr::symbols::placeholder_N{}, left,
      series);
  auto series_right =
      OpVoid{};  // expr::recreate_series(expr::symbols::i_<I0, 0>{} !=
                 // expr::val<arg_N>, expr_right, series);
  return series_left + series_right;
}

template <typename E, typename... Ts, int... I0s, int... P0s, typename B,
          typename C, typename E0, int I0, int P0, typename G0, typename... Ls,
          typename... Rs>
auto apply_operators_sum_case(
    SymbolicSum<E, Substitution<Ts...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B,
                C> const& series,
    OpExpression<E0> const& applied,
    SymbolicDerivative<DynamicVariable<G0>> const& symbol,
    symphas::lib::types_list<SymbolicTernaryCase<
        expr::symbols::i_<I0, P0>, DynamicVariable<G0>, Ls, Rs>...>) {
  constexpr int ind_N = symphas::lib::index_of_value<int, I0, I0s...>;

  using v_type = expr::symbols::v_<expr::symbols::i_<I0, P0>>;

  auto arg =
      get_arg(std::get<size_t(ind_N)>(series.substitution), symbol.index);

  auto left_case = expr::transform::swap_grid<v_type>(
      SymbolicCase(expr::symbols::i_<I0, P0>{} = DynamicVariable<G0>{}),
      *static_cast<E0 const*>(&applied), arg);
  // auto right_case = expr::transform::swap_grid<v_type>
  //	(SymbolicCase(expr::symbols::i_<I0, P0>{} != DynamicVariable<G0>{}),
  //*static_cast<E0 const*>(&applied), arg);

  auto left =
      expr::transform::swap_grid<SymbolicCaseSwap<>>(left_case, OpVoid{});
  // auto right = expr::transform::swap_grid<SymbolicCaseSwap<>>(right_case,
  // OpVoid{});

  auto series_left = expr::recreate_series(
      expr::symbols::i_<I0, P0>{} = expr::symbols::placeholder_N{}, left,
      series);
  auto series_right =
      OpVoid{};  // expr::recreate_series(expr::symbols::i_<I0, 0>{} !=
                 // expr::val<arg_N>, expr_right, series);
  return expr::transform::swap_grid<expr::symbols::placeholder_N_symbol>(
      series_left + series_right, symbol.index);
}

template <typename E, typename... Ts, int... I0s, int... P0s, typename B,
          typename C, typename E0, int I0, int P0, typename G0, typename... Ls,
          typename... Rs>
auto apply_operators_sum_case(
    SymbolicSum<E, Substitution<Ts...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B,
                C> const& series,
    OpExpression<E0> const& applied,
    SymbolicFunctionalDerivative<G0> const& symbol,
    symphas::lib::types_list<
        SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, Ls, Rs>...>
        cases) {
  return apply_operators_sum_case(series, *static_cast<E0 const*>(&applied),
                                  SymbolicDerivative<G0>{}, cases);
}

template <typename E, typename... Ts, int... I0s, int... P0s, typename B,
          typename C, typename E0, int I0, int P0, typename G0, typename... Ls,
          typename... Rs>
auto apply_operators_sum_case(
    SymbolicSum<E, Substitution<Ts...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B,
                C> const& series,
    OpExpression<E0> const& applied,
    SymbolicFunctionalDerivative<DynamicVariable<G0>> const& symbol,
    symphas::lib::types_list<SymbolicTernaryCase<
        expr::symbols::i_<I0, P0>, DynamicVariable<G0>, Ls, Rs>...>
        cases) {
  return apply_operators_sum_case(
      series, *static_cast<E0 const*>(&applied),
      SymbolicDerivative<DynamicVariable<G0>>(symbol.index), cases);
}

template <typename E, typename... Ts, int... I0s, int... P0s, typename B,
          typename C, typename E0, int I0, int P0, size_t Z0, typename G0,
          typename... Ls, typename... Rs>
auto apply_operators_sum_case(
    SymbolicSum<E, Substitution<Ts...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B,
                C> const& series,
    OpExpression<E0> const& applied,
    SymbolicFunctionalDerivative<Variable<Z0, G0>> const& symbol,
    symphas::lib::types_list<SymbolicTernaryCase<expr::symbols::i_<I0, P0>,
                                                 Variable<Z0, G0>, Ls, Rs>...>
        cases) {
  return apply_operators_sum_case(series, *static_cast<E0 const*>(&applied),
                                  SymbolicDerivative<Variable<Z0, G0>>{},
                                  cases);
}

template <typename E, typename T, typename Seq, typename B, typename C,
          typename E0, typename G>
auto apply_operators_sum(SymbolicSum<E, T, Seq, B, C> const& series,
                         OpExpression<E0> const& applied,
                         SymbolicDerivative<G> const& symbol,
                         symphas::lib::types_list<>) {
  return *static_cast<E0 const*>(&applied);
}

// template<typename V, typename E, typename... Ts,
//	int... I0s, int... P0s, typename A, typename B, typename C, typename...
// Rest> auto apply_operators_sum_recurse( 	OpSum<V, E,
//		Substitution<SymbolicDataArray<Ts>...>,
//		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
// C> const& sum, 	symphas::lib::types_list<Rest...>);

template <typename E, typename T, typename Seq, typename B, typename C,
          typename E0, typename G, typename... Cs, typename... Rest>
auto apply_operators_sum(
    SymbolicSum<E, T, Seq, B, C> const& series, OpExpression<E0> const& applied,
    SymbolicDerivative<G> const& symbol,
    symphas::lib::types_list<symphas::lib::types_list<Cs...>, Rest...>) {
  return (apply_operators_sum_case(series, *static_cast<E0 const*>(&applied),
                                   symbol, symphas::lib::types_list<Cs...>{}) +
          ... +
          apply_operators_sum_case(series, *static_cast<E0 const*>(&applied),
                                   symbol, Rest{}));
}

template <typename E, typename T, typename Seq, typename B, typename C,
          typename E0, typename G>
auto apply_operators_sum(SymbolicSum<E, T, Seq, B, C> const& series,
                         OpExpression<E0> const& applied,
                         SymbolicDerivative<G> const& symbol) {
  using ops_t = typename filter_for_case<expr::op_types_t<E0>>::type;
  return apply_operators_sum(series, *static_cast<E0 const*>(&applied), symbol,
                             ops_t{});
}

// template<typename V, typename E, typename... Ts,
//	int... I0s, int... P0s, typename A, typename B, typename C, typename...
// Rest> auto apply_operators_sum_recurse( 	OpSum<V, E,
//		Substitution<SymbolicDataArray<Ts>...>,
//		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
// C> const& sum, 	symphas::lib::types_list<Rest...>)
//{
//	auto expr = expr::coeff(sum) * sum.data.substitute_placeholders(sum.f);
//	return apply_operators_sum(sum.data, expr,
// symphas::lib::types_list<Rest...>{});
// }

template <typename... Vs, size_t O, typename V, typename V0, typename E,
          typename... Ts, int... I0s, int... P0s, typename A, typename B,
          typename C, typename... GGs, typename... Gs, size_t... Ns>
auto apply_operators_sum_multiple(
    std::index_sequence<O>, V const& coeff,
    OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
          symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
          C> const& sum,
    std::tuple<SymbolicDerivative<GGs>...> const& symbols,
    std::tuple<Gs...> const& args, std::index_sequence<Ns...>) {
  auto expr = sum.data.substitute_placeholders(sum.f);
  return coeff * expr::coeff(sum) *
         (expr::transform::swap_grid<Vs>(
              apply_operators_sum(sum.data,
                                  apply_operators(expr::make_derivative<O, GGs>(
                                      expr, std::get<Ns>(symbols))),
                                  std::get<Ns>(symbols)),
              std::get<Ns>(args)) +
          ...);
}

template <size_t O, typename V, typename V0, typename E, typename... Ts,
          int... I0s, int... P0s, typename A, typename B, typename... Vs,
          size_t Z, typename GG>
auto apply_operators_sum(
    std::index_sequence<O>, V const& coeff,
    OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
          symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
          symphas::lib::types_list<Vs...>> const& sum,
    SymbolicDerivative<Variable<Z, GG>> const& symbol) {
  auto expr = sum.data.substitute_placeholders(sum.f);
  return coeff * expr::coeff(sum) *
         apply_operators_sum(
             sum.data,
             apply_operators(
                 expr::make_derivative<O, Variable<Z, GG>>(expr, symbol)),
             symbol);
}

template <size_t O, typename V, typename V0, typename E, typename... Ts,
          int... I0s, int... P0s, typename A, typename B, typename... Vs,
          typename GG>
auto apply_operators_sum(
    std::index_sequence<O>, V const& coeff,
    OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
          symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
          symphas::lib::types_list<Vs...>> const& sum,
    SymbolicDerivative<DynamicVariable<GG>> const& symbol) {
  auto expr = sum.data.substitute_placeholders(sum.f);
  return coeff * expr::coeff(sum) *
         apply_operators_sum(
             sum.data,
             apply_operators(
                 expr::make_derivative<O, DynamicVariable<GG>>(expr, symbol)),
             symbol);
}

template <typename T>
struct is_variational_trait {
  static const bool value = false;
};

template <typename T>
struct is_variational_trait<expr::variational_t<T>> {
  static const bool value = true;
};

template <
    size_t O, typename V, typename V0, typename E, typename... Ts, int... I0s,
    int... P0s, typename A, typename B, typename... Vs, typename GG,
    std::enable_if_t<(!expr::is_id_variable<GG> &&
                      !expr::is_functional_derivative<SymbolicDerivative<GG>>),
                     int> = 0>
auto apply_operators_sum(
    std::index_sequence<O>, V const& coeff,
    OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
          symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
          symphas::lib::types_list<Vs...>> const& sum,
    SymbolicDerivative<GG> const& symbol) {
  if constexpr (expr::is_expression<GG>) {
    using vs_seq = vars_in_ops_t<op_types_t<GG>>;
    constexpr size_t L = vs_seq::size();

    // if multiple variables are in the expression with respect to differentiate
    if constexpr (L > 1) {
      return OpVoid{};
    } else {
      constexpr size_t D = expr::grid_dim<GG>::value;
      if constexpr (vs_seq::size() > 0) {
        constexpr size_t Z0 = symphas::lib::seq_index_value<0, vs_seq>::value;

        auto dvs =
            std::make_tuple(SymbolicDerivative{expr::transform::swap_grid<Z0>(
                GG{}, GridSymbol<
                          expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>,
                          D>{})}...);
        auto args = get_args<Z0>(sum.data.substitution,
                                 std::make_index_sequence<sizeof...(Ts)>{});
        return apply_operators_sum_multiple<
            expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>...>(
            std::index_sequence<O>{}, coeff, sum, dvs, args,
            std::make_index_sequence<sizeof...(Ts)>{});
      } else {
        auto dvs =
            std::make_tuple(SymbolicDerivative{expr::transform::swap_grid(
                op_types_t<GG>{}, GG{},
                GridSymbol<
                    expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>,
                    D>{})}...);
        // auto dvs = symphas::lib::types_list<
        //	decltype(expr::transform::swap_grid(op_types_t<GG>{},
        // std::declval<GG>(),
        // GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>,
        // D>{}))... > {};
        auto args = get_args(sum.data.substitution,
                             std::make_index_sequence<sizeof...(Ts)>{});
        return apply_operators_sum_multiple<
            expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>...>(
            std::index_sequence<O>{}, coeff, sum, dvs, args,
            std::make_index_sequence<sizeof...(Ts)>{});
      }
    }
  } else {
    auto expr = sum.data.substitute_placeholders(sum.f);
    return coeff * expr::coeff(sum) *
           apply_operators_sum(
               sum.data,
               apply_operators(expr::make_derivative<O, GG>(expr, symbol)),
               symbol);
  }
}

template <size_t O, typename V, typename V0, typename E, typename... Ts,
          int... I0s, int... P0s, typename A, typename B, typename... Vs,
          size_t Z, typename GG>
auto apply_operators_sum(
    std::index_sequence<O>, V const& coeff,
    OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
          symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
          symphas::lib::types_list<Vs...>> const& sum,
    SymbolicFunctionalDerivative<Variable<Z, GG>> const& symbol) {
  constexpr size_t D = expr::grid_dim<GG>::value;
  auto symbols = std::make_tuple(
      SymbolicFunctionalDerivative<GridSymbol<
          expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>, D>>()...);
  auto args = get_args<Z>(sum.data.substitution,
                          std::make_index_sequence<sizeof...(Ts)>{});

  return apply_operators_sum_multiple<
      expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>...>(
      std::index_sequence<O>{}, coeff, sum, symbols, args,
      std::make_index_sequence<sizeof...(Ts)>{});
}

template <size_t O, typename V, typename V0, typename E, typename... Ts,
          int... I0s, int... P0s, typename A, typename B, typename... Vs,
          typename GG>
auto apply_operators_sum(
    std::index_sequence<O>, V const& coeff,
    OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
          symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
          symphas::lib::types_list<Vs...>> const& sum,
    SymbolicFunctionalDerivative<DynamicVariable<GG>> const& symbol) {
  constexpr size_t D = expr::grid_dim<GG>::value;
  auto symbols = std::make_tuple(
      SymbolicFunctionalDerivative<
          GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>, D>>(
          symbol.index)...);
  auto args = get_args(sum.data.substitution, symbol.index,
                       std::make_index_sequence<sizeof...(Ts)>{});

  auto result = apply_operators_sum_multiple<
      expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>...>(
      std::index_sequence<O>{}, coeff, sum, symbols, args,
      std::make_index_sequence<sizeof...(Ts)>{});
  return expr::transform::swap_grid<
      expr::symbols::placeholder_N_symbol_<0>,
      OpCoeffSwap<expr::symbols::placeholder_N_symbol_<0>>>(
      result, symbol.index, symbol.index);
}

}  // namespace

template <size_t O, typename V, typename V0, typename E, typename... Ts,
          int... I0s, int... P0s, typename A, typename B, typename... Vs,
          typename GG>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
                       symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
                       A, B, symphas::lib::types_list<Vs...>>,
                 SymbolicDerivative<GG>> const& e) {
  return apply_operators_sum(std::index_sequence<O>{}, expr::coeff(e),
                             expr::get_enclosed_expression(e), e.solver);
}

template <size_t O, typename V, size_t Z, int I0, int P0, typename V0,
          typename E, typename... Ts, int... I0s, int... P0s, typename A,
          typename B, typename... Vs>
auto apply_operators(
    OpDerivative<std::index_sequence<O>, V,
                 OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
                       symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
                       A, B, symphas::lib::types_list<Vs...>>,
                 SymbolicDerivative<expr::symbols::v_id_type<
                     expr::symbols::i_<I0, P0>>>> const& e) {
  auto sum = expr::get_enclosed_expression(e);
  return expr::coeff(e) *
         apply_operators_sum(
             sum.data,
             apply_operators(
                 expr::make_derivative<
                     O, expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>>(
                     sum.data.e)),
             e.solver);
}

template <typename V, typename E, typename... Ts, int... I0s, int... P0s,
          typename A, typename B, typename... Vs>
auto apply_operators(
    OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>,
          symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B,
          symphas::lib::types_list<Vs...>> const& sum) {
  auto expr = apply_operators(sum.data.substitute_placeholders(sum.f));
  return expr::coeff(sum) * expr::recreate_series(expr, sum.data);
}

}  // namespace expr

namespace symphas::internal {

template <typename E1, typename E2, size_t... Rs, size_t R = sizeof...(Rs)>
auto dot_tensor_components(OpExpression<E1> const& a, OpExpression<E2> const& b,
                           std::index_sequence<Rs...>) {
  return (((expr::make_row_vector<Rs, R>() * (*static_cast<E1 const*>(&a))) *
           (expr::make_row_vector<Rs, R>() * (*static_cast<E2 const*>(&b)))) +
          ...);
}

template <typename E1, typename E2, size_t... Rs, size_t R = sizeof...(Rs)>
auto mul_tensor_components_rc(OpExpression<E1> const& a,
                              OpExpression<E2> const& b,
                              std::index_sequence<Rs...>) {
  return ((((*static_cast<E1 const*>(&a)) * expr::make_column_vector<Rs, R>()) *
           (expr::make_row_vector<Rs, R>() * (*static_cast<E2 const*>(&b)))) +
          ...);
}

template <size_t R, size_t P, typename E1, typename E2, size_t... Rs,
          size_t... Ps>
auto _mul_tensor_components_cr(OpExpression<E1> const& a,
                               OpExpression<E2> const& b,
                               std::index_sequence<Rs...>,
                               std::index_sequence<Ps...>) {
  return (
      (expr::make_tensor<Rs, Ps, R, P>() *
       ((expr::make_row_vector<Rs, R>() * (*static_cast<E1 const*>(&a))) *
        ((*static_cast<E2 const*>(&b)) * expr::make_column_vector<Ps, P>()))) +
      ...);
}

template <typename E1, typename E2, size_t... Rs, size_t R = sizeof...(Rs)>
auto mul_tensor_components_cr(OpExpression<E1> const& a,
                              OpExpression<E2> const& b,
                              std::index_sequence<Rs...>) {
  return _mul_tensor_components_cr<R, R>(
      *static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b),
      symphas::lib::seq_join_t<
          std::index_sequence<>,
          type_ignore_index<Rs, std::index_sequence<Rs...>>...>{},
      symphas::lib::seq_join_t<
          std::index_sequence<>,
          symphas::lib::seq_repeating_value_t<R, size_t, Rs>...>{});
}

template <size_t R0, size_t R, size_t P0, size_t P, typename E1, typename E2,
          size_t... Qs, size_t Q = sizeof...(Qs)>
auto get_dot_product_at(OpExpression<E1> const& a, OpExpression<E2> const& b,
                        std::index_sequence<Qs...>) {
  return (
      (((expr::make_row_vector<R0, R>() * (*static_cast<E1 const*>(&a))) *
        expr::make_column_vector<Qs, Q>()) *
       (expr::make_row_vector<Qs, Q>() *
        ((*static_cast<E2 const*>(&b)) * expr::make_column_vector<P0, P>()))) +
      ...);
}

template <size_t R, size_t P, typename E1, typename E2, size_t... Rs,
          size_t... Qs, size_t... Ps, size_t Q = sizeof...(Qs)>
auto _mul_tensor_components(OpExpression<E1> const& a,
                            OpExpression<E2> const& b,
                            std::index_sequence<Rs...>,
                            std::index_sequence<Qs...>,
                            std::index_sequence<Ps...>) {
  return ((expr::make_tensor<Rs, Ps, R, P>() *
           get_dot_product_at<Rs, R, Ps, P>(*static_cast<E1 const*>(&a),
                                            *static_cast<E2 const*>(&b),
                                            std::index_sequence<Qs...>{})) +
          ...);
}

template <typename E1, typename E2, size_t... Rs, size_t... Qs, size_t... Ps,
          size_t R = sizeof...(Rs), size_t Q = sizeof...(Qs),
          size_t P = sizeof...(Ps)>
auto mul_tensor_components(OpExpression<E1> const& a, OpExpression<E2> const& b,
                           std::index_sequence<Rs...>,
                           std::index_sequence<Qs...>,
                           std::index_sequence<Ps...>) {
  return _mul_tensor_components<R, P>(
      *static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b),
      symphas::lib::seq_join_t<
          std::index_sequence<>,
          type_ignore_index<Ps, std::index_sequence<Rs...>>...>{},
      std::index_sequence<Qs...>{},
      symphas::lib::seq_join_t<
          std::index_sequence<>,
          symphas::lib::seq_repeating_value_t<R, size_t, Ps>...>{});
}

template <typename E1, typename E2, size_t R1 = expr::eval_type<E1>::rank,
          size_t R2 = expr::eval_type<E2>::rank,
          size_t Q1 = expr::eval_type<E1>::template rank_<1>,
          size_t Q2 = expr::eval_type<E2>::template rank_<1>>
auto dot(OpExpression<E1> const& a, OpExpression<E2> const& b) {
  // in the case when the expressions don't have tensors, but are simply scalar
  if constexpr (R1 == 0 && R2 == 0 && Q1 == 0 && Q2 == 0) {
    return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
  }
  // in the special case when 1D tensors are being multiplied
  else if constexpr (R1 == 1 && Q2 == 1 && Q1 == 1 && R1 == 1) {
    return (symphas::internal::tensor_cancel{} *
            (*static_cast<E1 const*>(&a))) *
           (symphas::internal::tensor_cancel{} * (*static_cast<E2 const*>(&b)));
  }
  // multiply a row type by a column type
  else if constexpr (R1 == 1 && Q2 == 1 && Q1 == R2 && R2 > 1) {
    return mul_tensor_components_rc(*static_cast<E1 const*>(&a),
                                    *static_cast<E2 const*>(&b),
                                    std::make_index_sequence<R2>{});
  }
  // multiply a column type by a row type
  else if constexpr (R2 == 1 && Q1 == 1 && Q2 == R1 && R1 > 1) {
    return mul_tensor_components_cr(*static_cast<E1 const*>(&a),
                                    *static_cast<E2 const*>(&b),
                                    std::make_index_sequence<R1>{});
  }
  // dot product of two column vector types
  else if constexpr (R1 == R2 && Q1 == 1 && Q2 == 1) {
    return dot_tensor_components(*static_cast<E1 const*>(&a),
                                 *static_cast<E2 const*>(&b),
                                 std::make_index_sequence<R1>{});
  }
  // matrix multiplication
  else if constexpr (Q1 == R2) {
    return mul_tensor_components(
        *static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b),
        std::make_index_sequence<R1>{}, std::make_index_sequence<Q1>{},
        std::make_index_sequence<Q2>{});
  }
  // multiplying incompatible types just gives their multiplication
  else {
    return expr::make_mul(*static_cast<E1 const*>(&a),
                          *static_cast<E2 const*>(&b));
  }
}

template <typename E1, typename E2>
auto dot(OpExpression<E1> const& a, OpOperator<E2> const& b) {
  return expr::make_mul(*static_cast<E1 const*>(&a),
                        *static_cast<E2 const*>(&b));
}

template <typename E1, typename E2>
auto dot(OpOperator<E1> const& a, OpExpression<E2> const& b) {
  return expr::make_mul(*static_cast<E1 const*>(&a),
                        *static_cast<E2 const*>(&b));
}

template <typename E1, typename E2>
auto dot(OpOperator<E1> const& a, OpOperator<E2> const& b) {
  return expr::make_mul(*static_cast<E1 const*>(&a),
                        *static_cast<E2 const*>(&b));
}

template <typename E1, typename E2>
auto dot(OpOperatorChain<OpIdentity, E1> const& a,
         OpOperatorChain<OpIdentity, E2> const& b) {
  return dot(a.g, b.g);
}

template <typename E1>
auto dot(OpExpression<E1> const& a, OpVoid) {
  return OpVoid{};
}

template <typename E2>
auto dot(OpVoid, OpExpression<E2> const& b) {
  return OpVoid{};
}

inline auto dot(OpVoid, OpVoid) { return OpVoid{}; }
}  // namespace symphas::internal

namespace expr {
template <typename E1, typename E2>
auto dot(OpExpression<E1> const& a, OpExpression<E2> const& b) {
  return symphas::internal::dot(*static_cast<E1 const*>(&a),
                                *static_cast<E2 const*>(&b));
}

template <typename E1, typename E2>
auto dot(OpExpression<E1> const& a, OpOperator<E2> const& b) {
  return symphas::internal::dot(*static_cast<E1 const*>(&a),
                                *static_cast<E2 const*>(&b));
}

template <typename E1, typename E2>
auto dot(OpOperator<E1> const& a, OpExpression<E2> const& b) {
  return symphas::internal::dot(*static_cast<E1 const*>(&a),
                                *static_cast<E2 const*>(&b));
}

template <typename E1, typename E2>
auto dot(OpOperator<E1> const& a, OpOperator<E2> const& b) {
  return symphas::internal::dot(*static_cast<E1 const*>(&a),
                                *static_cast<E2 const*>(&b));
}

template <typename T, size_t N, size_t D>
auto transpose(OpTensor<T, N, D> const& tensor) {
  return expr::make_tensor<0, D, 1, N>(T(tensor));
}

template <typename T, size_t N0, size_t N1, size_t D0, size_t D1>
auto transpose(OpTensor<T, N0, N1, D0, D1> const& tensor) {
  return expr::make_tensor<N1, N0, D1, D0>(T(tensor));
}

template <size_t R1, size_t R2, size_t P, typename E, size_t... Qs>
auto compute_transpose_row(OpExpression<E> const& e,
                           std::index_sequence<Qs...>) {
  return (
      (expr::make_tensor<Qs, P, R2, R1>() *
       (expr::make_row_vector<P, R1>() *
        (*static_cast<E const*>(&e) * expr::make_column_vector<Qs, R2>()))) +
      ...);
}

template <size_t R1, size_t R2, typename E, size_t... Ps, size_t... Qs>
auto compute_transpose(OpExpression<E> const& e, std::index_sequence<Ps...>,
                       std::index_sequence<Qs...>) {
  return (compute_transpose_row<R1, R2, Ps>(*static_cast<E const*>(&e),
                                            std::index_sequence<Qs...>{}) +
          ...);
}

template <typename E, size_t R1 = expr::eval_type<E>::rank,
          size_t R2 = expr::eval_type<E>::template rank_<1>>
auto transpose(OpExpression<E> const& e) {
  if constexpr ((R1 == 0 && R2 == 0) || (R1 == 1 && R2 == 1)) {
    return *static_cast<E const*>(&e);
  } else {
    return compute_transpose<R1, R2>(*static_cast<E const*>(&e),
                                     std::make_index_sequence<R1>{},
                                     std::make_index_sequence<R2>{});
  }
}
}  // namespace expr

namespace symphas::internal {

template <size_t Z, typename G, typename... Gs>
auto filter_variables(std::tuple<Variable<Z, G>, Gs...> const& data_list);
template <size_t Z0, typename G0>
auto sort_variables(std::tuple<Variable<Z0, G0>> const& data_list);

inline auto filter_variables(std::tuple<> const& data_list) {
  return std::make_tuple();
}

template <typename G, typename... Gs>
auto filter_variables(std::tuple<G, Gs...> const& data_list) {
  return filter_variables(symphas::lib::get_tuple_ge<1>(data_list));
}

template <size_t Z, typename G, typename... Gs>
auto filter_variables(std::tuple<Variable<Z, G>, Gs...> const& data_list) {
  return std::tuple_cat(
      std::make_tuple(std::get<0>(data_list)),
      filter_variables(symphas::lib::get_tuple_ge<1>(data_list)));
}

inline auto sort_variables(std::tuple<> const& data_list) {
  return std::make_tuple();
}

template <size_t Z0, typename G0>
auto sort_variables(std::tuple<Variable<Z0, G0>> const& data_list) {
  return data_list;
}

template <size_t Z0, typename G0>
auto sort_variables(Variable<Z0, G0> const& data0,
                    std::tuple<> const& data_list) {
  return std::make_tuple(data0);
}

template <size_t Z0, size_t Z1, typename G0, typename G1, typename... Gs>
auto sort_variables(Variable<Z0, G0> const& data0,
                    std::tuple<Variable<Z1, G1>, Gs...> const& data_list) {
  if constexpr (Z0 > Z1) {
    return std::tuple_cat(std::make_tuple(std::get<0>(data_list), data0),
                          symphas::lib::get_tuple_ge<1>(data_list));
  } else {
    return std::tuple_cat(std::make_tuple(data0, std::get<0>(data_list)),
                          symphas::lib::get_tuple_ge<1>(data_list));
  }
}

template <size_t Z0, size_t Z1, typename G0, typename G1, typename... Gs>
auto sort_variables(
    std::tuple<Variable<Z0, G0>, Variable<Z1, G1>, Gs...> const& data_list) {
  if constexpr (Z0 > Z1) {
    return sort_variables(
        std::get<1>(data_list),
        sort_variables(
            std::get<0>(data_list),
            sort_variables(symphas::lib::get_tuple_ge<2>(data_list))));
  } else {
    return sort_variables(
        std::get<0>(data_list),
        sort_variables(
            std::get<1>(data_list),
            sort_variables(symphas::lib::get_tuple_ge<2>(data_list))));
  }
}

template <typename... Gs>
auto index_variables(std::tuple<Gs...> const& data_list,
                     std::index_sequence<>) {
  return std::make_tuple();
}

template <size_t I0, size_t... Is>
auto index_variables(std::tuple<> const& data_list,
                     std::index_sequence<I0, Is...>) {
  return std::make_tuple(OpVoid{}, type_ignore_index<Is, OpVoid>{}...);
}

template <size_t Z, Axis ax, typename G>
auto remove_component(Variable<Z, VectorComponent<ax, G>> const& data0) {
  return Variable<Z, G>(*static_cast<G const*>(&data0));
}

template <size_t Z, typename G>
auto remove_component(Variable<Z, G> const& data0) {
  return data0;
}

template <size_t Z0, typename G0, typename... Gs, size_t I0, size_t... Is>
auto index_variables(std::tuple<Variable<Z0, G0>, Gs...> const& data_list,
                     std::index_sequence<I0, Is...>) {
  if constexpr (Z0 == I0) {
    return std::tuple_cat(
        std::make_tuple(remove_component(std::get<0>(data_list))),
        index_variables(symphas::lib::get_tuple_ge<1>(data_list),
                        std::index_sequence<Is...>{}));
  } else if constexpr (Z0 < I0) {
    return std::tuple_cat(
        std::make_tuple(OpVoid{}),
        index_variables(symphas::lib::get_tuple_ge<1>(data_list),
                        std::index_sequence<Is...>{}));
  } else {
    return std::tuple_cat(
        std::make_tuple(OpVoid{}),
        index_variables(data_list, std::index_sequence<Is...>{}));
  }
}

template <size_t Z0, typename G0, size_t... Zs, typename... Gs>
auto index_variables(
    std::tuple<Variable<Z0, G0>, Variable<Zs, Gs>...> const& data_list) {
  constexpr size_t Zm =
      symphas::lib::seq_index_value<sizeof...(Zs),
                                    std::index_sequence<Z0, Zs...>>::value;
  return index_variables(data_list, std::make_index_sequence<Zm + 1>{});
}

inline auto index_variables(std::tuple<> const& data_list) {
  return std::make_tuple();
}
}  // namespace symphas::internal

namespace expr {
//! Returns a list of all the variables from the expression, sorted.
/*!
 * A list of all the variables from the OpTerm elements in the expression is
 * aggregrated, duplicates are removed and then the list is sorted by variable
 * ID, i.e. `Z` is the ID in `Variable<Z, G>`. The variables are placed in a
 * tuple according to their ID, and OpVoid is placed where no ID exists.
 *
 * \param e The expression in which to search variables.
 */
template <typename E>
auto get_indexed_variable_list(OpExpression<E> const& e) {
  return symphas::internal::index_variables(
      symphas::internal::sort_variables(symphas::internal::filter_variables(
          expr::data_list(*static_cast<E const*>(&e)))));
}
}  // namespace expr

/*
 *
 * Multiplication between two expressions that are vector type
 * is assumed to be the dot product.
 *
 ******************************************************************************/

template <typename E1, typename E2,
          typename std::enable_if_t<(expr::eval_type<E1>::rank > 0 &&
                                     expr::eval_type<E2>::rank > 0),
                                    int> = 0>
auto operator*(OpExpression<E1> const& a, OpExpression<E2> const& b) {
  return expr::dot(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}

template <
    typename V0, typename... G0s, expr::exp_key_t... X0s, typename V1,
    typename W0, typename... G1s, expr::exp_key_t... X1s, typename W1,
    typename Dd, typename E, typename Sp,
    typename = std::enable_if_t<
        (expr::is_combinable<G0s...> &&
         symphas::lib::filter_seq_t<
             std::integer_sequence<expr::exp_key_t, X0s...>,
             std::integer_sequence<expr::exp_key_t, X1s...>>::size() == 0 &&
         symphas::lib::filter_seq_t<
             std::integer_sequence<expr::exp_key_t, X1s...>,
             std::integer_sequence<expr::exp_key_t, X0s...>>::size() == 0 &&
         ((symphas::lib::index_of_type<G0s, G1s...> >= 0) && ...) &&
         ((symphas::lib::index_of_type<G1s, G0s...> >= 0) && ...)),
        int>>
auto operator+(OpBinaryMul<OpTerms<V0, Term<G0s, X0s>...>,
                           OpDerivative<Dd, V1, E, Sp>> const& a,
               OpBinaryMul<OpTerms<W0, Term<G1s, X1s>...>,
                           OpDerivative<Dd, W1, E, Sp>> const& b) {
  return (expr::coeff(b.a) * expr::coeff(b.b) +
          expr::coeff(a.a) * expr::coeff(a.b)) *
         (OpTerms<OpIdentity, Term<G0s, X0s>...>(OpIdentity{},
                                                 expr::terms_after_first(a.a)) *
          expr::make_derivative<Dd>(expr::get_enclosed_expression(a.b),
                                    a.b.solver));
}

template <
    typename V0, typename... G0s, expr::exp_key_t... X0s, typename V1,
    typename W0, typename... G1s, expr::exp_key_t... X1s, typename W1,
    typename Dd, typename E, typename Sp,
    typename = std::enable_if_t<
        (expr::is_combinable<G0s...> &&
         symphas::lib::filter_seq_t<
             std::integer_sequence<expr::exp_key_t, X0s...>,
             std::integer_sequence<expr::exp_key_t, X1s...>>::size() == 0 &&
         symphas::lib::filter_seq_t<
             std::integer_sequence<expr::exp_key_t, X1s...>,
             std::integer_sequence<expr::exp_key_t, X0s...>>::size() == 0 &&
         ((symphas::lib::index_of_type<G0s, G1s...> >= 0) && ...) &&
         ((symphas::lib::index_of_type<G1s, G0s...> >= 0) && ...)),
        int>>
auto operator-(OpBinaryMul<OpTerms<V0, Term<G0s, X0s>...>,
                           OpDerivative<Dd, V1, E, Sp>> const& a,
               OpBinaryMul<OpTerms<W0, Term<G1s, X1s>...>,
                           OpDerivative<Dd, W1, E, Sp>> const& b) {
  return (expr::coeff(b.a) * expr::coeff(b.b) -
          expr::coeff(a.a) * expr::coeff(a.b)) *
         (OpTerms<OpIdentity, Term<G0s, X0s>...>(OpIdentity{},
                                                 expr::terms_after_first(a.a)) *
          expr::make_derivative<Dd>(expr::get_enclosed_expression(a.b),
                                    a.b.solver));
}

namespace expr {

template <typename A, typename B>
struct numeric_range_dist {
  using type = sub_result_t<B, A>;
};

template <size_t N, size_t D, size_t N0, size_t D0, bool first_neg,
          bool second_neg>
struct numeric_range_count_impl;

template <size_t N, size_t D, size_t N0, size_t D0>
struct numeric_range_count_impl<N, D, N0, D0, false, false> {
  using type = sub_result_t<OpFractionLiteral<N0, D0>, OpFractionLiteral<N, D>>;
};

template <size_t N, size_t D, size_t N0, size_t D0>
struct numeric_range_count_impl<N, D, N0, D0, false, true> {
  using type =
      sub_result_t<OpNegFractionLiteral<N0, D0>, OpFractionLiteral<N, D>>;
};

template <size_t N, size_t D, size_t N0, size_t D0>
struct numeric_range_count_impl<N, D, N0, D0, true, false> {
  using type =
      sub_result_t<OpFractionLiteral<N0, D0>, OpNegFractionLiteral<N, D>>;
};

template <size_t N, size_t D, size_t N0, size_t D0>
struct numeric_range_count_impl<N, D, N0, D0, true, true> {
  using type =
      sub_result_t<OpNegFractionLiteral<N0, D0>, OpNegFractionLiteral<N, D>>;
};

template <typename A>
struct numeric_range_value {
  constexpr static int value = A{}.eval();
};

template <int N, size_t D, int N0, size_t D0>
struct numeric_range_count {
  using type = typename numeric_range_count_impl < (N < 0) ? size_t(-N)
                                                           : size_t(N),
        D, (N0 < 0) ? size_t(-N0) : size_t(N0), D0, N < 0, N0<0>::type;
  constexpr static int value = numeric_range_value<type>::value;
};

//! Represents a numeric range of values.
/*!
 * Represents a numeric range of values.
 *
 * \tparam N The numerator of the starting value.
 * \tparam D The denominator of the starting value.
 * \tparam N0 The numerator of the ending value.
 * \tparam D0 The denominator of the ending value.
 * \tparam C The number (and direction) of steps.
 */
template <int N, size_t D, int N0, size_t D0, int C>
struct numeric_range_state {};

//! Represents a constant interval.
/*!
 * Represents a constant interval.
 *
 * \tparam N The numerator of the starting value.
 * \tparam D The denominator of the starting value.
 * \tparam N0 The numerator of the ending value.
 * \tparam D0 The denominator of the ending value.
 */
template <int N, size_t D, int N0, size_t D0>
struct numeric_interval_state {
  template <typename A, typename B>
  constexpr numeric_interval_state(A, B) {}

  constexpr auto operator&(OpIdentity) const {
    return numeric_range_state<N, D, N0, D0, 1>{};
  }

  constexpr auto operator&(OpNegIdentity) const {
    return numeric_range_state<N, D, N0, D0, -1>{};
  }

  template <size_t C>
  constexpr auto operator&(OpFractionLiteral<C, 1>) const {
    return numeric_range_state<N, D, N0, D0, int(C)>{};
  }

  template <size_t C>
  constexpr auto operator&(OpNegFractionLiteral<C, 1>) const {
    return numeric_range_state<N, D, N0, D0, -int(C)>{};
  }

  constexpr auto operator|(OpIdentity) const {
    constexpr int C = numeric_range_count<N, D, N0, D0>::value;
    return numeric_range_state<N, D, N0, D0, C>{};
  }

  constexpr auto operator|(OpNegIdentity) const {
    constexpr int C = numeric_range_count<N, D, N0, D0>::value;
    return numeric_range_state<N, D, N0, D0, C>{};
  }

  template <size_t dN, size_t dD>
  constexpr auto operator|(OpFractionLiteral<dN, dD>) const {
    constexpr int C =
        div_result_t<typename numeric_range_count<N, D, N0, D0>::type,
                     OpFractionLiteral<dN, dD>>{}
            .eval();
    return numeric_range_state<N, D, N0, D0, C>{};
  }

  template <size_t dN, size_t dD>
  constexpr auto operator|(OpNegFractionLiteral<dN, dD>) const {
    constexpr int C =
        div_result_t<typename numeric_range_count<N, D, N0, D0>::type,
                     OpFractionLiteral<dN, dD>>{}
            .eval();
    return numeric_range_state<N, D, N0, D0, C>{};
  }
};

numeric_interval_state(OpIdentity, OpIdentity)
    -> numeric_interval_state<1, 1, 1, 1>;
numeric_interval_state(OpIdentity, OpNegIdentity)
    -> numeric_interval_state<1, 1, -1, 1>;
numeric_interval_state(OpIdentity, OpVoid)
    -> numeric_interval_state<1, 1, 0, 1>;
template <size_t N, size_t D>
numeric_interval_state(OpIdentity, OpFractionLiteral<N, D>)
    -> numeric_interval_state<1, 1, int(N), D>;
template <size_t N, size_t D>
numeric_interval_state(OpIdentity, OpNegFractionLiteral<N, D>)
    -> numeric_interval_state<1, 1, -int(N), D>;

numeric_interval_state(OpNegIdentity, OpIdentity)
    -> numeric_interval_state<-1, 1, 1, 1>;
numeric_interval_state(OpNegIdentity, OpNegIdentity)
    -> numeric_interval_state<-1, 1, -1, 1>;
numeric_interval_state(OpNegIdentity, OpVoid)
    -> numeric_interval_state<-1, 1, 0, 1>;
template <size_t N, size_t D>
numeric_interval_state(OpNegIdentity, OpFractionLiteral<N, D>)
    -> numeric_interval_state<-1, 1, int(N), D>;
template <size_t N, size_t D>
numeric_interval_state(OpNegIdentity, OpNegFractionLiteral<N, D>)
    -> numeric_interval_state<-1, 1, -int(N), D>;

numeric_interval_state(OpVoid, OpIdentity)
    -> numeric_interval_state<0, 1, 1, 1>;
numeric_interval_state(OpVoid, OpNegIdentity)
    -> numeric_interval_state<0, 1, -1, 1>;
numeric_interval_state(OpVoid, OpVoid) -> numeric_interval_state<0, 1, 0, 1>;
template <size_t N, size_t D>
numeric_interval_state(OpVoid, OpFractionLiteral<N, D>)
    -> numeric_interval_state<0, 1, int(N), D>;
template <size_t N, size_t D>
numeric_interval_state(OpVoid, OpNegFractionLiteral<N, D>)
    -> numeric_interval_state<0, 1, -int(N), D>;

template <size_t N, size_t D>
numeric_interval_state(OpFractionLiteral<N, D>, OpIdentity)
    -> numeric_interval_state<int(N), D, 1, 1>;
template <size_t N, size_t D>
numeric_interval_state(OpFractionLiteral<N, D>, OpNegIdentity)
    -> numeric_interval_state<int(N), D, -1, 1>;
template <size_t N, size_t D>
numeric_interval_state(OpFractionLiteral<N, D>, OpVoid)
    -> numeric_interval_state<int(N), D, 0, 1>;
template <size_t N, size_t D, size_t N0, size_t D0>
numeric_interval_state(OpFractionLiteral<N, D>, OpFractionLiteral<N0, D0>)
    -> numeric_interval_state<int(N), D, int(N0), D0>;
template <size_t N, size_t D, size_t N0, size_t D0>
numeric_interval_state(OpFractionLiteral<N, D>, OpNegFractionLiteral<N0, D0>)
    -> numeric_interval_state<int(N), D, -int(N0), D0>;

template <size_t N, size_t D>
numeric_interval_state(OpNegFractionLiteral<N, D>, OpIdentity)
    -> numeric_interval_state<-int(N), D, 1, 1>;
template <size_t N, size_t D>
numeric_interval_state(OpNegFractionLiteral<N, D>, OpNegIdentity)
    -> numeric_interval_state<-int(N), D, -1, 1>;
template <size_t N, size_t D>
numeric_interval_state(OpNegFractionLiteral<N, D>, OpVoid)
    -> numeric_interval_state<-int(N), D, 0, 1>;
template <size_t N, size_t D, size_t N0, size_t D0>
numeric_interval_state(OpNegFractionLiteral<N, D>, OpFractionLiteral<N0, D0>)
    -> numeric_interval_state<-int(N), D, int(N0), D0>;
template <size_t N, size_t D, size_t N0, size_t D0>
numeric_interval_state(OpNegFractionLiteral<N, D>, OpNegFractionLiteral<N0, D0>)
    -> numeric_interval_state<-int(N), D, -int(N0), D0>;

template <typename A, typename B>
constexpr auto numeric_range(A, B) {
  return numeric_interval_state(A{}, B{});
}

//! Represents the beginning of a constant interval.
/*!
 * Represents the beginning of a constant interval.
 *
 * \tparam N The numerator of the starting value.
 * \tparam D The denominator of the starting value.
 */
template <int N, size_t D>
struct numeric_range_state_start {
  template <typename A>
  constexpr numeric_range_state_start(A) {}

  template <typename B>
  auto operator>(B) const {
    if constexpr (N < 0) {
      return numeric_range(-expr::make_fraction<size_t(N), D>(), B{});
    } else {
      return numeric_range(expr::make_fraction<size_t(N), D>(), B{});
    }
  }
};

numeric_range_state_start(OpIdentity) -> numeric_range_state_start<1, 1>;
numeric_range_state_start(OpVoid) -> numeric_range_state_start<0, 1>;
numeric_range_state_start(OpNegIdentity) -> numeric_range_state_start<-1, 1>;
template <size_t N, size_t D>
numeric_range_state_start(OpFractionLiteral<N, D>)
    -> numeric_range_state_start<int(N), D>;
template <size_t N, size_t D>
numeric_range_state_start(OpNegFractionLiteral<N, D>)
    -> numeric_range_state_start<-int(N), D>;

template <typename A>
constexpr auto numeric_range(A) {
  return numeric_range_state_start(A{});
}
}  // namespace expr

constexpr auto OpVoid::operator--(int) const {
  return expr::numeric_range(OpVoid{});
}

constexpr auto OpIdentity::operator--(int) const {
  return expr::numeric_range(OpIdentity{});
}

constexpr auto OpNegIdentity::operator--(int) const {
  return expr::numeric_range(OpNegIdentity{});
}

namespace symphas::internal {
template <typename G, typename V, typename E, typename T, typename working_grid>
auto to_optimized(OpIntegral<V, E, T> const& e,
                  OpTerm<OpIdentity, working_grid*> const& data) {
  return expr::transform::swap_expression<G>(e, data);
}
}  // namespace symphas::internal

template <typename V, typename E, typename T, typename G>
struct OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>>
    : OpExpression<OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>>> {
  using this_type = OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>>;
  using working_grid = expr::storage_type_t<G>;

  using expr_t = std::invoke_result_t<
      decltype(&symphas::internal::to_optimized<G, V, E, T, working_grid>),
      OpIntegral<V, E, T>, OpTerm<OpIdentity, working_grid*>>;

  OpOptimized() : working{0}, e{}, term{} {}

  OpOptimized(OpBinaryMul<OpIntegral<V, E, T>, G> const& e)
      : working{expr::data_dimensions(e.b)},
        e{symphas::internal::to_optimized<G>(e.a, expr::make_term(&working))},
        term{e.b} {
    symphas::internal::update_temporary_grid(working, term);
    expr::prune::update(this->e);
  }

  OpOptimized(OpOptimized const& other) : OpOptimized(other.get_expression()) {}
  OpOptimized(OpOptimized&& other) noexcept : OpOptimized() {
    swap(*this, other);
  }
  OpOptimized& operator=(OpOptimized other) {
    swap(*this, other);
    return *this;
  }

  friend void swap(OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>>& first,
                   OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>>& second) {
    using std::swap;
    swap(first.working, second.working);
    swap(first.e, second.e);
    swap(first.term, second.term);
  }

  auto eval(iter_type n) const {
    return e.eval(n) * expr::BaseData<working_grid>::get(working, n);
  }

  template <typename eval_handler_type, typename... condition_ts>
  void update(eval_handler_type const& eval_handler,
              symphas::lib::types_list<condition_ts...>) {
    auto region = expr::iterable_domain(term);
    if (!grid::is_fully_overlapping(grid::get_iterable_domain(working),
                                    region)) {
      symphas::internal::update_temporary_grid(working, term);
    }
    eval_handler.result(term, working, region);
    expr::prune::update<condition_ts...>(e);
  }

  template <typename eval_handler_type>
  void update(eval_handler_type const& eval_handler) {
    update(eval_handler, symphas::lib::types_list<>{});
  }

  auto operator-() const;

#ifdef PRINTABLE_EQUATIONS

  size_t print(FILE* out) const { return get_expression().print(out); }

  size_t print(char* out) const { return get_expression().print(out); }

  size_t print_length() const { return get_expression().print_length(); }
#endif

  auto get_expression() const {
    return expr::transform::swap_grid<OpTerm<OpIdentity, working_grid*>>(e,
                                                                         term) *
           term;
  }

  working_grid working;
  expr_t e;
  G term;
};

template <typename E>
OpOptimized(E) -> OpOptimized<E>;

namespace expr::transform {
template <typename E>
auto optimize(OpEvaluable<E> const& e) {
  return *static_cast<E const*>(&e);
}

template <typename V, typename E, typename T, typename G,
          typename = std::enable_if_t<
              expr::satisfies<E, expr::contains_matching_anywhere<G>>, int>>
auto optimize(OpBinaryMul<OpIntegral<V, E, T>, G> const& e) {
  return OpOptimized(e);
}

template <typename... Es, size_t... Is>
auto optimize(OpAdd<Es...> const& e, std::index_sequence<Is...>) {
  return (optimize(expr::get<Is>(e)) + ...);
}

template <typename... Es>
auto optimize(OpAdd<Es...> const& e) {
  return optimize(e, std::make_index_sequence<sizeof...(Es)>{});
}
}  // namespace expr::transform

template <typename V, typename E, typename T, typename G>
auto OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>>::operator-() const {
  return expr::transform::optimize(-get_expression());
}
