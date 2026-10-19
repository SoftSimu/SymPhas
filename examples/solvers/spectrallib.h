
/* ***************************************************************************
 * This file is part of the SymPhas package, containing a framework for
 * implementing solvers for phase-field problems with compile-time symbolic
 * algebra.
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
 * This file supports the functionality of the semi-implicit Fourier
 * spectral solver (solversp.h)
 *
 * ***************************************************************************
 */

#pragma once

#include "expressiontypeincludes.h"
#include "gridpair.h"
#include "indexseqhelpers.h"
#include "spslib.h"

#include <functional>
#include <memory>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

namespace solver_sp {

namespace {

// **************************************************************************************

// swap for the primary variables (in model.system<N>())
template <Axis ax, size_t Z, size_t Q0, size_t N, typename T, typename E>
auto swap_var_apply_axes(MultiBlock<N, T> const& grid,
                         OpExpression<E> const& expr, const len_type* dims) {
  auto data = expr::as_grid_data<N>(grid, dims);
#ifdef PRINTABLE_EQUATIONS
  auto v = expr::get_variable<Q0>(*static_cast<E const*>(&expr));
  auto name = expr::get_fourier_name(expr::get_op_name(v));
  auto replace = expr::as_component<ax>(
      expr::as_variable<Z>(NamedData(std::move(data), name)));
#else
  auto replace = expr::as_component<ax>(expr::as_variable<Z>(std::move(data)));
#endif
  auto e = expr::transform::swap_grid<ax, Q0>(*static_cast<const E*>(&expr),
                                              replace);
  constexpr size_t I = (ax == Axis::X)   ? 0
                       : (ax == Axis::Y) ? 1
                       : (ax == Axis::Z) ? 2
                                         : -1;

  if constexpr (I + 1 < N) {
    constexpr Axis bx = (I == 0) ? Axis::Y : Axis::Z;
    return swap_var_apply_axes<bx, Z, Q0>(grid, e, dims);
  } else {
    return e;
  }
}

// swap for the primary variables (in model.system<N>())
template <size_t Z, size_t Q0, size_t N, typename T, typename E>
auto swap_var_apply(MultiBlock<N, T> const& grid, OpExpression<E> const& expr,
                    const len_type* dims) {
  auto data = expr::as_grid_data<N>(grid, dims);
#ifdef PRINTABLE_EQUATIONS
  auto name = expr::get_fourier_name(
      expr::get_op_name(expr::get_variable<Q0>(*static_cast<E const*>(&expr))));
  auto replace = expr::as_variable<Z>(NamedData(std::move(data), name));
#else
  auto replace = expr::as_variable<Z>(std::move(data));
#endif
  auto e = swap_var_apply_axes<Axis::X, Z, Q0>(
      grid, *static_cast<const E*>(&expr), dims);
  return expr::transform::swap_grid<Q0>(e, replace);
}

// swap for the primary variables (in model.system<N>())
template <size_t Z, size_t Q0, typename T, typename E>
auto swap_var_apply(T* const& grid, OpExpression<E> const& expr,
                    const len_type* dims) {
  constexpr size_t D = expr::grid_dim<E>::value;
  auto data = expr::as_grid_data<D>(grid, dims);

#ifdef PRINTABLE_EQUATIONS
  auto name = expr::get_fourier_name(
      expr::get_op_name(expr::get_variable<Q0>(*static_cast<E const*>(&expr))));
  auto replace = expr::as_variable<Z>(NamedData(std::move(data), name));
#else
  auto replace = expr::as_variable<Z>(std::move(data));
#endif
  return expr::transform::swap_grid<Q0>(*static_cast<const E*>(&expr), replace);
}

template <typename... S, typename E>
auto swap_var_apply(std::tuple<S...> const&, E&& expr, std::index_sequence<>) {
  return std::forward<E>(expr);
}

template <typename T, size_t D>
auto set_system_dimensions(System<T, D> const& s, len_type (&dims)[D]) {
  for (iter_type i = 0; i < D; ++i) {
    dims[i] = s.dims[i];
  }
}

// swap for the primary variables (in model.system<N>())
template <typename... S, typename E, size_t Q0, size_t... Qs>
auto swap_var_apply(std::tuple<S...> const& systems,
                    OpExpression<E> const& expr,
                    std::index_sequence<Q0, Qs...>) {
  constexpr size_t D = expr::grid_dim<E>::value;
  len_type dims[D];
  set_system_dimensions(std::get<Q0>(systems), dims);

  auto e = swap_var_apply<sizeof...(S) + Q0, Q0>(
      std::get<Q0>(systems).frame_t, *static_cast<E const*>(&expr), dims);

  if constexpr (sizeof...(Qs) == 0) {
    return e;
  } else {
    return swap_var_apply(systems, e, std::index_sequence<Qs...>{});
  }
}

}  // namespace

template <typename E>
auto drop_hcts(OpExpression<E> const& e) {
  return (*static_cast<E const*>(&e));
}

template <typename E>
auto drop_hcts(OpMap<symphas::internal::HCTS, OpIdentity, E> const& e) {
  return expr::get_enclosed_expression(e);
}

template <typename T, typename E>
auto sthc_apply_on_scalar(OpExpression<E> const& e) {
  using rt = typename symphas::internal::real_space_type<T>::type;
  if constexpr (std::is_same<rt, scalar_t>::value) {
    return expr::sthc(*static_cast<E const*>(&e));
  } else {
    return (*static_cast<E const*>(&e));
  }
}

template <typename T, typename E>
auto hcts_apply_on_scalar(OpExpression<E> const& e) {
  using rt = typename symphas::internal::real_space_type<T>::type;
  if constexpr (std::is_same<rt, scalar_t>::value) {
    return expr::hcts(*static_cast<E const*>(&e));
  } else {
    return (*static_cast<E const*>(&e));
  }
}

template <typename V, expr::NoiseType nt, typename T, size_t D, typename E,
          typename... Ts>
auto hcts_apply_on_scalar(OpSymbolicEval<V, NoiseData<nt, T, D, Grid>,
                                         SymbolicFunction<E, Ts...>> const& e) {
  return e;
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <size_t D, typename E, size_t R = expr::eval_type<E>::rank,
          size_t... Rs>
auto form_A_op(OpExpression<E> const& l_op, double dt, len_type const* dims,
               std::index_sequence<Rs...>) {
  auto A_expr =
      ((expr::make_column_vector<Rs, R>() *
        expr::exp(expr::make_literal(dt) * (expr::make_row_vector<Rs, R>() *
                                            (*static_cast<E const*>(&l_op))))) +
       ...);
  return A_expr;
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <size_t D, typename E, size_t R = expr::eval_type<E>::rank,
          typename std::enable_if_t<(R == 0), int> = 0>
auto form_A_op(OpExpression<E> const& l_op, double dt, len_type const* dims) {
  auto A_expr =
      expr::exp(expr::make_literal(dt) * (*static_cast<E const*>(&l_op)));
  return A_expr;
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <size_t D, typename E, size_t R = expr::eval_type<E>::rank,
          typename std::enable_if_t<(R > 0), int> = 0>
auto form_A_op(OpExpression<E> const& l_op, double dt, len_type const* dims) {
  return form_A_op<D>(*static_cast<E const*>(&l_op), dt, dims,
                      std::make_index_sequence<R>{});
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <size_t D>
auto form_A_op(OpVoid, double dt, len_type const* dims) {
  return OpIdentity{};
}

namespace {
//! (exp(dt L) - 1) / L, including its limit dt as L -> 0.
/*!
 * Evaluating (exp(z) - 1) / L directly gives 0/0 where L = 0 and loses
 * digits to cancellation where |z| = |dt L| is small, so a series is used
 * there.
 */
template <typename T>
T etd_b_value(T L, double dt) {
  T z = dt * L;
  if (std::abs(z) < 1e-4) {
    return dt * (T(1) + z * (1. / 2 + z * (1. / 6 + z * (1. / 24))));
  }
  return (std::exp(z) - T(1)) / L;
}

//! Evaluates the \f$B\f$ operator of a scalar linear operator into a grid.
template <size_t D, typename E>
auto make_B_term(OpExpression<E> const& l_op, double dt,
                 len_type const* dims) {
  using T = expr::eval_type_t<E>;
  Grid<T, D> B(dims);
  expr::result(*static_cast<E const*>(&l_op), B);
  for (iter_type n = 0; n < B.len; ++n) {
    B[n] = etd_b_value(B[n], dt);
  }
  return expr::make_term(NamedData(std::move(B), "(exp(L dt) - 1)/L"));
}
}  // namespace

//! Given the linear expression, returns the \f$B\f$ operator.
template <size_t D, typename E, size_t R = expr::eval_type<E>::rank,
          size_t... Rs>
auto form_B_op(OpExpression<E> const& l_op, double dt, len_type const* dims,
               std::index_sequence<Rs...>) {
  return ((expr::make_column_vector<Rs, R>() *
           make_B_term<D>(expr::make_row_vector<Rs, R>() *
                              (*static_cast<E const*>(&l_op)),
                          dt, dims)) +
          ...);
}

//! Given the linear expression, returns the \f$B\f$ operator.
template <size_t D, typename E, size_t R = expr::eval_type<E>::rank,
          typename std::enable_if_t<(R == 0), int> = 0>
auto form_B_op(OpExpression<E> const& l_op, double dt, len_type const* dims) {
  return make_B_term<D>(*static_cast<E const*>(&l_op), dt, dims);
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <size_t D, typename E, size_t R = expr::eval_type<E>::rank,
          typename std::enable_if_t<(R > 0), int> = 0>
auto form_B_op(OpExpression<E> const& l_op, double dt, len_type const* dims) {
  return form_B_op<D>(*static_cast<E const*>(&l_op), dt, dims,
                      std::make_index_sequence<R>{});
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <size_t D>
auto form_B_op(OpVoid, double dt, len_type const* dims) {
  return expr::make_literal(dt);
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T, typename A, size_t R = expr::eval_type<A>::rank,
          size_t... Rs>
auto get_A_term(OpExpression<A> const& A_expression,
                std::index_sequence<Rs...>);

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T, typename T0>
auto get_A_term(OpLiteral<T0> const& A_expression) {
  return A_expression;
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T>
inline auto get_A_term(OpIdentity) {
  return OpIdentity{};
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T>
inline auto get_A_term(OpVoid) {
  return OpVoid{};
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T, typename A, size_t R = expr::eval_type<A>::rank,
          typename std::enable_if_t<(R == 0), int> = 0>
auto get_A_term(OpExpression<A> const& A_expression) {
  auto a = expr::transform::to_grid(
      sthc_apply_on_scalar<T>(*static_cast<A const*>(&A_expression)));
  auto A_term = expr::make_term(
      NamedData(std::move(a), *static_cast<A const*>(&A_expression)));
  return A_term;
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T, typename A, size_t R = expr::eval_type<A>::rank,
          typename std::enable_if_t<(R > 0), int> = 0>
auto get_A_term(OpExpression<A> const& A_expression) {
  return get_A_term<T>(*static_cast<A const*>(&A_expression),
                       std::make_index_sequence<R>{});
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T, typename A, size_t R, size_t... Rs>
auto get_A_term(OpExpression<A> const& A_expression,
                std::index_sequence<Rs...>) {
  return ((expr::make_tensor<Rs, Rs, R, R>() *
           // The scalar get_A_term applies the half-complex arrangement (sthc)
           // to each component; applying it here as well would arrange twice.
           get_A_term<T>(expr::make_row_vector<Rs, R>() *
                         (*static_cast<A const*>(&A_expression)))) +
          ...);
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T, typename T0>
auto get_B_term(OpLiteral<T0> const& B_expression) {
  return B_expression;
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T, typename B, size_t R = expr::eval_type<B>::rank,
          size_t... Rs>
auto get_B_term(OpExpression<B> const& B_expression,
                std::index_sequence<Rs...>);

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T, typename B, size_t R = expr::eval_type<B>::rank,
          typename std::enable_if_t<(R == 0), int> = 0>
auto get_B_term(OpExpression<B> const& B_expression) {
  return expr::make_term(
      NamedData(expr::transform::to_grid(*static_cast<B const*>(&B_expression)),
                (*static_cast<B const*>(&B_expression))));
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T, typename B, size_t R = expr::eval_type<B>::rank,
          typename std::enable_if_t<(R > 0), int> = 0>
auto get_B_term(OpExpression<B> const& B_expression) {
  return get_B_term<T>(*static_cast<B const*>(&B_expression),
                       std::make_index_sequence<R>{});
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T>
inline auto get_B_term(OpIdentity) {
  return OpIdentity{};
}

//! Given the linear expression, returns the \f$A\f$ operator.
template <typename T, typename B, size_t R, size_t... Rs>
auto get_B_term(OpExpression<B> const& B_expression,
                std::index_sequence<Rs...>) {
  return ((expr::make_tensor<Rs, Rs, R, R>() *
           get_B_term<T>(expr::make_row_vector<Rs, R>() *
                         *static_cast<B const*>(&B_expression))) +
          ...);
}

//! Swap real for Fourier variables in the given expression.
/*!
 * Take a list of all the systems appearing in the phase field problem,
 * then for each index associated variable in the expression (which is the
 * real space order parameter), make the equivalent expression instead with
 *
 */
template <typename... S, typename E>
decltype(auto) swap_var(std::tuple<S...> const& systems,
                        OpExpression<E> const& expr) {
  return swap_var_apply(systems, *static_cast<const E*>(&expr),
                        expr::vars<E>::get_ids());
}

//! Apply the convolution theorem to construct the nonlinear scheme.
/*!
 * Apply the convolution theorem to construct the nonlinear scheme.
 */
template <size_t Zn, typename T1, typename T2, typename... S, typename E1,
          typename E2>
auto anti_convolve(std::tuple<S...> const&, OpExpression<E1> const& a,
                   OpExpression<E2> const& b) {
  auto ae = expr::make_inv_fourier_map<T1>(*static_cast<E1 const*>(&a));
  auto be = expr::make_inv_fourier_map<T2>(*static_cast<E2 const*>(&b));
  return expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename S2,
          typename V2, typename E1, typename E2>
auto anti_convolve(std::tuple<S...> const&, OpExpression<E1> const& a,
                   OpFourierTransform<S2, V2, E2> const& b) {
  auto e = expr::get_enclosed_expression(b);
  return expr::coeff(b) *
         expr::make_fourier_map(
             expr::make_inv_fourier_map<T1>(*static_cast<E1 const*>(&a)) * e);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename S1,
          typename V1, typename E1, typename E2>
auto anti_convolve(std::tuple<S...> const&,
                   OpFourierTransform<S1, V1, E1> const& a,
                   OpExpression<E2> const& b) {
  auto e = expr::get_enclosed_expression(a);
  return expr::coeff(a) *
         expr::make_fourier_map(
             e * expr::make_inv_fourier_map<T2>(*static_cast<E2 const*>(&b)));
}

template <size_t Zn, typename T1, typename T2, typename... S, typename S1,
          typename V1, typename E1, typename S2, typename V2, typename E2>
auto anti_convolve(std::tuple<S...> const&,
                   OpFourierTransform<S1, V1, E1> const& a,
                   OpFourierTransform<S2, V2, E2> const& b) {
  auto ae = expr::get_enclosed_expression(a);
  auto be = expr::get_enclosed_expression(b);
  return (expr::coeff(a) * expr::coeff(b)) * expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename T,
          size_t Z, typename G, typename E2,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables,
                   OpTerm<T, Variable<Z, G>> const& a,
                   OpExpression<E2> const& b) {
  auto ae = expr::transform::swap_grid<Z>(a, std::get<Z - Zn>(variables));
  auto be = expr::make_inv_fourier_map<T2>(*static_cast<E2 const*>(&b));
  return expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename T,
          size_t Z, typename G, typename S2, typename V2, typename E2,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables,
                   OpTerm<T, Variable<Z, G>> const& a,
                   OpFourierTransform<S2, V2, E2> const& b) {
  auto ae = expr::transform::swap_grid<Z>(a, std::get<Z - Zn>(variables));
  auto be = expr::get_enclosed_expression(b);
  return expr::coeff(b) * expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename E1,
          typename T, size_t Z, typename G, size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables, OpExpression<E1> const& a,
                   OpTerm<T, Variable<Z, G>> const& b) {
  auto ae = expr::make_inv_fourier_map<T1>(*static_cast<E1 const*>(&a));
  auto be = expr::transform::swap_grid<Z>(b, std::get<Z - Zn>(variables));
  return expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename S1,
          typename V1, typename E1, typename T, size_t Z, typename G,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables,
                   OpFourierTransform<S1, V1, E1> const& a,
                   OpTerm<T, Variable<Z, G>> const& b) {
  auto ae = expr::get_enclosed_expression(a);
  auto be = expr::transform::swap_grid<Z>(b, std::get<Z - Zn>(variables));
  return expr::coeff(a) * expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename T,
          size_t Z, Axis ax, typename G, typename E2,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables,
                   OpTerm<T, Variable<Z, VectorComponent<ax, G>>> const& a,
                   OpExpression<E2> const& b) {
  auto ae = expr::transform::swap_grid<ax, Z>(
      a, expr::as_component<ax>(std::get<Z - Zn>(variables)));
  auto be = expr::make_inv_fourier_map<T2>(*static_cast<E2 const*>(&b));
  return expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename T,
          size_t Z, Axis ax, typename G, typename S2, typename V2, typename E2,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables,
                   OpTerm<T, Variable<Z, VectorComponent<ax, G>>> const& a,
                   OpFourierTransform<S2, V2, E2> const& b) {
  auto ae = expr::transform::swap_grid<ax, Z>(
      a, expr::as_component<ax>(std::get<Z - Zn>(variables)));
  auto be = expr::get_enclosed_expression(b);
  return expr::coeff(b) * expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename E1,
          typename T, size_t Z, Axis ax, typename G,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables, OpExpression<E1> const& a,
                   OpTerm<T, Variable<Z, VectorComponent<ax, G>>> const& b) {
  auto ae = expr::make_inv_fourier_map<T1>(*static_cast<E1 const*>(&a));
  auto be = expr::transform::swap_grid<ax, Z>(
      b, expr::as_component<ax>(std::get<Z - Zn>(variables)));
  return expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename S1,
          typename V1, typename E1, typename T, size_t Z, Axis ax, typename G,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables,
                   OpFourierTransform<S1, V1, E1> const& a,
                   OpTerm<T, Variable<Z, VectorComponent<ax, G>>> const& b) {
  auto ae = expr::get_enclosed_expression(a);
  auto be = expr::transform::swap_grid<ax, Z>(
      b, expr::as_component<ax>(std::get<Z - Zn>(variables)));
  return expr::coeff(a) * expr::make_fourier_map(ae * be);
}

// where there is hcts

template <typename T, typename G>
using OpTermHCTS = OpMap<symphas::internal::HCTS, OpIdentity, OpTerm<T, G>>;

template <size_t Zn, typename T1, typename T2, typename... S, typename T,
          size_t Z, typename G, typename E2,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables,
                   OpTermHCTS<T, Variable<Z, G>> const& a,
                   OpExpression<E2> const& b) {
  auto ae = expr::transform::swap_grid<Z>(a.e, std::get<Z - Zn>(variables));
  auto be = expr::make_inv_fourier_map<T2>(*static_cast<E2 const*>(&b));
  return expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename T,
          size_t Z, typename G, typename S2, typename V2, typename E2,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables,
                   OpTermHCTS<T, Variable<Z, G>> const& a,
                   OpFourierTransform<S2, V2, E2> const& b) {
  auto ae = expr::transform::swap_grid<Z>(a.e, std::get<Z - Zn>(variables));
  auto be = expr::get_enclosed_expression(b);
  return expr::coeff(b) * expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename E1,
          typename T, size_t Z, typename G, size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables, OpExpression<E1> const& a,
                   OpTermHCTS<T, Variable<Z, G>> const& b) {
  auto ae = expr::make_inv_fourier_map<T1>(*static_cast<E1 const*>(&a));
  auto be = expr::transform::swap_grid<Z>(b.e, std::get<Z - Zn>(variables));
  return expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename S1,
          typename V1, typename E1, typename T, size_t Z, typename G,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables,
                   OpFourierTransform<S1, V1, E1> const& a,
                   OpTermHCTS<T, Variable<Z, G>> const& b) {
  auto ae = expr::get_enclosed_expression(a);
  auto be = expr::transform::swap_grid<Z>(b.e, std::get<Z - Zn>(variables));
  return expr::coeff(a) * expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename T,
          size_t Z, Axis ax, typename G, typename E2,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables,
                   OpTermHCTS<T, Variable<Z, VectorComponent<ax, G>>> const& a,
                   OpExpression<E2> const& b) {
  auto ae = expr::transform::swap_grid<ax, Z>(
      a.e, expr::as_component<ax>(std::get<Z - Zn>(variables)));
  auto be = expr::make_inv_fourier_map<T2>(*static_cast<E2 const*>(&b));
  return expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename T,
          size_t Z, Axis ax, typename G, typename S2, typename V2, typename E2,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(std::tuple<S...> const& variables,
                   OpTermHCTS<T, Variable<Z, VectorComponent<ax, G>>> const& a,
                   OpFourierTransform<S2, V2, E2> const& b) {
  auto ae = expr::transform::swap_grid<ax, Z>(
      a.e, expr::as_component<ax>(std::get<Z - Zn>(variables)));
  auto be = expr::get_enclosed_expression(b);
  return expr::coeff(b) * expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename E1,
          typename T, size_t Z, Axis ax, typename G,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(
    std::tuple<S...> const& variables, OpExpression<E1> const& a,
    OpTermHCTS<T, Variable<Z, VectorComponent<ax, G>>> const& b) {
  auto ae = expr::make_inv_fourier_map<T1>(*static_cast<E1 const*>(&a));
  auto be = expr::transform::swap_grid<ax, Z>(
      b.e, expr::as_component<ax>(std::get<Z - Zn>(variables)));
  return expr::make_fourier_map(ae * be);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename S1,
          typename V1, typename E1, typename T, size_t Z, Axis ax, typename G,
          size_t D = expr::grid_dim<G>::value>
auto anti_convolve(
    std::tuple<S...> const& variables, OpFourierTransform<S1, V1, E1> const& a,
    OpTermHCTS<T, Variable<Z, VectorComponent<ax, G>>> const& b) {
  auto ae = expr::get_enclosed_expression(a);
  auto be = expr::transform::swap_grid<ax, Z>(
      b.e, expr::as_component<ax>(std::get<Z - Zn>(variables)));
  return expr::coeff(a) * expr::make_fourier_map(ae * be);
}

// end of terms

template <size_t Zn, typename T1, typename T2, typename... S, typename... As,
          typename E1, size_t... Is>
auto anti_convolve(std::tuple<S...> const& variables, E1 const& a,
                   OpAdd<As...> const& b, std::index_sequence<Is...>) {
  return (anti_convolve<Zn, T1, T2>(variables, a, expr::get<Is>(b)) + ...);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename... As,
          typename E1>
auto anti_convolve(std::tuple<S...> const& variables, E1 const& a,
                   OpAdd<As...> const& b) {
  return anti_convolve<Zn, T1, T2>(variables, a, b,
                                   std::make_index_sequence<sizeof...(As)>{});
}

template <size_t Zn, typename T1, typename T2, typename... S, typename... As,
          typename E2, size_t... Is>
auto anti_convolve(std::tuple<S...> const& variables, OpAdd<As...> const& a,
                   E2 const& b, std::index_sequence<Is...>) {
  return (anti_convolve<Zn, T1, T2>(variables, expr::get<Is>(a), b) + ...);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename... As,
          typename E2>
auto anti_convolve(std::tuple<S...> const& variables, OpAdd<As...> const& a,
                   E2 const& b) {
  return anti_convolve<Zn, T1, T2>(variables, a, b,
                                   std::make_index_sequence<sizeof...(As)>{});
}

template <size_t Zn, typename T1, typename T2, typename... S, typename... As,
          typename... Bs, size_t... Is, size_t... Js>
auto anti_convolve(std::tuple<S...> const& variables, OpAdd<As...> const& a,
                   OpAdd<Bs...> const& b, std::index_sequence<Is...>,
                   std::index_sequence<Js...>) {
  return (
      anti_convolve<Zn, T1, T2>(variables, expr::get<Is>(a), expr::get<Js>(b)) +
      ...);
}

template <size_t Zn, typename T1, typename T2, typename... S, typename... As,
          typename... Bs>
auto anti_convolve(std::tuple<S...> const& variables, OpAdd<As...> const& a,
                   OpAdd<Bs...> const& b) {
  using a_seq = std::make_index_sequence<sizeof...(As)>;
  using b_seq = std::make_index_sequence<sizeof...(Bs)>;
  using a_seq_repeat =
      symphas::lib::repeat_each_element_t<sizeof...(Bs), a_seq>;
  using b_seq_repeat =
      symphas::lib::repeat_all_elements_t<sizeof...(As), b_seq>;
  return anti_convolve<Zn, T1, T2>(variables, a, b, a_seq_repeat{},
                                   b_seq_repeat{});
}

//! Convert terms of an expression to Fourier space.
/*!
 * The terms which cannot be converted are directly returned. For example,
 * constants and variables are directly returned. Only operators,
 * multiplications and convolutions are converted.
 *
 * This is similar to expr::transform::to_ft(), but it ignores constants
 * and variables. It assumes that variables have been substituted for the
 * Fourier space equivalent.
 */

template <size_t Z0, size_t D, typename... S, typename B, typename Dd,
          typename V, typename E, typename Sp>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpDerivative<Dd, V, E, Sp> const& e, double const* h,
                         const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename B, size_t O, typename V,
          typename Sp>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpOperatorDerivative<O, V, Sp> const& e,
                         double const* h, const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename B, typename A1,
          typename A2, typename E>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpCombination<A1, A2, E> const& e, double const* h,
                         const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename B, typename A1,
          typename A2, typename E>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpChain<A1, A2, E> const& e, double const* h,
                         const len_type* dims);
template <size_t Z0, size_t D, size_t Z, typename... S, typename B, typename T,
          typename G, typename = std::enable_if_t<(Z < sizeof...(S)), int>>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpTerm<T, Variable<Z, G>> const& e, double const* h,
                         const len_type* dims);
// template<size_t Z0, size_t D, expr::NoiseType nt, typename... S, typename B,
// typename T, typename V> auto construct_nonlinear(std::tuple<S...> const&
// systems, OpExpression<B> const& bop, 	OpTerm<V, NoiseData<nt, T, D>>
// const& e,
// double const* h, const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename B, typename... Es>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop, OpAdd<Es...> const& e,
                         double const* h, const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename B, expr::exp_key_t X,
          typename V, typename E>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop, OpPow<X, V, E> const& e,
                         double const* h, const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename B, typename E1,
          typename E2>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpBinaryMul<E1, E2> const& e, double const* h,
                         const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename B, typename E1,
          typename E2>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpBinaryDiv<E1, E2> const& e, double const* h,
                         const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename B, typename E2>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpBinaryDiv<OpIdentity, E2> const& e, double const* h,
                         const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename B, typename E>
auto construct_nonlinear(
    std::tuple<S...> const& systems, OpExpression<B> const& bop,
    OpMap<symphas::internal::PoissonSolver, OpIdentity, E> const& e,
    double const* h, const len_type* dims);

//! Inverse transforms shared by all the terms of one equation.
/*!
 * Several terms of an equation often need the same field in real space, such
 * as dn/dx. Each distinct inverse transform is made once, stored here, and
 * every term refers to its result; the transforms are updated once per step,
 * before the equation is evaluated (SpectralData::update).
 *
 * Transforms are identified by the type of their Fourier-space expression.
 * The type only determines the value when the expression holds no run-time
 * coefficients (variables, wave vectors and type-level constants only), so
 * expressions with an OpLiteral are not shared.
 */
struct InverseMapCache {
  std::vector<std::shared_ptr<void>> maps;
  std::vector<std::function<void()>> updaters;
  std::unordered_map<std::string, void*> index;

  void update() {
    for (auto& u : updaters) u();
  }

  //! The cache of the equation being constructed, or nullptr.
  static InverseMapCache*& current() {
    static InverseMapCache* c = nullptr;
    return c;
  }

  //! A stored transform; exposes the grid that holds its result.
  template <typename M>
  struct Entry : M {
    Entry(M const& m) : M(m) {}
    auto& grid() { return this->data; }
  };

  template <typename M>
  Entry<M>& get(M const& map) {
    std::string key = typeid(M).name();
    auto it = index.find(key);
    if (it != index.end()) {
      return *static_cast<Entry<M>*>(it->second);
    }
    auto m = std::make_shared<Entry<M>>(map);
    Entry<M>* ptr = m.get();
    maps.push_back(m);
    index[key] = ptr;
    updaters.push_back([ptr, first = true]() mutable {
      M& base = *ptr;
      if (first) {
        base.allocate();
        first = false;
      }
      expr::prune::update(base);
    });
    return *ptr;
  }
};

//! Whether a type holds run-time values that its type does not determine.
/*!
 * Run-time coefficients (OpLiteral) and data grids (e.g. a stored
 * wave-vector polynomial) are such values. Variables are identified by
 * their index, and wave vectors by their axis and order (all transforms of
 * one equation are on the same grid), so their types determine them.
 */
template <typename T>
struct holds_literal : std::false_type {};
template <typename T>
struct holds_literal<OpLiteral<T>> : std::true_type {};
template <typename T>
struct holds_literal<NamedData<T>> : std::true_type {};
template <size_t O, size_t D>
struct holds_literal<NamedData<WaveVectorData<O, D>>> : std::false_type {};
template <Axis ax, size_t O, size_t D>
struct holds_literal<NamedData<WaveVectorDataAxis<ax, O, D>>> : std::false_type {};
template <template <typename...> typename C, typename... Ts>
struct holds_literal<C<Ts...>>
    : std::integral_constant<bool, (holds_literal<Ts>::value || ...)> {};
template <typename V, size_t... Ns>
struct holds_literal<OpTensor<V, Ns...>> : holds_literal<V> {};
template <size_t Z, typename G>
struct holds_literal<Variable<Z, G>> : std::false_type {};
template <typename G, expr::exp_key_t X>
struct holds_literal<Term<G, X>> : holds_literal<G> {};
template <Axis ax, typename G>
struct holds_literal<VectorComponent<ax, G>> : holds_literal<G> {};

template <typename T, typename V, typename E, size_t... Rs>
auto shared_inv_fourier_map_components(V const& value, E const& e,
                                       std::index_sequence<Rs...>);

//! An inverse transform, shared through the current InverseMapCache.
/*!
 * Same result as expr::make_inv_fourier_map(value, e). A vector expression
 * is split into components, as make_inv_fourier_map does, and each scalar
 * transform without run-time coefficients is shared.
 */
template <typename T, typename V, typename E>
auto shared_inv_fourier_map(V const& value, OpExpression<E> const& e) {
  constexpr size_t R = expr::eval_type<E>::rank;
  using rt = typename symphas::internal::real_space_type<T>::type;
  if constexpr (R > 0) {
    return shared_inv_fourier_map_components<T>(
        value, *static_cast<E const*>(&e), std::make_index_sequence<R>{});
  } else if constexpr (std::is_same<rt, scalar_t>::value &&
                       !holds_literal<E>::value) {
    auto map = expr::make_inv_fourier_map<T>(OpIdentity{}, *static_cast<E const*>(&e));
    using M = decltype(map);
    InverseMapCache* cache = InverseMapCache::current();
    if (cache == nullptr) {
      fprintf(SYMPHAS_ERR, "spectral solver: no transform cache while forming the scheme\n");
      exit(1);
    }
    auto& shared = cache->get(map);
    constexpr size_t D = expr::grid_dim<E>::value;
    return value * expr::make_term(NamedData(
                       std::ref(static_cast<Grid<scalar_t, D>&>(shared.grid())),
                       static_cast<M const&>(shared)));
  } else {
    return expr::make_inv_fourier_map<T>(value, *static_cast<E const*>(&e));
  }
}

template <typename T, typename V>
auto shared_inv_fourier_map(V const& value, OpVoid) {
  return OpVoid{};
}

template <typename T, typename V, typename E, size_t... Rs>
auto shared_inv_fourier_map_components(V const& value, E const& e,
                                       std::index_sequence<Rs...>) {
  constexpr size_t R = sizeof...(Rs);
  return (shared_inv_fourier_map<T>(expr::make_column_vector<Rs, R>() * value,
                                    expr::make_row_vector<Rs, R>() * e) +
          ...);
}

//! Real-space form of a nonlinear expression (see the definitions below). All
//! overloads are declared here, so that they see each other when nested.
template <size_t Z0, size_t D, typename... S, typename E>
auto real_space_form(std::tuple<S...> const& systems, OpExpression<E> const& e,
                     double const* h, const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename V, typename... Ts>
auto real_space_form(std::tuple<S...> const& systems,
                     OpTerms<V, Ts...> const& e, double const* h,
                     const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename E>
auto real_space_form(std::tuple<S...> const& systems, OpOperator<E> const& e,
                     double const* h, const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename E>
auto real_space_form(std::tuple<S...> const& systems,
                     OpOperatorChain<OpIdentity, E> const& e, double const* h,
                     const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename E1, typename E2>
auto real_space_form(std::tuple<S...> const& systems,
                     OpBinaryMul<E1, E2> const& e, double const* h,
                     const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename... Es>
auto real_space_form(std::tuple<S...> const& systems, OpAdd<Es...> const& e,
                     double const* h, const len_type* dims);
template <size_t Z0, size_t D, typename... S, expr::exp_key_t X, typename V,
          typename E>
auto real_space_form(std::tuple<S...> const& systems, OpPow<X, V, E> const& e,
                     double const* h, const len_type* dims);
template <size_t Z0, size_t D, typename... S, typename Dd, typename V,
          typename E, typename Sp>
auto real_space_form(std::tuple<S...> const& systems,
                     OpDerivative<Dd, V, E, Sp> const& e, double const* h,
                     const len_type* dims);

template <typename S>
// using variable_type = decltype(std::declval<S>()[0]);
using variable_type = typename grid::value_type_of<S>::type;

namespace {
// Some expressions (e.g. applied derivatives) only declare their rank and
// report expr::symbols::Symbol as their evaluation type, so the real space
// type is reconstructed from the rank.
template <typename E, typename T = expr::eval_type_t<E>,
          size_t R = expr::eval_type<E>::rank>
struct real_space_eval_type {
  using type = T;
};

template <typename E, size_t R>
struct real_space_eval_type<E, expr::symbols::Symbol, R> {
  using type = std::conditional_t<(R == 0), scalar_t, vector_t<R>>;
};

template <typename E>
using real_space_eval_t = typename real_space_eval_type<E>::type;
}  // namespace

//! Fourier transform of a spatially uniform value of 1.
/*!
 * The transform is unnormalized (the inverse divides by the number of
 * points N), so a uniform value of 1 transforms to N at k = 0 and zero at
 * every other wave vector. Index 0 is k = 0 in both the full and half-complex
 * arrangements.
 */
template <size_t D>
auto uniform_ft(const len_type* dims) {
  Grid<scalar_t, D> delta(dims);
  std::fill(delta.values, delta.values + delta.len, scalar_t{});
  delta.values[0] = static_cast<scalar_t>(delta.len);
  return expr::make_term(NamedData(std::move(delta), "N.delta_k"));
}

template <size_t Z0, size_t D, typename... S, typename B, typename E>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop, OpExpression<E> const& e,
                         double const* h, const len_type* dims) {
  using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;

  auto B_term = get_B_term<vt>(*static_cast<B const*>(&bop));
  if constexpr (expr::grid_dim<E>::value == 0) {
    // A term without spatial data (a constant such as an applied field, or a
    // factor that only depends on time) is uniform in space. to_ft() returns
    // it unchanged, which has no data to be evaluated over and contributed
    // nothing.
    return B_term * (*static_cast<E const*>(&e) * uniform_ft<D>(dims));
  } else {
    auto ft = expr::transform::to_ft<D>(*static_cast<E const*>(&e), h, dims);
    return B_term * ft;
  }
}

template <size_t Z0, size_t D, typename... S, typename B, typename Dd,
          typename V, typename E, typename Sp>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpDerivative<Dd, V, E, Sp> const& e, double const* h,
                         const len_type* dims) {
  using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;

  auto [op, en] = expr::split::separate_operator(e);
  auto B_expression =
      (*static_cast<B const*>(&bop)) * expr::transform::to_ft<D>(op, h, dims);
  auto B_term = get_B_term<vt>(B_expression);

  auto ft = construct_nonlinear<Z0, D>(systems, OpIdentity{}, en, h, dims);
  return B_term * ft;
}

template <size_t Z0, size_t D, typename... S, typename B, typename V,
          typename E1, typename E2>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpConvolution<V, E1, E2> const& e, double const* h,
                         const len_type* dims) {
  return expr::coeff(e) *
         construct_nonlinear<Z0, D>(systems, *static_cast<B const*>(&bop), e.a,
                                    h, dims) *
         construct_nonlinear<Z0, D>(systems, OpIdentity{},
                                    *static_cast<B const*>(&bop), e.b, h, dims);
}

template <size_t Z0, size_t D, typename... S, typename B, typename V,
          template <typename, size_t> typename grid_type, typename E>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpConvolution<V, GaussianSmoothing<D, grid_type>, E> const& e,
                         double const* h, const len_type* dims) {
  return expr::coeff(e) * expr::transform::to_ft<D>(e.smoother, h, dims) *
         construct_nonlinear<Z0, D>(systems, *static_cast<B const*>(&bop),
                                    expr::get_enclosed_expression(e), h, dims);
}

template <size_t Z0, size_t D, typename... S, typename B, size_t O, typename V,
          typename Sp>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpOperatorDerivative<O, V, Sp> const& e,
                         double const* h, const len_type* dims) {
  using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;

  auto [op, en] = expr::split::separate_operator(e);
  auto B_expression =
      (*static_cast<B const*>(&bop)) * expr::transform::to_ft<D>(op, h, dims);
  auto B_term = get_B_term<vt>(B_expression);

  return B_term *
         construct_nonlinear<Z0, D>(systems, OpIdentity{}, en, h, dims);
}

template <size_t Z0, size_t D, typename... S, typename B, typename A1,
          typename A2, typename E>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpCombination<A1, A2, E> const& e, double const* h,
                         const len_type* dims) {
  using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;

  auto [op, en] = expr::split::separate_operator(e);
  auto B_expression =
      (*static_cast<B const*>(&bop)) * expr::transform::to_ft<D>(op, h, dims);
  auto B_term = get_B_term<vt>(B_expression);

  return B_term *
         construct_nonlinear<Z0, D>(systems, OpIdentity{}, en, h, dims);
}

template <size_t Z0, size_t D, typename... S, typename B, typename A1,
          typename A2, typename E>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpChain<A1, A2, E> const& e, double const* h,
                         const len_type* dims) {
  using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;

  auto [op, en] = expr::split::separate_operator(e);
  auto B_expression =
      (*static_cast<B const*>(&bop)) * expr::transform::to_ft<D>(op, h, dims);
  auto B_term = get_B_term<vt>(B_expression);

  return B_term *
         construct_nonlinear<Z0, D>(systems, OpIdentity{}, en, h, dims);
}

template <size_t Z0, size_t D, size_t Z, typename... S, typename B, typename T,
          typename G, typename>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpTerm<T, Variable<Z, G>> const& e, double const* h,
                         const len_type* dims) {
  using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;

  auto B_term = get_B_term<vt>(*static_cast<B const*>(&bop));
  return B_term *
         solver_sp::hcts_apply_on_scalar<vt>(solver_sp::swap_var(systems, e));
}

template <size_t Z0, size_t D, typename... S, typename B, typename E>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpOperatorChain<OpIdentity, E> const& e,
                         double const* h, const len_type* dims) {
  return construct_nonlinear<Z0, D>(systems, *static_cast<B const*>(&bop), e.g,
                                    h, dims);
}
/*
 * Transform grouping.
 *
 * The nonlinear part of an equation is a sum of terms, and the transform is
 * linear, so all terms that share one Fourier-space multiplier can be summed
 * in real space and transformed together:
 *  - pointwise products of fields (n^3, n^2 M, ...) are summed and transformed
 *    once;
 *  - derivatives of nonlinear expressions, such as lap(n^2) or
 *    lap((M.grad n) div M), are grouped by derivative, and each group is
 *    transformed once and multiplied by the derivative's wave-vector factor;
 *  - the other terms (constants, terms that are linear in Fourier space such
 *    as the Poisson solve) are constructed one by one as before.
 * Inside a product, each factor that is not a field (a derivative, a sum
 * with derivatives, ...) is brought to real space with one inverse transform
 * (real_space_form), instead of transforming the product of two factors at
 * a time (anti_convolve), which also transformed every intermediate product
 * forward and back.
 */
namespace sp_group {

template <typename Sp>
struct is_symbolic_solver : std::false_type {};
template <typename G>
struct is_symbolic_solver<SymbolicDerivative<G>> : std::true_type {};

//! Whether evaluating E needs a pointwise product of fields in real space.
template <typename E>
struct needs_real_space : std::false_type {};

template <typename V, typename... Ts>
struct needs_real_space<OpTerms<V, Ts...>>
    : std::integral_constant<bool, (sizeof...(Ts) > 0) &&
                                       !expr::is_linear<OpTerms<V, Ts...>>::value> {};

template <typename E1, typename E2>
struct needs_real_space<OpBinaryMul<E1, E2>>
    : std::integral_constant<bool, expr::is_coeff<E1>   ? needs_real_space<E2>::value
                                   : expr::is_coeff<E2> ? needs_real_space<E1>::value
                                                        : true> {};

template <expr::exp_key_t X, typename V, typename E>
struct needs_real_space<OpPow<X, V, E>> : std::true_type {};

template <typename... Es>
struct needs_real_space<OpAdd<Es...>>
    : std::integral_constant<bool, (needs_real_space<Es>::value || ...)> {};

//! The derivative a term is grouped under, or void if it is not a spatial
//! derivative of a nonlinear expression with a scalar coefficient. Even-order derivative<ax, O> are
//! powers of the laplacian, whatever their axis tag, so they share one key.
template <typename E>
struct deriv_key {
  using type = void;
};

template <typename Dd, typename V, typename E, typename Sp>
struct deriv_key<OpDerivative<Dd, V, E, Sp>> {
  static const bool even_laplacian =
      !Dd::is_directional && !Dd::is_mixed && (Dd::order % 2 == 0);
  using norm_t =
      std::conditional_t<even_laplacian,
                         typename Solver<Sp>::template derivative<Axis::X, Dd::order>,
                         Dd>;
  // A tensor coefficient (e.g. a basis vector) is kept outside: moved into
  // the enclosed product it would be pushed into its factors.
  using type = std::conditional_t<!is_symbolic_solver<Sp>::value &&
                                      needs_real_space<E>::value &&
                                      !expr::is_tensor<V>,
                                  norm_t, void>;
};

template <typename E>
using deriv_key_t = typename deriv_key<E>::type;

//! 0: other term, 1: pointwise product, 2: derivative of a nonlinear expression
template <typename E>
constexpr int kind = !std::is_void<deriv_key_t<E>>::value ? 2
                     : needs_real_space<E>::value         ? 1
                                                          : 0;

template <typename... Ks>
struct key_list {};

template <typename L, typename K>
struct key_add;

template <typename... Ks, typename K>
struct key_add<key_list<Ks...>, K> {
  using type = std::conditional_t<std::is_void<K>::value ||
                                      (std::is_same<K, Ks>::value || ...),
                                  key_list<Ks...>, key_list<Ks..., K>>;
};

template <typename L, typename... Es>
struct keys_of {
  using type = L;
};

template <typename L, typename E0, typename... Es>
struct keys_of<L, E0, Es...> {
  using type =
      typename keys_of<typename key_add<L, deriv_key_t<E0>>::type, Es...>::type;
};

//! The solver of the first term grouped under K.
template <typename K, typename... Es, size_t I0, size_t... Is>
auto const& group_solver(OpAdd<Es...> const& e, std::index_sequence<I0, Is...>) {
  using E0 = std::decay_t<decltype(expr::get<I0>(e))>;
  if constexpr (std::is_same<deriv_key_t<E0>, K>::value) {
    return expr::get<I0>(e).solver;
  } else {
    return group_solver<K>(e, std::index_sequence<Is...>{});
  }
}

//! The enclosed expression, with its coefficient, of a term grouped under K.
template <typename K, typename E>
auto group_part(E const& term) {
  if constexpr (std::is_same<deriv_key_t<E>, K>::value) {
    return expr::coeff(term) * expr::get_enclosed_expression(term);
  } else {
    return OpVoid{};
  }
}

}  // namespace sp_group

//! Real-space form of a nonlinear expression.
/*!
 * Returns an expression that is evaluated pointwise in real space. Fields
 * and constants are used as they are; sums and products are formed from the
 * real-space forms of their parts; every other factor (a derivative, or a
 * sum with derivatives) is constructed in Fourier space and brought back
 * with one inverse transform, as anti_convolve did for each factor.
 */
template <size_t Z0, size_t D, typename... S, typename E>
auto real_space_form(std::tuple<S...> const& systems, OpExpression<E> const& e,
                     double const* h, const len_type* dims) {
  // Constants and expressions without grid data (functions of the time t)
  // have the same value at every point: they are used as they are.
  if constexpr (expr::is_coeff<E> || expr::grid_dim<E>::value == 0) {
    return *static_cast<E const*>(&e);
  } else {
    return shared_inv_fourier_map<real_space_eval_t<E>>(OpIdentity{},
        construct_nonlinear<Z0, D>(systems, OpIdentity{},
                                   *static_cast<E const*>(&e), h, dims));
  }
}

template <size_t Z0, size_t D, typename... S, typename V, typename... Ts>
auto real_space_form(std::tuple<S...> const& systems,
                     OpTerms<V, Ts...> const& e, double const* h,
                     const len_type* dims) {
  return e;
}

//! An operator used as a factor, such as the row vector of a dot product.
template <size_t Z0, size_t D, typename... S, typename E>
auto real_space_form(std::tuple<S...> const& systems, OpOperator<E> const& e,
                     double const* h, const len_type* dims) {
  return shared_inv_fourier_map<real_space_eval_t<E>>(OpIdentity{},
      construct_nonlinear<Z0, D>(systems, OpIdentity{},
                                 *static_cast<E const*>(&e), h, dims));
}

template <size_t Z0, size_t D, typename... S, typename E1, typename E2>
auto real_space_form(std::tuple<S...> const& systems,
                     OpBinaryMul<E1, E2> const& e, double const* h,
                     const len_type* dims) {
  if constexpr (expr::grid_dim<OpBinaryMul<E1, E2>>::value == 0) {
    return e;  // no grid data (e.g. a function of the time t)
  } else {
    return expr::make_mul(real_space_form<Z0, D>(systems, e.a, h, dims),
                          real_space_form<Z0, D>(systems, e.b, h, dims));
  }
}

//! A power of an expression (e.g. (M.grad n)^5): the power of the real-space
//! form of its base, evaluated pointwise.
template <size_t Z0, size_t D, typename... S, expr::exp_key_t X, typename V,
          typename E>
auto real_space_form(std::tuple<S...> const& systems, OpPow<X, V, E> const& e,
                     double const* h, const len_type* dims) {
  if constexpr (expr::grid_dim<OpPow<X, V, E>>::value == 0) {
    return e;  // no grid data (e.g. a function of the time t)
  } else {
    return expr::make_pow<X>(e.value,
                             real_space_form<Z0, D>(systems, e.e, h, dims));
  }
}

//! A derivative: its coefficient is applied in real space, so that the
//! transform of the derivative itself can be shared with other terms.
template <size_t Z0, size_t D, typename... S, typename Dd, typename V,
          typename E, typename Sp>
auto real_space_form(std::tuple<S...> const& systems,
                     OpDerivative<Dd, V, E, Sp> const& e, double const* h,
                     const len_type* dims) {
  using vt = OpDerivative<Dd, V, E, Sp>;
  if constexpr (std::is_same<V, OpIdentity>::value) {
    return shared_inv_fourier_map<real_space_eval_t<vt>>(
        OpIdentity{}, construct_nonlinear<Z0, D>(systems, OpIdentity{}, e, h, dims));
  } else {
    return expr::coeff(e) *
           real_space_form<Z0, D>(
               systems,
               expr::make_derivative<Dd>(expr::get_enclosed_expression(e), e.solver),
               h, dims);
  }
}

//! The row vector of a dot product: when it holds products (the gradient
//! of a product expanded by the product rule) or is a vector of derivatives
//! (a gradient), it is kept and its content is formed in real space, so that
//! products are not transformed forward and back and the derivatives can be
//! shared with other terms.
template <size_t Z0, size_t D, typename... S, typename E>
auto real_space_form(std::tuple<S...> const& systems,
                     OpOperatorChain<OpIdentity, E> const& e, double const* h,
                     const len_type* dims) {
  if constexpr (sp_group::needs_real_space<E>::value ||
                (expr::eval_type<E>::rank > 0)) {
    auto content = real_space_form<Z0, D>(systems, e.g, h, dims);
    return OpOperatorChain<OpIdentity, decltype(content)>(OpIdentity{}, content);
  } else {
    return shared_inv_fourier_map<real_space_eval_t<OpOperatorChain<OpIdentity, E>>>(OpIdentity{},
        construct_nonlinear<Z0, D>(systems, OpIdentity{}, e, h, dims));
  }
}

namespace {
template <size_t Z0, size_t D, typename... S, typename... Es, size_t... Is>
auto real_space_form_adds(std::tuple<S...> const& systems,
                          OpAdd<Es...> const& e, double const* h,
                          const len_type* dims, std::index_sequence<Is...>) {
  return (OpVoid{} + ... +
          real_space_form<Z0, D>(systems, expr::get<Is>(e), h, dims));
}
}  // namespace

template <size_t Z0, size_t D, typename... S, typename... Es>
auto real_space_form(std::tuple<S...> const& systems, OpAdd<Es...> const& e,
                     double const* h, const len_type* dims) {
  // A sum without grid data (e.g. a function of the time t) is the same
  // everywhere and is used as it is. A sum with products is formed term by
  // term, and so is a vector-valued sum (a gradient, whose components are
  // single derivatives, which can be shared with other terms).
  if constexpr (expr::grid_dim<OpAdd<Es...>>::value == 0) {
    return e;
  } else if constexpr (sp_group::needs_real_space<OpAdd<Es...>>::value ||
                       (expr::eval_type<OpAdd<Es...>>::rank > 0)) {
    return real_space_form_adds<Z0, D>(
        systems, e, h, dims, std::make_index_sequence<sizeof...(Es)>{});
  } else {  // e.g. a sum of derivatives (div M): one inverse transform of the sum
    return shared_inv_fourier_map<real_space_eval_t<OpAdd<Es...>>>(OpIdentity{},
        construct_nonlinear<Z0, D>(systems, OpIdentity{}, e, h, dims));
  }
}

namespace {
template <size_t Z0, size_t D, typename... S, typename B, typename E>
auto construct_other(std::tuple<S...> const& systems,
                     OpExpression<B> const& bop, E const& term, double const* h,
                     const len_type* dims) {
  if constexpr (sp_group::kind<E> == 0) {
    using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;
    auto B_term = solver_sp::get_B_term<vt>(*static_cast<B const*>(&bop));
    return B_term *
           construct_nonlinear<Z0, D>(systems, OpIdentity{}, term, h, dims);
  } else {
    return OpVoid{};
  }
}

template <size_t Z0, size_t D, typename... S, typename E>
auto pointwise_part(std::tuple<S...> const& systems, E const& term,
                    double const* h, const len_type* dims) {
  if constexpr (sp_group::kind<E> == 1) {
    return real_space_form<Z0, D>(systems, term, h, dims);
  } else {
    return OpVoid{};
  }
}

template <typename K, size_t Z0, size_t D, typename... S, typename B,
          typename... Es, size_t... Is>
auto construct_group(std::tuple<S...> const& systems,
                     OpExpression<B> const& bop, OpAdd<Es...> const& e,
                     double const* h, const len_type* dims,
                     std::index_sequence<Is...> seq) {
  auto inner = (OpVoid{} + ... + sp_group::group_part<K>(expr::get<Is>(e)));
  auto const& solver = sp_group::group_solver<K>(e, seq);
  return construct_nonlinear<Z0, D>(systems, *static_cast<B const*>(&bop),
                                    expr::make_derivative<K>(inner, solver), h,
                                    dims);
}

template <size_t Z0, size_t D, typename... S, typename B, typename... Es,
          size_t... Is, typename... Ks>
auto construct_nonlinear_adds(std::tuple<S...> const& systems,
                              OpExpression<B> const& bop, OpAdd<Es...> const& e,
                              double const* h, const len_type* dims,
                              std::index_sequence<Is...> seq,
                              sp_group::key_list<Ks...>) {
  using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;
  B const& b = *static_cast<B const*>(&bop);

  auto others = (OpVoid{} + ... +
                 construct_other<Z0, D>(systems, b, expr::get<Is>(e), h, dims));
  auto pointwise = (OpVoid{} + ... + pointwise_part<Z0, D>(
                                         systems, expr::get<Is>(e), h, dims));
  auto groups = (OpVoid{} + ... +
                 construct_group<Ks, Z0, D>(systems, b, e, h, dims, seq));

  if constexpr (std::is_same<decltype(pointwise), OpVoid>::value) {
    return others + groups;
  } else {
    auto B_term = solver_sp::get_B_term<vt>(b);
    return others + B_term * expr::make_fourier_map(pointwise) + groups;
  }
}
}  // namespace

template <size_t Z0, size_t D, typename... S, typename B, typename... Es>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop, OpAdd<Es...> const& e,
                         double const* h, const len_type* dims) {
  return construct_nonlinear_adds<Z0, D>(
      systems, *static_cast<B const*>(&bop), e, h, dims,
      std::make_index_sequence<sizeof...(Es)>{},
      typename sp_group::keys_of<sp_group::key_list<>, Es...>::type{});
}

template <size_t Z0, size_t D, typename... S, typename B, typename E1,
          typename E2>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpBinaryMul<E1, E2> const& e, double const* h,
                         const len_type* dims) {
  using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;

  auto B_term = get_B_term<vt>(*static_cast<B const*>(&bop));

  // A constant factor (e.g. a basis vector) is unchanged by the transform.
  if constexpr (expr::is_coeff<E1>) {
    return B_term *
           (e.a * construct_nonlinear<Z0, D>(systems, OpIdentity{}, e.b, h,
                                             dims));
  } else if constexpr (expr::is_coeff<E2>) {
    return B_term * (construct_nonlinear<Z0, D>(systems, OpIdentity{}, e.a,
                                                h, dims) *
                     e.b);
  } else {
    // One forward transform of the product evaluated in real space.
    return B_term * expr::make_fourier_map(real_space_form<Z0, D>(systems, e, h, dims));
  }
}

//! A power of an expression: one forward transform of the power evaluated
//! in real space, as for a product.
template <size_t Z0, size_t D, typename... S, typename B, expr::exp_key_t X,
          typename V, typename E>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop, OpPow<X, V, E> const& e,
                         double const* h, const len_type* dims) {
  using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;
  auto B_term = get_B_term<vt>(*static_cast<B const*>(&bop));
  return B_term *
         expr::make_fourier_map(real_space_form<Z0, D>(systems, e, h, dims));
}

template <size_t Z0, size_t D, typename... S, typename B, typename E1,
          typename E2>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpBinaryDiv<E1, E2> const& e, double const* h,
                         const len_type* dims) {
  using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;

  auto B_term = get_B_term<vt>(*static_cast<B const*>(&bop));

  using T1 = real_space_eval_t<E1>;
  using T2 = real_space_eval_t<E2>;

  auto variables = expr::get_indexed_variable_list(e);

  auto ee = anti_convolve<sizeof...(S), T1, T2>(
      variables,
      (construct_nonlinear<Z0, D>(systems, OpIdentity{}, e.a, h, dims)),
      (construct_nonlinear<Z0, D>(systems, OpIdentity{}, expr::inverse(e.b), h,
                                  dims)));

  return B_term * ee;
}

template <size_t Z0, size_t D, typename... S, typename B, typename E2>
auto construct_nonlinear(std::tuple<S...> const& systems,
                         OpExpression<B> const& bop,
                         OpBinaryDiv<OpIdentity, E2> const& e, double const* h,
                         const len_type* dims) {
  using vt = variable_type<symphas::lib::type_at_index<Z0, S...>>;

  auto B_term = get_B_term<vt>(*static_cast<B const*>(&bop));
  auto ft = expr::transform::to_ft<D>(e, h, dims);
  return B_term * ft;
}

namespace {
//! Fourier space kernel of PoissonSolver, mu0 / |k|^2 (0 at k = 0).
/*!
 * PoissonSolver(E) returns A solving lap(A) = -mu0 E, so in Fourier space
 * A(k) = mu0 E(k) / |k|^2. The k = 0 mode is set to 0, as in the real space
 * solver. The wavenumbers are those of the transformed system data.
 */
template <size_t D>
auto poisson_kernel(double const* h, const len_type* dims) {
  WaveVectorData<2, D> k2(dims, h);  // evaluates to -|k|^2
  Grid<scalar_t, D> kernel(dims);
  for (iter_type n = 0; n < kernel.len; ++n) {
    scalar_t v = k2[n];
    kernel[n] = (n == 0 || v == 0)
                    ? 0.
                    : -symphas::internal::POISSON_MU0 / v;
  }
  return expr::make_term(NamedData(std::move(kernel), "mu0/k^2"));
}
}  // namespace

//! Solve the Poisson equation in Fourier space.
/*!
 * Instead of evaluating the source in real space and solving with an FFT
 * (which needs real space derivatives), the transform of the source is
 * constructed like any other nonlinear term and multiplied by mu0 / |k|^2.
 */
template <size_t Z0, size_t D, typename... S, typename B, typename E>
auto construct_nonlinear(
    std::tuple<S...> const& systems, OpExpression<B> const& bop,
    OpMap<symphas::internal::PoissonSolver, OpIdentity, E> const& e,
    double const* h, const len_type* dims) {
  return poisson_kernel<D>(h, dims) *
         construct_nonlinear<Z0, D>(systems, *static_cast<B const*>(&bop),
                                    expr::get_enclosed_expression(e), h, dims);
}

// **************************************************************************************

//! Construct a linear operator for the semi-implicit spectral method.
/*!
 * Expressions will be parsed in order to construct a linear operator for the
 * spectral solver. Expressions which are not used in the construction of the
 * spectral operator are pruned away from the terms used in the operator.
 * Therefore, the result is a pair, the first element containing the terms
 * acceptable for finalizing the operator, and the remaining terms are not used
 * in the operator. Typically, these are considered 'nonlinear' in the subject
 * variable.
 *
 * The default behaviour for expressions is to return them as the second
 * element in the pair, as it is not included in the operator.
 *
 * This is only compatible with linear expressions; applying this function to a
 * nonlinear will return an inconsistent expression with missing primary
 * variable in linear parts and transformed to Fourier space in nonlinear parts.
 *
 * Moreover, it is required to be a function of only a single variable;
 * otherwise the variables will be combined and the result will not be
 * consistent with a spectral linear operator.
 *
 * This algorithm follows an almost identical process as to_ft with different
 * termination and some other behaviour (such as multiplication), which includes
 * sorting terms which do not belong in the linear part.
 */
template <size_t Z, typename E>
auto get_l_op(OpExpression<E> const& e, double const* h);
template <size_t Z, typename... Es>
auto get_l_op(OpAdd<Es...> const& e, double const* h);
template <size_t Z, typename A1, typename A2, typename E>
auto get_l_op(OpChain<A1, A2, E> const& e, double const* h);
template <size_t Z, typename A1, typename A2, typename E>
auto get_l_op(OpCombination<A1, A2, E> const& e, double const* h);
template <size_t Z, typename Dd, typename V, typename E, typename Sp>
auto get_l_op(OpDerivative<Dd, V, E, Sp> const& e, double const* h);
// template <size_t Z, typename V, typename E1, typename E2>
// auto get_l_op(OpConvolution<V, E1, E2> const& e, double const* h);
//  template <size_t Z, typename V, size_t D, typename E>
//  auto get_l_op(OpConvolution<V, GaussianSmoothing<D>, E> const& e,
//                double const* h);
//  template <size_t Z, typename V, size_t D, typename G>
//  auto get_l_op(
//      OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>> const& e,
//      double const* h);

template <size_t Z, typename E>
constexpr bool l_op_compatible = expr::vars<E>::template only_id<Z>();

template <size_t Z, typename E>
auto get_l_op(OpExpression<E> const& e, double const*) {
  return std::make_pair(OpVoid{}, *static_cast<E const*>(&e));
}

template <size_t R0, typename G, typename T,
          size_t R = expr::eval_type<OpTerm<OpIdentity, G>>::rank>
auto column_vector_from_variable(OpTerm<T, G> const& e) {
  if constexpr (R == 0) {
    return OpIdentity{};
  } else {
    auto c = expr::make_column_vector<R0, R>();
    if constexpr (R0 + 1 < R) {
      return c + column_vector_from_variable<R0 + 1>(e);
    } else {
      return c;
    }
  }
}

template <size_t Z, typename T, typename G,
          typename std::enable_if_t<l_op_compatible<Z, G>, int> = 0>
auto get_l_op(OpTerm<T, G> const& e, double const*) {
  return std::make_pair(expr::coeff(e) * column_vector_from_variable<0>(e),
                        OpVoid{});
}

template <size_t Z, typename T, Axis ax, typename G,
          typename std::enable_if_t<l_op_compatible<Z, G>, int> = 0>
auto get_l_op(OpTerm<T, VectorComponent<ax, G>> const& e, double const*) {
  return std::make_pair(expr::coeff(e), OpVoid{});
}

namespace {
template <size_t Z, typename... Es, size_t... Is>
auto get_l_op_adds(OpAdd<Es...> const& e, double const* h,
                   std::index_sequence<Is...>) {
  auto&& a = std::make_tuple(get_l_op<Z>(expr::get<Is>(e), h)...);
  return std::make_pair(expr::add_all(std::get<Is>(a).first...),
                        expr::add_all(std::get<Is>(a).second...));
}
}  // namespace

template <size_t Z, typename... Es>
auto get_l_op(OpAdd<Es...> const& e, double const* h) {
  return get_l_op_adds<Z>(e, h, std::make_index_sequence<sizeof...(Es)>{});
}

template <size_t Z, typename A1, typename A2, typename E>
auto get_l_op(OpChain<A1, A2, E> const& e, double const* h) {
  constexpr size_t D = expr::grid_dim<E>::dimension;
  len_type dims[D];
  expr::fill_data_dimensions(e, dims);

  auto op_ft = expr::transform::to_ft<expr::grid_dim<E>::dimension>(
      e.combination, h, dims);
  auto&& [expr_l, non_op] = get_l_op<Z>(expr::get_enclosed_expression(e), h);

  return std::make_pair(op_ft * expr_l, e.combination(non_op));
}

template <size_t Z, typename A1, typename A2, typename E>
auto get_l_op(OpCombination<A1, A2, E> const& e, double const* h) {
  constexpr size_t D = expr::grid_dim<E>::dimension;
  len_type dims[D];
  expr::fill_data_dimensions(e, dims);

  auto op_ft = expr::transform::to_ft<expr::grid_dim<E>::dimension>(
      e.combination, h, dims);
  auto&& [expr_l, non_op] = get_l_op<Z>(expr::get_enclosed_expression(e), h);

  return std::make_pair(op_ft * expr_l, e.combination(non_op));
}

template <size_t Z, typename Dd, typename V, typename E, typename Sp>
auto get_l_op(OpDerivative<Dd, V, E, Sp> const& e, double const* h) {
  constexpr size_t D = expr::grid_dim<E>::dimension;
  len_type dims[D];
  expr::fill_data_dimensions(e, dims);

  constexpr size_t order = OpDerivative<Dd, V, E, Sp>::order;
  // constexpr Axis axis = OpDerivative<Dd, V, E, Sp>::axis;

  static_assert(order > 0);

  auto [op, en] = expr::split::separate_operator(e);
  auto&& [expr_l, non_op] = get_l_op<Z>(en, h);

  return std::make_pair(expr::transform::to_ft<D>(op, h, dims) * expr_l,
                        expr::make_derivative<Dd>(non_op, e.solver));
}
//
// template <size_t Z, typename V, size_t D, typename E>
// auto get_l_op(OpConvolution<V, GaussianSmoothing<D>, E> const& e,
//              double const* h) {
//  len_type dims[D];
//  expr::fill_data_dimensions(e, dims);
//  auto op_ft = expr::transform::to_ft<D>(e.a, h, dims);
//  auto&& [expr_l, non_op] = get_l_op<Z>(e.b, h);
//  return std::make_pair(op_ft * expr_l, expr::make_convolution(e.a, non_op));
//}

// template <size_t Z, typename V, size_t D, typename G>
// auto get_l_op(
//     OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>> const& e,
//     double const* h) {
//   len_type dims[D];
//   expr::fill_data_dimensions(e, dims);
//   auto op_ft = expr::transform::to_ft<D>(e.smoother, h, dims);
//   auto&& [expr_l, non_op] = get_l_op<Z>(expr::get_enclosed_expression(e), h);
//   return std::make_pair(op_ft * expr_l,
//                         expr::make_convolution(e.smoother, non_op));
// }
//
// template <size_t Z, typename V, typename E1, typename E2>
// auto get_l_op(OpConvolution<V, E1, E2> const& e, double const* h) {
//   if constexpr (expr::is_linear<E1>::value) {
//     constexpr size_t D = expr::grid_dim<E1>::dimension;
//     len_type dims[D];
//     expr::fill_data_dimensions(e, dims);
//     auto op_ft = expr::transform::to_ft<D>(e.a, h, dims);
//     auto&& [expr_l, non_op] = get_l_op<Z>(e.b, h);
//     return std::make_pair(op_ft * expr_l, expr::make_convolution(e.a,
//     non_op));
//   } else {
//     return std::make_pair(OpVoid{}, e);
//   }
// }

//! Construct a linear operator for the semi-implicit spectral method.
/*!
 * Expressions will be parsed in order to construct a linear operator for the
 * spectral solver. Expressions which are not used in the construction of the
 * spectral operator are pruned away from the terms used in the operator.
 * Therefore, the result is a pair, the first element containing the terms
 * acceptable for finalizing the operator, and the remaining terms are not used
 * in the operator. Typically, these are considered 'nonlinear' in the subject
 * variable.
 *
 * The default behaviour for expressions is to return them as the second
 * element in the pair, as it is not included in the operator.
 *
 * This is only compatible with linear expressions; applying this function to a
 * nonlinear will return an inconsistent expression with missing primary
 * variable in linear parts and transformed to Fourier space in nonlinear parts.
 *
 * Moreover, it is required to be a function of only a single variable;
 * otherwise the variables will be combined and the result will not be
 * consistent with a spectral linear operator.
 *
 * This algorithm follows an almost identical process as to_ft with different
 * termination and some other behaviour (such as multiplication), which includes
 * sorting terms which do not belong in the linear part.
 */
template <typename E>
auto get_l_op(OpExpression<E> const& e, double const* h);
template <typename... Es>
auto get_l_op(OpAdd<Es...> const& e, double const* h);
template <typename A1, typename A2, typename E>
auto get_l_op(OpChain<A1, A2, E> const& e, double const* h);
template <typename A1, typename A2, typename E>
auto get_l_op(OpCombination<A1, A2, E> const& e, double const* h);
template <typename Dd, typename V, typename E, typename Sp>
auto get_l_op(OpDerivative<Dd, V, E, Sp> const& e, double const* h);
// template <typename V, typename E1, typename E2>
// auto get_l_op(OpConvolution<V, E1, E2> const& e, double const* h);
//  template <typename V, size_t D, typename E>
//   auto get_l_op(OpConvolution<V, GaussianSmoothing<D>, E> const& e,
//                 double const* h);
//   template <typename V, size_t D, typename G>
//   auto get_l_op(
//       OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>> const& e,
//       double const* h);

template <typename E>
auto get_l_op(OpExpression<E> const& e, double const*) {
  return std::make_pair(OpVoid{}, *static_cast<E const*>(&e));
}

template <typename T, typename G>
auto get_l_op(OpTerm<T, G> const& e, double const*) {
  return std::make_pair(expr::coeff(e) * column_vector_from_variable<0>(e),
                        OpVoid{});
}

template <typename T, Axis ax, typename G>
auto get_l_op(OpTerm<T, VectorComponent<ax, G>> const& e, double const*) {
  return std::make_pair(expr::coeff(e), OpVoid{});
}

namespace {
template <typename... Es, size_t... Is>
auto get_l_op_adds(OpAdd<Es...> const& e, double const* h,
                   std::index_sequence<Is...>) {
  auto&& a = std::make_tuple(get_l_op(expr::get<Is>(e), h)...);
  return std::make_pair(expr::add_all(std::get<Is>(a).first...),
                        expr::add_all(std::get<Is>(a).second...));
}
}  // namespace

template <typename... Es>
auto get_l_op(OpAdd<Es...> const& e, double const* h) {
  return get_l_op_adds(e, h, std::make_index_sequence<sizeof...(Es)>{});
}

template <typename A1, typename A2, typename E>
auto get_l_op(OpChain<A1, A2, E> const& e, double const* h) {
  constexpr size_t D = expr::grid_dim<E>::dimension;
  len_type dims[D];
  expr::fill_data_dimensions(e, dims);

  auto op_ft = expr::transform::to_ft<expr::grid_dim<E>::dimension>(
      e.combination, h, dims);
  auto&& [expr_l, non_op] = get_l_op(expr::get_enclosed_expression(e), h);

  return std::make_pair(op_ft * expr_l, e.combination(non_op));
}

template <typename A1, typename A2, typename E>
auto get_l_op(OpCombination<A1, A2, E> const& e, double const* h) {
  constexpr size_t D = expr::grid_dim<E>::dimension;
  len_type dims[D];
  expr::fill_data_dimensions(e, dims);

  auto op_ft = expr::transform::to_ft<expr::grid_dim<E>::dimension>(
      e.combination, h, dims);
  auto&& [expr_l, non_op] = get_l_op(expr::get_enclosed_expression(e), h);

  return std::make_pair(op_ft * expr_l, e.combination(non_op));
}

template <typename Dd, typename V, typename E, typename Sp>
auto get_l_op(OpDerivative<Dd, V, E, Sp> const& e, double const* h) {
  constexpr size_t D = expr::grid_dim<E>::dimension;
  len_type dims[D];
  expr::fill_data_dimensions(e, dims);

  constexpr size_t order = OpDerivative<Dd, V, E, Sp>::order;
  // constexpr Axis axis = OpDerivative<Dd, V, E, Sp>::axis;

  static_assert(order > 0);

  auto [op, en] = expr::split::separate_operator(e);
  auto&& [expr_l, non_op] = get_l_op(en, h);

  return std::make_pair(expr::transform::to_ft<D>(op, h, dims) * expr_l,
                        expr::make_derivative<Dd>(non_op, e.solver));
}

// template <typename V, size_t D, typename E>
// auto get_l_op(OpConvolution<V, GaussianSmoothing<D>, E> const& e,
//               double const* h) {
//   len_type dims[D];
//   expr::fill_data_dimensions(e, dims);
//   auto op_ft = expr::transform::to_ft<D>(e.a, h, dims);
//   auto&& [expr_l, non_op] = get_l_op(e.b, h);
//   return std::make_pair(op_ft * expr_l, expr::make_convolution(e.a, non_op));
// }
//
// template <typename V, size_t D, typename G>
// auto get_l_op(
//     OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>> const& e,
//     double const* h) {
//   len_type dims[D];
//   expr::fill_data_dimensions(e, dims);
//   auto op_ft = expr::transform::to_ft<D>(e.smoother, h, dims);
//   auto&& [expr_l, non_op] = get_l_op(expr::get_enclosed_expression(e), h);
//   return std::make_pair(op_ft * expr_l,
//                         expr::make_convolution(e.smoother, non_op));
// }
//
// template <typename V, typename E1, typename E2>
// auto get_l_op(OpConvolution<V, E1, E2> const& e, double const* h) {
//  if constexpr (expr::is_linear<E1>::value) {
//    constexpr size_t D = expr::grid_dim<E1>::dimension;
//    len_type dims[D];
//    expr::fill_data_dimensions(e, dims);
//    auto op_ft = expr::transform::to_ft<D>(e.a, h, dims);
//    auto&& [expr_l, non_op] = get_l_op(e.b, h);
//    return std::make_pair(op_ft * expr_l, expr::make_convolution(e.a,
//    non_op));
//  } else {
//    return std::make_pair(OpVoid{}, e);
//  }
//}

// **************************************************************************************

}  // namespace solver_sp

template <typename E>
struct SpectralData {
  E scheme;
  //! Inverse transforms shared by the terms of the scheme.
  std::shared_ptr<solver_sp::InverseMapCache> cache;

  SpectralData(E const& scheme,
               std::shared_ptr<solver_sp::InverseMapCache> cache = nullptr)
      : scheme{scheme}, cache{cache} {}

  void update() {
    if (cache) cache->update();
    expr::prune::update(scheme);
  }
};
