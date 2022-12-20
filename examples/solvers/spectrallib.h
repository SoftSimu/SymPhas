
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

#include "spslib.h"
#include "spslibfftw.h"

#include "definitions.h"
#include "indexseqhelpers.h"
#include "gridpair.h"

#include "expressionexponentials.h"
#include "expressiontransforms.h"
#include "expressionlib.h"
#include "expressionrules.h"

namespace solver_sp
{

	namespace
	{

		// **************************************************************************************


		// swap for the primary variables (in model.system<N>())
		template<Axis ax, size_t Z, size_t Q0, size_t N, typename T, typename E>
		auto swap_var_apply_axes(MultiBlock<N, T> const& grid, OpExpression<E> const& expr, const len_type* dims)
		{
			auto data = expr::as_grid_data<N>(grid, dims);
#ifdef PRINTABLE_EQUATIONS
			auto v = expr::get_variable<Q0>(*static_cast<E const*>(&expr));
			auto name = expr::get_fourier_name(expr::get_op_name(v));
			auto replace = expr::as_component<ax>(expr::as_variable<Z>(NamedData(std::move(data), name)));
#else
			auto replace = expr::as_component<ax>(expr::as_variable<Z>(std::move(data)));
#endif
			auto e = expr::transform::swap_grid<ax, Q0>(*static_cast<const E*>(&expr), replace);
			constexpr size_t I = (ax == Axis::X) ? 0 : (ax == Axis::Y) ? 1 : (ax == Axis::Z) ? 2 : -1;

			if constexpr (I + 1 < N)
			{
				constexpr Axis bx = (I == 0) ? Axis::Y : Axis::Z;
				return swap_var_apply_axes<bx, Z, Q0>(grid, e, dims);
			}
			else
			{
				return e;
			}
		}

		// swap for the primary variables (in model.system<N>())
		template<size_t Z, size_t Q0, size_t N, typename T, typename E>
		auto swap_var_apply(MultiBlock<N, T> const& grid, OpExpression<E> const& expr, const len_type* dims)
		{
			auto data = expr::as_grid_data<N>(grid, dims);
#ifdef PRINTABLE_EQUATIONS
			auto name = expr::get_fourier_name(expr::get_op_name(expr::get_variable<Q0>(*static_cast<E const*>(&expr))));
			auto replace = expr::as_variable<Z>(NamedData(std::move(data), name));
#else
			auto replace = expr::as_variable<Z>(std::move(data));
#endif
			auto e = swap_var_apply_axes<Axis::X, Z, Q0>(grid, *static_cast<const E*>(&expr), dims);
			return expr::transform::swap_grid<Q0>(e, replace);
		}

		// swap for the primary variables (in model.system<N>())
		template<size_t Z, size_t Q0, typename T, typename E>
		auto swap_var_apply(T* const& grid, OpExpression<E> const& expr, const len_type* dims)
		{
			constexpr size_t D = expr::grid_dim<E>::value;
			auto data = expr::as_grid_data<D>(grid, dims);

#ifdef PRINTABLE_EQUATIONS
			auto name = expr::get_fourier_name(expr::get_op_name(expr::get_variable<Q0>(*static_cast<E const*>(&expr))));
			auto replace = expr::as_variable<Z>(NamedData(std::move(data), name));
#else
			auto replace = expr::as_variable<Z>(std::move(data));
#endif
			return expr::transform::swap_grid<Q0>(*static_cast<const E*>(&expr), replace);
		}

		template<typename... S, typename E>
		auto swap_var_apply(std::tuple<S...> const&, E&& expr, std::index_sequence<>)
		{
			return std::forward<E>(expr);
		}

		template<typename T, size_t D>
		auto set_system_dimensions(System<T, D> const& s, len_type(&dims)[D])
		{
			for (iter_type i = 0; i < D; ++i)
			{
				dims[i] = s.dims[i];
			}
		}

		// swap for the primary variables (in model.system<N>())
		template<typename... S, typename E, size_t Q0, size_t... Qs>
		auto swap_var_apply(std::tuple<S...> const& systems, OpExpression<E> const& expr, std::index_sequence<Q0, Qs...>)
		{
			constexpr size_t D = expr::grid_dim<E>::value;
			len_type dims[D];
			set_system_dimensions(std::get<Q0>(systems), dims);

 			auto e = swap_var_apply<sizeof...(S) + Q0, Q0>(
				std::get<Q0>(systems).frame_t, 
				*static_cast<E const*>(&expr),
				dims);

			if constexpr (sizeof...(Qs) == 0)
			{
				return e;
			}
			else
			{
				return swap_var_apply(systems, e,
					std::index_sequence<Qs...>{});
			}
		}


	}

	template<typename E>
	auto drop_hcts(OpExpression<E> const& e)
	{
		return (*static_cast<E const*>(&e));
	}

	template<typename E>
	auto drop_hcts(OpMap<symphas::internal::HCTS, OpIdentity, E> const& e)
	{
		return expr::get_enclosed_expression(e);
	}

	template<typename T, typename E>
	auto sthc_apply_on_scalar(OpExpression<E> const& e)
	{
		using rt = typename symphas::internal::real_space_type<T>::type;
		if constexpr (std::is_same<rt, scalar_t>::value)
		{
			return expr::sthc(*static_cast<E const*>(&e));
		}
		else
		{
			return (*static_cast<E const*>(&e));
		}
	}

	template<typename T, typename E>
	auto hcts_apply_on_scalar(OpExpression<E> const& e)
	{
		using rt = typename symphas::internal::real_space_type<T>::type;
		if constexpr (std::is_same<rt, scalar_t>::value)
		{
			return expr::hcts(*static_cast<E const*>(&e));
		}
		else
		{
			return (*static_cast<E const*>(&e));
		}
	}



	//! Given the linear expression, returns the \f$A\f$ operator.
	template<size_t D, typename E, size_t R = expr::eval_type<E>::rank, size_t... Rs>
	auto form_A_op(OpExpression<E> const& l_op, double dt, len_type const* dims, std::index_sequence<Rs...>)
	{
		auto A_expr =
			((expr::make_column_vector<Rs, R>() * expr::exp(
				expr::make_literal(dt) * (expr::make_row_vector<Rs, R>() * (*static_cast<E const*>(&l_op)))))
				+ ...);
		return A_expr;
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<size_t D, typename E, size_t R = expr::eval_type<E>::rank, typename std::enable_if_t<(R == 0), int> = 0>
	auto form_A_op(OpExpression<E> const& l_op, double dt, len_type const* dims)
	{
		auto A_expr = expr::exp(expr::make_literal(dt) * (*static_cast<E const*>(&l_op)));
		return A_expr;
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<size_t D, typename E, size_t R = expr::eval_type<E>::rank, typename std::enable_if_t<(R > 0), int> = 0>
		auto form_A_op(OpExpression<E> const& l_op, double dt, len_type const* dims)
	{
		return form_A_op<D>(*static_cast<E const*>(&l_op), dt, dims, std::make_index_sequence<R>{});
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<size_t D>
	auto form_A_op(OpVoid, double dt, len_type const* dims)
	{
		return OpIdentity{};
	}


	//! Given the linear expression, returns the \f$A\f$ operator.
	template<size_t D, typename E, size_t R = expr::eval_type<E>::rank, size_t... Rs>
	auto form_B_op(OpExpression<E> const& l_op, double dt, len_type const* dims, std::index_sequence<Rs...>)
	{
		auto A_expr = form_A_op<D>(*static_cast<E const*>(&l_op), dt, dims);
		auto B_expr = (
			(expr::make_column_vector<Rs, R>() * ((expr::make_row_vector<Rs, R>() * A_expr - expr::symbols::one) / (expr::make_row_vector<Rs, R>() * (*static_cast<E const*>(&l_op)))))
			+ ...);
		return B_expr;
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<size_t D, typename E, size_t R = expr::eval_type<E>::rank, typename std::enable_if_t<(R == 0), int> = 0>
	auto form_B_op(OpExpression<E> const& l_op, double dt, len_type const* dims)
	{
		auto A_expr = form_A_op<D>(*static_cast<E const*>(&l_op), dt, dims);
		auto B_expr = (A_expr - expr::symbols::one) / *static_cast<E const*>(&l_op);
		return B_expr;
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<size_t D, typename E, size_t R = expr::eval_type<E>::rank, typename std::enable_if_t<(R > 0), int> = 0>
	auto form_B_op(OpExpression<E> const& l_op, double dt, len_type const* dims)
	{
		return form_B_op<D>(*static_cast<E const*>(&l_op), dt, dims, std::make_index_sequence<R>{});
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<size_t D>
	auto form_B_op(OpVoid, double dt, len_type const* dims)
	{
		return expr::make_literal(dt);
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<typename T, typename A, size_t R = expr::eval_type<A>::rank, size_t... Rs>
	auto get_A_term(OpExpression<A> const& A_expression, std::index_sequence<Rs...>);

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<typename T, typename A, size_t R = expr::eval_type<A>::rank, typename std::enable_if_t<(R == 0), int> = 0>
	auto get_A_term(OpExpression<A> const& A_expression)
	{
		auto a = expr::transform::to_grid(sthc_apply_on_scalar<T>(*static_cast<A const*>(&A_expression)));
		auto A_term = expr::make_term(NamedData(std::move(a), *static_cast<A const*>(&A_expression)));
		return A_term;
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<typename T, typename A, size_t R = expr::eval_type<A>::rank, typename std::enable_if_t<(R > 0), int> = 0>
	auto get_A_term(OpExpression<A> const& A_expression)
	{
		return get_A_term<T>(*static_cast<A const*>(&A_expression), std::make_index_sequence<R>{});
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<typename T>
	inline auto get_A_term(OpIdentity)
	{
		return OpIdentity{};
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<typename T, typename A, size_t R, size_t... Rs>
	auto get_A_term(OpExpression<A> const& A_expression, std::index_sequence<Rs...>)
	{
		return ((expr::make_tensor<Rs, Rs, R, R>() * get_A_term<T>(
			expr::make_row_vector<Rs, R>() * sthc_apply_on_scalar<T>(*static_cast<A const*>(&A_expression))))
			+ ...);
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<typename T, typename T0>
	auto get_B_term(OpLiteral<T0> const& B_expression)
	{
		return B_expression;
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<typename T, typename B, size_t R = expr::eval_type<B>::rank, size_t... Rs>
	auto get_B_term(OpExpression<B> const& B_expression, std::index_sequence<Rs...>);

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<typename T, typename B, size_t R = expr::eval_type<B>::rank, typename std::enable_if_t<(R == 0), int> = 0>
	auto get_B_term(OpExpression<B> const& B_expression)
	{
		return expr::make_term(NamedData(expr::transform::to_grid(sthc_apply_on_scalar<T>(*static_cast<B const*>(&B_expression))), (*static_cast<B const*>(&B_expression))));
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<typename T, typename B, size_t R = expr::eval_type<B>::rank, typename std::enable_if_t<(R > 0), int> = 0>
	auto get_B_term(OpExpression<B> const& B_expression)
	{
		return get_B_term<T>(*static_cast<B const*>(&B_expression), std::make_index_sequence<R>{});
	}
	
	//! Given the linear expression, returns the \f$A\f$ operator.
	template<typename T>
	inline auto get_B_term(OpIdentity)
	{
		return OpIdentity{};
	}

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<typename T, typename B, size_t R, size_t... Rs>
	auto get_B_term(OpExpression<B> const& B_expression, std::index_sequence<Rs...>)
	{
		return ((expr::make_tensor<Rs, Rs, R, R>() * get_B_term<T>(
			expr::make_row_vector<Rs, R>() * sthc_apply_on_scalar<T>(*static_cast<B const*>(&B_expression))))
			+ ...);
	}

	//! Swap real for Fourier variables in the given expression.
	/*!
	 * Take a list of all the systems appearing in the phase field problem,
	 * then for each index associated variable in the expression (which is the
	 * real space order parameter), make the equivalent expression instead with
	 *
	 */
	template<typename... S, typename E>
	decltype(auto) swap_var(std::tuple<S...> const& systems, OpExpression<E> const& expr)
	{
		return swap_var_apply(systems, *static_cast<const E*>(&expr), expr::vars<E>::get_ids());
	}

	//! Apply the convolution theorem to construct the nonlinear scheme.
	/*!
	 * Apply the convolution theorem to construct the nonlinear scheme.
	 */
	template<size_t Zn, typename T1, typename T2, typename... S, typename E1, typename E2>
	auto anti_convolve(std::tuple<S...> const&, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		auto ae = expr::make_inv_fourier_map<T1>(hcts_apply_on_scalar<T1>(*static_cast<E1 const*>(&a)));
		auto be = expr::make_inv_fourier_map<T2>(hcts_apply_on_scalar<T2>(*static_cast<E2 const*>(&b)));
		return expr::make_fourier_map(ae * be);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename S2, typename V2, typename E1, typename E2>
	auto anti_convolve(std::tuple<S...> const&, OpExpression<E1> const& a, OpFourierTransform<S2, V2, E2> const& b)
	{
		auto e = expr::get_enclosed_expression(b);
		return expr::coeff(b) * expr::make_fourier_map(expr::make_inv_fourier_map<T1>(hcts_apply_on_scalar<T1>(*static_cast<E1 const*>(&a))) * e);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename S1, typename V1, typename E1, typename E2>
	auto anti_convolve(std::tuple<S...> const&, OpFourierTransform<S1, V1, E1> const& a, OpExpression<E2> const& b)
	{
		auto e = expr::get_enclosed_expression(a);
		return expr::coeff(a) * expr::make_fourier_map(e * expr::make_inv_fourier_map<T2>(hcts_apply_on_scalar<T2>(*static_cast<E2 const*>(&b))));
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename S1, typename V1, typename E1, typename S2, typename V2, typename E2>
	auto anti_convolve(std::tuple<S...> const&, OpFourierTransform<S1, V1, E1> const& a, OpFourierTransform<S2, V2, E2> const& b)
	{
		auto ae = expr::get_enclosed_expression(a);
		auto be = expr::get_enclosed_expression(a);
		return (expr::coeff(a) * expr::coeff(b)) * expr::make_fourier_map(ae * be);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename T, size_t Z, typename G, typename E2, size_t D = expr::grid_dim<G>::value>
	auto anti_convolve(std::tuple<S...> const& variables, OpTerm<T, Variable<Z, G>> const& a, OpExpression<E2> const& b)
	{
		auto ae = expr::transform::swap_grid<Z>(a, std::get<Z - Zn>(variables));
		auto be = expr::make_inv_fourier_map<T2>(hcts_apply_on_scalar<T2>(*static_cast<E2 const*>(&b)));
		return expr::make_fourier_map(ae * be);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename T, size_t Z, typename G, typename S2, typename V2, typename E2, size_t D = expr::grid_dim<G>::value>
	auto anti_convolve(std::tuple<S...> const& variables, OpTerm<T, Variable<Z, G>> const& a, OpFourierTransform<S2, V2, E2> const& b)
	{
		auto ae = expr::transform::swap_grid<Z>(a, std::get<Z - Zn>(variables));
		auto be = expr::get_enclosed_expression(b);
		return expr::coeff(b) * expr::make_fourier_map(ae * be);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename E1, typename T, size_t Z, typename G, size_t D = expr::grid_dim<G>::value>
	auto anti_convolve(std::tuple<S...> const& variables, OpExpression<E1> const& a, OpTerm<T, Variable<Z, G>> const& b)
	{
		auto ae = expr::make_inv_fourier_map<T1>(hcts_apply_on_scalar<T1>(*static_cast<E1 const*>(&a)));
		auto be = expr::transform::swap_grid<Z>(b, std::get<Z - Zn>(variables));
		return expr::make_fourier_map(ae * be);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename S1, typename V1, typename E1, typename T, size_t Z, typename G, size_t D = expr::grid_dim<G>::value>
	auto anti_convolve(std::tuple<S...> const& variables, OpFourierTransform<S1, V1, E1> const& a, OpTerm<T, Variable<Z, G>> const& b)
	{
		auto ae = expr::get_enclosed_expression(a);
		auto be = expr::transform::swap_grid<Z>(b, std::get<Z - Zn>(variables));
		return expr::coeff(a) * expr::make_fourier_map(ae * be);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename T, size_t Z, Axis ax, typename G, typename E2, size_t D = expr::grid_dim<G>::value>
	auto anti_convolve(std::tuple<S...> const& variables, OpTerm<T, Variable<Z, VectorComponent<ax, G>>> const& a, OpExpression<E2> const& b)
	{
		auto ae = expr::transform::swap_grid<ax, Z>(a, expr::as_component<ax>(std::get<Z - Zn>(variables)));
		auto be = expr::make_inv_fourier_map<T2>(hcts_apply_on_scalar<T2>(*static_cast<E2 const*>(&b)));
		return expr::make_fourier_map(ae * be);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename T, size_t Z, Axis ax, typename G, typename S2, typename V2, typename E2, size_t D = expr::grid_dim<G>::value>
	auto anti_convolve(std::tuple<S...> const& variables, OpTerm<T, Variable<Z, VectorComponent<ax, G>>> const& a, OpFourierTransform<S2, V2, E2> const& b)
	{
		auto ae = expr::transform::swap_grid<ax, Z>(a, expr::as_component<ax>(std::get<Z - Zn>(variables)));
		auto be = expr::get_enclosed_expression(b);
		return expr::coeff(b) * expr::make_fourier_map(ae * be);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename E1, typename T, size_t Z, Axis ax, typename G, size_t D = expr::grid_dim<G>::value>
	auto anti_convolve(std::tuple<S...> const& variables, OpExpression<E1> const& a, OpTerm<T, Variable<Z, VectorComponent<ax, G>>> const& b)
	{
		auto ae = expr::make_inv_fourier_map<T1>(hcts_apply_on_scalar<T1>(*static_cast<E1 const*>(&a)));
		auto be = expr::transform::swap_grid<ax, Z>(b, expr::as_component<ax>(std::get<Z - Zn>(variables)));
		return expr::make_fourier_map(ae * be);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename S1, typename V1, typename E1, typename T, size_t Z, Axis ax, typename G, size_t D = expr::grid_dim<G>::value>
	auto anti_convolve(std::tuple<S...> const& variables, OpFourierTransform<S1, V1, E1> const& a, OpTerm<T, Variable<Z, VectorComponent<ax, G>>> const& b)
	{
		auto ae = expr::get_enclosed_expression(a);
		auto be = expr::transform::swap_grid<ax, Z>(b, expr::as_component<ax>(std::get<Z - Zn>(variables)));
		return expr::coeff(a) * expr::make_fourier_map(ae * be);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename... As, typename E1, size_t... Is>
	auto anti_convolve(std::tuple<S...> const& variables, E1 const& a, OpAdd<As...> const& b, std::index_sequence<Is...>)
	{
		return (anti_convolve<Zn, T1, T2>(variables, a, expr::get<Is>(b)) + ...);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename... As, typename E1>
	auto anti_convolve(std::tuple<S...> const& variables, E1 const& a, OpAdd<As...> const& b)
	{
		return anti_convolve<Zn, T1, T2>(variables, a, b, std::make_index_sequence<sizeof...(As)>{});
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename... As, typename E2, size_t... Is>
	auto anti_convolve(std::tuple<S...> const& variables, OpAdd<As...> const& a, E2 const& b, std::index_sequence<Is...>)
	{
		return (anti_convolve<Zn, T1, T2>(variables, expr::get<Is>(a), b) +  ...);
	}

	template<size_t Zn, typename T1, typename T2, typename... S, typename... As, typename E2>
	auto anti_convolve(std::tuple<S...> const& variables, OpAdd<As...> const& a, E2 const& b)
	{
		return anti_convolve<Zn, T1, T2>(variables, a, b, std::make_index_sequence<sizeof...(As)>{});
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

	template<size_t Z0, size_t D, typename... S, typename B, typename Dd, typename V, typename E, typename Sp>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
		OpFuncDerivative<Dd, V, E, Sp> const& e, double const* h, const len_type* dims);
	template<size_t Z0, size_t D, typename... S, typename B, size_t O, typename V, typename Sp>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
		OpOperatorDerivative<O, V, Sp> const& e, double const* h, const len_type* dims);
	template<size_t Z0, size_t D, typename... S, typename B, typename A1, typename A2, typename E>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
		OpCombination<A1, A2, E> const& e, double const* h, const len_type* dims);
	template<size_t Z0, size_t D, typename... S, typename B, typename A1, typename A2, typename E>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
		OpChain<A1, A2, E> const& e, double const* h, const len_type* dims);
	template<size_t Z0, size_t D, size_t Z, typename... S, typename B, typename T, typename G,
		typename = std::enable_if_t<(Z < sizeof...(S)), int>>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
		OpTerm<T, Variable<Z, G>> const& e, double const* h, const len_type* dims);
	//template<size_t Z0, size_t D, NoiseType nt, typename... S, typename B, typename T, typename V>
	//auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
	//	OpTerm<V, NoiseData<nt, T, D>> const& e, double const* h, const len_type* dims);
	template<size_t Z0, size_t D, typename... S, typename B, typename... Es>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
		OpAdd<Es...> const& e, double const* h, const len_type* dims);
	template<size_t Z0, size_t D, typename... S, typename B, typename E1, typename E2>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
		OpBinaryMul<E1, E2> const& e, double const* h, const len_type* dims);
	template<size_t Z0, size_t D, typename... S, typename B, typename E1, typename E2>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
		OpBinaryDiv<E1, E2> const& e, double const* h, const len_type* dims);

	template<typename S>
	//using variable_type = decltype(std::declval<S>()[0]);
	using variable_type = typename grid::value_type_of<S>::type;

	template<size_t Z0, size_t D, typename... S, typename B, typename E>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
		OpExpression<E> const& e, double const* h, const len_type* dims)
	{
		using vt = variable_type<std::tuple_element_t<Z0, std::tuple<S...>>>;

		auto B_term = get_B_term<vt>(*static_cast<B const*>(&bop));
		auto ft = drop_hcts(expr::transform::to_ft<D>(*static_cast<E const*>(&e), h, dims));
		return B_term * ft;
	}

	template<size_t Z0, size_t D, typename... S, typename B, typename Dd, typename V, typename E, typename Sp>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
		OpFuncDerivative<Dd, V, E, Sp> const& e, double const* h, const len_type* dims)
	{
		using vt = variable_type<std::tuple_element_t<Z0, std::tuple<S...>>>;

		constexpr Axis axis = OpFuncDerivative<Dd, V, E, Sp>::axis;
		constexpr size_t order = OpFuncDerivative<Dd, V, E, Sp>::order;

		auto [op, en] = expr::split::separate_operator(e);
		auto B_expression = (*static_cast<B const*>(&bop)) * expr::transform::to_ft<D>(op, h, dims);
		auto B_term = get_B_term<vt>(B_expression);
		
		auto ft = construct_nonlinear<Z0, D>(systems, OpIdentity{}, en, h, dims);
		return B_term * ft;
	}

	template<size_t Z0, size_t D, typename... S, typename B, typename V, typename E1, typename E2>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop, 
		OpFuncConvolution<V, E1, E2> const& e, double const* h, const len_type* dims)
	{
		return expr::coeff(e) * construct_nonlinear<Z0, D>(systems, *static_cast<B const*>(&bop), e.a, h, dims)
			* construct_nonlinear<Z0, D>(systems, OpIdentity{}, *static_cast<B const*>(&bop), e.b, h, dims);
	}

	template<size_t Z0, size_t D, typename... S, typename B, typename V, typename E>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop, 
		OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e, double const* h, const len_type* dims)
	{
		return expr::coeff(e) * expr::transform::to_ft<D>(e.smoother, h, dims)
			* construct_nonlinear<Z0, D>(systems, *static_cast<B const*>(&bop), expr::get_enclosed_expression(e), h, dims);
	}

	template<size_t Z0, size_t D, typename... S, typename B, size_t O, typename V, typename Sp>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop, 
		OpOperatorDerivative<O, V, Sp> const& e, double const* h, const len_type* dims)
	{
		using vt = variable_type<std::tuple_element_t<Z0, std::tuple<S...>>>;

		auto [op, en] = expr::split::separate_operator(e);
		auto B_expression = (*static_cast<B const*>(&bop)) * expr::transform::to_ft<D>(op, h, dims);
		auto B_term = get_B_term<vt>(B_expression);

		return B_term * construct_nonlinear<Z0, D>(systems, OpIdentity{}, en, h, dims);
	}

	template<size_t Z0, size_t D, typename... S, typename B, typename A1, typename A2, typename E>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop, 
		OpCombination<A1, A2, E> const& e, double const* h, const len_type* dims)
	{
		using vt = variable_type<std::tuple_element_t<Z0, std::tuple<S...>>>;

		auto [op, en] = expr::split::separate_operator(e);
		auto B_expression = (*static_cast<B const*>(&bop)) * expr::transform::to_ft<D>(op, h, dims);
		auto B_term = get_B_term<vt>(B_expression);
		
		return B_term * construct_nonlinear<Z0, D>(systems, OpIdentity{}, en, h, dims);
	}

	template<size_t Z0, size_t D, typename... S, typename B, typename A1, typename A2, typename E>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop,
		OpChain<A1, A2, E> const& e, double const* h, const len_type* dims)
	{
		using vt = variable_type<std::tuple_element_t<Z0, std::tuple<S...>>>;

		auto [op, en] = expr::split::separate_operator(e);
		auto B_expression = (*static_cast<B const*>(&bop)) * expr::transform::to_ft<D>(op, h, dims);
		auto B_term = get_B_term<vt>(B_expression);

		return B_term * construct_nonlinear<Z0, D>(systems, OpIdentity{}, en, h, dims);
	}

	template<size_t Z0, size_t D, size_t Z, typename... S, typename B, typename T, typename G, typename>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop, 
		OpTerm<T, Variable<Z, G>> const& e, double const* h, const len_type* dims)
	{
		using vt = variable_type<std::tuple_element_t<Z0, std::tuple<S...>>>;

		auto B_term = get_B_term<vt>(*static_cast<B const*>(&bop));
		return B_term * solver_sp::swap_var(systems, e);
	}

	namespace
	{
		template<size_t Z0, size_t D, typename... S, typename B, typename... Es, size_t... Is>
		auto construct_nonlinear_adds(std::tuple<S...> const& systems, OpExpression<B> const& bop, 
			OpAdd<Es...> const& e, double const* h, const len_type* dims, std::index_sequence<Is...>)
		{
			using vt = variable_type<std::tuple_element_t<Z0, std::tuple<S...>>>;
			auto B_term = solver_sp::get_B_term<vt>(*static_cast<B const*>(&bop));

			return expr::make_add(B_term * construct_nonlinear<Z0, D>(systems, OpIdentity{}, expr::get<Is>(e), h, dims)...);
		}
	}

	template<size_t Z0, size_t D, typename... S, typename B, typename... Es>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop, 
		OpAdd<Es...> const& e, double const* h, const len_type* dims)
	{
		return construct_nonlinear_adds<Z0, D>(systems, *static_cast<B const*>(&bop), e, h, dims, std::make_index_sequence<sizeof...(Es)>{});
	}



	template<size_t Z0, size_t D, typename... S, typename B, typename E1, typename E2>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop, 
		OpBinaryMul<E1, E2> const& e, double const* h, const len_type* dims)
	{
		using vt = variable_type<std::tuple_element_t<Z0, std::tuple<S...>>>;

		auto B_term = get_B_term<vt>(*static_cast<B const*>(&bop));

		using T1 = expr::eval_type_t<E1>;
		using T2 = expr::eval_type_t<E2>;

		auto variables = expr::get_indexed_variable_list(e);

		auto ee = sthc_apply_on_scalar<vt>(anti_convolve<sizeof...(S), T1, T2>(variables,
			(construct_nonlinear<Z0, D>(systems, OpIdentity{}, e.a, h, dims)),
			(construct_nonlinear<Z0, D>(systems, OpIdentity{}, e.b, h, dims))));

		return B_term * ee;
	}

	template<size_t Z0, size_t D, typename... S, typename B, typename E1, typename E2>
	auto construct_nonlinear(std::tuple<S...> const& systems, OpExpression<B> const& bop, 
		OpBinaryDiv<E1, E2> const& e, double const* h, const len_type* dims)
	{
		using vt = variable_type<std::tuple_element_t<Z0, std::tuple<S...>>>;

		auto B_term = get_B_term<vt>(*static_cast<B const*>(&bop));

		using T1 = expr::eval_type_t<E1>;
		using T2 = expr::eval_type_t<E2>;

		auto variables = expr::get_indexed_variable_list(e);

		auto ee = sthc_apply_on_scalar<vt>(anti_convolve<sizeof...(S), T1, T2>(variables,
			(construct_nonlinear<Z0, D>(systems, OpIdentity{}, e.a, h, dims)),
			(construct_nonlinear<Z0, D>(systems, OpIdentity{}, expr::inverse(e.b), h, dims))));

		return B_term * ee;
	}

	// **************************************************************************************

	//! Construct a linear operator for the semi-implicit spectral method.
	/*! 
	 * Expressions will be parsed in order to construct a linear operator for the spectral
	 * solver. Expressions which are not used in the construction of the spectral operator
	 * are pruned away from the terms used in the operator. Therefore, the result is
	 * a pair, the first element containing the terms acceptable for finalizing
	 * the operator, and the remaining terms are not used in the operator. Typically,
	 * these are considered 'nonlinear' in the subject variable.
	 * 
	 * The default behaviour for expressions is to return them as the second
	 * element in the pair, as it is not included in the operator.
	 *
	 * This is only compatible with linear expressions; applying this function to a nonlinear
	 * will return an inconsistent expression with missing primary variable in linear parts and
	 * transformed to Fourier space in nonlinear parts.
	 * 
	 * Moreover, it is required to be a function of only a single variable; otherwise the variables
	 * will be combined and the result will not be consistent with a spectral linear operator.
	 *
	 * This algorithm follows an almost identical process as to_ft with different termination
	 * and some other behaviour (such as multiplication), which includes sorting terms which
	 * do not belong in the linear part.
	 */
	template<size_t Z, typename E>
	auto get_l_op(OpExpression<E> const& e, double const* h);
	template<size_t Z, typename... Es>
	auto get_l_op(OpAdd<Es...> const& e, double const* h);
	template<size_t Z, typename A1, typename A2, typename E>
	auto get_l_op(OpChain<A1, A2, E> const& e, double const* h);
	template<size_t Z, typename A1, typename A2, typename E>
	auto get_l_op(OpCombination<A1, A2, E> const& e, double const* h);
	template<size_t Z, typename Dd, typename V, typename E, typename Sp>
	auto get_l_op(OpFuncDerivative<Dd, V, E, Sp> const& e, double const* h);
	template<size_t Z, typename V, typename E1, typename E2>
	auto get_l_op(OpFuncConvolution<V, E1, E2> const& e, double const* h);
	template<size_t Z, typename V, size_t D, typename E>
	auto get_l_op(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e, double const* h);

	template<size_t Z, typename E>
	constexpr bool l_op_compatible = expr::vars<E>::template only_id<Z>();

	template<size_t Z, typename E>
	auto get_l_op(OpExpression<E> const& e, double const*)
	{
		return std::make_pair(OpVoid{}, *static_cast<E const*>(&e));
	}

	template<size_t R0, typename G, typename T, size_t R = expr::eval_type<OpTerm<OpIdentity, G>>::rank>
	auto column_vector_from_variable(OpTerm<T, G> const& e)
	{
		if constexpr (R == 0)
		{
			return OpIdentity{};
		}
		else
		{
			auto c = expr::make_column_vector<R0, R>();
			if constexpr (R0 + 1 < R)
			{
				return c + column_vector_from_variable<R0 + 1>(e);
			}
			else
			{
				return c;
			}
		}
	}

	template<size_t Z, typename T, typename G, 
		typename std::enable_if_t<l_op_compatible<Z, G>, int> = 0>
	auto get_l_op(OpTerm<T, G> const& e, double const*)
	{
		return std::make_pair(expr::coeff(e) * column_vector_from_variable<0>(e), OpVoid{});
	}

	//template<typename G, typename T, size_t R = expr::eval_type<OpTerm<T, VectorComponent<Axis::X, G>>>::rank>
	//auto column_vector_from_axis(OpTerm<T, VectorComponent<Axis::X, G>>)
	//{
	//	return expr::make_row_vector<0, R>();
	//}

	//template<typename G, typename T, size_t R = expr::eval_type<OpTerm<T, VectorComponent<Axis::X, G>>>::rank>
	//auto column_vector_from_axis(OpTerm<T, VectorComponent<Axis::Y, G>>)
	//{
	//	return expr::make_row_vector<1, R>();
	//}

	//template<typename G, typename T, size_t R = expr::eval_type<OpTerm<T, VectorComponent<Axis::X, G>>>::rank>
	//auto column_vector_from_axis(OpTerm<T, VectorComponent<Axis::Z, G>>)
	//{
	//	return expr::make_row_vector<2, R>();
	//}

	template<size_t Z, typename T, Axis ax, typename G,
		typename std::enable_if_t<l_op_compatible<Z, G>, int> = 0>
	auto get_l_op(OpTerm<T, VectorComponent<ax, G>> const& e, double const*)
	{
		// * column_vector_from_axis(e)
		return std::make_pair(expr::coeff(e), OpVoid{});
	}

	namespace
	{
		template<size_t Z, typename... Es, size_t... Is>
		auto get_l_op_adds(OpAdd<Es...> const& e, double const* h, std::index_sequence<Is...>)
		{
			auto&& a = std::make_tuple(get_l_op<Z>(expr::get<Is>(e), h)...);
			return std::make_pair(expr::add_all(std::get<Is>(a).first...), expr::add_all(std::get<Is>(a).second...));
		}
	}

	template<size_t Z, typename... Es>
	auto get_l_op(OpAdd<Es...> const& e, double const* h)
	{
		return get_l_op_adds<Z>(e, h, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<size_t Z, typename A1, typename A2, typename E>
	auto get_l_op(OpChain<A1, A2, E> const& e, double const* h)
	{
		constexpr size_t D = expr::grid_dim<E>::dimension;
		len_type dims[D];
		expr::fill_data_dimensions(e, dims);

		auto op_ft = expr::transform::to_ft<expr::grid_dim<E>::dimension>(e.combination, h, dims);
		auto&& [expr_l, non_op] = get_l_op<Z>(expr::get_enclosed_expression(e), h);

		return std::make_pair(op_ft * expr_l, e.combination(non_op));
	}

	template<size_t Z, typename A1, typename A2, typename E>
	auto get_l_op(OpCombination<A1, A2, E> const& e, double const* h)
	{
		constexpr size_t D = expr::grid_dim<E>::dimension;
		len_type dims[D];
		expr::fill_data_dimensions(e, dims);

		auto op_ft = expr::transform::to_ft<expr::grid_dim<E>::dimension>(e.combination, h, dims);
		auto&& [expr_l, non_op] = get_l_op<Z>(expr::get_enclosed_expression(e), h);

		return std::make_pair(op_ft * expr_l, e.combination(non_op));
	}

	template<size_t Z, typename Dd, typename V, typename E, typename Sp>
	auto get_l_op(OpFuncDerivative<Dd, V, E, Sp> const& e, double const* h)
	{
		constexpr size_t D = expr::grid_dim<E>::dimension;
		len_type dims[D];
		expr::fill_data_dimensions(e, dims);

		constexpr size_t order = OpFuncDerivative<Dd, V, E, Sp>::order;
		constexpr Axis axis = OpFuncDerivative<Dd, V, E, Sp>::axis;

		static_assert(order > 0);

		auto [op, en] = expr::split::separate_operator(e);
		auto&& [expr_l, non_op] = get_l_op<Z>(en, h);

		return std::make_pair(
			expr::transform::to_ft<D>(op, h, dims) * expr_l,
			expr::make_derivative<Dd>(non_op, e.solver));
	}





	// **************************************************************************************


}


template<typename E>
struct SpectralData
{
	E scheme;

	SpectralData(E const& scheme) : scheme{ scheme } {}

	void update()
	{
		expr::prune::update(scheme);
	}
};


