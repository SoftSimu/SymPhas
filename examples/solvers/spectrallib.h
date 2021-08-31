
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

	template<size_t D, typename E>
	auto form_A_op(OpExpression<E> const& l_op, double dt, len_type const*);



	namespace
	{

		// **************************************************************************************

		/*
		 * given the (differential) operator, converts it into the expression form required
		 * to determine the B operator
		 */

		template<size_t D, typename A1, typename A2, typename L>
		auto convert_nl_derivative(OpOperatorChain<A1, A2> const& e, 
			OpExpression<L> const& l_op, double dt, double const* h, const len_type* dims);
		template<size_t D, size_t O, typename Sp, typename V, typename L>
		auto convert_nl_derivative(OpOperatorDerivative<O, Sp, V> const&, 
			OpExpression<L> const& l_op, double dt, double const* h, const len_type*);
		template<size_t D, typename T, typename L>
		auto convert_nl_derivative(OpLiteral<T> const& a, OpExpression<L> const& l_op, double dt, double const* h, const len_type* dims);



		template<size_t D, typename L, typename std::enable_if_t<!std::is_same<L, OpVoid>::value, int> = 0>
		auto convert_nl_derivative(OpIdentity const,
			OpExpression<L> const& l_op, double dt, double const*, const len_type* dims)
		{
			return (form_A_op<D>(*static_cast<L const*>(&l_op), dt, dims) - OpIdentity{}) / ((*static_cast<L const*>(&l_op)));
		}

		template<size_t D, typename L, typename std::enable_if_t<std::is_same<L, OpVoid>::value, int> = 0>
		auto convert_nl_derivative(OpIdentity const,
			OpExpression<L> const&, double dt, double const*, const len_type*)
		{
			return expr::make_literal(dt);
		}

		template<size_t D, typename T, typename L>
		auto convert_nl_derivative(OpLiteral<T> const& a,
			OpExpression<L> const& l_op, double dt, double const* h, const len_type* dims)
		{
			return a.value * convert_nl_derivative<D>(OpIdentity{}, *static_cast<L const*>(&l_op), dt, h, dims);
		}

		template<size_t D, typename A1, typename A2, typename L>
		auto convert_nl_derivative(OpOperatorCombination<A1, A2> const& e, 
			OpExpression<L> const& l_op, double dt, double const* h, const len_type* dims)
		{
			return convert_nl_derivative<D>(e.f, *static_cast<L const*>(&l_op), dt, h, dims) 
				+ convert_nl_derivative<D>(e.g, *static_cast<L const*>(&l_op), dt, h, dims);
		}

		template<size_t D, typename A1, typename A2, typename L>
		auto convert_nl_derivative(OpOperatorChain<A1, A2> const& e, 
			OpExpression<L> const& l_op, double dt, double const* h, const len_type* dims)
		{
			return convert_nl_derivative<D>(e.f, *static_cast<L const*>(&l_op), dt, h, dims) 
				* convert_nl_derivative<D>(e.g, *static_cast<L const*>(&l_op), dt, h, dims);
		}

		template<size_t D, size_t O, typename Sp, typename V, typename L>
		auto convert_nl_derivative(OpOperatorDerivative<O, Sp, V> const& e, 
			OpExpression<L> const& l_op, double dt, double const* h, const len_type* dims)
		{
			auto kk = expr::transform::to_ft<D>(e, h, dims);
			return kk * convert_nl_derivative<D>(OpIdentity{}, *static_cast<L const*>(&l_op), dt, h, dims);
		}




		/*
		 * takes the operators that have been factored from the nonlinear parts, and for each member
		 * of the tuple, it will compute the appropriate operator
		 * this also requires the original linear operator, which is passed as the second parameter
		 */
		template<size_t D, typename D_op, typename L>
		auto get_nl_op(OpExpression<D_op> const& d_op, OpExpression<L> const& l_op, double dt, double const* h, len_type const* dims)
		{
			return convert_nl_derivative<D>(*static_cast<D_op const*>(&d_op), *static_cast<L const*>(&l_op), dt, h, dims);
		}

		template<size_t D, typename D_op>
		auto get_nl_op(OpExpression<D_op> const&, OpVoid const, double dt, double const*, len_type const*)
		{
			return expr::make_literal(dt);
		}



		template<size_t D, typename... E_ops, typename L, size_t... Is>
		auto get_nl_op(std::tuple<E_ops...> const& nls, L&& l_op, 
			double dt, double const* h, len_type const* dims, std::index_sequence<Is...>)
		{
			return std::make_tuple(get_nl_op<D>(std::get<Is>(nls), std::forward<L>(l_op), dt, h, dims)...);
		}


		// **************************************************************************************


		template<typename... S, typename E>
		auto swap_var_apply(std::tuple<S...> const&, E&& expr, std::index_sequence<>)
		{
			return std::forward<E>(expr);
		}

		// swap for the primary variables (in model.system<N>())
		template<typename... S, typename E, size_t Q0, size_t... Qs>
		auto swap_var_apply(std::tuple<S...> const& systems, OpExpression<E> const& expr, std::index_sequence<Q0, Qs...>)
		{
#ifdef PRINTABLE_EQUATIONS
			auto replace = NamedData(
				reinterpret_cast<complex_t*>(std::get<Q0>(systems).frame_t),
				std::string(SYEX_FT_OF_OP_FMT_A)
				+ std::string(expr::get_op_name(expr::property::get_data_variable<Q0>(expr)))
				+ std::string(SYEX_FT_OF_OP_FMT_B)
			);
#else
			auto replace = reinterpret_cast<complex_t*>(std::get<Q0>(systems).frame_t);

#endif

			if constexpr (sizeof...(Qs) == 0)
			{
				return expr::transform::swap_grid<Q0>(*static_cast<const E*>(&expr), replace);
			}
			else
			{
				return swap_var_apply(
					systems, 
					expr::transform::swap_grid<Q0>(*static_cast<const E*>(&expr), replace), 
					std::index_sequence<Qs...>{});
			}
		}



		// **************************************************************************************


		/*
		 * this object is necessary for the underlying memory management
		 * it applies an existing object as a using declaration by substituting one of the parameters
		 */

		template<typename T, size_t D>
		using B_working = FourierGrid<T, complex_t, D>;



		/* construct the working tuple
		 */

		template<typename T_src, size_t D, typename E>
		auto make_working_tuple(OpExpression<E>& B_expression, const len_type* dims)
		{
			auto B = *static_cast<E const*>(&B_expression);

#ifdef PRINTABLE_EQUATIONS
			char* name = new char[B.print_length() + 1];
			B.print(name);
			auto B_op = expr::make_op(NamedData(expr::transform::to_fftw_grid(B), name));
			delete[] name;
#else
			auto B_op = expr::make_op(expr::transform::to_fftw_grid(B));
#endif
			
			return std::make_pair(B_op, B_working<T_src, D>(dims));
		}

		template<typename T_src, size_t D, typename T>
		auto make_working_tuple(OpLiteral<T> const& e, const len_type* dims)
		{
			return std::make_pair(e, B_working<T_src, D>(dims));
		}


		template<typename T_src, size_t D, typename... Ts, size_t... Is>
		decltype(auto) pair_B_with_working(std::tuple<Ts...>& t, const len_type* dims, std::index_sequence<Is...>)
		{
			return std::make_tuple(make_working_tuple<T_src, D>(std::get<Is>(t), dims)...);
		}

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

	template<size_t D, typename A1, typename A2>
	auto convert_to_ft(OpOperatorCombination<A1, A2> const& e, double const* h, const len_type* dims);
	template<size_t D, typename A1, typename A2>
	auto convert_to_ft(OpOperatorChain<A1, A2> const& e, double const* h, const len_type* dims);
	template<size_t D, typename E1, typename E2>
	auto convert_to_ft(OpBinarySub<E1, E2> const& e, double const* h, const len_type* dims);
	template<size_t D, typename E1, typename E2>
	auto convert_to_ft(OpBinaryAdd<E1, E2> const& e, double const* h, const len_type* dims);
	template<size_t D, typename E1, typename E2>
	auto convert_to_ft(OpBinaryMul<E1, E2> const& e, double const* h, const len_type* dims);
	template<size_t D, typename A1, typename A2, typename E>
	auto convert_to_ft(OpChain<A1, A2, E> const& e, double const* h, const len_type* dims);
	template<size_t D, typename A1, typename A2, typename E>
	auto convert_to_ft(OpCombination<A1, A2, E> const& e, double const* h, const len_type* dims);

	//

	template<size_t D, typename E>
	auto convert_to_ft(OpExpression<E> const& e, double const*, const len_type*)
	{
		return *static_cast<E const*>(&e);
	}

	template<size_t D, typename V, typename E1, typename E2>
	auto convert_to_ft(OpFuncConvolution<V, E1, E2> const& e, double const* h, const len_type* dims)
	{
		return expr::make_literal(e.value) * convert_to_ft<D>(e.a, h, dims) * convert_to_ft<D>(e.b, h, dims);
	}

	template<size_t D, typename V, typename E>
	auto convert_to_ft(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e, double const* h, const len_type* dims)
	{
		return expr::make_literal(e.value) * expr::transform::to_ft<D>(e.smoother, h, dims) * convert_to_ft<D>(expr::compound_get::expr(e), h, dims);
	}

	template<size_t D, size_t O, typename V, typename Sp>
	auto convert_to_ft(OpOperatorDerivative<O, V, Sp> const& e, double const* h, const len_type* dims)
	{
		return expr::make_literal(e.value) * expr::transform::to_ft<D>(e, h, dims);
	}

	template<size_t D, typename A1, typename A2>
	auto convert_to_ft(OpOperatorCombination<A1, A2> const& e, double const* h, const len_type* dims)
	{
		return convert_to_ft<D>(e.f, h, dims) + convert_to_ft<D>(e.g, h, dims);
	}

	template<size_t D, typename A1, typename A2, typename E>
	auto convert_to_ft(OpCombination<A1, A2, E> const& e, double const* h, const len_type* dims)
	{
		return convert_to_ft<D>(e.combination, h, dims) * convert_to_ft<D>(expr::compound_get::expr(e), h, dims);
	}

	template<size_t D, typename A1, typename A2>
	auto convert_to_ft(OpOperatorChain<A1, A2> const& e, double const* h, const len_type* dims)
	{
		return convert_to_ft<D>(e.f, h, dims) * convert_to_ft<D>(e.g, h, dims);
	}

	template<size_t D, typename A1, typename A2, typename E>
	auto convert_to_ft(OpChain<A1, A2, E> const& e, double const* h, const len_type* dims)
	{
		return convert_to_ft<D>(e.combination, h, dims) * convert_to_ft<D>(expr::compound_get::expr(e), h, dims);
	}

	template<size_t D, typename E1, typename E2>
	auto convert_to_ft(OpBinaryAdd<E1, E2> const& e, double const* h, const len_type* dims)
	{
		return convert_to_ft<D>(e.a, h, dims) + convert_to_ft<D>(e.b, h, dims);
	}

	template<size_t D, typename E1, typename E2>
	auto convert_to_ft(OpBinarySub<E1, E2> const& e, double const* h, const len_type* dims)
	{
		return convert_to_ft<D>(e.a, h, dims) - convert_to_ft<D>(e.b, h, dims);
	}

	template<size_t D, typename E1, typename E2>
	auto convert_to_ft(OpBinaryMul<E1, E2> const& e, double const*, const len_type*)
	{
		return expr::make_convolution(e.a, e.b);
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
	template<size_t Z, typename E1, typename E2>
	auto get_l_op(OpBinaryAdd<E1, E2> const& e, double const* h);
	template<size_t Z, typename E1, typename E2>
	auto get_l_op(OpBinarySub<E1, E2> const& e, double const* h);
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
	constexpr bool l_op_compatible = expr::property::vars<E>::template only_id<Z>();

	template<size_t Z, typename E>
	auto get_l_op(OpExpression<E> const& e, double const*)
	{
		return std::make_pair(OpVoid{}, *static_cast<E const*>(e));
	}

	template<size_t Z, typename T>
	auto get_l_op(OpLiteral<T> const e, double const*)
	{
		return std::make_pair(e, OpVoid{});
	}

	template<size_t Z>
	inline auto get_l_op(OpNegIdentity const, double const*)
	{
		return std::make_pair(OpNegIdentity{}, OpVoid{});
	}

	template<size_t Z>
	inline auto get_l_op(OpIdentity const, double const*)
	{
		return std::make_pair(OpIdentity{}, OpVoid{});
	}

	template<size_t Z>
	inline auto get_l_op(OpVoid const, double const*)
	{
		return std::make_pair(OpVoid{}, OpVoid{});
	}

	template<size_t Z, typename T, typename G, 
		typename std::enable_if_t<l_op_compatible<Z, G>, int> = 0>
	auto get_l_op(OpLVariable<T, G> const& e, double const*)
	{
		return std::make_pair(expr::make_literal(e.value), OpVoid{});
	}

	template<size_t Z, typename T, typename G,
		typename std::enable_if_t<!l_op_compatible<Z, G>, int> = 0>
		auto get_l_op(OpLVariable<T, G> const& e, double const*)
	{
		return std::make_pair(OpVoid{}, e);
	}

	template<size_t Z, typename V, typename E1, typename E2>
	auto get_l_op(OpFuncConvolution<V, E1, E2> const& e, double const*)
	{
		return std::make_pair(e, OpVoid{});// expr::make_literal(e.value)* l_op_a* l_op_b;
	}

	template<size_t Z, typename V, size_t D, typename E>
	auto get_l_op(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e, double const*)
	{
		return std::make_pair(OpVoid{}, e);
	}

	template<size_t Z, typename E1, typename E2>
	auto get_l_op(OpBinaryAdd<E1, E2> const& e, double const* h)
	{
		auto&& [l_op_a, non_op_a] = get_l_op<Z>(e.a, h);
		auto&& [l_op_b, non_op_b] = get_l_op<Z>(e.b, h);
		return std::make_pair(l_op_a + l_op_b, non_op_a + non_op_b);
	}

	template<size_t Z, typename E1, typename E2>
	auto get_l_op(OpBinarySub<E1, E2> const& e, double const* h)
	{
		auto&& [l_op_a, non_op_a] = get_l_op<Z>(e.a, h);
		auto&& [l_op_b, non_op_b] = get_l_op<Z>(e.b, h);
		return std::make_pair(l_op_a - l_op_b, non_op_a - non_op_b);
	}

	template<size_t Z, typename A1, typename A2, typename E>
	auto get_l_op(OpChain<A1, A2, E> const& e, double const* h)
	{
		const len_type* dims = expr::property::data_dimensions(expr::compound_get::template expr(e));

		auto chain_l = convert_to_ft<expr::grid_dim<E>::dimension>(e.combination, h, dims);
		auto&& [expr_l, non_op] = get_l_op<Z>(expr::compound_get::template expr(e), h);

		return std::make_pair(chain_l * expr_l, non_op);
	}

	template<size_t Z, typename A1, typename A2, typename E>
	auto get_l_op(OpCombination<A1, A2, E> const& e, double const* h)
	{
		const len_type* dims = expr::property::data_dimensions(expr::compound_get::template expr(e));

		auto combination_l = convert_to_ft<expr::grid_dim<E>::dimension>(e.combination, h, dims);
		auto&& [expr_l, non_op] = get_l_op<Z>(expr::compound_get::template expr(e), h);

		return std::make_pair(combination_l * expr_l, non_op);
	}

	template<size_t Z, typename Dd, typename V, typename E, typename Sp>
	auto get_l_op(OpFuncDerivative<Dd, V, E, Sp> const& e, double const* h)
	{
		const len_type* dims = expr::property::data_dimensions(e);

		auto deriv_l = convert_to_ft<expr::grid_dim<E>::dimension>(expr::make_operator_derivative<Dd::order>(e.solver), h, dims);
		auto&& [expr_l, non_op] = get_l_op<Z>(expr::compound_get::template expr(e), h);

		return std::make_pair(expr::make_literal(e.value) * deriv_l * expr_l, non_op);
	}





	// **************************************************************************************

	template<size_t D, typename... Ps, size_t... Is>
	decltype(auto) get_nonlinear_exprs(std::tuple<Ps...> const& npl, double const* h, const len_type* dims, std::index_sequence<Is...>)
	{
		return std::make_tuple(std::make_pair(convert_to_ft<D>(std::get<0>(std::get<Is>(npl)), h, dims), std::get<1>(std::get<Is>(npl)))...);
	}

	/* given a tuple where each element is a pair, where that pair consists of the expression
	 * which is nonlinear in Z and the other expression which is purely nonlinear, then it will
	 * only modify the expressions which are nonlinear in Z by converting them to Fourier space
	 */
	template<size_t D, typename... Ps>
	decltype(auto) get_nonlinear_exprs(std::tuple<Ps...> const& npl, double const* h, const len_type* dims)
	{
		return get_nonlinear_exprs<D>(npl, h, dims, std::make_index_sequence<sizeof...(Ps)>{});
	}





	// **************************************************************************************


	//! Convert the tuple of given expressions into nonlinear operators. 
	/*!
	 * Requires a tuple of all the given operators in the same order as the 
	 * derivative (operators) are generated; order matters so that the operators 
	 * can be associated with the correct expression ultimately. It will map the given 
	 * expression to the correct function so it can return an operator for the 
	 * given nonlinear part.
	 */

	template<size_t D, typename... E_ops, typename L>
	decltype(auto) get_nl_op(std::tuple<E_ops...> const& nls, L&& l_op, double dt, double const* h, len_type const* dims)
	{
		return get_nl_op<D>(nls, std::forward<L>(l_op), dt, h, dims, std::make_index_sequence<sizeof...(E_ops)>{});
	}

	// **************************************************************************************


	/* initiate the variable swapping
	 */

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
		return swap_var_apply(systems, *static_cast<const E*>(&expr), expr::property::vars<E>::get_ids());
	}



	// **************************************************************************************




	/* 
	 * Form the operators:
	 * the respective operators are formed using the functions from the spectral 
	 * solver namespace, and additionally pre-evaluated over its entire domain.
	 ***************************************************************************/

	//! Given the linear expression, returns the \f$A\f$ operator.
	template<size_t D, typename E>
	auto form_A_op(OpExpression<E> const& l_op, double dt, len_type const*)
	{
		auto A_expr = expr::exp(expr::make_literal(dt) * (*static_cast<E const*>(&l_op)));
		return A_expr;
	}

	//! Given a list of nonlinear expressions, returns the \f$B\f$ operators.	
	template<size_t D, typename T, typename L>
	decltype(auto) form_B_ops(T&& nl_tuple, L&& l_op, double dt, double const* h, len_type const* dims)
	{
		return solver_sp::get_nl_op<D>(std::forward<T>(nl_tuple), std::forward<L>(l_op), dt, h, dims);
	}


	// **************************************************************************************



	//! Return an object for managing the nonlinear expression.
	/*!
	 * Combines the given tuple of \f$B\f$ operators into a new tuple containing  
	 * the \f$B\f$ operator associated with a grid that can compute the in 
	 * place Fourier transform of an array of data of the prescribed dimensions.
	 * 
	 * \param t The tuple of \f$B\f$ operators.
	 * \param dims The dimensions of the system.
	 */
	template<typename T_src, size_t D, typename... Ts>
	decltype(auto) pair_B_with_working(std::tuple<Ts...>& t, const len_type* dims)
	{
		return pair_B_with_working<T_src, D>(t, dims, std::make_index_sequence<sizeof...(Ts)>{});
	}


	//! Creates the expression for the nonlinear part of the spectral form.
	/*!
	 * Given the object used for the \f$B\f$ operator and the pair of the
	 * nonlinear parts, construct the nonlinear expression of the spectral form.
	 * This overload considers that both parts of the pair are nonzero. In this
	 * case, the first part represents an expression which has had Fourier terms
	 * substituted, and the second part represents real space; therefore the
	 * second element is instead incorporated through the data element of the
	 * \f$B\f$ operator object.
	 * 
	 * \param B_pack The object containing the evaluated \f$B\f$ operator as
	 * well as data which can perform a Fourier transform of real data.
	 * \param nls The pair of nonlinear expressions.
	 */
	template<typename B, typename F, typename LNL, typename NL>
	auto expr_nl(std::pair<B, F>& B_pack, std::pair<LNL, NL> const& nls, size_t = 0)
	{
		auto& [l_nl, nl] = nls;
		auto& [B_op, b] = B_pack;

#ifdef PRINTABLE_EQUATIONS
		char* nlstr = new char[nl.print_length() + SYEX_FT_OF_EXPR_FMT_LEN + 1];
		
		size_t n = sprintf(nlstr, SYEX_FT_OF_EXPR_FMT_A);
		n += nl.print(nlstr + n);
		n += sprintf(nlstr + n, SYEX_FT_OF_EXPR_FMT_B);

		std::string s{ nlstr };
		delete[] nlstr;

		return expr::make_mul(B_op, (expr::make_op(NamedData(b.values, s)) + l_nl));
#else
		return expr::make_mul(B_op, (expr::make_op(b.values) + l_nl));
#endif
	}

	//! Creates the expression for the nonlinear part of the spectral form.
	/*!
	 * Given the object used for the \f$B\f$ operator and the pair of the
	 * nonlinear parts, construct the nonlinear expression of the spectral form.
	 * This overload considers only the fully nonlinear part is nonzero. In this
	 * case, the second part represents real space, so the second element is 
	 * instead incorporated through the data element of the \f$B\f$ operator 
	 * object.
	 * 
	 * \param B_pack The object containing the evaluated \f$B\f$ operator as
	 * well as data which can perform a Fourier transform of real data.
	 * \param nls The pair of nonlinear expressions.
	 */
	template<typename B, typename F, typename NL>
	auto expr_nl(std::pair<B, F>& B_pack, std::pair<OpVoid, NL> const& nls, size_t = 0)
	{
		auto& nl = std::get<1>(nls);
		auto& [B_op, b] = B_pack;

#ifdef PRINTABLE_EQUATIONS
		char* nlstr = new char[nl.print_length() + SYEX_FT_OF_EXPR_FMT_LEN + 1];

		size_t n = sprintf(nlstr, SYEX_FT_OF_EXPR_FMT_A);
		n += nl.print(nlstr + n);
		n += sprintf(nlstr + n, SYEX_FT_OF_EXPR_FMT_B);

		std::string s{ nlstr };
		delete[] nlstr;

		return expr::make_mul(B_op, expr::make_op(NamedData(b.values, s)));
#else
		return expr::make_mul(B_op, expr::make_op(b.values));

#endif

	}

	//! Creates the expression for the nonlinear part of the spectral form.
	/*!
	 * Given the object used for the \f$B\f$ operator and the pair of the
	 * nonlinear parts, construct the nonlinear expression of the spectral form.
	 * This overload considers only the first part of the pair is nonzero. That
	 * is, the first part represents an expression which has had Fourier terms
	 * substituted, and can be directly multiplied by the \f$B\f$ operator.
	 * 
	 * \param B_pack The object containing the evaluated \f$B\f$ operator as
	 * well as data which can perform a Fourier transform of real data.
	 * \param nls The pair of nonlinear expressions.
	 */
	template<typename B, typename F, typename LNL>
	auto expr_nl(std::pair<B, F>& B_pack, std::pair<LNL, OpVoid> const& nls, size_t = 0)
	{
		auto& l_nl = std::get<0>(nls);
		return expr::make_mul(std::get<0>(B_pack), l_nl);
	}

	template<typename F, typename LNL>
	auto expr_nl(std::pair<OpIdentity, F>&, std::pair<LNL, OpVoid> const& nls, size_t = 0)
	{
		return std::get<0>(nls);
	}


	template<typename... Bs, typename... LNLs, typename... NLs, size_t... Is>
	auto expr_nls(std::tuple<Bs...>& Bs_tuple, 
		std::tuple<std::pair<LNLs, NLs>...>& nls_tuple, 
		std::index_sequence<Is...>)
	{
		return ((OpVoid{} + ... + expr_nl(std::get<Is>(Bs_tuple), std::get<Is>(nls_tuple), Is)));
	}


	//! Update the nonlinear expression of the spectral form.
	/*!
	 * Given the nonlinear pair of expressions, where both expressions are
	 * is nonzero, evaluate the second and store the result into the data object
	 * managing the \f$B\f$ operator.
	 *
	 * \param B_pack The object containing the evaluated \f$B\f$ operator and
	 * a grid into which the nonlinear expression is evaluated to.
	 * \param nls The pair of nonlinear expressions.
	 */
	template<typename B, typename F, typename LNL, typename NL>
	auto update_nl(std::pair<B, F>& B_pack, std::pair<LNL, NL>& nls)
	{
		auto& b = std::get<1>(B_pack);
		expr::result(std::get<1>(nls), b.values_src_cast(), expr::property::data_len(std::get<1>(nls)));
		b.update();
	}

	//! Update the nonlinear expression of the spectral form.
	/*!
	 * Updating the nonlinear expression only needs to be performed on the
	 * second element of the nonlinear expression pair, which is the nonlinear
	 * expression representing the real space. Thus this function doesn't do
	 * anything.
	 */
	template<typename B_type, typename pair_type, typename LNL>
	auto update_nl(std::pair<B_type, pair_type>&, std::pair<LNL, OpVoid>&) {}

	template<typename pair_type, typename NL>
	auto update_nl(std::pair<double, pair_type>&, std::pair<OpVoid, NL>&) {}
	template<typename pair_type, typename LNL, typename NL>
	auto update_nl(std::pair<double, pair_type>&, std::pair<LNL, NL>&) {}



	template<typename... Bs, typename... LNLs, typename... NLs, size_t... Is>
	auto update_nls(std::tuple<Bs...>& Bs_tuple, 
		std::tuple<std::pair<LNLs, NLs>...>& nls_tuple, 
		std::index_sequence<Is...>)
	{
		return ((update_nl(std::get<Is>(Bs_tuple), std::get<Is>(nls_tuple)), ...));
	}

	//! Creates a list of the nonlinear expressions in the spectral form.
	template<typename... Bs, typename... LNLs, typename... NLs>
	auto expr_nls(std::tuple<Bs...>& Bs_tuple, 
		std::tuple<std::pair<LNLs, NLs>...>& nls_tuple)
	{
		constexpr size_t LEN = std::tuple_size<std::tuple<std::pair<LNLs, NLs>...>>::value;
		return expr_nls(Bs_tuple, nls_tuple, std::make_index_sequence<LEN>{});
	}

	//! Updates the expressions in the spectral form.
	template<typename... Bs, typename... LNLs, typename... NLs>
	auto update_nls(std::tuple<Bs...>& Bs_tuple, std::tuple<std::pair<LNLs, NLs>...>& nls_tuple)
	{
		constexpr size_t LEN = std::tuple_size<std::tuple<std::pair<LNLs, NLs>...>>::value;
		return update_nls(Bs_tuple, nls_tuple, std::make_index_sequence<LEN>{});
	}


}





template<typename At, typename Bs, typename NLs>
struct SpectralData
{

	At A_op;
	Bs Bs_tuple;
	NLs nls_tuple;
	const complex_t* values;
	std::string name;

protected:

	auto make_evolution();

public:

	using evolution_type = typename std::invoke_result_t<decltype(&SpectralData<At, Bs, NLs>::make_evolution), SpectralData>;
	evolution_type evolution_equation;


	SpectralData(At const& A_op, Bs const& Bs_tuple, NLs const& nls_tuple, const complex_t* values, std::string name = "?") :
		A_op{ A_op }, Bs_tuple{ Bs_tuple }, nls_tuple{ nls_tuple }, values{ values }, name{ name },
		evolution_equation{ make_evolution() } {}
	SpectralData(SpectralData<At, Bs, NLs> const& other) : SpectralData(other.A_op, other.Bs_tuple, other.nls_tuple, other.values, other.name) {}
	SpectralData(SpectralData<At, Bs, NLs>&& other) noexcept : SpectralData(other.A_op, other.Bs_tuple, other.nls_tuple, other.values, other.name) {}


	SpectralData<At, Bs, NLs>& operator=(SpectralData<At, Bs, NLs> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(SpectralData<At, Bs, NLs>& first, SpectralData<At, Bs, NLs>& second)
	{
		using std::swap;
		swap(first.A_op, second.A_op);
		swap(first.Bs_tuple, second.Bs_tuple);
		swap(first.nls_tuple, second.nls_tuple);
		swap(first.values, second.values);
		swap(first.name, second.name);
	}

	void update()
	{
		solver_sp::update_nls(Bs_tuple, nls_tuple);
	}
};


template<typename At, typename Bs, typename NLs>
auto SpectralData<At, Bs, NLs>::make_evolution()
{
	return A_op 
		* expr::make_op(NamedData(values, std::string(SYEX_FT_OF_OP_FMT_A) + name + std::string(SYEX_FT_OF_OP_FMT_B))) 
		+ solver_sp::expr_nls(Bs_tuple, nls_tuple);
}



