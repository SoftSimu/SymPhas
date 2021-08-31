
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


#include <type_traits>

#include "expressionconvolution.h"
#include "expressionderivatives.h"


namespace expr
{
	//! Prune an expression to update the state of all nested terms.
	/*! 
	 * Specifies algorithms for different expression types that will update
	 * the underlying data before the expression can be evaluated.
	 */
	namespace prune {}
}

namespace expr::prune
{

	namespace
	{
		template<typename E>
		void _update(E&) {}

		template<typename V, typename E>
		void _update(OpExponential<V, E>& e);
		template<typename V, typename T, typename G>
		void _update(OpExponential<V, OpLVariable<T, G>>&);
		template<typename V, typename E1, typename E2>
		void _update(OpFuncConvolution<V, E1, E2>& e);
		template<typename V, typename E, typename F, typename... Args>
		void _update(OpFunc<V, E, F, Args...>& e);
		template<typename V, size_t D, typename E>
		void _update(OpFuncConvolution<V, GaussianSmoothing<D>, E>& e);
		template<typename V, size_t D, typename G>
		void _update(OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>>& e);
		template<typename Dd, typename V, typename G, typename Sp>
		void _update(OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp>&);
		template<typename Dd, typename V, typename E, typename Sp>
		void _update(OpFuncDerivative<Dd, V, E, Sp>& e);
		template<typename A1, typename A2>
		void _update(OpBinaryAdd<A1, A2>& e);
		template<typename A1, typename A2>
		void _update(OpBinarySub<A1, A2>& e);
		template<typename A1, typename A2>
		void _update(OpBinaryMul<A1, A2>& e);
		template<typename A1, typename A2>
		void _update(OpBinaryDiv<A1, A2>& e);
		template<typename A1, typename A2, typename E>
		void _update(OpCombination<A1, A2, E>& e);
		template<typename A1, typename A2, typename E>
		void _update(OpChain<A1, A2, E>& e);
		template<typename G, typename V, typename E>
		void _update(OpMap<G, V, E>& e);



		template<typename V, typename E>
		void _update(OpExponential<V, E>& e)
		{
			_update(expr::compound_get::template expr(e));
			e.update();
		}

		template<typename V, typename T, typename G>
		void _update(OpExponential<V, OpLVariable<T, G>>&) {}

		template<typename V, typename E1, typename E2>
		void _update(OpFuncConvolution<V, E1, E2>& e)
		{
			_update(e.a);
			_update(e.b);
			e.update();
		}

		template<typename V, typename E, typename F, typename... Args>
		void _update(OpFunc<V, E, F, Args...>& e)
		{
			_update(expr::compound_get::template expr(e));
		}

		template<typename V, size_t D, typename E>
		void _update(OpFuncConvolution<V, GaussianSmoothing<D>, E>& e)
		{
			_update(expr::compound_get::template expr(e));
			e.update();
		}

		template<typename V, size_t D, typename G>
		void _update(OpFuncConvolution<V, GaussianSmoothing<D>, OpLVariable<OpIdentity, G>>& e)
		{
			e.update();
		}


		/* derivative pruning
		 */

		template<typename Dd, typename V, typename G, typename Sp>
		void _update(OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp>&) {}

		template<typename Dd, typename V, typename E, typename Sp>
		void _update(OpFuncDerivative<Dd, V, E, Sp>& e)

		{
			_update(expr::compound_get::template expr(e));
			e.update();
		}


		/* binary op pruning
		 */

		template<typename A1, typename A2>
		void _update(OpBinaryAdd<A1, A2>& e)
		{
			_update(e.a);
			_update(e.b);
		}

		template<typename A1, typename A2>
		void _update(OpBinarySub<A1, A2>& e)
		{
			_update(e.a);
			_update(e.b);
		}

		template<typename A1, typename A2>
		void _update(OpBinaryMul<A1, A2>& e)
		{
			_update(e.a);
			_update(e.b);
		}

		template<typename A1, typename A2>
		void _update(OpBinaryDiv<A1, A2>& e)
		{
			_update(e.a);
			_update(e.b);
		}

		/* operator pruning
		 */

		template<typename A1, typename A2, typename E>
		void _update(OpCombination<A1, A2, E>& e)
		{
			_update(expr::compound_get::template expr(e));
			e.update();
		}
		template<typename A1, typename A2, typename E>
		void _update(OpChain<A1, A2, E>& e)
		{
			_update(expr::compound_get::template expr(e));
			_update(e.inner);
			_update(e.outer);
		}

		template<typename G, typename V, typename E>
		void _update(OpMap<G, V, E>& e)
		{
			_update(expr::compound_get::template expr(e));
			e.update();
		}
	}

	//! Update underlying the given expression.
	/*!
	 * For expressions which store intermediate data, such as derivatives, this
	 * data must be updated before the expression can be evaluated. This will
	 * traverse the expression tree and perform all necessary updating.
	 * 
	 * \param e The expression that is updated.
	 */
	template<typename E>
	void update(OpExpression<E>& e)
	{
		_update(*static_cast<E*>(&e));
	}

}


// ******************************************************************************************


namespace expr
{



	namespace
	{
		/* returns whether a derivative is present in the expression by
		 * checking that the derivative index is greater than the index associated
		 * with a constant
		 */

		template<typename E>
		struct expr_has_deriv
		{
			static const bool value = expr::property::derivative_index<2, E>::value > 1;
		};

	}

	//! Distribute derivatives applied to linear combinations.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<typename E>
	auto distribute_operators(OpExpression<E> const& e)
	{
		return *static_cast<E const*>(&e);
	}

	//! Distribute derivatives applied to linear combinations.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The operator which is distributed.
	 */
	template<typename E>
	auto distribute_operators(OpOperator<E> const& e)
	{
		return *static_cast<E const*>(&e);
	}


	//! Distribute derivatives applied to linear combinations.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<typename A1, typename A2, typename E>
	auto distribute_operators(OpChain<A1, A2, E> const& e);

	//! Distribute derivatives applied to linear combinations.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<typename A1, typename A2, typename E>
	auto distribute_operators(OpCombination<A1, A2, E> const& e);

	//! Distribute derivatives applied to linear combinations.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<typename Dd, typename V, typename E, typename Sp,
		typename Enable = typename std::enable_if<(expr_has_deriv<E>::value&& is_combination<E>::value)>::type>
		auto distribute_operators(OpFuncDerivative<Dd, V, E, Sp> const& e);
	template<typename Dd, typename V, typename E1, typename E2, typename Sp>
	auto distribute_operators(OpFuncDerivative<Dd, V, OpBinaryAdd<E1, E2>, Sp> const& e);
	template<typename Dd, typename V, typename E1, typename E2, typename Sp>
	auto distribute_operators(OpFuncDerivative<Dd, V, OpBinarySub<E1, E2>, Sp> const& e);

	//! Distribute derivatives applied to linear combinations.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<typename E1, typename E2>
	auto distribute_operators(OpBinaryAdd<E1, E2> const& e);

	//! Distribute derivatives applied to linear combinations.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<typename E1, typename E2>
	auto distribute_operators(OpBinarySub<E1, E2> const& e);

	//! Distribute derivatives applied to linear combinations.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<typename E1, typename E2>
	auto distribute_operators(OpBinaryMul<E1, E2> const& e);


	namespace
	{
		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<expr_has_deriv<E>::value, int> = 0>
			auto distribute_operators_chain(OpChain<A1, A2, E> const& e)
		{
			return distribute_operators(expr::distribute_operator(e.combination, expr::compound_get::template expr(e)));
		}

		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<!expr_has_deriv<E>::value, int> = 0>
			auto distribute_operators_chain(OpChain<A1, A2, E> const& e)
		{
			return expr::expand_operator(e.combination, expr::compound_get::template expr(e));
		}

	}

	template<typename A1, typename A2, typename E>
	auto distribute_operators(OpChain<A1, A2, E> const& e)
	{
		return distribute_operators_chain(e);
	}

	namespace
	{

		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<expr_has_deriv<E>::value, int> = 0>
			auto distribute_operators_combination(OpCombination<A1, A2, E> const& e)
		{
			return distribute_operators(expr::distribute_operator(e.combination, expr::compound_get::template expr(e)));
		}

		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<!expr_has_deriv<E>::value, int> = 0>
			auto distribute_operators_combination(OpCombination<A1, A2, E> const& e)
		{
			return expr::expand_operator(e.combination, expr::compound_get::template expr(e));
		}
	}

	template<typename A1, typename A2, typename E>
	auto distribute_operators(OpCombination<A1, A2, E> const& e)
	{
		return distribute_operators_combination(e);
	}


	template<typename Dd, typename V, typename E, typename Sp, typename Enable>
	auto distribute_operators(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return distribute_operators(expr::make_operator_derivative<Dd::order>(e.value, e.solver)
			* expr::compound_get::template expr(e));
	}

	template<typename Dd, typename V, typename E1, typename E2, typename Sp>
	auto distribute_operators(OpFuncDerivative<Dd, V, OpBinaryAdd<E1, E2>, Sp> const& e)
	{
		auto&& add = expr::compound_get::template expr(e);
		auto&& d = expr::make_operator_derivative<Dd::order>(e.value, e.solver);
		return distribute_operators(d * add.a + d * add.b);
	}

	template<typename Dd, typename V, typename E1, typename E2, typename Sp>
	auto distribute_operators(OpFuncDerivative<Dd, V, OpBinarySub<E1, E2>, Sp> const& e)
	{
		auto&& sub = expr::compound_get::template expr(e);
		auto&& d = expr::make_operator_derivative<Dd::order>(e.value, e.solver);
		return distribute_operators(d * sub.a - d * sub.b);
	}

	template<typename E1, typename E2>
	auto distribute_operators(OpBinaryAdd<E1, E2> const& e)
	{
		return distribute_operators(e.a) + distribute_operators(e.b);
	}

	template<typename E1, typename E2>
	auto distribute_operators(OpBinarySub<E1, E2> const& e)
	{
		return distribute_operators(e.a) - distribute_operators(e.b);
	}

	template<typename E1, typename E2>
	auto distribute_operators(OpBinaryMul<E1, E2> const& e)
	{
		return distribute_operators(e.a) * distribute_operators(e.b);
	}

	template<typename E1, typename E2>
	auto distribute_operators(OpBinaryDiv<E1, E2> const& e)
	{
		return distribute_operators(e.a) / distribute_operators(e.b);
	}

}

// ******************************************************************************************


namespace expr
{
	//! Defines elements which transform an expression into another one.
	/*!
	 * Implements the functions that transform an expression of a given
	 * form into another form. The functions typically constructing another 
	 * expression based on some rules.
	 */
	namespace transform {}
}

namespace symphas::internal
{
	//! Constructs a grid of the prescribed dimension.
	/*!
	 * A grid is constructed using the given dimension, and the expression
	 * evaluation type is used as the grid underlying type.
	 * 
	 * The expression is evaluated into the grid.
	 * 
	 * \tparam D The dimension of the grid.
	 */
	template<size_t D>
	struct construct_grid_of_dimension
	{
		//! Create a grid and populate it with the evaluated expression.
		/*!
		 * The grid is created using the dimensions of the underlying
		 * expression data, and then populated by evaluating the expression
		 * and storing the result in the grid.
		 * 
		 * \param e The expression with which to populate the grid.
		 */
		template<typename E, typename T = typename expr::eval_type<E>::type>
		auto operator()(OpExpression<E>& e)
		{
			expr::prune::update(e);
			Grid<T, D> result(expr::property::data_dimensions(e));
			expr::result(e, result.values, result.len);
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
	template<>
	struct construct_grid_of_dimension<0>
	{
		//! Create a grid and populate it with the evaluated expression.
		/*!
		 * The grid is created using the dimensions of the underlying
		 * expression data, and then populated by evaluating the expression
		 * and storing the result in the grid.
		 *
		 * \param e The expression with which to populate the grid.
		 */
		template<typename E, typename T = typename expr::eval_type<E>::type>
		auto operator()(OpExpression<E>& e)
		{
			expr::prune::update(e);
			Block<T> result(expr::property::data_len(e));
			expr::result(e, result.values, result.len);
			return result;
		}
	};
}


namespace expr::transform
{

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
	template<typename E, size_t D = expr::grid_dim<E>::value>
	auto to_grid(OpExpression<E>& e)
	{
		return symphas::internal::construct_grid_of_dimension<D>{}(e);
	}


	//! Evaluating a literal as a grid returns just the value of the literal.
	/*!
	 * Evaluating a literal as a grid returns just the value of the literal.
	 * See expr::transform::to_grid(OpExpressioN<E>&)
	 */
	template<typename T>
	auto to_grid(OpLiteral<T> const e)
	{
		return e.eval();
	}

	//! Evaluating an identity as a grid returns just the value of the literal.
	/*!
	 * Evaluating an identity as a grid returns just the value of the literal.
	 * See expr::transform::to_grid(OpExpressioN<E>&)
	 */
	inline auto to_grid(OpIdentity const e)
	{
		return e.eval();
	}

	//! Evaluating an identity as a grid returns just the value of the literal.
	/*!
	 * Evaluating an identity as a grid returns just the value of the literal.
	 * See expr::transform::to_grid(OpExpressioN<E>&)
	 */
	inline auto to_grid(OpNegIdentity const e)
	{
		return e.eval();
	}

	//! Evaluating an identity as a grid returns just the value of the literal.
	/*!
	 * Evaluating an identity as a grid returns just the value of the literal.
	 * See expr::transform::to_grid(OpExpressioN<E>&)
	 */
	inline auto to_grid(OpVoid const e)
	{
		return e.eval();
	}

	namespace
	{

		//! A list of expressions is evaluated into a list of grids.
		template<typename... Es, size_t... Is>
		decltype(auto) to_grid(std::tuple<Es...>& e, std::index_sequence<Is...>)
		{
			return std::make_tuple(expr::transform::to_grid(std::get<Is>(e))...);
		}
	}


	//! A list of expressions is evaluated into a list of grids.
	/*!
	 * Converts the expressions into grids, the order of the returned grid list
	 * is consistent with the list of expressions. Constructs the list by
	 * applying expr::to_grid(OpExpression<E>&) to each element of the list.
	 *
	 * \param es The list of expressions to evaluate into grids.
	 */
	template<typename... Es>
	decltype(auto) to_grid(std::tuple<Es...>& es)
	{
		return expr::transform::to_grid(es, std::make_index_sequence<sizeof...(Es)>{});
	}


	//! Create a grid where values are adjusted to account for FFTW routines.
	/*!
	 * Converts the expression into a grid with indexing consistent with the 
	 * the FFTW routine for that type. That is, 
	 * expected. The index of the values is adjusted based on what is expected
	 * in order to perform FFTW functions.
	 *
	 * \tparam T_src The known source type of the underlying data or evaluation
	 * type of the expression, which is considered when determining what index
	 * scheme to use for evaluation.
	 * \tparam D The dimension of the grid to create.
	 * \tparam T The result type of the expression.
	 */
	template<typename E, 
		typename T = typename expr::eval_type<E>::type, 
		size_t D = expr::grid_dim<E>::value,
		typename std::enable_if_t<expr::has_state<E>::value, int> = 0>
	Grid<T, D> to_fftw_grid(OpExpression<E>& e)
	{
		expr::prune::update(*static_cast<E*>(&e));

		const len_type* dims = expr::property::data_dimensions(e);

		len_type fftw_dims[D];
		fftw_dims[0] = dims[0] / 2 + 1;
		std::copy(dims + 1, dims + D, fftw_dims + 1);

		Grid<T, D> grid(dims);
		symphas::dft::iterate_rc<T, D>(grid.values, *static_cast<E const*>(&e), dims);

		return grid;
	}

	template<typename E, 
		typename T = typename expr::eval_type<E>::type, 
		size_t D = expr::grid_dim<E>::value,
		typename std::enable_if_t<!expr::has_state<E>::value, int> = 0>
		Grid<T, D> to_fftw_grid(OpExpression<E> const& e)
	{
		const len_type* dims = expr::property::data_dimensions(e);

		len_type fftw_dims[D];
		fftw_dims[0] = dims[0] / 2 + 1;
		std::copy(dims + 1, dims + D, fftw_dims + 1);

		Grid<T, D> grid(dims);
		symphas::dft::iterate_rc<T, D>(grid.values, *static_cast<E const*>(&e), dims);

		return grid;
	}

	//! Turning a constant to a grid returns the constant.
	/*!
	 * Turning a constant to a grid returns the constant.
	 * 
	 * \param e The constant to turn into a grid.
	 */
	template<typename T>
	auto to_fftw_grid(OpLiteral<T> const& e)
	{
		return e.eval();
	}

	//! Turning a constant to a grid returns the constant.
	/*!
	 * Turning a constant to a grid returns the constant.
	 *
	 * \param e The constant to turn into a grid.
	 */
	inline auto to_fftw_grid(OpIdentity const& e)
	{
		return e.eval();
	}

	//! Turning a constant to a grid returns the constant.
	/*!
	 * Turning a constant to a grid returns the constant.
	 *
	 * \param e The constant to turn into a grid.
	 */
	inline auto to_fftw_grid(OpNegIdentity const& e)
	{
		return e.eval();
	}

	//! Turning a zero to a grid returns zero.
	/*!
	 * Turning a zero to a grid returns the zero. Zero is the additive
	 * identity.
	 *
	 * \param e The zero.
	 */
	inline auto to_fftw_grid(OpVoid const& e)
	{
		return e.eval();
	}

	namespace
	{
		template<typename... Es, size_t... Is>
		decltype(auto) to_fftw_grid(std::tuple<Es...>& e, std::index_sequence<Is...>)
		{
			return std::make_tuple(expr::transform::to_fftw_grid(std::get<Is>(e))...);
		}
	}

	//! Create a list of grids initialized to the result of the expressions.
	/*!
	 * Converts the operator into a grid with indexing consistent with the type
	 * expected. Constructs the list by
	 * applying expr::to_fftw_grid(OpExpression<E>&) to each element of the list.
	 *
	 * \param es The list of expressions to evaluate into grids.
	 */
	template<typename... Es>
	decltype(auto) to_fftw_grid(std::tuple<Es...>& es)
	{
		return expr::transform::to_fftw_grid(es, std::make_index_sequence<sizeof...(Es)>{});
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
	template<size_t D, typename A1, typename A2>
	auto to_ft(OpOperatorCombination<A1, A2> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename A1, typename A2>
	auto to_ft(OpOperatorChain<A1, A2> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename G, typename V, typename E>
	auto to_ft(OpMap<G, V, E> const& e, double const* h, const len_type* dims);

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
	template<size_t D, size_t O, typename V, typename Sp>
	auto to_ft(OpOperatorDerivative<O, V, Sp> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename V, typename E1, typename E2>
	auto to_ft(OpFuncConvolution<V, E1, E2> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename V, typename E>
	auto to_ft(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename E1, typename E2>
	auto to_ft(OpBinaryAdd<E1, E2> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename E1, typename E2>
	auto to_ft(OpBinarySub<E1, E2> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename A1, typename A2, typename E>
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
	template<size_t D, typename A1, typename A2, typename E>
	auto to_ft(OpCombination<A1, A2, E> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename Dd, typename V, typename E, typename Sp>
	auto to_ft(OpFuncDerivative<Dd, V, E, Sp> const& e, double const* h, const len_type* dims);

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
	template<size_t D>
	auto to_ft(GaussianSmoothing<D> const& e, double const*, const len_type* dims);



	template<size_t D, typename E>
	auto to_ft(OpExpression<E> const& e, double const* h, const len_type*)
	{
		return expr::make_fourier_map(*static_cast<E const*>(&e));
	}

	template<size_t D>
	auto to_ft(OpVoid const, double const* h, const len_type*)
	{
		return OpVoid{};
	}

	template<size_t D>
	auto to_ft(OpIdentity const, double const* h, const len_type*)
	{
		return OpIdentity{};
	}

	template<size_t D>
	auto to_ft(OpNegIdentity const, double const* h, const len_type* dims)
	{
		return -to_ft<D>(OpIdentity{}, h, dims);
	}

	template<size_t D, typename T>
	auto to_ft(OpLiteral<T> const e, double const* h, const len_type* dims)
	{
		return e.value * to_ft<D>(OpIdentity{}, h, dims);
	}

	template<size_t D, typename T, typename A2>
	auto to_ft(OpOperatorCombination<OpLiteral<T>, A2> const& e, double const* h, const len_type* dims)
	{
		return expr::make_literal(e.g) + to_ft<D>(e.g, h, dims);
	}

	template<size_t D, typename A1, typename T>
	auto to_ft(OpOperatorCombination<A1, OpLiteral<T>> const& e, double const* h, const len_type* dims)
	{
		return to_ft<D>(e.f, h, dims) + expr::make_literal(e.g);
	}

	template<size_t D, typename A2>
	auto to_ft(OpOperatorCombination<OpIdentity, A2> const& e, double const* h, const len_type* dims)
	{
		return e.f + to_ft<D>(e.g, h, dims);
	}

	template<size_t D, typename A1>
	auto to_ft(OpOperatorCombination<A1, OpIdentity> const& e, double const* h, const len_type* dims)
	{
		return to_ft<D>(e.f, h, dims) + e.g;
	}

	template<size_t D, typename A2>
	auto to_ft(OpOperatorCombination<OpNegIdentity, A2> const& e, double const* h, const len_type* dims)
	{
		return e.f + to_ft<D>(e.g, h, dims);
	}

	template<size_t D, typename A1>
	auto to_ft(OpOperatorCombination<A1, OpNegIdentity> const& e, double const* h, const len_type* dims)
	{
		return to_ft<D>(e.f, h, dims) + e.g;
	}

	template<size_t D, typename A1, typename A2>
	auto to_ft(OpOperatorCombination<A1, A2> const& e, double const* h, const len_type* dims)
	{
		return to_ft<D>(e.f, h, dims) + to_ft<D>(e.g, h, dims);
	}

	template<size_t D, typename A1, typename A2>
	auto to_ft(OpOperatorChain<A1, A2> const& e, double const* h, const len_type* dims)
	{
		return to_ft<D>(e.f, h, dims) * to_ft<D>(e.g, h, dims);
	}

	template<size_t D, typename G, typename V, typename E>
	auto to_ft(OpMap<G, V, E> const& e, double const*, const len_type* dims)
	{
		return e.value * expr::compound_get::template expr(e);
	}


	template<size_t D, size_t O, typename V, typename Sp>
	auto to_ft(OpOperatorDerivative<O, V, Sp> const& e, double const* h, const len_type* dims)
	{
		return expr::make_literal(e.value * (((O / 2) % 2 == 0) ? 1 : -1)) * expr::make_op(K_Grid<O, D>(dims, h));
	}

	template<size_t D, typename V, typename E1, typename E2>
	auto to_ft(OpFuncConvolution<V, E1, E2> const& e, double const* h, const len_type* dims)
	{
		return expr::make_literal(e.value) * to_ft<D>(e.a, h, dims) * to_ft<D>(e.b, h, dims);
	}

	template<size_t D>
	auto to_ft(GaussianSmoothing<D> const& e, double const*, const len_type*)
	{
#ifdef PRINTABLE_EQUATIONS
		char* gname = new char[e.print_length() + STR_ARR_LEN(SYEX_FT_OF_OP_FMT_A SYEX_FT_OF_OP_FMT_B)];
		size_t n = sprintf(gname, SYEX_FT_OF_OP_FMT_A);
		n += e.print(gname + n);
		n =+ sprintf(gname + n, SYEX_FT_OF_OP_FMT_B);
		auto op = expr::make_op(NamedData(to_fftw_grid(e), gname));

		delete[] gname;
		return op;
#else
		return expr::make_op(to_fftw_grid(e));
#endif
	}

	template<size_t D, typename V, typename E>
	auto to_ft(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e, double const* h, const len_type* dims)
	{
		return expr::make_literal(e.value) * to_ft<D>(e.smoother, h, dims) * to_ft<D>(expr::compound_get::template expr(e), h, dims);
	}

	template<size_t D, typename E1, typename E2>
	auto to_ft(OpBinaryAdd<E1, E2> const& e, double const* h, const len_type* dims)
	{
		return to_ft<D>(e.a, h, dims) + to_ft<D>(e.b, h, dims);
	}

	template<size_t D, typename E1, typename E2>
	auto to_ft(OpBinarySub<E1, E2> const& e, double const* h, const len_type* dims)
	{
		return to_ft<D>(e.a, h, dims) - to_ft<D>(e.b, h, dims);
	}

	/* conversion of operators into fourier space
	 */

	template<size_t D, typename A1, typename A2, typename E>
	auto to_ft(OpChain<A1, A2, E> const& e, double const* h, const len_type* dims)
	{
		return to_ft<D>(e.combination, h, dims) * to_ft<D>(expr::compound_get::template expr(e), h, dims);
	}

	template<size_t D, typename A1, typename A2, typename E>
	auto to_ft(OpCombination<A1, A2, E> const& e, double const* h, const len_type* dims)
	{
		return to_ft<D>(e.combination, h, dims) * to_ft<D>(expr::compound_get::template expr(e), h, dims);
	}

	template<size_t D, typename Dd, typename V, typename E, typename Sp>
	auto to_ft(OpFuncDerivative<Dd, V, E, Sp> const& e, double const* h, const len_type* dims)
	{
		return expr::make_literal(e.value) 
			* to_ft<D>(expr::make_operator_derivative<Dd::order>(e.solver), h, dims) 
			* * to_ft<D>(expr::compound_get::template expr(e), h, dims);
	}


	// **************************************************************************************


	/* algorithm that takes the id of a variable and swaps all instances of the variable a
	 * different grid
	 * accounts for both being given a grid and being given a usual variable; doesn't handle
	 * swapping an expression in place of a variable
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
	template<size_t Z, typename E, typename G_F>
	decltype(auto) swap_grid(E const& e, G_F&& g);

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
	template<size_t Z, typename G, typename G_F>
	decltype(auto) swap_grid(Variable<Z, G> const&, G_F&& g);

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
	template<size_t Z, typename T, typename G, typename G_F>
	decltype(auto) swap_grid(OpLVariable<T, Variable<Z, G>> const& v, G_F&& g);

	//! Swap a data term in the expression.
	/*!
	 * Implementation of a successful search, where the given variable term
	 * associated with the prescribed index will be switched with the given
	 * expression.
	 *
	 * \param v The term which is swapped.
	 * \param g The expression that replaces the variable.
	 *
	 * \tparam Z The index of the variable to change, which matches the index
	 * of the variable associated with the variable term.
	 */
	template<size_t Z, typename T, typename G, typename E>
	auto swap_grid(OpLVariable<T, Variable<Z, G>> v, OpExpression<E> const& e);

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
	template<size_t Z, typename T, typename... Gs, typename G_F>
	auto swap_grid(OpNLVariable<T, Gs...> const& e, G_F&& g);

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
	template<size_t Z, typename E1, typename E2, typename G_F>
	auto swap_grid(OpBinaryAdd<E1, E2> const& e, G_F&& g);

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
	template<size_t Z, typename E1, typename E2, typename G_F>
	auto swap_grid(OpBinarySub<E1, E2> const& e, G_F&& g);

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
	template<size_t Z, typename E1, typename E2, typename G_F>
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
	template<size_t Z, typename V, typename E1, typename E2, typename G_F>
	auto swap_grid(OpFuncConvolution<V, E1, E2> const& e, G_F&& g);

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
	template<size_t Z, typename V, size_t D, typename E, typename G_F>
	auto swap_grid(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e, G_F&& g);

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
	template<size_t Z, typename Dd, typename V, typename E, typename Sp, typename G_F>
	auto swap_grid(OpFuncDerivative<Dd, V, E, Sp> const& e, G_F&& g);

	//

	namespace
	{

		/* swapping in the tuple
		 */
		template<size_t Z, typename... Gs, typename G_F, size_t... Is>
		auto _swap_grid_tuple(std::tuple<Gs...> const& ts, G_F&& g, std::index_sequence<Is...>)
		{
			return std::make_tuple(swap_grid<Z>(std::get<Is>(ts), std::forward<G_F>(g))...);
		}
	}


	template<size_t Z, typename E, typename G_F>
	decltype(auto) swap_grid(E const& e, G_F&&)
	{
		return e;
	}


	template<size_t Z, typename G, typename G_F>
	decltype(auto) swap_grid(Variable<Z, G> const&, G_F&& g)
	{
		return expr::as_variable<Z>(std::forward<G_F>(g));
	}

	/* opvariable swapping (lowest level in the recursion)
	 */

	template<size_t Z, typename T, typename G, typename G_F>
	decltype(auto) swap_grid(OpLVariable<T, Variable<Z, G>> const& v, G_F&& g)
	{
		return expr::make_op<Z>(v.value, std::forward<G_F>(g));
	}

	template<size_t Z, typename T, typename G, typename G_F>
	decltype(auto) swap_grid(OpLVariable<T, Variable<Z, NamedData<G>>> const& v, G_F&& g)
	{
		return expr::make_op<Z>(v.value, std::forward<G_F>(g));
	}

	// swap with an expression (including another oplvariable)
	template<size_t Z, typename T, typename G, typename E>
	auto swap_grid(OpLVariable<T, Variable<Z, G>> const& v, OpExpression<E> const& e)
	{
		return v.value * (*static_cast<const E*>(&e));
	}

	template<size_t Z, typename T, typename... Gs, typename G_F>
	auto swap_grid(OpNLVariable<T, Gs...> const& e, G_F&& g)
	{
		return OpNLVariable(e.value, _swap_grid_tuple<Z>(e.datas, std::forward<G_F>(g), std::make_index_sequence<sizeof...(Gs)>{}));
	}


	/* recursive swapping
	 */

	template<size_t Z, typename E1, typename E2, typename G_F>
	auto swap_grid(OpBinaryAdd<E1, E2> const& e, G_F&& g)
	{
		return swap_grid<Z>(e.a, std::forward<G_F>(g)) + swap_grid<Z>(e.b, std::forward<G_F>(g));
	}

	template<size_t Z, typename E1, typename E2, typename G_F>
	auto swap_grid(OpBinarySub<E1, E2> const& e, G_F&& g)
	{
		return swap_grid<Z>(e.a, std::forward<G_F>(g)) - swap_grid<Z>(e.b, std::forward<G_F>(g));
	}

	template<size_t Z, typename E1, typename E2, typename G_F>
	auto swap_grid(OpBinaryMul<E1, E2> const& e, G_F&& g)
	{
		return swap_grid<Z>(e.a, std::forward<G_F>(g)) * swap_grid<Z>(e.b, std::forward<G_F>(g));
	}

	template<size_t Z, typename V, typename E1, typename E2, typename G_F>
	auto swap_grid(OpFuncConvolution<V, E1, E2> const& e, G_F&& g)
	{
		return expr::make_convolution(
			e.value,
			swap_grid<Z>(e.a, std::forward<G_F>(g)), 
			swap_grid<Z>(e.b, std::forward<G_F>(g)));
	}

	template<size_t Z, typename V, size_t D, typename E, typename G_F>
	auto swap_grid(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e, G_F&& g)
	{
		return expr::make_convolution(
			e.value, 
			swap_grid<Z>(expr::compound_get::template expr(e), std::forward<G_F>(g)), 
			e.smoother);
	}


	template<size_t Z, typename A1, typename A2, typename E, typename G_F>
	auto swap_grid(OpChain<A1, A2, E> const& e, G_F&& g)
	{
		return e.combination * swap_grid<Z>(expr::compound_get::expr(e), std::forward<G_F>(g));
	}

	template<size_t Z, typename A1, typename A2, typename E, typename G_F>
	auto swap_grid(OpCombination<A1, A2, E> const& e, G_F&& g)
	{
		return e.combination * swap_grid<Z>(expr::compound_get::expr(e), std::forward<G_F>(g));
	}

	/* a derivative is treated slightly differently; an operator needs to be applied because the
	 * the swapped expression could result in a grid of a different type
	 */

	template<size_t Z, typename Dd, typename V, typename E, typename Sp, typename G_F>
	auto swap_grid(OpFuncDerivative<Dd, V, E, Sp> const& e, G_F&& g)
	{
		return expr::make_operator_derivative<Dd::order>(e.value, e.solver)
			* swap_grid<Z>(expr::compound_get::template expr(e), std::forward<G_F>(g));
	}

}



namespace expr
{

	//! Defines ways to split functions based on criteria.
	/*! 
	 * Implements all functions that will split an expression based
	 * on certain criteria.
	 */
	namespace split {}
}


namespace expr::split
{

	namespace
	{

		/* given the derivative index, retrieve the lowest order
		 */

		template<size_t I>
		struct min_order_from_index
		{
		protected:

			static constexpr size_t get_value()
			{
				iter_type OD = 0;
				while (true)
				{
					if (I & (1ull << OD))
					{
						return OD;
					}
					++OD;
				}
			}

		public:
			static const size_t value = get_value();
		};




		// a compile time check that the order of the derivative equals O
		template<size_t O, typename deriv_type>
		constexpr bool _check_equality_derivative()
		{
			return O == deriv_type::order;
		}





		/* convenience functions to pack an expression into either side of a pair
		 */

		template<typename E>
		auto pack_left(OpExpression<E> const& e)
		{
			return std::make_pair(*static_cast<const E*>(&e), OpVoid{});
		}

		template<typename E>
		auto pack_right(OpExpression<E> const& e)
		{
			return std::make_pair(OpVoid{}, *static_cast<const E*>(&e));
		}

	}

	//! Determine the smallest derivative order of the expression.
	/*!
	 * Determine the smallest derivative order of the expression.
	 * 
	 * \tparam E The expression type from which to get the minimum derivative
	 * order.
	 */
	template<typename E>
	struct min_derivative_order
	{
		static const size_t value = min_order_from_index<expr::property::derivative_index<0, E>::value>::value;
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
	template<typename Dd, typename V, typename E, typename Sp>
	auto by_linear(OpFuncDerivative<Dd, V, E, Sp> const& e);

	//! Split all linear variables from nonlinear variables.
	/*!
	 * Using the available predicates to determine linearity/nonlinearity of
	 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
	 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
	 * for linearity. Similarly, the `NL` part is nonlinear.
	 *
	 * \param e The expression which is split.
	 */
	template<typename A1, typename A2, typename E>
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
	template<typename A1, typename A2, typename E>
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
	template<typename E>
	auto by_linear(OpExpression<E> const& e);

	//! Split all linear variables from nonlinear variables.
	/*!
	 * Using the available predicates to determine linearity/nonlinearity of
	 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
	 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
	 * for linearity. Similarly, the `NL` part is nonlinear.
	 *
	 * \param e The expression which is split.
	 */
	template<typename T, typename G>
	auto by_linear(OpLVariable<T, G> const& e);

	//! Split all linear variables from nonlinear variables.
	/*!
	 * Using the available predicates to determine linearity/nonlinearity of
	 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
	 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
	 * for linearity. Similarly, the `NL` part is nonlinear.
	 *
	 * \param e The expression which is split.
	 */
	template<typename V, typename E1, typename E2>
	auto by_linear(OpFuncConvolution<V, E1, E2> const& e);

	//! Split all linear variables from nonlinear variables.
	/*!
	 * Using the available predicates to determine linearity/nonlinearity of
	 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
	 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
	 * for linearity. Similarly, the `NL` part is nonlinear.
	 *
	 * \param e The expression which is split.
	 */
	template<typename E1, typename E2>
	auto by_linear(OpBinaryAdd<E1, E2> const& e);

	//! Split all linear variables from nonlinear variables.
	/*!
	 * Using the available predicates to determine linearity/nonlinearity of
	 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
	 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
	 * for linearity. Similarly, the `NL` part is nonlinear.
	 *
	 * \param e The expression which is split.
	 */
	template<typename E1, typename E2>
	auto by_linear(OpBinarySub<E1, E2> const& e);
	

	template<typename E>
	auto by_linear(E e)
	{
		return pack_right(e);
	}

	// derivative

	namespace
	{

		template<typename Dd, typename V, typename E, typename T,
			typename std::enable_if_t<(expression_predicate<E>::linear), int> = 0>
			auto by_linear_derivative(OpFuncDerivative<Dd, V, E, T> const& e)
		{
			return pack_left(e);
		}

		template<typename Dd, typename V, typename E, typename T,
			typename std::enable_if_t<(!expression_predicate<E>::linear && expression_predicate<E>::combination), int> = 0>
			auto by_linear_derivative(OpFuncDerivative<Dd, V, E, T> const& e)
		{
			constexpr size_t order = OpFuncDerivative<Dd, V, E, T>::order;
			return expr::split::by_linear(
				expr::distribute_operators(
					expr::make_operator_derivative<order>(e.value, e.solver) * expr::compound_get::template expr(e)));
		}

		template<typename Dd, typename V, typename E, typename T,
			typename std::enable_if_t<(expression_predicate<E>::nonlinear), int> = 0>
			auto by_linear_derivative(OpFuncDerivative<Dd, V, E, T> const& e)
		{
			return pack_right(e);
		}
	}

	template<typename Dd, typename V, typename E, typename Sp>
	auto by_linear(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return by_linear_derivative(e);
	}


	// op combination

	namespace
	{
		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<(expression_predicate<OpCombination<A1, A2, E>>::linear), int> = 0>
			auto by_linear_combination(OpCombination<A1, A2, E> const& e)
		{
			return pack_left(e);
		}

		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<
			(!expression_predicate<OpCombination<A1, A2, E>>::linear
				&& (!expression_predicate<E>::nonlinear&& expression_predicate<E>::combination)), int> = 0>
			auto by_linear_combination(OpCombination<A1, A2, E> const& e)
		{
			return expr::split::by_linear(expr::distribute_operator(e.combination, expr::compound_get::template expr(e)));
		}

		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<expression_predicate<E>::nonlinear, int> = 0>
			auto by_linear_combination(OpCombination<A1, A2, E> const& e)
		{
			return pack_right(e);
		}
	}

	template<typename A1, typename A2, typename E>
	auto by_linear(OpCombination<A1, A2, E> const& e)
	{
		return by_linear_combination(e);
	}


	// op chain

	namespace
	{
		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<(expression_predicate<OpChain<A1, A2, E>>::linear), int> = 0>
			auto by_linear_chain(OpChain<A1, A2, E> const& e)
		{
			return pack_left(e);
		}

		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<
			(!expression_predicate<OpChain<A1, A2, E>>::linear
				&& expression_predicate<E>::combination), int> = 0>
			auto by_linear_chain(OpChain<A1, A2, E> const& e)
		{
			return expr::split::by_linear(expr::distribute_operator(e.combination, expr::compound_get::template expr(e)));
		}

		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<(expression_predicate<OpChain<A1, A2, E>>::nonlinear), int> = 0>
			auto by_linear_chain(OpChain<A1, A2, E> const& e)
		{
			return pack_right(e);
		}
	}

	template<typename A1, typename A2, typename E>
	auto by_linear(OpChain<A1, A2, E> const& e)
	{
		return by_linear_chain(e);
	}


	// convolution

	namespace
	{

		template<typename V, typename E1, typename E2,
			typename std::enable_if_t<expression_predicate<OpFuncConvolution<V, E1, E2>>::linear, int> = 0>
			auto by_linear_convolution(OpFuncConvolution<V, E1, E2> const& e)
		{
			return pack_left(e);
		}

		template<typename V, typename E1, typename E2,
			typename std::enable_if_t<!expression_predicate<OpFuncConvolution<V, E1, E2>>::linear, int> = 0>
			auto by_linear_convolution(OpFuncConvolution<V, E1, E2> const& e)
		{
			return pack_right(e);
		}

	}

	template<typename V, typename E1, typename E2>
	auto by_linear(OpFuncConvolution<V, E1, E2> const& e)
	{
		return by_linear_convolution(e);
	}

	// handling general expressions (recursion termination) and binary operations (recursive)

	template<typename T, typename G>
	auto by_linear(OpLVariable<T, G> const& e)
	{
		return pack_left(e);
	}



	template<typename E1, typename E2>
	auto by_linear(OpBinaryAdd<E1, E2> const& e)
	{
		auto a = by_linear(e.a);
		auto b = by_linear(e.b);
		return std::make_pair(std::get<0>(a) + std::get<0>(b), std::get<1>(a) + std::get<1>(b));
	}

	template<typename E1, typename E2>
	auto by_linear(OpBinarySub<E1, E2> const& e)
	{
		auto a = by_linear(e.a);
		auto b = by_linear(e.b);
		return std::make_pair(std::get<0>(a) - std::get<0>(b), std::get<1>(a) - std::get<1>(b));
	}

	template<typename E>
	auto by_linear(OpExpression<E> const& e)
	{
		return by_linear(*static_cast<E const*>(&e));
	}



	// **************************************************************************************



	 // there are several different cases a compound operator might have might happen:
	 // 1. E is a function of only V(Z)
	 // 2. E is a function of V(Z) and others, and contains a term that's only of V(Z)
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
	template<size_t Z, typename E1, typename E2>
	auto separate_var(OpBinaryAdd<E1, E2> const& e);

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
	template<size_t Z, typename E1, typename E2>
	auto separate_var(OpBinarySub<E1, E2> const& e);

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
	template<size_t Z, typename E1, typename E2>
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
	template<size_t Z, typename Dd, typename V, typename E, typename Sp>
	auto separate_var(OpFuncDerivative<Dd, V, E, Sp> const& e);

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
	template<size_t Z, typename Dd, typename V, typename G, typename Sp>
	auto separate_var(OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp> const& e);

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
	template<size_t Z, typename A1, typename A2, typename E>
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
	template<size_t Z, typename A1, typename A2, typename E>
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
	template<size_t Z, typename V, typename E1, typename E2>
	auto separate_var(OpFuncConvolution<V, E1, E2> const& e);


	//

	template<size_t Z, typename E>
	auto separate_var(E const& e)
	{
		return pack_right(e);
	}
	
	template<size_t Z, typename T, typename G>
	auto separate_var(OpLVariable<T, Variable<Z, G>> const& e)
	{
		return pack_left(e);
	}

	namespace
	{
		/*
		 * Used in avoiding "if constexpr" constructs.
		 */

		template<size_t Z, typename E1, typename E2>
		constexpr bool svc_pred = (expr::property::vars<E1>::template only_id<Z>() && expr::property::vars<E2>::template only_id<Z>());
		template<size_t Z, typename E>
		constexpr bool svcg_pred = (expr::property::vars<E>::template only_id<Z>());
		template<size_t Z, typename E1, typename E2>
		constexpr bool svm_pred = (expr::property::vars<E1>::template only_id<Z>() && expr::property::vars<E2>::template only_id<Z>());
		template<size_t Z, typename E>
		constexpr bool svd_pred_1 = (expr::property::vars<E>::template only_id<Z>());
		template<size_t Z, typename E>
		constexpr bool svd_pred_2 = (expression_predicate<E>::combination && expr::property::vars<E>::template has_id<Z>());
		template<size_t Z, typename E>
		constexpr bool svcc_pred_1 = (expr::property::vars<E>::template only_id<Z>());
		template<size_t Z, typename E>
		constexpr bool svcc_pred_2 = (expression_predicate<E>::combination && expr::property::vars<E>::template has_id<Z>());

	}

	namespace
	{

		template<size_t Z, typename V, typename E1, typename E2,
			typename std::enable_if_t<svc_pred<Z, E1, E2>, int> = 0>
			auto separate_var_convolution(OpFuncConvolution<V, E1, E2> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, typename V, typename E1, typename E2,
			typename std::enable_if_t<!svc_pred<Z, E1, E2>, int> = 0>
			auto separate_var_convolution(OpFuncConvolution<V, E1, E2> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename V, typename E1, typename E2>
	auto separate_var(OpFuncConvolution<V, E1, E2> const& e)
	{
		return separate_var_convolution<Z>(e);
	}


	namespace
	{
		template<size_t Z, typename V, size_t D, typename E,
			typename std::enable_if_t<svcg_pred<Z, E>, int> = 0>
			auto separate_var_convolution_g(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, typename V, size_t D, typename E,
			typename std::enable_if_t<!svcg_pred<Z, E>, int> = 0>
			auto separate_var_convolution_g(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename V, size_t D, typename E>
	auto separate_var(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e)
	{
		return separate_var_convolution_g<Z>(e);
	}



	template<size_t Z, typename E1, typename E2>
	auto separate_var(OpBinaryAdd<E1, E2> const& e)
	{
		auto a = separate_var<Z>(e.a);
		auto b = separate_var<Z>(e.b);
		return std::make_pair(std::get<0>(a) + std::get<0>(b), std::get<1>(a) + std::get<1>(b));
	}

	template<size_t Z, typename E1, typename E2>
	auto separate_var(OpBinarySub<E1, E2> const& e)
	{
		auto a = separate_var<Z>(e.a);
		auto b = separate_var<Z>(e.b);
		return std::make_pair(std::get<0>(a) - std::get<0>(b), std::get<1>(a) - std::get<1>(b));
	}

	namespace
	{
		template<size_t Z, typename E1, typename E2,
			typename std::enable_if_t<svm_pred<Z, E1, E2>, int> = 0>
			auto separate_var_mul(OpBinaryMul<E1, E2> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, typename E1, typename E2,
			typename std::enable_if_t<!svm_pred<Z, E1, E2>, int> = 0>
			auto separate_var_mul(OpBinaryMul<E1, E2> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename E1, typename E2>
	auto separate_var(OpBinaryMul<E1, E2> const& e)
	{
		return separate_var_mul<Z>(e);
	}


	namespace
	{
		template<size_t Z, typename Dd, typename V, typename E, typename T,
			typename std::enable_if_t<svd_pred_1<Z, E>, int> = 0>
			auto separate_var_derivative(OpFuncDerivative<Dd, V, E, T> const& e)
		{
			return pack_left(e);
		}


		template<size_t Z, typename Dd, typename V, typename E, typename T,
			typename std::enable_if_t<(!svd_pred_1<Z, E>&& svd_pred_2<Z, E>), int> = 0>
			auto separate_var_derivative(OpFuncDerivative<Dd, V, E, T> const& e)
		{
			constexpr size_t order = OpFuncDerivative<Dd, V, E, T>::order;
			return separate_var<Z>(expr::make_operator_derivative<order>(e.value, e.solver) * expr::compound_get::template expr(e));

		}

		template<size_t Z, typename Dd, typename V, typename E, typename T,
			typename std::enable_if_t<!(svd_pred_1<Z, E> || svd_pred_2<Z, E>), int> = 0>
			auto separate_var_derivative(OpFuncDerivative<Dd, V, E, T> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename Dd, typename V, typename E, typename Sp>
	auto separate_var(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return separate_var_derivative<Z>(e);
	}



	namespace
	{
		template<size_t Z, typename Dd, typename V, typename G, typename T,
			typename std::enable_if_t<svd_pred_1<Z, G>, int> = 0>
			auto separate_var_derivative_lop(OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, T> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, typename Dd, typename V, typename G, typename T,
			typename std::enable_if_t<!svd_pred_1<Z, G>, int> = 0>
			auto separate_var_derivative_lop(OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, T> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename Dd, typename V, typename G, typename Sp>
	auto separate_var(OpFuncDerivative<Dd, V, OpLVariable<OpIdentity, G>, Sp> const& e)
	{
		return separate_var_derivative_lop<Z>(e);
	}

	namespace
	{
		template<size_t Z, typename A1, typename A2, typename E,
			typename std::enable_if_t<svcc_pred_1<Z, E>, int> = 0>
			auto separate_var_chain(OpChain<A1, A2, E> const& e)
		{
			return pack_left(e);

		}

		template<size_t Z, typename A1, typename A2, typename E,
			typename std::enable_if_t<(!svcc_pred_1<Z, E>&& svcc_pred_2<Z, E>), int> = 0>
			auto separate_var_chain(OpChain<A1, A2, E> const& e)
		{
			return separate_var<Z>(expr::distribute_operator(e.combination, expr::compound_get::template expr(e)));
		}
	}

	template<size_t Z, typename A1, typename A2, typename E,
		typename std::enable_if_t<!(svcc_pred_1<Z, E> || svcc_pred_2<Z, E>), int> = 0>
	auto separate_var_chain(OpChain<A1, A2, E> const& e)
	{
		return pack_right(e);
	}

	template<size_t Z, typename A1, typename A2, typename E>
	auto separate_var(OpChain<A1, A2, E> const& e)
	{
		return separate_var_chain<Z>(e);
	}


	namespace
	{
		template<size_t Z, typename A1, typename A2, typename E,
			typename std::enable_if_t<svcc_pred_1<Z, E>, int> = 0>
			auto separate_var_combination(OpCombination<A1, A2, E> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, typename A1, typename A2, typename E,
			typename std::enable_if_t<(!svcc_pred_1<Z, E>&& svcc_pred_2<Z, E>), int> = 0>
			auto separate_var_combination(OpCombination<A1, A2, E> const& e)
		{
			return separate_var<Z>(expr::distribute_operator(e.combination, expr::compound_get::template expr(e)));
		}

		template<size_t Z, typename A1, typename A2, typename E,
			typename std::enable_if_t<!(svcc_pred_1<Z, E> || svcc_pred_2<Z, E>), int> = 0>
			auto separate_var_combination(OpCombination<A1, A2, E> const& e)
		{
			return pack_right(e);
		}
	}


	template<size_t Z, typename A1, typename A2, typename E>
	auto separate_var(OpCombination<A1, A2, E> const& e)
	{
		return separate_var_combination<Z>(e);
	}


	// **************************************************************************************

	// factors the given derivative order (note that this is a WHOLE order) and put the operator 
	// in the first part of a pair when derivatives are factored, which will always be returned with 
	// OpIdentity coefficient

	// note that a combination operator derivative can't be factored from a linear combination;
	// inherent limitation of type-deterministic symbolic algebra (same limitation as for factoring)
	// hence why only the whole derivative may be factored


	// if the derivative order is 0 or if the expression has a derivative index which is too
	// small, then return a pair with the identity in place of the factor

	//! Factor out the derivatives from the expression.
	/*!
	 * Factors out the given derivative order. Derivatives of higher orders will
	 * be correspondingly "divided", so that the result which is a pair containing
	 * the expressions `D` and `E` 
	 * will recover the final expression by multiplying `D` to `E`. That is,
	 * `D` will be an expression of only derivative operators, so that the given
	 * expression would be recovered by applying `D` to `E`.
	 *
	 * The derivative order should be less than the minimum order, otherwise
	 * the factorization would be undefined, in which case the original
	 * expression is returned as `E`.
	 * 
	 * \param e The expression which is factored.
	 * 
	 * \tparam OD The order of the derivative that is factored from the expression.
	 */
	template<size_t OD, typename A1, typename A2, typename Enable = typename std::enable_if<min_derivative_order<OpOperatorChain<A1, A2>>::value >= OD>::type>
	auto factor_deriv(OpOperatorCombination<A1, A2> const& e);

	//! Factor out the derivatives from the expression.
	/*!
	 * Factors out the given derivative order. Derivatives of higher orders will
	 * be correspondingly "divided", so that the result which is a pair containing
	 * the expressions `D` and `E`
	 * will recover the final expression by multiplying `D` to `E`. That is,
	 * `D` will be an expression of only derivative operators, so that the given
	 * expression would be recovered by applying `D` to `E`.
	 *
	 * The derivative order should be less than the minimum order, otherwise
	 * the factorization would be undefined, in which case the original
	 * expression is returned as `E`.
	 *
	 * \param e The expression which is factored.
	 *
	 * \tparam OD The order of the derivative that is factored from the expression.
	 */
	template<size_t OD, typename A1, typename A2, typename Enable = typename std::enable_if<min_derivative_order<A1>::value >= OD>::type>
	auto factor_deriv(OpOperatorChain<A1, A2> const& e);

	//! Factor out the derivatives from the expression.
	/*!
	 * Factors out the given derivative order. Derivatives of higher orders will
	 * be correspondingly "divided", so that the result which is a pair containing
	 * the expressions `D` and `E`
	 * will recover the final expression by multiplying `D` to `E`. That is,
	 * `D` will be an expression of only derivative operators, so that the given
	 * expression would be recovered by applying `D` to `E`.
	 *
	 * The derivative order should be less than the minimum order, otherwise
	 * the factorization would be undefined, in which case the original
	 * expression is returned as `E`.
	 *
	 * \param e The expression which is factored.
	 *
	 * \tparam OD The order of the derivative that is factored from the expression.
	 */
	template<size_t OD, typename A1, typename A2, typename E>
	auto factor_deriv(OpCombination<A1, A2, E> const& e);

	//! Factor out the derivatives from the expression.
	/*!
	 * Factors out the given derivative order. Derivatives of higher orders will
	 * be correspondingly "divided", so that the result which is a pair containing
	 * the expressions `D` and `E`
	 * will recover the final expression by multiplying `D` to `E`. That is,
	 * `D` will be an expression of only derivative operators, so that the given
	 * expression would be recovered by applying `D` to `E`.
	 *
	 * The derivative order should be less than the minimum order, otherwise
	 * the factorization would be undefined, in which case the original
	 * expression is returned as `E`.
	 *
	 * \param e The expression which is factored.
	 *
	 * \tparam OD The order of the derivative that is factored from the expression.
	 */
	template<size_t OD, typename A1, typename A2, typename E>
	auto factor_deriv(OpChain<A1, A2, E> const& e);

	//! Factor out the derivatives from the expression.
	/*!
	 * Factors out the given derivative order. Derivatives of higher orders will
	 * be correspondingly "divided", so that the result which is a pair containing
	 * the expressions `D` and `E`
	 * will recover the final expression by multiplying `D` to `E`. That is,
	 * `D` will be an expression of only derivative operators, so that the given
	 * expression would be recovered by applying `D` to `E`.
	 *
	 * The derivative order should be less than the minimum order, otherwise
	 * the factorization would be undefined, in which case the original
	 * expression is returned as `E`.
	 *
	 * \param e The expression which is factored.
	 *
	 * \tparam OD The order of the derivative that is factored from the expression.
	 */
	template<size_t OD, typename E1, typename E2, typename Enable = typename std::enable_if<(min_derivative_order<E1>::value >= OD && min_derivative_order<E2>::value >= OD)>::type>
	auto factor_deriv(OpBinaryAdd<E1, E2> const& e);

	//! Factor out the derivatives from the expression.
	/*!
	 * Factors out the given derivative order. Derivatives of higher orders will
	 * be correspondingly "divided", so that the result which is a pair containing
	 * the expressions `D` and `E`
	 * will recover the final expression by multiplying `D` to `E`. That is,
	 * `D` will be an expression of only derivative operators, so that the given
	 * expression would be recovered by applying `D` to `E`.
	 *
	 * The derivative order should be less than the minimum order, otherwise
	 * the factorization would be undefined, in which case the original
	 * expression is returned as `E`.
	 *
	 * \param e The expression which is factored.
	 *
	 * \tparam OD The order of the derivative that is factored from the expression.
	 */
	template<size_t OD, typename E1, typename E2, typename Enable = typename std::enable_if<(min_derivative_order<E1>::value >= OD && min_derivative_order<E2>::value >= OD)>::type>
	auto factor_deriv(OpBinarySub<E1, E2> const& e);


	template<size_t OD, size_t O, typename V, typename Sp, typename std::enable_if<(O == OD), int>::type = 0>
	auto factor_deriv(OpOperatorDerivative<O, V, Sp> const& e)
	{
		return std::make_pair(expr::make_operator_derivative<OD>(e.solver), expr::make_literal(e.value));
	}

	template<size_t OD, size_t O, typename V, typename Sp, typename std::enable_if<(O > OD), int>::type = 0>
	auto factor_deriv(OpOperatorDerivative<O, V, Sp> const& e)
	{
		return std::make_pair(expr::make_operator_derivative<OD>(e.solver), expr::make_operator_derivative<O - OD>(e.value, e.solver));
	}

	template<size_t OD, typename Dd, typename V, typename E, typename Sp, 
		size_t O = OpFuncDerivative<Dd, V, E, Sp>::order,
		typename std::enable_if<(O == OD), int>::type = 0>
	auto factor_deriv(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return std::make_pair(expr::make_operator_derivative<OD>(e.solver), expr::make_literal(e.value) * expr::compound_get::template expr(e));
	}

	template<size_t OD, typename Dd, typename V, typename E, typename Sp, 
		size_t O = OpFuncDerivative<Dd, V, E, Sp>::order,
		typename std::enable_if_t<(O != OD), int>::type = 0>
	auto factor_deriv(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return std::make_pair(expr::make_operator_derivative<OD>(e.solver), expr::make_operator_derivative<O - OD>(e.value, e.solver) * expr::compound_get::template expr(e));
	}


	//! Factor out the derivatives from the expression.
	/*!
	 * Factors out the given derivative order. Derivatives of higher orders will
	 * be correspondingly "divided", so that the result which is a pair containing
	 * the expressions `D` and `E`
	 * will recover the final expression by multiplying `D` to `E`. That is,
	 * `D` will be an expression of only derivative operators, so that the given
	 * expression would be recovered by applying `D` to `E`.
	 *
	 * The derivative order should be less than the minimum order, otherwise
	 * the factorization would be undefined, in which case the original
	 * expression is returned as `E`.
	 *
	 * \param e The expression which is factored.
	 *
	 * \tparam OD The order of the derivative that is factored from the expression.
	 */
	template<size_t OD, typename E, typename std::enable_if<(OD == 0 || min_derivative_order<E>::value < OD), int>::type = 0>
	auto factor_deriv(OpExpression<E> const& e)
	{
		return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
	}


	template<size_t OD, typename A1, typename A2, typename Enable>
	auto factor_deriv(OpOperatorCombination<A1, A2> const& e)
	{
		auto a = factor_deriv<OD>(e.f);
		auto b = factor_deriv<OD>(e.g);
		return std::make_pair(a.first, a.second + b.second);
	}

	template<size_t OD, typename A1, typename A2, typename Enable>
	auto factor_deriv(OpOperatorChain<A1, A2> const& e)
	{
		auto a = factor_deriv<OD>(e.f);
		return std::make_pair(a.first, a.second * e.g);
	}

	template<size_t OD, typename A1, typename A2, typename E>
	auto factor_deriv(OpCombination<A1, A2, E> const& e)
	{
		auto fac = factor_deriv<OD>(e.combination);
		return std::make_pair(fac.first, fac.second * expr::compound_get::template expr(e));
	}

	template<size_t OD, typename A1, typename A2, typename E>
	auto factor_deriv(OpChain<A1, A2, E> const& e)
	{
		auto fac = factor_deriv<OD>(e.combination);
		return std::make_pair(fac.first, fac.second * expr::compound_get::template expr(e));
	}

	template<size_t OD, typename E1, typename E2, typename Enable>
	auto factor_deriv(OpBinaryAdd<E1, E2> const& e)
	{
		auto a = factor_deriv<OD>(e.a);
		auto b = factor_deriv<OD>(e.b);
		return std::make_pair(a.first, a.second + b.second);
	}


	template<size_t OD, typename E1, typename E2, typename Enable>
	auto factor_deriv(OpBinarySub<E1, E2> const& e)
	{
		auto a = factor_deriv<OD>(e.a);
		auto b = factor_deriv<OD>(e.b);
		return std::make_pair(a.first, a.second - b.second);
	}



	//! Factor out the derivatives from the expression.
	/*!
	 * Factors out the given derivative order. Derivatives of higher orders will
	 * be correspondingly "divided", so that the result which is a pair containing
	 * the expressions `D` and `E`
	 * will recover the final expression by multiplying `D` to `E`. That is,
	 * `D` will be an expression of only derivative operators, so that the given
	 * expression would be recovered by applying `D` to `E`.
	 *
	 * The derivative order must be less than the minimum order, otherwise
	 * the factorization would be undefined.
	 *
	 * \param e The expression which is factored.
	 *
	 * \tparam OD The order of the derivative that is factored from the expression.
	 */
	template<typename E>
	auto factor_deriv(OpExpression<E> const& e)
	{
		return factor_deriv<min_derivative_order<E>::value>(*static_cast<E const*>(&e));
	}


	// **************************************************************************************


	namespace
	{

		template<size_t O, typename E, typename std::enable_if_t<(O == 0 && min_derivative_order<E>::value == 0), int> = 0>
		auto separate_deriv(OpExpression<E> const& e)
		{
			return pack_left(*static_cast<E const*>(&e));
		}

		template<size_t O, typename E, typename std::enable_if_t<!(O == 0 && min_derivative_order<E>::value == 0), int> = 0>
		auto separate_deriv(OpExpression<E> const& e)
		{
			return pack_right(*static_cast<E const*>(&e));
		}

		template<size_t O, typename Dd, typename V, typename E, typename Sp, 
			typename std::enable_if<_check_equality_derivative<O, OpFuncDerivative<Dd, V, E, Sp>>(), int>::type = 0>
		auto separate_deriv(OpFuncDerivative<Dd, V, E, Sp> const& e)
		{
			return pack_left(e);
		}

		template<size_t O, typename A1, typename A2, typename E,
			typename std::enable_if<min_derivative_order<OpOperatorCombination<A1, A2>>::value == O, int>::type = 0>
		auto separate_deriv(OpCombination<A1, A2, E> const& e)
		{
			return pack_left(e);
		}

		template<size_t O, typename A1, typename A2, typename E, 
			typename std::enable_if<min_derivative_order<OpOperatorChain<A1, A2>>::value == O, int>::type = 0>
		auto separate_deriv(OpChain<A1, A2, E> const& e)
		{
			return pack_left(e);
		}

		template<size_t O, typename E1, typename E2>
		auto separate_deriv(OpBinaryAdd<E1, E2> const& e)
		{
			auto [d_i_a, rest_a] = separate_deriv<O>(e.a);
			auto [d_i_b, rest_b] = separate_deriv<O>(e.b);
			return std::make_pair(d_i_a + d_i_b, rest_a + rest_b);
		}

		template<size_t O, typename E1, typename E2>
		auto separate_deriv(OpBinarySub<E1, E2> const& e)
		{
			auto [d_i_a, rest_a] = separate_deriv<O>(e.a);
			auto [d_i_b, rest_b] = separate_deriv<O>(e.b);
			return std::make_pair(d_i_a - d_i_b, rest_a - rest_b);
		}



		template<size_t O, typename E>
		auto _separate_deriv(OpExpression<E> const& e);
		
		template<size_t O, size_t next_order, typename Dd, typename R, typename std::enable_if_t<(next_order > O), int> = 0>
		auto _separate_deriv_return(Dd&& d_i, R&& rest)
		{
			return std::tuple_cat(
				std::make_tuple(factor_deriv<O>(std::forward<Dd>(d_i))), 
				_separate_deriv<next_order>(std::forward<R>(rest)));
		}

		template<size_t O, size_t next_order, typename Dd, typename R, typename std::enable_if_t<(next_order <= O), int> = 0>
		auto _separate_deriv_return(Dd&& d_i, R&&)
		{
			return std::make_tuple(factor_deriv<O>(std::forward<Dd>(d_i)));
		}

		template<size_t O, typename E>
		auto _separate_deriv(OpExpression<E> const& e)
		{
			auto [d_i, rest] = separate_deriv<O>(*static_cast<E const*>(&e));
			static_assert(!std::is_same<decltype(d_i), OpVoid>::value);
			constexpr size_t next_order = min_derivative_order<decltype(rest)>::value;
			return _separate_deriv_return<O, next_order>(d_i, rest);
		}


	}

	//! Separate an expression by derivative.
	/*! 
	 * Nothing can be separated from the zero value, and this will
	 * immediately return an empty list.
	 */
	inline auto separate_deriv(OpVoid)
	{
		return std::make_tuple();
	}


	//! Separate an expression by derivative.
	/*!
	 * Separate an expression by the derivative, starting at the derivative with
	 * the lowest index. The entry point will always assume there's nothing lower
	 * than the order given.
	 *
	 * Each operator corresponding to the derivative index will be put as the
	 * first member of a pair, with the expression that it operates on as the
	 * second member of that pair; subsequently all the operators are
	 * recursively determined and each pair is put together as part of a tuple.
	 *
	 * Concisely, this algorithm groups expressions of different operators and
	 * creates a list of pairs of `D` and `E`. The original expression is
	 * recovered by taking `D * E` of all pairs and adding the results together.
	 *
	 * The derivative searching is **not** a greedy algorithm, meaning constituents
	 * of an operator chain will be pruned off from the rest of the chain. This
	 * means that if there is a term such that
	 * some derivative `d2` is applied to another derivative `d4`, and the
	 * algorithm is searching for `d2`, it would accept this term, hence
	 * the `E` part of the pair will still contain a term with a derivative.
	 * 
	 * \param e The expression that is split.
	 */
	template<typename E>
	auto separate_deriv(OpExpression<E> const& e)
	{
		return _separate_deriv<min_derivative_order<E>::value>(*static_cast<E const*>(&e));
	}


	// **************************************************************************************

	namespace
	{

		// factors the expression by the given variable
		// the oplvariables and opnlvariables with matching type C are factored
		// N is the number of types of C to factor out

		template<typename C>
		struct divide_variable
		{
			template<typename T, typename G, typename std::enable_if<(expr::factor_count<C, G>::value >= 0), int>::type = 0>
			auto operator()(OpLVariable<T, G> const& e)
			{
				return std::make_pair(
					expr::make_op(typename expr::original_data_type<C>::type(expr::BaseData<G>::get(e.data))), 
					e / expr::make_op(typename expr::original_data_type<C>::type(expr::BaseData<G>::get(e.data))));
			}
		};

		template<size_t Z>
		struct divide_variable<std::index_sequence<Z>>
		{
			template<typename T, typename G>
			auto operator()(OpLVariable<T, Variable<Z, G>> const& e)
			{
				return std::make_pair(OpLVariable<OpIdentity, Variable<Z, G>>(e.data), expr::make_literal(e.value));
			}
		};


		template<size_t N, typename C, typename E1, typename E2, typename Enable = typename std::enable_if<(expr::factor_count<C, OpBinaryAdd<E1, E2>>::value >= N)>::type>
		auto _factor(OpBinaryAdd<E1, E2> const& e);
		template<size_t N, typename C, typename E1, typename E2, typename Enable = typename std::enable_if<(expr::factor_count<C, OpBinarySub<E1, E2>>::value >= N)>::type>
		auto _factor(OpBinarySub<E1, E2> const& e);
		template<size_t N, typename C, typename E1, typename E2, typename Enable = typename std::enable_if<(expr::factor_count<C, OpBinaryMul<E1, E2>>::value >= N)>::type>
		auto _factor(OpBinaryMul<E1, E2> const& e);
		template<size_t N, typename C, typename E1, typename E2, typename Enable = typename std::enable_if<(expr::factor_count<C, OpBinaryDiv<E1, E2>>::value >= N)>::type>
		auto _factor(OpBinaryDiv<E1, E2> const& e);


		template<size_t N, typename C, typename E, typename std::enable_if<(N == 0 || expr::factor_count<C, E>::value < N), int>::type = 0>
		auto _factor(OpExpression<E> const& e)
		{
			return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
		}

		template<size_t N, typename C, typename T, typename G, typename std::enable_if<(N > 0 && expr::factor_count<C, G>::value >= N), int>::type = 0>
		auto _factor(OpLVariable<T, G> const& e)
		{
			return divide_variable<C>{}(e);
		}


		template<size_t N, typename C, typename T, typename G0, typename G1,
			typename std::enable_if<(expr::factor_count<C, G0>::value == 1), int>::type = 0>
			auto _factor_nl_ge2(OpNLVariable<T, G0, G1> const& e)
		{
			auto a = _factor<1, C>(OpLVariable<OpIdentity, G0>(std::get<0>(e.datas)));
			auto b = _factor<N - 1, C>(OpLVariable<OpIdentity, G1>(std::get<1>(e.datas)));
			return std::make_pair(a.first * b.first, expr::make_literal(e.value) * a.second * b.second);
		}

		template<size_t N, typename C, typename T, typename G0, typename G1,
			typename std::enable_if<(expr::factor_count<C, G0>::value != 1), int>::type = 0>
			auto _factor_nl_ge2(OpNLVariable<T, G0, G1> const& e)
		{
			auto b = _factor<1, C>(OpLVariable<OpIdentity, G1>(std::get<1>(e.datas)));
			return std::make_pair(b.first, expr::make_literal(e.value) * OpLVariable<OpIdentity, G0>(std::get<0>(e.datas)) * b.second);
		}

		template<size_t N, typename C, typename T, typename G0, typename G1, 
			typename std::enable_if<(N >= 2), int>::type = 0>
		auto _factor(OpNLVariable<T, G0, G1> const& e)
		{
			return _factor_nl_ge2<N, C>(e);
		}



		template<size_t N, typename C, typename T, typename G0, typename G1,
			typename std::enable_if<(expr::factor_count<C, G0>::value == 1), int>::type = 0>
			auto _factor_nl_lt2(OpNLVariable<T, G0, G1> const& e)
		{
			auto a = _factor<1, C>(OpLVariable<OpIdentity, G0>(std::get<0>(e.datas)));
			return std::make_pair(a.first, expr::make_literal(e.value) * OpLVariable<OpIdentity, G1>(std::get<1>(e.datas)) * a.second);
		}

		template<size_t N, typename C, typename T, typename G0, typename G1,
			typename std::enable_if<(expr::factor_count<C, G0>::value != 1), int>::type = 0>
			auto _factor_nl_lt2(OpNLVariable<T, G0, G1> const& e)
		{
			auto b = _factor<1, C>(OpLVariable<OpIdentity, G1>(std::get<1>(e.datas)));
			return std::make_pair(b.first, expr::make_literal(e.value) * OpLVariable<OpIdentity, G0>(std::get<0>(e.datas)) * b.second);
		}

		template<size_t N, typename C, typename T, typename G0, typename G1, 
			typename std::enable_if<(N < 2), int>::type = 0>
		auto _factor(OpNLVariable<T, G0, G1> const& e)
		{
			return _factor_nl_lt2<N, C>(e);
		}


		template<size_t N, typename C, typename T, typename G0, typename G1, typename... Gs, 
			typename std::enable_if<(expr::factor_count<C, G0>::value == 1), int>::type = 0>
		auto _factor_nl(OpNLVariable<T, G0, G1, Gs...> const& e)
		{
			auto a = _factor<1, C>(OpLVariable<OpIdentity, G0>(std::get<0>(e.datas)));
			auto b = _factor<N - 1, C>(OpNLVariable<OpIdentity, G1, Gs...>(symphas::lib::get_tuple_ge<1>(e.datas)));
			return std::make_pair(a.first * b.first, expr::make_literal(e.value) * a.second * b.second);
		}
		
		template<size_t N, typename C, typename T, typename G0, typename G1, typename... Gs, 
			typename std::enable_if<(expr::factor_count<C, G0>::value != 1), int>::type = 0>
		auto _factor_nl(OpNLVariable<T, G0, G1, Gs...> const& e)
		{
			auto b = _factor<N, C>(OpNLVariable<OpIdentity, G1, Gs...>(symphas::lib::get_tuple_ge<1>(e.datas)));
			return std::make_pair(b.first, expr::make_literal(e.value) * OpLVariable<OpIdentity, G0>(std::get<0>(e.datas)) * b.second);
		}

		template<size_t N, typename C, typename T, typename G0, typename G1, typename... Gs, typename std::enable_if<(N > 0), int>::type = 0>
		auto _factor(OpNLVariable<T, G0, G1, Gs...> const& e)
		{
			return _factor_nl<N, C>(e);
		}

		template<size_t N, typename C, typename E1, typename E2, typename Enable>
		auto _factor(OpBinaryAdd<E1, E2> const& e)
		{
			auto a = _factor<N, C>(e.a);
			auto b = _factor<N, C>(e.b);
			return std::make_pair(a.first, a.second + b.second);
		}

		template<size_t N, typename C, typename E1, typename E2, typename Enable>
		auto _factor(OpBinarySub<E1, E2> const& e)
		{
			auto a = _factor<N, C>(e.a);
			auto b = _factor<N, C>(e.b);
			return std::make_pair(a.first, a.second - b.second);
		}




		template<size_t N, typename C, typename E1, typename E2, 
			size_t AN = expr::factor_count<C, E1>::value, typename std::enable_if_t<(N - AN > 0), int> = 0>
		auto _factor_sift(OpBinaryMul<E1, E2> const& e)
		{
			auto a = _factor<AN, C>(e.a);
			auto b = _factor<N - AN, C>(e.b);
			return std::make_pair(a.first * b.first, a.second * b.second);
		}
		template<size_t N, typename C, typename E1, typename E2,
			size_t AN = expr::factor_count<C, E1>::value, typename std::enable_if_t<(N - AN <= 0), int> = 0>
			auto _factor_sift(OpBinaryMul<E1, E2> const& e)
		{
			auto a = _factor<N, C>(e.a);
			return std::make_pair(a.first, a.second * e.b);
		}

		template<size_t N, typename C, typename E1, typename E2, typename Enable>
		auto _factor(OpBinaryMul<E1, E2> const& e)
		{
			/* this algorithm takes into account the edge case that the same type of variable
			 * is separated by a strict multiplication, although in general, a multiplication rule
			 * will put them together into an nlvariable
			 */
			return _factor_sift<N, C>(e);
		}


		template<size_t N, typename C, typename E1, typename E2, typename Enable>
		auto _factor(OpBinaryDiv<E1, E2> const& e)
		{
			auto a = _factor<N, C>(e.a);
			return std::make_pair(a.first, a.second / e.b);
		}



		template<typename C>
		struct is_variable_data_factor
		{
			static const bool value = false;
		};

		template<size_t Z>
		struct is_variable_data_factor<expr::variable_data_factor<Z>>
		{
			static const bool value = true;
		};

		// automatically selects the largest possible number of datas to factor
		template<typename C0, typename E, typename std::enable_if<(expr::grid_can_combine<C0>::value || is_variable_data_factor<C0>::value), int>::type = 0>
		auto _factor(OpExpression<E> const& e)
		{
			return _factor<expr::factor_count<C0, E>::value, C0>(*static_cast<E const*>(&e));
		}

		template<typename C0, typename E, typename std::enable_if<(!expr::grid_can_combine<C0>::value && !is_variable_data_factor<C0>::value), int>::type = 0>
		auto _factor(OpExpression<E> const& e)
		{
			return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
		}

		template<typename C0, typename... Cs>
		struct factor_pred
		{
			static const bool value = ((std::is_same<C0, Cs>::value || ...));
		};

		template<typename C0>
		struct factor_pred<C0>
		{
			static const bool value = false;
		};
	}


	//! Factor an expression by the given terms.
	/*!
	 * Factor an expression by the given terms which are represented by
	 * the given types, not passed explicitly.
	 *
	 * The types to factor by must all be unique. This produces a maximally
	 * factored expression where all the types have been potentially factored.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam C0 The first term to factor.
	 * \tparam Cs... The remaining terms to factor.
	 */
	template<typename C0, typename... Cs, typename E, typename std::enable_if_t<(sizeof...(Cs) == 0 && !factor_pred<C0>::value), int> = 0>
	auto factor(OpExpression<E> const& e)
	{
		auto a = _factor<C0>(*static_cast<const E*>(&e));
		return std::make_pair(a.first, a.second);
	}

	//! Factor an expression by the given terms.
	/*!
	 * Factor an expression by the given terms which are represented by
	 * the given types, not passed explicitly.
	 *
	 * The types to factor by must all be unique. This produces a maximally
	 * factored expression where all the types have been potentially factored.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam C0 The first term to factor.
	 * \tparam Cs... The remaining terms to factor.
	 */
	template<typename C0, typename... Cs, typename E, typename std::enable_if_t<(sizeof...(Cs) > 0 && !factor_pred<C0, Cs...>::value), int> = 0>
	auto factor(OpExpression<E> const& e)
	{
		auto a = _factor<C0>(*static_cast<const E*>(&e));
		auto b = factor<Cs...>(a.second);
		return std::make_pair(a.first * b.first, b.second);
	}

	//! Factor an expression by the given terms the given number of times.
	/*!
	 * Factor an expression by the given terms which are represented by
	 * the given types, not passed explicitly. Each term is factored out
	 * the given number of times.
	 *
	 * The types to factor by must all be unique. 
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam N The number of times to factor each type out.
	 * \tparam C0 The first term to factor.
	 * \tparam Cs... The remaining terms to factor.
	 */
	template<size_t N, typename C0, typename... Cs, typename E, typename std::enable_if_t<(sizeof...(Cs) == 0 && !factor_pred<C0>::value), int> = 0>
	auto factor(OpExpression<E> const& e)
	{
		constexpr size_t min_order = fixed_min<expr::factor_count<C0, E>::value, N>;
		auto a = _factor<min_order, C0>(*static_cast<const E*>(&e));
		return std::make_pair(a.first, a.second);
	}

	//! Factor an expression by the given terms the given number of times.
	/*!
	 * Factor an expression by the given terms which are represented by
	 * the given types, not passed explicitly. Each term is factored out
	 * the given number of times.
	 *
	 * The types to factor by must all be unique.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam N The number of times to factor each type out.
	 * \tparam C0 The first term to factor.
	 * \tparam Cs... The remaining terms to factor.
	 */
	template<size_t N, typename C0, typename... Cs, typename E, typename std::enable_if_t<(sizeof...(Cs) > 0 && !factor_pred<C0, Cs...>::value), int> = 0>
	auto factor(OpExpression<E> const& e)
	{
		constexpr size_t min_order = fixed_min<expr::factor_count<C0, E>::value, N>;
		auto a = _factor<min_order, C0>(*static_cast<const E*>(&e));
		auto b = factor<N, Cs...>(a.second);
		return std::make_pair(a.first * b.first, b.second);
	}


	//! Make a list of the result of each factorization.
	/*!
	 * For each given term, the original given expression is factored. The
	 * results are concatenated into a list.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam C0 The first term to factor.
	 * \tparam Cs... The remaining terms to factor.
	 */
	template<typename C0, typename... Cs, typename E, typename std::enable_if_t<(sizeof...(Cs) == 0 && !factor_pred<C0>::value), int> = 0>
	auto factor_tuple(OpExpression<E> const& e)
	{
		return std::make_tuple(factor<C0>(*static_cast<const E*>(&e)));
	}

	//! Make a list of the result of each factorization.
	/*!
	 * For each given term, the original given expression is factored. The
	 * results are concatenated into a list.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam C0 The first term to factor.
	 * \tparam Cs... The remaining terms to factor.
	 */
	template<typename C0, typename... Cs, typename E, typename std::enable_if_t<(sizeof...(Cs) > 0 && !factor_pred<C0, Cs...>::value), int> = 0>
	auto factor_tuple(OpExpression<E> const& e)
	{
		auto a = factor<C0>(*static_cast<const E*>(&e));
		auto b = factor_tuple<Cs...>(e);
		return std::tuple_cat(a, b);
	}

	//! Keep factoring the given expression once by all terms.
	/*!
	 * The expression is continuously factored by all terms once, and the
	 * new factored expression is returned as the second element of a pair,
	 * where the first element is the product of all the terms
	 * that could be factored.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam C0 The first term to factor.
	 * \tparam Cs... The remaining terms to factor.
	 */
	template<typename C0, typename... Cs, typename E, typename std::enable_if_t<(sizeof...(Cs) == 0), int> = 0>
	auto factor_list(OpExpression<E> const& e)
	{
		constexpr size_t min_order = fixed_min<expr::factor_count<C0, E>::value, 1>;
		auto a = _factor<min_order, C0>(*static_cast<const E*>(&e));
		return std::make_pair(a.first, a.second);
	}

	//! Keep factoring the given expression once by all terms.
	/*!
	 * The expression is continuously factored by all terms once, and the
	 * new factored expression is returned as the second element of a pair,
	 * where the first element is the product of all the terms
	 * that could be factored.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam C0 The first term to factor.
	 * \tparam Cs... The remaining terms to factor.
	 */
	template<typename C0, typename... Cs, typename E, typename std::enable_if_t<(sizeof...(Cs) > 0), int> = 0>
	auto factor_list(OpExpression<E> const& e)
	{
		constexpr size_t min_order = fixed_min<expr::factor_count<C0, E>::value, 1>;

		auto a = _factor<min_order, C0>(*static_cast<const E*>(&e));
		auto b = factor_list<Cs...>(a.second);
		return std::make_pair(a.first * b.first, b.second);
	}



	//! Keep factoring the given expression once by all terms.
	/*!
	 * The expression is continuously factored by all variable indices once,
	 * and the new factored expression is returned as the second element of a pair,
	 * where the first element is the product of all the terms
	 * that could be factored.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam Z0 The first variable index to factor.
	 * \tparam Cs... The remaining variable indices to factor.
	 */
	template<size_t Z0, size_t... Zs, typename E, typename std::enable_if_t<(sizeof...(Zs) == 0), int> = 0>
	auto factor_list(OpExpression<E> const& e)
	{
		auto a = _factor<expr::variable_data_factor<Z0>>(*static_cast<const E*>(&e));
		return std::make_pair(a.first, a.second);
	}


	//! Keep factoring the given expression once by all terms.
	/*!
	 * The expression is continuously factored by all variable indices once,
	 * and the new factored expression is returned as the second element of a pair,
	 * where the first element is the product of all the terms
	 * that could be factored.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam Z0 The first variable index to factor.
	 * \tparam Cs... The remaining variable indices to factor.
	 */
	template<size_t Z0, size_t... Zs, typename E, typename std::enable_if_t<(sizeof...(Zs) > 0), int> = 0>
	auto factor_list(OpExpression<E> const& e)
	{
		auto a = _factor<expr::variable_data_factor<Z0>>(*static_cast<const E*>(&e));
		auto b = factor_list<Zs...>(a.second);
		return std::make_pair(a.first * b.first, b.second);
	}


}



namespace expr
{



	//! Applies division using the provided types as factors.
	/*!
	 * In the general case, no factoring is performed, but when some factors are
	 * passed (and present), then the division will factor those provided elements.
	 * The given types are provided in a tuple and have to be factors. If no
	 * tuple is provided, a regular division will be returned.
	 */
	template<typename G>
	struct divide_with_factors
	{

		//! Returns the division term between two expressions.
		template<typename E1, typename E2>
		auto operator()(OpExpression<E1> const& a, OpExpression<E2> const& b)
		{
			return expr::make_div(*static_cast<E1 const*>(a), *static_cast<E2 const*>(b));
		}
	};

	//! Specialization of expr::divide_with_factors with no factors given.
	template<>
	struct divide_with_factors<std::tuple<>>
	{
		template<typename E1, typename E2>
		auto operator()(OpExpression<E1> const& a, OpExpression<E2> const& b)
		{
			return expr::make_div(a, b);
		}
	};

	//! Specialization of expr::divide_with_factors a list of factors given.
	template<size_t N, typename G0, typename... Ts>
	struct divide_with_factors<std::tuple<std::pair<std::index_sequence<N>, G0>, Ts...>>
	{
		template<typename E1, typename E2>
		auto operator()(OpExpression<E1> const& a, OpExpression<E2> const& b)
		{
			auto f = expr::split::factor<N, G0>(a);
			auto g = expr::split::factor<N, G0>(b);
			return divide_with_factors<std::tuple<Ts...>>::get(f.second, g.second);
		}
	};


}


//! The division operator is overloaded to apply factoring in general.
template<typename E1, typename E2, typename std::enable_if<expr::factor_list_all<E1, E2>::value, int>::type = 0>
auto operator/(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::divide_with_factors<typename expr::factor_list_all<E1, E2>::type>{}(a, b);
}

//! The division operator is overloaded to apply factoring in general.
template<typename E1, typename E2, typename std::enable_if<!expr::factor_list_all<E1, E2>::value, int>::type = 0>
auto operator/(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_div(a, b);
}





template<typename E, typename T, typename G, typename std::enable_if<(expr::factor_count<G, E>::value > 0), int>::type = 0>
auto operator/(OpExpression<E> const& a, OpLVariable<T, G> const& b)
{
	auto f = expr::split::factor<1, G>(a);
	return f.second / expr::make_literal(b.value);
}

template<typename E, typename T, typename G, typename std::enable_if<(expr::factor_count<G, E>::value > 0), int>::type = 0>
auto operator/(OpLVariable<T, G> const& a, OpExpression<E> const& b)
{
	auto f = expr::split::factor<1, G>(b);
	return expr::make_literal(a.value) / f.second;
}

template<typename E, typename T, typename... Gs, 
	size_t N = ((expr::factor_count<Gs, E>::value + ...)), typename std::enable_if<(N > 0), int>::type = 0>
auto operator/(OpExpression<E> const& a, OpNLVariable<T, Gs...> const& b)
{
	auto f = expr::split::factor_list<Gs...>(a);
	return f.second / (b / f.first);
}

template<typename E, typename T, typename... Gs, 
	size_t N = ((expr::factor_count<Gs, E>::value + ...)), typename std::enable_if<(N > 0), int>::type = 0>
auto operator/(OpNLVariable<T, Gs...> const& a, OpExpression<E> const& b)
{
	auto f = expr::split::factor_list<Gs...>(b);
	return (a / f.first) / f.second;
}


template<typename T1, typename T2, typename... G1s, typename G2, 
	size_t N = ((expr::factor_count<G2, G1s>::value + ...)), typename std::enable_if<(N > 0), int>::type = 0>
auto operator/(OpNLVariable<T1, G1s...> const& a, OpLVariable<T2, G2> const&)
{
	auto f = expr::split::factor<1, G2>(a);
	return f.second;
}

template<typename T1, typename T2, typename G1, typename... G2s, 
	size_t N = ((expr::factor_count<G1, G2s>::value + ...)), typename std::enable_if<(N > 0), int>::type = 0>
auto operator/(OpLVariable<T1, G1> const&, OpNLVariable<T2, G2s...> const& b)
{
	auto f = expr::split::factor<1, G1>(b);
	return expr::inverse(f.second);
}


template<typename T1, typename T2, typename... G1s, typename... G2s, 
	size_t N = ((expr::factor_count_list<G1s, G2s...>::value + ...)), typename std::enable_if<(N > 0), int>::type = 0>
auto operator/(OpNLVariable<T1, G1s...> const& a, OpNLVariable<T2, G2s...> const& b)
{
	using factor_list = typename expr::grid_types<OpNLVariable<T2, G2s...>>::type;
	auto f = expr::split::factor_list<G2s...>(a);
	auto g = expr::split::factor_list<G1s...>(b);
	return f.second / g.second;
}










