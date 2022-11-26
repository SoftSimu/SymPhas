
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

#include "expressionrules.h"
#include "expressionconvolution.h"
#include "expressionderivatives.h"
#include "expressionexponentials.h"
#include "expressiontypeincludes.h"
#include "expressionfunctions.h"
#include "expressionproperties.h"

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
		void _update(OpExpression<E>& e) {}

		template<typename V, typename E>
		void _update(OpExponential<V, E>& e);
		template<typename V, typename E1, typename E2>
		void _update(OpFuncConvolution<V, E1, E2>& e);
		template<typename V, typename E, typename F, typename... Args>
		void _update(OpFunc<V, E, F, Args...>& e);
		template<typename V, size_t D, typename E>
		void _update(OpFuncConvolution<V, GaussianSmoothing<D>, E>& e);
		template<typename V, size_t D, typename G>
		void _update(OpFuncConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>& e);
		template<typename Dd, typename V, typename G, typename Sp>
		void _update(OpFuncDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>&);
		template<typename Dd, typename V, typename E, typename Sp>
		void _update(OpFuncDerivative<Dd, V, E, Sp>& e);
		template<typename... Es>
		void _update(OpAdd<Es...>& e);
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
			_update(expr::get_enclosed_expression(e));
		}

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
			_update(expr::get_enclosed_expression(e));
		}

		template<typename V, size_t D, typename E>
		void _update(OpFuncConvolution<V, GaussianSmoothing<D>, E>& e)
		{
			_update(expr::get_enclosed_expression(e));
			e.update();
		}

		template<typename V, size_t D, typename G>
		void _update(OpFuncConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>& e)
		{
			e.update();
		}


		/* derivative pruning
		 */

		template<size_t O, typename V, typename E, typename G>
		void _update(OpFuncDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G>>& e) {}

		template<typename Dd, typename V, typename G, typename Sp>
		void _update(OpFuncDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>&) {}

		template<typename Dd, typename V, typename E, typename Sp>
		void _update(OpFuncDerivative<Dd, V, E, Sp>& e)
		{
			_update(expr::get_enclosed_expression(e));
			e.update();
		}


		/* binary op pruning
		 */

		template<typename... Es, size_t... Is>
		void _update(OpAdd<Es...>& e, std::index_sequence<Is...>)
		{
			(_update(expr::get<Is>(e)), ...);
		}

		template<typename... Es>
		void _update(OpAdd<Es...>& e)
		{
			_update(e, std::make_index_sequence<sizeof...(Es)>{});
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
			if constexpr (expr::has_state<E>::value)
			{
				_update(expr::get_enclosed_expression(e));
			}
			e.update();
		}
		template<typename A1, typename A2, typename E>
		void _update(OpChain<A1, A2, E>& e)
		{
			if constexpr (expr::has_state<E>::value)
			{
				_update(expr::get_enclosed_expression(e));
			}
			e.update();
		}

		template<typename G, typename V, typename E>
		void _update(OpMap<G, V, E>& e)
		{
			_update(expr::get_enclosed_expression(e));
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
		protected:

			static const size_t deriv_order_1 = derivative_index_raw<1>::value;
			static const size_t deriv_order_0 = derivative_index_raw<0>::value;

		public:

			static const bool value = expr::derivative_index<deriv_order_1, E>::value > deriv_order_0;
		};

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
	template<typename E>
	auto apply_operators(OpExpression<E> const& e)
	{
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
	template<typename E>
	auto apply_operators(OpOperator<E> const& e)
	{
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
	template<typename A1, typename A2>
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
	template<typename A1, typename A2>
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
	template<typename A1, typename A2, typename E>
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
	template<typename A1, typename A2, typename E>
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
	template<typename Dd, typename V, typename E, typename Sp,
		typename = typename std::enable_if_t<expr_has_deriv<E>::value, int>>
	auto apply_operators(OpFuncDerivative<Dd, V, E, Sp> const& e);

	//! Distribute operators so they are applied to individual expressions.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<typename Dd1, typename Dd2, typename V1, typename V2, typename E, typename Sp>
	auto apply_operators(OpFuncDerivative<Dd1, V1, OpFuncDerivative<Dd2, V2, E, Sp>, Sp> const& e);

	//! Distribute operators so they are applied to individual expressions.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<typename Dd, typename V, typename... Es, typename Sp>
	auto apply_operators(OpFuncDerivative<Dd, V, OpAdd<Es...>, Sp> const& e);

	//! Distribute operators so they are applied to individual expressions.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<typename... Es>
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
	template<typename E1, typename E2>
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
	template<typename E1, typename E2>
	auto apply_operators(OpBinaryDiv<E1, E2> const& e);


	/// TODO: Add implementation for OpExponential


	 //! Implementation of the product rule.
	template<size_t O, typename G, typename V, typename E1, typename E2, typename std::enable_if_t<(O > 0), int> = 0>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpBinaryDiv<E1, E2>, SymbolicDerivative<G>> const& e);

	//! Implementation of the division rule.
	template<size_t O, typename G, typename V, typename E1, typename E2, typename std::enable_if_t<(O > 0), int> = 0>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpBinaryMul<E1, E2>, SymbolicDerivative<G>> const& e);
	
	template<size_t O, typename V, typename V1, typename... Gs, expr::exp_key_t... Xs, typename G0,
		size_t N = expr::factor_count_list<G0, Term<Gs, Xs>...>::value, 
		typename = std::enable_if_t<(O > 0 && N >= O && !expr::is_expression<G0>), int>>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpTerms<V1, Term<Gs, Xs>...>, SymbolicDerivative<G0>> const& e);
	
	template<size_t O, typename V, typename V1, typename... Gs, expr::exp_key_t... Xs, typename G0,
		size_t N = expr::factor_count_list<G0, Term<Gs, Xs>...>::value, typename = std::enable_if_t<(O > 0 && N >= O), int>>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpTerms<V1, Term<Gs, Xs>...>, SymbolicDerivative<OpTerm<OpIdentity, G0>>> const& e);

	template<size_t O, typename V0, auto f, typename V1, typename E, typename G0>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V0, OpFuncApply<f, V1, E>, SymbolicDerivative<G0>> const& e);

	template<size_t O, typename V, typename... Es, typename G0>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpAdd<Es...>, SymbolicDerivative<G0>> const& e);

	template<size_t O, typename G, typename V, typename E, typename std::enable_if_t<!expr_has_deriv<E>::value, int> = 0>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G>> const& e)
	{
		return OpVoid{};
	}


	template<size_t O, typename V, typename V1, typename E, typename Dd, typename Sp, 
		typename = std::enable_if_t<!std::is_same<V1, OpIdentity>::value, int>>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpFuncDerivative<Dd, V1, E, Sp>, SymbolicDerivative<OpFuncDerivative<Dd, OpIdentity, E, Sp>>> const& e)
	{
		auto&& expr = expr::get_enclosed_expression(e);
		return expr::make_literal(e.value * expr.value);
	}

	template<size_t O, typename V, typename V1, typename E, typename Dd, typename Sp, typename G0, 
		typename = std::enable_if_t<!std::is_same<G0, OpFuncDerivative<Dd, V1, E, Sp>>::value, int>>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpFuncDerivative<Dd, V1, E, Sp>, SymbolicDerivative<G0>> const& e)
	{
		return OpVoid{};
	}

	template<size_t O, typename V, typename V1, typename E, typename Dd, typename G00, typename G01,
		typename = std::enable_if_t<!std::is_same<G01, OpFuncDerivative<Dd, V1, E, SymbolicDerivative<G00>>>::value, int>>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpFuncDerivative<Dd, V1, E, SymbolicDerivative<G00>>, SymbolicDerivative<G01>> const& e)
	{
		auto&& expr = expr::get_enclosed_expression(e);
		auto&& d = expr::make_operator_derivative<O, G01>(e.value);
		return apply_operators(d * apply_operators(expr));
	}

	template<typename G, typename V, typename E>
	auto apply_operators(OpFuncDerivative<std::index_sequence<0>, V, E, SymbolicDerivative<G>> const& e)
	{
		return expr::make_literal(e.value) * expr::get_enclosed_expression(e);
	}


	template<typename V, typename E>
	auto apply_operators(OpFuncDerivative<std::index_sequence<1>, V, E, SymbolicDerivative<E>> const& e)
	{
		return expr::make_literal(e.value);
	}

	//! Implementation of the quotient rule.
	template<size_t O, typename G, typename V, typename E1, typename E2, typename std::enable_if_t<(O > 0), int>>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpBinaryDiv<E1, E2>, SymbolicDerivative<G>> const& e)
	{
		auto&& expr = expr::get_enclosed_expression(e);
		auto lhs = apply_operators(expr::make_derivative<O, G>(expr.a)) * expr.b;
		auto rhs = expr.a * apply_operators(expr::make_derivative<O, G>(expr.b));
		return expr::make_literal(e.value) * (lhs - rhs) / (expr.b * expr.b);
	}

	//! Implementation of the product rule.
	template<size_t O, typename G, typename V, typename E1, typename E2, typename std::enable_if_t<(O > 0), int>>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpBinaryMul<E1, E2>, SymbolicDerivative<G>> const& e)
	{
		auto&& expr = expr::get_enclosed_expression(e);
		auto lhs = apply_operators(expr::make_derivative<O, G>(expr.a)) * expr.b;
		auto rhs = expr.a * apply_operators(expr::make_derivative<O, G>(expr.b));
		return expr::make_literal(e.value) * (lhs + rhs);
	}

	namespace
	{
		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<expr_has_deriv<E>::value, int> = 0>
		auto apply_operators_chain(OpChain<A1, A2, E> const& e)
		{
			return apply_operators(
				expr::distribute_operator(e.combination.f, apply_operators(
					expr::distribute_operator(e.combination.g, 
						apply_operators(expr::get_enclosed_expression(e))))
				));
		}

		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<!expr_has_deriv<E>::value, int> = 0>
		auto apply_operators_chain(OpChain<A1, A2, E> const& e)
		{
			return apply_operators(expr::expand_operator(e.combination, expr::get_enclosed_expression(e)));
		}

		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<expr_has_deriv<E>::value, int> = 0>
		auto apply_operators_combination(OpCombination<A1, A2, E> const& e)
		{
			return apply_operators(expr::distribute_operator(
				e.combination.f, 
				apply_operators(expr::get_enclosed_expression(e))))
			+ apply_operators(expr::distribute_operator(
				e.combination.g, 
				apply_operators(expr::get_enclosed_expression(e))));
		}

		template<typename A1, typename A2, typename E,
			typename std::enable_if_t<!expr_has_deriv<E>::value, int> = 0>
		auto apply_operators_combination(OpCombination<A1, A2, E> const& e)
		{
			return apply_operators(expr::expand_operator(e.combination, expr::get_enclosed_expression(e)));
		}

		template<typename... Es, size_t... Is>
		auto apply_operators_adds(OpAdd<Es...> const& e, std::index_sequence<Is...>)
		{
			return (apply_operators(expr::get<Is>(e)) + ...);
		}

		template<typename Dd, typename Sp, typename... Es, size_t... Is>
		auto apply_operators_adds(Sp const& solver, OpAdd<Es...> const& e, std::index_sequence<Is...>)
		{
			return (apply_operators(expr::make_derivative<Dd>(apply_operators(expr::get<Is>(e)), solver)) + ...);
		}

		template<size_t O, typename G0, typename... Es, size_t... Is>
		auto apply_operators_adds(SymbolicDerivative<G0>, OpAdd<Es...> const& e, std::index_sequence<Is...>)
		{
			return (apply_operators(expr::make_derivative<O, G0>(apply_operators(expr::get<Is>(e)))) + ...);
		}
	}

	template<typename A1, typename A2>
	auto apply_operators(OpOperatorChain<A1, A2> const& e)
	{
		return apply_operators(apply_operators(e.f)(apply_operators(e.g)));
	}

	template<typename A1, typename A2>
	auto apply_operators(OpOperatorCombination<A1, A2> const& e)
	{
		return apply_operators(apply_operators(e.f) + apply_operators(e.g));
	}

	template<typename A1, typename A2, typename E>
	auto apply_operators(OpChain<A1, A2, E> const& e)
	{
		return apply_operators_chain(e);
	}

	template<typename A1, typename A2, typename E>
	auto apply_operators(OpCombination<A1, A2, E> const& e)
	{
		return apply_operators_combination(e);
	}

	namespace
	{

		template<size_t O, typename Sp, Axis... axs>
		auto break_up_derivative(solver_op_type<Sp> solver, symphas::lib::axis_list<axs...>)
		{
			return (expr::make_operator_directional_derivative<axs, O>(solver) + ...);
		}
	}

	template<typename Dd1, typename Dd2, typename V1, typename V2, typename E, typename Sp>
	auto apply_operators(OpFuncDerivative<Dd1, V1, OpFuncDerivative<Dd2, V2, E, Sp>, Sp> const& e)
	{
		using d1 = OpFuncDerivative<Dd1, V1, OpFuncDerivative<Dd2, V2, E, Sp>, Sp>;
		using d2 = OpFuncDerivative<Dd2, V2, E, Sp>;

		constexpr size_t D = expr::grid_dim<d1>::value;

		constexpr size_t order1 = d1::order;
		constexpr size_t order2 = d2::order;

		auto enclosed1 = expr::get_enclosed_expression(e);
		auto enclosed2 = expr::get_enclosed_expression(enclosed1);
		auto enclosed = apply_operators(expr::coeff(e) * expr::coeff(enclosed1) * enclosed2);

		// if they are both directional (a derivative of 1 is always directional)
		if constexpr (Dd1::is_directional && Dd2::is_directional)
		{
			if constexpr (Dd1::axis != Dd2::axis)
			{
				return expr::make_operator_mixed_derivative<Dd1::axis, order1, Dd2::axis, order2>(e.solver)(enclosed);
			}
			else
			{
				return expr::make_operator_directional_derivative<Dd1::axis, order1 + order2>(e.solver)(enclosed);
			}
		}
		else if constexpr (Dd1::is_directional)
		{
			if constexpr (order2 % 2 == 1)
			{
				return (expr::make_operator_directional_derivative<Dd1::axis, order1>(e.solver) *
					expr::make_operator_directional_derivative<Dd2::axis, 1>(e.solver) *
					break_up_derivative<order2 - 1>(e.solver, symphas::lib::make_axis_list<D>()))(enclosed);
			}
			else // order 2 even
			{
				if constexpr (order1 == 1)
				{
					return expr::make_derivative<Dd1::axis, order1 + order2>(enclosed, e.solver);
				}
				else
				{
					return (expr::make_operator_directional_derivative<Dd1::axis, order1>(e.solver) *
						break_up_derivative<order2>(e.solver, symphas::lib::make_axis_list<D>()))(enclosed);
				}
			}
		}
		else if constexpr (Dd2::is_directional)
		{
			if constexpr (order1 % 2 == 1)
			{
				return (expr::make_operator_directional_derivative<Dd2::axis, order2>(e.solver) *
					expr::make_operator_directional_derivative<Dd1::axis, 1>(e.solver) *
					break_up_derivative<order1 - 1>(e.solver, symphas::lib::make_axis_list<D>()))(enclosed);
			}
			else // order 2 even
			{
				if constexpr (order2 == 1)
				{
					return expr::make_derivative<Dd2::axis, order2 + order1>(enclosed, e.solver);
				}
				else
				{
					return (expr::make_operator_directional_derivative<Dd2::axis, order2>(e.solver) *
						break_up_derivative<order1>(e.solver, symphas::lib::make_axis_list<D>()))(enclosed);
				}
			}
		}
		else 
		{
			if constexpr (order1 % 2 == 0 && order2 % 2 == 0)
			{
				return expr::make_derivative<Dd1::axis, order1 + order2>(enclosed, e.solver);
			}
			else if constexpr (order1 % 2 == 1)
			{
				return expr::make_derivative<Dd1::axis, order1 + order2>(enclosed, e.solver);
			}
			else if constexpr (order2 % 2 == 1)
			{
				return expr::make_derivative<Dd2::axis, order1 + order2>(enclosed, e.solver);
			}
			else
			{
				return (expr::make_operator_directional_derivative<Dd1::axis, 1>(e.solver) *
					expr::make_operator_directional_derivative<Dd2::axis, 1>(e.solver) *
					break_up_derivative<order1 + order2 - 2>(e.solver, symphas::lib::make_axis_list<D>()))(enclosed);
			}
		}
	}

	template<typename Dd, typename V, typename E, typename Sp, typename>
	auto apply_operators(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		return apply_operators(expr::make_derivative<Dd>(apply_operators(expr::coeff(e) * expr::get_enclosed_expression(e)), e.solver));
	}

	template<typename Dd, typename V, typename... Es, typename Sp>
	auto apply_operators(OpFuncDerivative<Dd, V, OpAdd<Es...>, Sp> const& e)
	{
		auto&& add = expr::get_enclosed_expression(e);
		return e.value * apply_operators_adds<Dd>(e.solver, add, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<size_t O, typename V, typename... Es, typename G0>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpAdd<Es...>, SymbolicDerivative<G0>> const& e)
	{
		auto&& add = expr::get_enclosed_expression(e);
		return e.value * apply_operators_adds<O>(SymbolicDerivative<G0>{}, add, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<typename... Es>
	auto apply_operators(OpAdd<Es...> const& e)
	{
		return apply_operators_adds(e, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<typename E1, typename E2>
	auto apply_operators(OpBinaryMul<E1, E2> const& e)
	{
		return apply_operators(e.a) * apply_operators(e.b);
	}

	template<typename E1, typename E2>
	auto apply_operators(OpBinaryDiv<E1, E2> const& e)
	{
		return apply_operators(e.a) / apply_operators(e.b);
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
		auto operator()(OpExpression<E> const& e)
		{
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
		auto operator()(OpExpression<E> const& e)
		{
			expr::prune::update(*const_cast<E*>(static_cast<E const*>(&e)));
			Block<T> result(expr::data_len(*static_cast<E const*>(&e)));
			expr::result(*static_cast<E const*>(&e), result.values, result.len);
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
	auto to_grid(OpExpression<E> const& e)
	{
		return symphas::internal::construct_grid_of_dimension<D>{}(*static_cast<E const*>(&e));
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

	//! Evaluating an identity or fraction as a grid returns just the value of the literal.
	/*!
	 * Evaluating an identity or fraction as a grid returns just the value of the literal.
	 * See expr::transform::to_grid(OpExpression<E>&)
	 */
	template<typename coeff_t, typename = std::enable_if_t<(expr::is_identity<coeff_t> || expr::is_fraction<coeff_t>), int>>
	auto to_grid(coeff_t)
	{
		return coeff_t{}.eval();
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


	////! Create a grid where values are adjusted to account for FFTW routines.
	///*!
	// * Converts the expression into a grid with indexing consistent with the 
	// * the FFTW routine for that type. That is, 
	// * expected. The index of the values is adjusted based on what is expected
	// * in order to perform FFTW functions.
	// *
	// * \tparam T_src The known source type of the underlying data or evaluation
	// * type of the expression, which is considered when determining what index
	// * scheme to use for evaluation.
	// * \tparam D The dimension of the grid to create.
	// * \tparam T The result type of the expression.
	// */
	//template<typename E, 
	//	typename T = typename expr::eval_type<E>::type, 
	//	size_t D = expr::grid_dim<E>::value>
	//Grid<T, D> to_fftw_grid(OpExpression<E> const& e)
	//{
	//	expr::prune::update(*const_cast<E*>(static_cast<E const*>(&e)));

	//	const len_type* dims = expr::data_dimensions(*static_cast<E const*>(&e));

	//	len_type fftw_dims[D];
	//	fftw_dims[0] = dims[0] / 2 + 1;
	//	std::copy(dims + 1, dims + D, fftw_dims + 1);

	//	Grid<T, D> grid(dims);
	//	symphas::dft::iterate_rc<T, D>(grid.values, *static_cast<E const*>(&e), dims);

	//	return grid;
	//}


	////! Turning a constant to a grid returns the constant.
	///*!
	// * Turning a constant to a grid returns the constant.
	// *
	// * \param e The constant to turn into a grid.
	// */
	//template<typename coeff_t, typename = std::enable_if_t<expr::is_coeff<coeff_t>, int>>
	//auto to_fftw_grid(coeff_t const& c)
	//{
	//	return c;
	//}

	//namespace
	//{
	//	template<typename... Es, size_t... Is>
	//	decltype(auto) to_fftw_grid(std::tuple<Es...>& e, std::index_sequence<Is...>)
	//	{
	//		return std::make_tuple(expr::transform::to_fftw_grid(std::get<Is>(e))...);
	//	}
	//}

	////! Create a list of grids initialized to the result of the expressions.
	///*!
	// * Converts the operator into a grid with indexing consistent with the type
	// * expected. Constructs the list by
	// * applying expr::to_fftw_grid(OpExpression<E>&) to each element of the list.
	// *
	// * \param es The list of expressions to evaluate into grids.
	// */
	//template<typename... Es>
	//decltype(auto) to_fftw_grid(std::tuple<Es...>& es)
	//{
	//	return expr::transform::to_fftw_grid(es, std::make_index_sequence<sizeof...(Es)>{});
	//}


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
	template<size_t D, typename T0, typename T1, typename V, typename E>
	auto to_ft(OpMap<MapGridFourier<T0, T1, D>, V, E> const& e, double const*, const len_type* dims);

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
	template<size_t D, Axis ax, size_t O, typename V, typename Sp>
	auto to_ft(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename... Es>
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
	template<size_t D, typename E1, typename E2>
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
	template<size_t D, typename E>
	auto to_ft(OpBinaryDiv<OpIdentity, E> const& e, double const* h, const len_type* dims);

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
	template<size_t D, size_t O, typename T>
	auto to_ft(OpTerm<T, k_grid_type<O, D>> const& e, double const* h, const len_type*)
	{
		//return expr::make_operator_derivative<O>(Solver<void>{});
	}

	template<size_t D, Axis ax, size_t O, typename T>
	auto to_ft(OpTerm<T, k_grid_axis_type<ax, O, D>> const& e, double const* h, const len_type*)
	{
		//return symphas::internal::nth_derivative_apply<ax, O, Sp>::template get(Solver<void>{});
	}

	template<size_t D, Axis ax, size_t O, typename T>
	auto to_ft(OpTerm<T, k_grid_component_type<ax, O, D>> const& e, double const* h, const len_type*)
	{
		//return symphas::internal::make_directional_operator_derivative<O>(Solver<void>{});
	}

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
		return e * to_ft<D>(OpIdentity{}, h, dims);
	}

	template<size_t D, size_t N, size_t D0>
	auto to_ft(OpFractionLiteral<N, D>, double const* h, const len_type* dims)
	{
		return OpFractionLiteral<N, D>{} * to_ft<D>(OpIdentity{}, h, dims);
	}

	template<size_t D, size_t N, size_t D0>
	auto to_ft(OpNegFractionLiteral<N, D>, double const* h, const len_type* dims)
	{
		return OpNegFractionLiteral<N, D>{} * to_ft<D>(OpIdentity{}, h, dims);
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

	template<size_t D, typename T0, typename T1, typename V, typename E>
	auto to_ft(OpMap<MapGridInverseFourier<T0, T1, D>, V, E> const& e, double const*, const len_type* dims)
	{
		return expr::coeff(e) * expr::get_enclosed_expression(e);
	}

	template<size_t D, size_t O, typename V, typename Sp>
	auto to_ft(OpOperatorDerivative<O, V, Sp> const& e, double const* h, const len_type* dims)
	{
		return expr::coeff(e) * expr::make_term(k_grid_type<O, D>(dims, h));
	}

	template<size_t D, Axis ax, size_t O, typename V, typename Sp>
	auto to_ft(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e, double const* h, const len_type* dims)
	{
		return expr::coeff(e) * expr::make_term(k_grid_axis_type<ax, O, D>(dims, h));
	}

	template<size_t D, typename V, typename E1, typename E2>
	auto to_ft(OpFuncConvolution<V, E1, E2> const& e, double const* h, const len_type* dims)
	{
		return expr::coeff(e) * to_ft<D>(e.a, h, dims) * to_ft<D>(e.b, h, dims);
	}

	template<size_t D>
	auto to_ft(GaussianSmoothing<D> const& e, double const*, const len_type*)
	{
#ifdef PRINTABLE_EQUATIONS
		char* gname = new char[e.print_length() + STR_ARR_LEN(SYEX_FT_OF_OP_FMT_A SYEX_FT_OF_OP_FMT_B)];
		size_t n = sprintf(gname, SYEX_FT_OF_OP_FMT_A);
		n += e.print(gname + n);
		n += sprintf(gname + n, SYEX_FT_OF_OP_FMT_B);
		auto op = expr::make_term(NamedData(to_grid(e), gname));

		delete[] gname;
		return op;
#else
		return expr::make_term(to_grid(e));
#endif
	}

	template<size_t D, typename V, typename E>
	auto to_ft(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e, double const* h, const len_type* dims)
	{
		return expr::coeff(e) * to_ft<D>(e.smoother, h, dims) * to_ft<D>(expr::get_enclosed_expression(e), h, dims);
	}

	namespace
	{
		template<size_t D, typename... Es, size_t... Is>
		auto to_ft_adds(OpAdd<Es...> const& e, double const* h, const len_type* dims, std::index_sequence<Is...>)
		{
			return (to_ft<D>(expr::get<Is>(e), h, dims) + ...);
		}
	}

	template<size_t D, typename... Es>
	auto to_ft(OpAdd<Es...> const& e, double const* h, const len_type* dims)
	{
		return to_ft_adds<D>(e, h, dims, std::make_index_sequence<sizeof...(Es)>{});
	}


	template<size_t D, typename E1, typename E2>
	auto to_ft(OpBinaryMul<E1, E2> const& e, double const* h, const len_type* dims)
	{
		return expr::make_convolution(to_ft<D>(e.a, h, dims), to_ft<D>(e.b, h, dims));
	}

	template<size_t D, typename E>
	auto to_ft(OpBinaryDiv<OpIdentity, E> const& e, double const* h, const len_type* dims)
	{
		return expr::make_fourier_map(e);
	}

	template<size_t D, typename E1, typename E2>
	auto to_ft(OpBinaryDiv<E1, E2> const& e, double const* h, const len_type* dims)
	{
		return expr::make_convolution(to_ft<D>(e.a, h, dims), to_ft<D>(expr::inverse(e.b), h, dims));
	}

	/* conversion of operators into fourier space
	 */

	template<size_t D, typename A1, typename A2, typename E>
	auto to_ft(OpChain<A1, A2, E> const& e, double const* h, const len_type* dims)
	{
		return to_ft<D>(e.combination, h, dims) * to_ft<D>(expr::get_enclosed_expression(e), h, dims);
	}

	template<size_t D, typename A1, typename A2, typename E>
	auto to_ft(OpCombination<A1, A2, E> const& e, double const* h, const len_type* dims)
	{
		return to_ft<D>(e.combination, h, dims) * to_ft<D>(expr::get_enclosed_expression(e), h, dims);
	}

	template<size_t D, typename Dd, typename V, typename E, typename Sp>
	auto to_ft(OpFuncDerivative<Dd, V, E, Sp> const& e, double const* h, const len_type* dims)
	{
		constexpr Axis axis = OpFuncDerivative<Dd, V, E, Sp>::axis;
		constexpr size_t order = OpFuncDerivative<Dd, V, E, Sp>::order;

		if constexpr (Dd::is_directional)
		{
			return expr::coeff(e) * expr::make_term(k_grid_axis_type<axis, order, D>(dims, h))
				* to_ft<D>(expr::get_enclosed_expression(e), h, dims);
		}
		else
		{
			if constexpr (order % 2 == 0)
			{
				return expr::coeff(e) * expr::make_term(k_grid_type<order, D>(dims, h))
					* to_ft<D>(expr::get_enclosed_expression(e), h, dims);
			}
			else
			{
				return expr::coeff(e) * expr::make_term(k_grid_component_type<axis, order, D>(dims, h))
					* to_ft<D>(expr::get_enclosed_expression(e), h, dims);
			}
		}

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
	decltype(auto) swap_grid(OpExpression<E> const& e, G_F&& g);

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
	template<size_t Z, typename V, typename... Gs, exp_key_t... Xs, typename G_F>
	decltype(auto) swap_grid(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g);

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
	template<size_t Z, typename... Es, typename G_F>
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
	template<size_t Z, typename E1, typename E2, typename G_F>
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
	template<size_t Z, typename A1, typename A2, typename G_F>
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
	template<size_t Z, typename A1, typename A2, typename E, typename G_F>
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
	template<size_t Z, typename A1, typename A2, typename G_F>
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
	template<size_t Z, typename A1, typename A2, typename E, typename G_F>
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
	template<size_t Z, typename V, typename E, typename F, typename Arg0, typename... Args, typename G_F>
	auto swap_grid(OpFunc<V, E, F, Arg0, Args...> const& e, G_F&& g);

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
	template<size_t Z, auto f, typename V, typename E, typename G_F>
	auto swap_grid(OpFuncApply<f, V, E> const& e, G_F&& g);


	namespace
	{
		template<typename V>
		auto sift_term(V value)
		{
			return value;
		}

		template<typename V, typename G0, exp_key_t X0, typename... G1s, exp_key_t... X1s>
		auto sift_term(V value, Term<G0, X0> const& term0, Term<G1s, X1s> const&... rest)
		{
			return OpTerms(value, term0, rest...);
		}

		template<typename V, typename... G1s, exp_key_t... X1s, size_t... Ns, size_t... Ms, typename G_F>
		auto pick_terms(OpTerms<V, Term<G1s, X1s>...> const& a,
			std::index_sequence<Ns...>, std::index_sequence<Ms...>, 
			std::integer_sequence<exp_key_t>, G_F&& g)
		{
			return a;
		}


		template<typename V, typename... G1s, exp_key_t... X1s, size_t... Ns, size_t... Ms, exp_key_t X, typename G_F>
		auto pick_terms(OpTerms<V, Term<G1s, X1s>...> const& a, 
			std::index_sequence<Ns...>, std::index_sequence<Ms...>, 
			std::integer_sequence<exp_key_t, X>, G_F&& g)
		{
			if constexpr (!expr::is_expression<G_F>)
			{
				return sift_term(expr::coeff(a), expr::get<Ns>(a)..., Term(std::forward<G_F>(g), expr::get<Ms>(a)...));
			}
			else
			{
				if constexpr (_Xk_t<X>::D == 1)
				{
					if constexpr (_Xk_t<X>::N == 1)
					{
						return sift_term(expr::coeff(a), expr::get<Ns>(a)...) 
							* std::forward<G_F>(g) 
							* sift_term(OpIdentity{}, expr::get<Ms>(a)...);
					}
					else
					{
						return sift_term(expr::coeff(a), expr::get<Ns>(a)...)
							* expr::pow<_Xk_t<X>::N>(std::forward<G_F>(g)) 
							* sift_term(OpIdentity{}, expr::get<Ms>(a)...);
					}
				}
				else
				{
					return sift_term(expr::coeff(a), expr::get<Ns>(a)...)
						* expr::exp(
							((_Xk_t<X>::sign) ? OpNegIdentity{} : OpIdentity{}) * expr::make_fraction<_Xk_t<X>::N, _Xk_t<X>::D>()
							* expr::log(std::forward<G_F>(g)))
						* sift_term(OpIdentity{}, expr::get<Ms>(a)...);
				}

			}
		}


		template<typename V, typename... Gs, exp_key_t... Xs, size_t... Is, bool... fs, typename G_F>
		decltype(auto) swap_terms(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<Is...>, std::integer_sequence<bool, fs...>, G_F&& g)
		{
			using ::symphas::lib::seq_join_t;
			
			using swap_seq_t = seq_join_t<
				std::index_sequence<>,
				std::conditional_t<
					fs,
					std::index_sequence<Is>,
					std::index_sequence<>>...
				>;

			if constexpr (swap_seq_t::size() > 0)
			{
				constexpr size_t swap_index = symphas::lib::seq_index_value<0, swap_seq_t>::value;

				using pick_terms1_t = symphas::lib::seq_add_t<
					std::make_index_sequence<swap_index>,
					symphas::lib::seq_repeating_value_t<swap_index, size_t, 1>>;
				using pick_terms2_t = symphas::lib::seq_add_t<
					std::make_index_sequence<sizeof...(Is) - swap_index - 1>,
					symphas::lib::seq_repeating_value_t<sizeof...(Is) - swap_index - 1, size_t, swap_index + 2>>;

				using pick_power_t = symphas::lib::type_at_index<swap_index, std::integer_sequence<exp_key_t, Xs>...>;

				return pick_terms(e, pick_terms1_t{}, pick_terms2_t{}, pick_power_t{}, std::forward<G_F>(g));
			}
			else
			{
				return e;
			}
		}

		template<size_t Z, typename... Es, typename G_F, size_t... Is>
		auto swap_grid_adds(OpAdd<Es...> const& e, G_F&& g, std::index_sequence<Is...>)
		{
			return (swap_grid<Z>(expr::get<Is>(e), std::forward<G_F>(g)) + ...);
		}
	}

	template<size_t Z, typename E, typename G_F>
	decltype(auto) swap_grid(OpExpression<E> const& e, G_F&&)
	{
		return *static_cast<E const*>(&e);
	}

	template<size_t Z, typename V, typename... Gs, exp_key_t... Xs, typename G_F>
	decltype(auto) swap_grid(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g)
	{
		using mask_t = std::integer_sequence<bool, (expr::factor_count<Variable<Z>, Gs>::value > 0)...>;
		return swap_terms(e, std::make_index_sequence<sizeof...(Gs)>{}, mask_t{}, std::forward<G_F>(g));
	}

	template<size_t Z, typename... Es, typename G_F>
	auto swap_grid(OpAdd<Es...> const& e, G_F&& g)
	{
		return swap_grid_adds<Z>(e, std::forward<G_F>(g), std::make_index_sequence<sizeof...(Es)>{});
	}

	template<size_t Z, typename E1, typename E2, typename G_F>
	auto swap_grid(OpBinaryMul<E1, E2> const& e, G_F&& g)
	{
		return swap_grid<Z>(e.a, std::forward<G_F>(g)) * swap_grid<Z>(e.b, std::forward<G_F>(g));
	}

	template<size_t Z, typename E1, typename E2, typename G_F>
	auto swap_grid(OpBinaryDiv<E1, E2> const& e, G_F&& g)
	{
		return swap_grid<Z>(e.a, std::forward<G_F>(g)) / swap_grid<Z>(e.b, std::forward<G_F>(g));
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
			swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), 
			e.smoother);
	}

	template<size_t Z, typename A1, typename A2, typename G_F>
	auto swap_grid(OpOperatorChain<A1, A2> const& e, G_F&& g)
	{
		return OpOperatorChain(swap_grid<Z>(e.f, std::forward<G_F>(g)), swap_grid<Z>(e.g, std::forward<G_F>(g)));
	}

	template<size_t Z, typename A1, typename A2, typename E, typename G_F>
	auto swap_grid(OpChain<A1, A2, E> const& e, G_F&& g)
	{
		return swap_grid<Z>(e.combination, std::forward<G_F>(g))(swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
	}

	template<size_t Z, typename A1, typename A2, typename G_F>
	auto swap_grid(OpOperatorCombination<A1, A2> const& e, G_F&& g)
	{
		return OpOperatorCombination(swap_grid<Z>(e.f, std::forward<G_F>(g)), swap_grid<Z>(e.g, std::forward<G_F>(g)));
	}

	template<size_t Z, typename A1, typename A2, typename E, typename G_F>
	auto swap_grid(OpCombination<A1, A2, E> const& e, G_F&& g)
	{
		return swap_grid<Z>(e.combination, std::forward<G_F>(g))(swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
	}

	template<size_t Z, typename Dd, typename V, typename E, typename Sp, typename G_F>
	auto swap_grid(OpFuncDerivative<Dd, V, E, Sp> const& e, G_F&& g)
	{
		constexpr size_t order = OpFuncDerivative<Dd, V, E, Sp>::order;
		constexpr Axis axis = OpFuncDerivative<Dd, V, E, Sp>::axis;
		return symphas::internal::nth_derivative_apply<axis, order, Sp>::template get(
			e.value, swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), e.solver);
	}

	template<size_t Z, typename V, typename E, typename F, typename Arg0, typename... Args, typename G_F>
	auto swap_grid(OpFunc<V, E, F, Arg0, Args...> const& e, G_F&& g)
	{
		return OpFunc(e.name, e.value, swap_grid<Z>(e.e, std::forward<G_F>(g)), e.f, e.tt);
	}

	template<size_t Z, auto f, typename V, typename E, typename G_F>
	auto swap_grid(OpFuncApply<f, V, E> const& e, G_F&& g)
	{
        auto eg = swap_grid<Z>(e.e, std::forward<G_F>(g));
		return OpFuncApply<f, V, decltype(eg)>(e.value, eg);
	}




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
	template<typename Sg, typename G_F>
	decltype(auto) swap_grid(OpExpression<Sg> const& e, G_F&& g)
	{
		return std::forward<G_F>(g);
	}

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
	template<typename Sg, typename E, typename G_F, typename = std::enable_if_t<!std::is_same<Sg, E>::value, int>>
	decltype(auto) swap_grid(OpExpression<E> const& e, G_F&& g)
	{
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
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename T, typename G, typename G_F, 
		typename std::enable_if_t<(expr::is_same_base<Sg, G> && !expr::is_expression<G_F>), int> = 0>
	decltype(auto) swap_grid(OpTerm<T, G> const& v, G_F&& g);

	//! Swap a data term in the expression.
	/*!
	 * Implementation of a successful search, where the given variable term
	 * associated with the prescribed type will be switched with the given
	 * expression.
	 *
	 * \param v The term which is swapped.
	 * \param g The element which will replace the variable.
	 *
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename V, typename... Gs, exp_key_t... Xs, typename G_F>
	decltype(auto) swap_grid(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g);

	//! Swap a data term in the expression.
	/*!
	 * Swaps the instance of the variable term which matches the data type
	 * for a different term or expression. If the term is not found, the
	 * expression is returned unchanged.
	 *
	 * \param e The expression to search for the term to swap.
	 * \param g The element which will replace the variable.
	 *
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename... Es, typename G_F>
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
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename E1, typename E2, typename G_F>
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
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename E1, typename E2, typename G_F>
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
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename A1, typename A2, typename G_F>
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
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename A1, typename A2, typename E, typename G_F>
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
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename A1, typename A2, typename G_F>
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
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename A1, typename A2, typename E, typename G_F>
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
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename V, typename E1, typename E2, typename G_F>
	auto swap_grid(OpFuncConvolution<V, E1, E2> const& e, G_F&& g);

	//! Swap a data term in the expression.
	/*!
	 * Swaps the instance of the variable term which matches the data type
	 * for a different term or expression. If the term is not found, the
	 * expression is returned unchanged.
	 *
	 * \param e The expression to search for the term to swap.
	 * \param g The element which will replace the variable.
	 *
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename V, size_t D, typename E, typename G_F>
	auto swap_grid(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e, G_F&& g);

	//! Swap a data term in the expression.
	/*!
	 * Swaps the instance of the variable term which matches the data type
	 * for a different term or expression. If the term is not found, the
	 * expression is returned unchanged.
	 *
	 * \param e The expression to search for the term to swap.
	 * \param g The element which will replace the variable.
	 *
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename Dd, typename V, typename E, typename Sp, typename G_F>
	auto swap_grid(OpFuncDerivative<Dd, V, E, Sp> const& e, G_F&& g);


	//! Swap a data term in the expression.
	/*!
	 * Swaps the instance of the variable term which matches the data type
	 * for a different term or expression. If the term is not found, the
	 * expression is returned unchanged.
	 *
	 * \param e The expression to search for the term to swap.
	 * \param g The element which will replace the variable.
	 *
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, typename V, typename E, typename F, typename Arg0, typename... Args, typename G_F>
	auto swap_grid(OpFunc<V, E, F, Arg0, Args...> const& e, G_F&& g);


	//! Swap a data term in the expression.
	/*!
	 * Swaps the instance of the variable term which matches the data type
	 * for a different term or expression. If the term is not found, the
	 * expression is returned unchanged.
	 *
	 * \param e The expression to search for the term to swap.
	 * \param g The element which will replace the variable.
	 *
	 * 	param Sg The type of the grid to match for the swap.
	 */
	template<typename Sg, auto f, typename V, typename E, typename G_F>
	auto swap_grid(OpFuncApply<f, V, E> const& e, G_F&& g);


	template<typename Sg, typename V, typename... Gs, exp_key_t... Xs, typename G_F>
	decltype(auto) swap_grid(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g)
	{
		using mask_t = std::integer_sequence<bool, (expr::factor_count<Sg, Gs>::value > 0)...>;
		return swap_terms(e, std::make_index_sequence<sizeof...(Gs)>{}, mask_t{}, std::forward<G_F>(g));
	}

	namespace
	{
		template<typename Sg, typename... Es, typename G_F, size_t... Is>
		auto swap_grid_adds(OpAdd<Es...> const& e, G_F&& g, std::index_sequence<Is...>)
		{
			return (swap_grid<Sg>(expr::get<Is>(e), std::forward<G_F>(g)) + ...);
		}
	}

	template<typename Sg, typename... Es, typename G_F>
	auto swap_grid(OpAdd<Es...> const& e, G_F&& g)
	{
		return swap_grid_adds<Sg>(e, std::forward<G_F>(g), std::make_index_sequence<sizeof...(Es)>{});
	}

	template<typename Sg, typename E1, typename E2, typename G_F>
	auto swap_grid(OpBinaryMul<E1, E2> const& e, G_F&& g)
	{
		return swap_grid<Sg>(e.a, std::forward<G_F>(g)) * swap_grid<Sg>(e.b, std::forward<G_F>(g));
	}

	template<typename Sg, typename E1, typename E2, typename G_F>
	auto swap_grid(OpBinaryDiv<E1, E2> const& e, G_F&& g)
	{
		return swap_grid<Sg>(e.a, std::forward<G_F>(g)) / swap_grid<Sg>(e.b, std::forward<G_F>(g));
	}

	template<typename Sg, typename V, typename E1, typename E2, typename G_F>
	auto swap_grid(OpFuncConvolution<V, E1, E2> const& e, G_F&& g)
	{
		return expr::make_convolution(
			e.value,
			swap_grid<Sg>(e.a, std::forward<G_F>(g)),
			swap_grid<Sg>(e.b, std::forward<G_F>(g)));
	}

	template<typename Sg, typename V, size_t D, typename E, typename G_F>
	auto swap_grid(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e, G_F&& g)
	{
		return expr::make_convolution(
			e.value,
			swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)),
			e.smoother);
	}

	template<typename Sg, typename A1, typename A2, typename G_F>
	auto swap_grid(OpOperatorChain<A1, A2> const& e, G_F&& g)
	{
		return OpOperatorChain(swap_grid<Sg>(e.f, std::forward<G_F>(g)), swap_grid<Sg>(e.g, std::forward<G_F>(g)));
	}

	template<typename Sg, typename A1, typename A2, typename E, typename G_F>
	auto swap_grid(OpChain<A1, A2, E> const& e, G_F&& g)
	{
		return e.combination * swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g));
	}

	template<typename Sg, typename A1, typename A2, typename G_F>
	auto swap_grid(OpOperatorCombination<A1, A2> const& e, G_F&& g)
	{
		return OpOperatorCombination(swap_grid<Sg>(e.f, std::forward<G_F>(g)), swap_grid<Sg>(e.g, std::forward<G_F>(g)));
	}

	template<typename Sg, typename A1, typename A2, typename E, typename G_F>
	auto swap_grid(OpCombination<A1, A2, E> const& e, G_F&& g)
	{
		return e.combination * swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g));
	}

	template<typename Sg, typename Dd, typename V, typename E, typename Sp, typename G_F>
	auto swap_grid(OpFuncDerivative<Dd, V, E, Sp> const& e, G_F&& g)
	{
		constexpr size_t order = OpFuncDerivative<Dd, V, E, Sp>::order;
		constexpr Axis axis = OpFuncDerivative<Dd, V, E, Sp>::axis;
		return symphas::internal::nth_derivative_apply<axis, order, Sp>::template get(
			e.value, swap_grid<Sp>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), e.solver);;
	}

	template<typename Sg, typename V, typename E, typename F, typename Arg0, typename... Args, typename G_F>
	auto swap_grid(OpFunc<V, E, F, Arg0, Args...> const& e, G_F&& g)
	{
		return OpFunc(e.name, e.value, swap_grid<Sg>(e.e, std::forward<G_F>(g)), e.f, e.tt);
	}

	template<typename Sg, auto f, typename V, typename E, typename G_F>
	auto swap_grid(OpFuncApply<f, V, E> const& e, G_F&& g)
	{
        auto eg = swap_grid<Sg>(e.e, std::forward<G_F>(g));
		return OpFuncApply<f, V, decltype(eg)>(e.value, eg);
	}



	template<Axis ax, size_t Z, typename E, typename G_F>
	auto swap_grid(OpExpression<E> const& e, G_F&& g)
	{
		return swap_grid<Variable<Z, VectorComponent<ax>>>(*static_cast<E const*>(&e), std::forward<G_F>(g));
	}

	template<Axis ax, typename Sg, typename E, typename G_F>
	auto swap_grid(OpExpression<E> const& e, G_F&& g)
	{
		return swap_grid<VectorComponent<ax, Sg>>(*static_cast<E const*>(&e), std::forward<G_F>(g));
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

		template<typename... E0s, typename... E1s>
		auto adds_expand_pair(std::pair<E0s, E1s> const& ...pairs)
		{
			return std::make_pair((pairs.first + ...), (pairs.second + ...));
		}

		template<typename E0, typename E1, typename... E0s, typename... E1s>
		auto adds_expand_pair_no_first(std::pair<E0, E1> const& pair0, std::pair<E0s, E1s> const& ...pairs)
		{
			return std::make_pair(pair0.first, pair0.second + (pairs.second + ...));
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
		static const size_t value = min_order_from_index<expr::derivative_index<0, E>::value>::value;
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
	template<typename... Es>
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
	template<typename A, typename B>
	auto by_linear(OpBinaryDiv<A, B> const& e);

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
			return expr::split::by_linear(expr::apply_operators(e));
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
				&& (!expression_predicate<E>::nonlinear && expression_predicate<E>::combination)), int> = 0>
		auto by_linear_combination(OpCombination<A1, A2, E> const& e)
		{
			return expr::split::by_linear(expr::apply_operators(e));
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
			return expr::split::by_linear(expr::apply_operators(e));
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
	auto by_linear(OpTerm<T, G> const& e)
	{
		return pack_left(e);
	}

	namespace
	{
		template<typename... Es, size_t... Is>
		auto by_linear_adds(OpAdd<Es...> const& e, std::index_sequence<Is...>)
		{
			return adds_expand_pair(by_linear(expr::get<Is>(e))...);
			//return std::make_pair((std::get<0>(std::get<Is>(a)) + ...), (std::get<1>(std::get<Is>(a)) + ...));
		}

		// tries to find if there's an expression in the numerator that is linear through the
		// denominator.
		// **************************************************************************************

		template<typename... E0s, typename... E1s, typename... Ts>
		auto filter_linear(std::tuple<E0s...> const& terms0, Ts const& ...ts);

		template<typename... E0s, size_t... Is>
		auto filter_linear(std::tuple<E0s...> const& terms, std::index_sequence<Is...>)
		{
			return std::make_tuple(std::get<Is>(terms)...);
		}

		template<typename... E0s, typename... E1s, typename... Ts, size_t... Is>
		auto filter_linear(std::tuple<E0s...> const& terms0, std::index_sequence<Is...>, 
			std::tuple<E1s...> const& terms1, Ts const& ...ts)
		{
			using filtered_indices_t = symphas::lib::seq_join_t<
				std::index_sequence<>,
				std::conditional_t<
					(symphas::lib::index_of_type<E0s, E1s...> >= 0 && !std::is_same<E0s, OpVoid>::value),
					std::index_sequence<Is>,
					std::index_sequence<>>...>;

			return filter_linear(filter_linear(terms0, filtered_indices_t{}), ts...);
		}

		template<typename... E0s, typename... E1s, typename... Ts>
		auto filter_linear(std::tuple<E0s...> const& terms0, Ts const& ...ts)
		{
			return filter_linear(terms0, std::make_index_sequence<sizeof...(E0s)>{}, ts...);
		}

		template<typename E>
		auto get_linear_term(OpExpression<E> const& e)
		{
			auto [l, _] = expr::split::by_linear(*static_cast<E const*>(&e));
			return std::make_tuple(l);
		}

		template<typename A, typename B>
		auto get_linear_term(OpBinaryDiv<A, B> const& e)
		{
			return std::make_tuple(OpVoid{});
		}


		// \p a are the terms from the numerator, and \p e is the term from the denominator which is being
		// checked that it can divide at least one of them.
		template<typename E, typename A>
		auto check_term_in_num(OpExpression<E> const& e, OpExpression<A> const& a, std::index_sequence<>)
		{
			auto div = (*static_cast<A const*>(&a)) / (*static_cast<E const*>(&e));
			return get_linear_term(div);
		}

		// \p a are the terms from the numerator, and \p e is the term from the denominator which is being
		// checked that it can divide at least one of them.
		template<typename E, typename... As, size_t... Is>
		auto check_term_in_num(OpExpression<E> const& e, OpAdd<As...> const& a, std::index_sequence<Is...>)
		{
			auto coeff = check_term_in_num(*static_cast<E const*>(&e), expr::get<sizeof...(Is)>(a), std::index_sequence<>{});
			auto rest = check_term_in_num(*static_cast<E const*>(&e), (expr::get<Is>(a) + ...), std::make_index_sequence<sizeof...(Is) - 1>{});

			if constexpr (std::is_same<std::tuple<OpVoid>, decltype(coeff)>::value)
			{
				return rest;
			}
			else
			{
				return std::tuple_cat(coeff, rest);
			}
		}

		// \p a are the terms from the numerator, and \p e is the term from the denominator which is being
		// checked that it can divide at least one of them.
		template<typename... Es, typename... As, size_t... Js>
		auto check_term_in_num(OpAdd<Es...> const& e, OpAdd<As...> const& a, std::index_sequence<Js...>)
		{
			auto coeffs = filter_linear(check_term_in_num(expr::get<Js>(e), a, std::make_index_sequence<sizeof...(As) - 1>{})...);
			return coeffs;
		}

		template<typename A, typename... Bs>
		auto by_linear_divs(OpBinaryDiv<A, OpAdd<Bs...>> const& e)
		{
			return std::tuple<OpVoid>{};
		}

		template<typename T0, typename... Ts, typename... As, typename B, size_t... Is>
		auto by_linear_divs(OpBinaryDiv<OpAdd<As...>, B> const& e)
		{
			return check_term_in_num(e.b, e.a, std::make_index_sequence<sizeof...(As)>{});
		}

		template<typename... As, typename... Bs>
		auto by_linear_divs(OpBinaryDiv<OpAdd<As...>, OpAdd<Bs...>> const& e)
		{
			return check_term_in_num(e.b, e.a, std::make_index_sequence<sizeof...(Bs)>{});
		}

		template<typename... Ls, size_t... Is>
		auto add_linear_terms(std::tuple<Ls...> const& linear_terms, std::index_sequence<Is...>)
		{
			return (std::get<Is>(linear_terms) + ...);
		}
	}

	template<typename... Es>
	auto by_linear(OpAdd<Es...> const& e)
	{
		return by_linear_adds(e, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<typename A, typename B>
	auto by_linear(OpBinaryDiv<A, B> const& e)
	{
		auto linear_terms = by_linear_divs(e);
		using seq_t = std::make_index_sequence<std::tuple_size<decltype(linear_terms)>::value>;
		
		auto l = add_linear_terms(linear_terms, seq_t{});
		return std::make_pair(l, e - l);
	}

	template<typename E>
	auto by_linear(OpExpression<E> const& e)
	{
		return pack_right(*static_cast<E const*>(&e));
	}

	// **************************************************************************************



	 // there are several different cases a compound operator might happen:
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
	template<size_t Z, typename... Es>
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
	template<size_t Z, typename E1, typename E2>
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
	auto separate_var(OpFuncDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e);

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
	template<size_t Z, typename V, size_t D, typename E>
	auto separate_var(OpFuncConvolution<V, GaussianSmoothing<D>, E> const& e);

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
	template<size_t Z, auto f, typename V, typename E>
	auto separate_var(OpFuncApply<f, V, E> const& e);


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
	template<size_t Z, typename V, typename E, typename F, typename Arg0, typename... Args>
	auto separate_var(OpFunc<V, E, F, Arg0, Args...> const& e);


	//

	template<size_t Z, typename E>
	auto separate_var(E const& e)
	{
		return pack_right(e);
	}
	
	template<size_t Z, typename T, typename G>
	auto separate_var(OpTerm<T, Variable<Z, G>> const& e)
	{
		return pack_left(e);
	}

	namespace
	{
		/*
		 * Used in avoiding "if constexpr" constructs.
		 */

		template<size_t Z, typename E1, typename E2>
		constexpr bool svc_pred = (expr::vars<E1>::template only_id<Z>() && expr::vars<E2>::template only_id<Z>());
		template<size_t Z, typename E>
		constexpr bool svcg_pred = (expr::vars<E>::template only_id<Z>());
		template<size_t Z, typename E1, typename E2>
		constexpr bool svm_pred = (expr::vars<E1>::template only_id<Z>() && expr::vars<E2>::template only_id<Z>());
		template<size_t Z, typename E1, typename E2>
		constexpr bool svd_pred = (!expr::vars<E1>::template only_id<Z>() && expr::vars<E2>::template only_id<Z>() && !svm_pred<Z, E1, E2>);
		template<size_t Z, typename E>
		constexpr bool svd_pred_1 = (expr::vars<E>::template only_id<Z>());
		template<size_t Z, typename E>
		constexpr bool svd_pred_2 = (expression_predicate<E>::combination && expr::vars<E>::template has_id<Z>());
		template<size_t Z, typename E>
		constexpr bool svcc_pred_1 = (expr::vars<E>::template only_id<Z>());
		template<size_t Z, typename E>
		constexpr bool svcc_pred_2 = (expression_predicate<E>::combination && expr::vars<E>::template has_id<Z>());

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

	namespace
	{
		template<size_t Z, typename... Es, size_t... Is>
		auto separate_var_adds(OpAdd<Es...> const& e, std::index_sequence<Is...>)
		{
			return adds_expand_pair(separate_var<Z>(expr::get<Is>(e))...);
			//return std::make_pair((std::get<0>(std::get<Is>(a)) + ...), (std::get<1>(std::get<Is>(a)) + ...));
		}
	}

	template<size_t Z, typename... Es>
	auto separate_var(OpAdd<Es...> const& e)
	{
		return separate_var_adds<Z>(e, std::make_index_sequence<sizeof...(Es)>{});
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
		template<size_t Z, typename E1, typename E2,
			typename std::enable_if_t<svd_pred<Z, E1, E2>, int> = 0>
		auto separate_var_div(OpBinaryDiv<E1, E2> const& e)
		{
			auto [a, b] = separate_var<Z>(e.a);
			return std::make_pair(a / e.b, b / e.b);
		}

		template<size_t Z, typename E1, typename E2,
			typename std::enable_if_t<!svm_pred<Z, E1, E2>, int> = 0>
		auto separate_var_div(OpBinaryDiv<E1, E2> const& e)
		{
			return pack_right(e);
		}

		template<size_t Z, typename E1, typename E2,
			typename std::enable_if_t<svm_pred<Z, E1, E2>, int> = 0>
		auto separate_var_div(OpBinaryDiv<E1, E2> const& e)
		{
			return pack_left(e);
		}
	}

	template<size_t Z, typename E1, typename E2>
	auto separate_var(OpBinaryDiv<E1, E2> const& e)
	{
		return separate_var_div<Z>(e);
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
			return separate_var<Z>(expr::apply_operators(e));

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
		auto separate_var_derivative_lop(OpFuncDerivative<Dd, V, OpTerm<OpIdentity, G>, T> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, typename Dd, typename V, typename G, typename T,
			typename std::enable_if_t<!svd_pred_1<Z, G>, int> = 0>
		auto separate_var_derivative_lop(OpFuncDerivative<Dd, V, OpTerm<OpIdentity, G>, T> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename Dd, typename V, typename G, typename Sp>
	auto separate_var(OpFuncDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e)
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
			return separate_var<Z>(expr::apply_operators(e));
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
			return separate_var<Z>(expr::apply_operators(e));
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


	namespace
	{
		template<size_t Z, auto f, typename V, typename E,
			typename std::enable_if_t<svcg_pred<Z, E>, int> = 0>
			auto separate_var_function_apply(OpFuncApply<f, V, E> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, auto f, typename V, typename E,
			typename std::enable_if_t<!svcg_pred<Z, E>, int> = 0>
			auto separate_var_function_apply(OpFuncApply<f, V, E> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, auto f, typename V, typename E>
	auto separate_var(OpFuncApply<f, V, E> const& e)
	{
		return separate_var_function_apply<Z>(e);
	}

	namespace
	{
		template<size_t Z, typename V, typename E, typename F, typename Arg0, typename... Args,
			typename std::enable_if_t<svcg_pred<Z, E>, int> = 0>
		auto separate_var_function(OpFunc<V, E, F, Arg0, Args...> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, typename V, typename E, typename F, typename Arg0, typename... Args,
			typename std::enable_if_t<!svcg_pred<Z, E>, int> = 0>
		auto separate_var_function(OpFunc<V, E, F, Arg0, Args...> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename V, typename E, typename F, typename Arg0, typename... Args>
	auto separate_var(OpFunc<V, E, F, Arg0, Args...> const& e)
	{
		return separate_var_function<Z>(e);
	}


	// **************************************************************************************

	//// factors out a derivative order and axis and puts the corresponding operator 
	//// in the first part of a pair when derivatives are factored, which will always be returned with 
	//// OpIdentity coefficient

	//// note that a combination operator derivative can't be factored from a linear combination;
	//// inherent limitation of type-deterministic symbolic algebra (same limitation as for factoring)
	//// hence why only the whole derivative may be factored


	//// if the derivative order is 0 or if the expression has a derivative index which is too
	//// small, then return a pair with the identity in place of the factor

	////! Factor out the derivatives from the expression by axis and order.
	///*!
	// * Factors out the given derivative order. Derivatives of higher orders will
	// * be correspondingly "divided", so that the result which is a pair containing
	// * the expressions `D` and `E` 
	// * will recover the final expression by multiplying `D` to `E`. That is,
	// * `D` will be an expression of only derivative operators, so that the given
	// * expression would be recovered by applying `D` to `E`.
	// *
	// * The derivative order should be less than the minimum order, otherwise
	// * the factorization would be undefined, in which case the original
	// * expression is returned as `E`.
	// * 
	// * \param e The expression which is factored.
	// * 
	// * \tparam OD The order of the derivative that is factored from the expression.
	// */
	//template<Axis ax, size_t OD, typename A1, typename A2, 
	//	typename Enable = typename std::enable_if<min_derivative_order<OpOperatorChain<A1, A2>>::value >= OD>::type>
	//auto factor_deriv(OpOperatorCombination<A1, A2> const& e);

	////! Factor out the derivatives from the expression.
	///*!
	// * Factors out derivatives based on the order and axis. Derivatives of higher orders 
	// * will be correspondingly "divided", so that the result which is a pair containing
	// * the expressions `D` and `E`
	// * will recover the final expression by multiplying `D` to `E`. That is,
	// * `D` will be an expression of only derivative operators, so that the given
	// * expression would be recovered by applying `D` to `E`.
	// *
	// * The derivative order should be less than the minimum order, otherwise
	// * the factorization would be undefined, in which case the original
	// * expression is returned as `E`.
	// *
	// * \param e The expression which is factored.
	// *
	// * \tparam OD The order of the derivative that is factored from the expression.
	// */
	//template<Axis ax, size_t OD, typename A1, typename A2,
	//	typename Enable = typename std::enable_if<min_derivative_order<A1>::value >= OD>::type>
	//auto factor_deriv(OpOperatorChain<A1, A2> const& e);

	////! Factor out the derivatives from the expression.
	///*!
	// * Factors out the given derivative order. Derivatives of higher orders will
	// * be correspondingly "divided", so that the result which is a pair containing
	// * the expressions `D` and `E`
	// * will recover the final expression by multiplying `D` to `E`. That is,
	// * `D` will be an expression of only derivative operators, so that the given
	// * expression would be recovered by applying `D` to `E`.
	// *
	// * The derivative order should be less than the minimum order, otherwise
	// * the factorization would be undefined, in which case the original
	// * expression is returned as `E`.
	// *
	// * \param e The expression which is factored.
	// *
	// * \tparam OD The order of the derivative that is factored from the expression.
	// */
	//template<Axis ax, size_t OD, typename A1, typename A2, typename E>
	//auto factor_deriv(OpCombination<A1, A2, E> const& e);

	////! Factor out the derivatives from the expression.
	///*!
	// * Factors out the given derivative order. Derivatives of higher orders will
	// * be correspondingly "divided", so that the result which is a pair containing
	// * the expressions `D` and `E`
	// * will recover the final expression by multiplying `D` to `E`. That is,
	// * `D` will be an expression of only derivative operators, so that the given
	// * expression would be recovered by applying `D` to `E`.
	// *
	// * The derivative order should be less than the minimum order, otherwise
	// * the factorization would be undefined, in which case the original
	// * expression is returned as `E`.
	// *
	// * \param e The expression which is factored.
	// *
	// * \tparam OD The order of the derivative that is factored from the expression.
	// */
	//template<Axis ax, size_t OD, typename A1, typename A2, typename E>
	//auto factor_deriv(OpChain<A1, A2, E> const& e);

	////! Factor out the derivatives from the expression.
	///*!
	// * Factors out the given derivative order. Derivatives of higher orders will
	// * be correspondingly "divided", so that the result which is a pair containing
	// * the expressions `D` and `E`
	// * will recover the final expression by multiplying `D` to `E`. That is,
	// * `D` will be an expression of only derivative operators, so that the given
	// * expression would be recovered by applying `D` to `E`.
	// *
	// * The derivative order should be less than the minimum order, otherwise
	// * the factorization would be undefined, in which case the original
	// * expression is returned as `E`.
	// *
	// * \param e The expression which is factored.
	// *
	// * \tparam OD The order of the derivative that is factored from the expression.
	// */
	//template<Axis ax, size_t OD, typename... Es,
	//	typename Enable = typename std::enable_if<((min_derivative_order<Es>::value >= OD) && ...)>::type>
	//auto factor_deriv(OpAdd<Es...> const& e);

	//template<Axis ax, size_t OD, size_t O, typename V, typename Sp, typename std::enable_if<(O == OD), int>::type = 0>
	//auto factor_deriv(OpOperatorDerivative<O, V, Sp> const& e)
	//{
	//	return std::make_pair(expr::make_operator_derivative<OD>(e.solver), expr::make_literal(e.value));
	//}

	//template<Axis ax, size_t OD, size_t O, typename V, typename Sp, typename std::enable_if<(O > OD), int>::type = 0>
	//auto factor_deriv(OpOperatorDerivative<O, V, Sp> const& e)
	//{
	//	return std::make_pair(expr::make_operator_derivative<OD>(e.solver), expr::make_operator_derivative<O - OD>(e.value, e.solver));
	//}

	//template<Axis ax, size_t OD, typename Dd, typename V, typename E, typename Sp,
	//	size_t O = OpFuncDerivative<Dd, V, E, Sp>::order,
	//	typename std::enable_if<(O == OD), int>::type = 0>
	//auto factor_deriv(OpFuncDerivative<Dd, V, E, Sp> const& e)
	//{
	//	if constexpr (O % 2 == 0)
	//	{
	//		return std::make_pair(expr::make_operator_derivative<OD>(e.solver), expr::make_literal(e.value) * expr::get_enclosed_expression(e));
	//	}
	//	else
	//	{
	//		return std::make_pair(expr::make_operator_derivative<OD>(e.solver), expr::make_literal(e.value) * expr::get_enclosed_expression(e));
	//	}
	//}

	//template<Axis ax, size_t OD, typename Dd, typename V, typename E, typename Sp,
	//	size_t O = OpFuncDerivative<Dd, V, E, Sp>::order,
	//	typename std::enable_if_t<(O != OD), int>::type = 0>
	//auto factor_deriv(OpFuncDerivative<Dd, V, E, Sp> const& e)
	//{
	//	return std::make_pair(expr::make_operator_derivative<OD>(e.solver), expr::make_operator_derivative<O - OD>(e.value, e.solver) * expr::get_enclosed_expression(e));
	//}


	////! Factor out the derivatives from the expression.
	///*!
	// * Factors out the given derivative order. Derivatives of higher orders will
	// * be correspondingly "divided", so that the result which is a pair containing
	// * the expressions `D` and `E`
	// * will recover the final expression by multiplying `D` to `E`. That is,
	// * `D` will be an expression of only derivative operators, so that the given
	// * expression would be recovered by applying `D` to `E`.
	// *
	// * The derivative order should be less than the minimum order, otherwise
	// * the factorization would be undefined, in which case the original
	// * expression is returned as `E`.
	// *
	// * \param e The expression which is factored.
	// *
	// * \tparam OD The order of the derivative that is factored from the expression.
	// */
	//template<Axis ax, size_t OD, typename E, typename std::enable_if<(OD == 0 || min_derivative_order<E>::value < OD), int>::type = 0>
	//auto factor_deriv(OpExpression<E> const& e)
	//{
	//	return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
	//}


	//template<Axis ax, size_t OD, typename A1, typename A2, typename Enable>
	//auto factor_deriv(OpOperatorCombination<A1, A2> const& e)
	//{
	//	auto a = factor_deriv<ax, OD>(e.f);
	//	auto b = factor_deriv<ax, OD>(e.g);
	//	return std::make_pair(a.first, a.second + b.second);
	//}

	//template<Axis ax, size_t OD, typename A1, typename A2, typename Enable>
	//auto factor_deriv(OpOperatorChain<A1, A2> const& e)
	//{
	//	auto a = factor_deriv<ax, OD>(e.f);
	//	return std::make_pair(a.first, a.second * e.g);
	//}

	//template<Axis ax, size_t OD, typename A1, typename A2, typename E>
	//auto factor_deriv(OpCombination<A1, A2, E> const& e)
	//{
	//	auto fac = factor_deriv<ax, OD>(e.combination);
	//	return std::make_pair(fac.first, fac.second * expr::get_enclosed_expression(e));
	//}

	//template<Axis ax, size_t OD, typename A1, typename A2, typename E>
	//auto factor_deriv(OpChain<A1, A2, E> const& e)
	//{
	//	auto fac = factor_deriv<ax, OD>(e.combination);
	//	return std::make_pair(fac.first, fac.second * expr::get_enclosed_expression(e));
	//}

	//namespace
	//{

	//	template<Axis ax, size_t OD, typename... Es, size_t... Is>
	//	auto factor_deriv_adds(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	//	{
	//		return adds_expand_pair_no_first(factor_deriv<ax, OD>(expr::get<Is>(e))...);
	//		//return std::make_pair(std::get<0>(a).first, (std::get<Is>(a).second + ...));
	//	}
	//}

	//template<Axis ax, size_t OD, typename... Es, typename Enable>
	//auto factor_deriv(OpAdd<Es...> const& e)
	//{
	//	return factor_deriv_adds<ax, OD>(e, std::make_index_sequence<sizeof...(Es)>{});
	//}

	////! Factor out the derivatives from the expression.
	///*!
	// * Factors out the given derivative order. Derivatives of higher orders will
	// * be correspondingly "divided", so that the result which is a pair containing
	// * the expressions `D` and `E`
	// * will recover the final expression by multiplying `D` to `E`. That is,
	// * `D` will be an expression of only derivative operators, so that the given
	// * expression would be recovered by applying `D` to `E`.
	// *
	// * The derivative order must be less than the minimum order, otherwise
	// * the factorization would be undefined.
	// *
	// * \param e The expression which is factored.
	// *
	// * \tparam OD The order of the derivative that is factored from the expression.
	// */
	//template<Axis ax, typename E>
	//auto factor_deriv(OpExpression<E> const& e)
	//{
	//	return factor_deriv<ax, min_derivative_order<E>::value>(*static_cast<E const*>(&e));
	//}


	// **************************************************************************************


	//! Separates an operator from the expression, if there is one applied to it.
	/*!
	 * The operator is returned in the first entry of the pair, and the expression it is
	 * applied to is returned in the second. The original expression that is passed to this
	 * function is recovered by applying the operator to the enclosed expression. For example,
	 * given an expression U which is obtained by applying an operator D to an enclosed expression
	 * E, then U = D(E). Note that in general, it is NOT the case that U = D * U, since this is
	 * taking the dot product.
	 * 
	 * The return value of this function will be (D, E), so that U is recovered by applying D 
	 * to E.
	 */
	template<typename E>
	auto separate_operator(OpExpression<E> const& e)
	{
		return pack_right(*static_cast<E const*>(&e));
	}

	//! Separates an operator from the expression, if there is one applied to it.
	/*!
	 * The operator is returned in the first entry of the pair, and the expression it is
	 * applied to is returned in the second. The original expression that is passed to this
	 * function is recovered by applying the operator to the enclosed expression. For example,
	 * given an expression U which is obtained by applying an operator D to an enclosed expression
	 * E, then U = D(E). Note that in general, it is NOT the case that U = D * U, since this is
	 * taking the dot product.
	 *
	 * The return value of this function will be (D, E), so that U is recovered by applying D
	 * to E.
	 */
	template<typename Dd, typename V, typename E, typename Sp, size_t R = expr::eval_type<E>::rank>
	auto separate_operator(OpFuncDerivative<Dd, V, E, Sp> const& e)
	{
		constexpr size_t order = OpFuncDerivative<Dd, V, E, Sp>::order;
		constexpr Axis axis = OpFuncDerivative<Dd, V, E, Sp>::axis;
		if constexpr (order % 2 == 1)
		{
			auto op_ax = expr::make_operator_directional_derivative<axis, 1>(e.solver);
			if constexpr (order == 1)
			{
				return std::make_pair(op_ax, expr::coeff(e) * expr::get_enclosed_expression(e));
			}
			else
			{
				auto op = op_ax * expr::make_operator_derivative<order - 1>(e.solver);
				return std::make_pair(op, expr::coeff(e) * expr::get_enclosed_expression(e));
			}
		}
		else
		{
			return std::make_pair(expr::make_operator_derivative<order>(e.solver), expr::coeff(e) * expr::get_enclosed_expression(e));
		}
	}


	//! Separates an operator from the expression, if there is one applied to it.
	/*!
	 * The operator is returned in the first entry of the pair, and the expression it is
	 * applied to is returned in the second. The original expression that is passed to this
	 * function is recovered by applying the operator to the enclosed expression. For example,
	 * given an expression U which is obtained by applying an operator D to an enclosed expression
	 * E, then U = D(E). Note that in general, it is NOT the case that U = D * U, since this is
	 * taking the dot product.
	 *
	 * The return value of this function will be (D, E), so that U is recovered by applying D
	 * to E.
	 */
	template<Axis ax, size_t O, typename A1, typename A2, typename E,
		typename std::enable_if<min_derivative_order<OpOperatorCombination<A1, A2>>::value == O, int>::type = 0>
	auto separate_operator(OpCombination<A1, A2, E> const& e)
	{
		return std::make_pair(e.combination, e.e);
	}

	//! Separates an operator from the expression, if there is one applied to it.
	/*!
	 * The operator is returned in the first entry of the pair, and the expression it is
	 * applied to is returned in the second. The original expression that is passed to this
	 * function is recovered by applying the operator to the enclosed expression. For example,
	 * given an expression U which is obtained by applying an operator D to an enclosed expression
	 * E, then U = D(E). Note that in general, it is NOT the case that U = D * U, since this is
	 * taking the dot product.
	 *
	 * The return value of this function will be (D, E), so that U is recovered by applying D
	 * to E.
	 */
	template<typename A1, typename A2, typename E>
	auto separate_operator(OpCombination<A1, A2, E> const& e)
	{
		return separate_operator<Axis::X, 0>(e);
	}

	//! Separates an operator from the expression, if there is one applied to it.
	/*!
	 * The operator is returned in the first entry of the pair, and the expression it is
	 * applied to is returned in the second. The original expression that is passed to this
	 * function is recovered by applying the operator to the enclosed expression. For example,
	 * given an expression U which is obtained by applying an operator D to an enclosed expression
	 * E, then U = D(E). Note that in general, it is NOT the case that U = D * U, since this is
	 * taking the dot product.
	 *
	 * The return value of this function will be (D, E), so that U is recovered by applying D
	 * to E.
	 */
	template<Axis ax, size_t O, typename A1, typename A2, typename E,
		typename std::enable_if<min_derivative_order<OpOperatorChain<A1, A2>>::value == O, int>::type = 0>
	auto separate_operator(OpChain<A1, A2, E> const& e)
	{
		return std::make_pair(e.combination, e.e);
	}

	//! Separates an operator from the expression, if there is one applied to it.
	/*!
	 * The operator is returned in the first entry of the pair, and the expression it is
	 * applied to is returned in the second. The original expression that is passed to this
	 * function is recovered by applying the operator to the enclosed expression. For example,
	 * given an expression U which is obtained by applying an operator D to an enclosed expression
	 * E, then U = D(E). Note that in general, it is NOT the case that U = D * U, since this is
	 * taking the dot product.
	 *
	 * The return value of this function will be (D, E), so that U is recovered by applying D
	 * to E.
	 */
	template<typename A1, typename A2, typename E>
	auto separate_operator(OpChain<A1, A2, E> const& e)
	{
		return separate_operator<Axis::X, 0>(e);
	}


	// **************************************************************************************


	namespace
	{


		// factors the expression by the given variable
		// the OpTerms and opnlvariables with matching type C are factored
		// N is the number of types of C to factor out

		template<typename C>
		struct divide_variable
		{
			template<typename T, typename G, typename std::enable_if<(expr::factor_count<C, G>::value > 0), int>::type = 0>
				auto operator()(OpTerm<T, G> const& e)
			{
				return std::make_pair(OpTerm<OpIdentity, G>(e.data), expr::make_literal(e.value));
			}

			template<typename T, typename G, typename std::enable_if<(expr::factor_count<C, G>::value == 0), int>::type = 0>
			auto operator()(OpTerm<T, G> const& e)
			{
				return std::make_pair(OpIdentity{}, e);
			}
		};

		template<size_t Z>
		struct divide_variable<Variable<Z>>
		{
			template<typename T, typename G>
			auto operator()(OpTerm<T, Variable<Z, G>> const& e)
			{
				return std::make_pair(OpTerm<OpIdentity, Variable<Z, G>>(e.data), expr::make_literal(e.value));
			}

			template<typename T, typename G>
			auto operator()(OpTerm<T, G> const& e)
			{
				return std::make_pair(OpIdentity{}, e);
			}

		};


		template<size_t N, typename C, typename... Es, 
			typename Enable = typename std::enable_if<(expr::factor_count<C, OpAdd<Es...>>::value >= N)>::type>
		auto _factor(OpAdd<Es...> const& e);
		template<size_t N, typename C, typename E1, typename E2, 
			typename Enable = typename std::enable_if<(expr::factor_count<C, OpBinaryMul<E1, E2>>::value >= N)>::type>
		auto _factor(OpBinaryMul<E1, E2> const& e);
		template<size_t N, typename C, typename E1, typename E2, 
			typename Enable = typename std::enable_if<(expr::factor_count<C, OpBinaryDiv<E1, E2>>::value >= N)>::type>
		auto _factor(OpBinaryDiv<E1, E2> const& e);
		template<size_t N, typename C, typename V, typename... Gs, exp_key_t... Xs>
		auto _factor(OpTerms<V, Term<Gs, Xs>...> const& e);



		template<size_t N, typename C, typename E,
			typename std::enable_if<(N == 0 || expr::factor_count<C, E>::value < N), int>::type = 0>
		auto _factor(OpExpression<E> const& e)
		{
			return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
		}


		// Helper functions for factoring.
		// **************************************************************************************

		template<size_t N, typename C, typename... Es, size_t... Is>
		auto _factor_adds(OpAdd<Es...> const& e, std::index_sequence<Is...>)
		{
			return adds_expand_pair_no_first(_factor<N, C>(expr::get<Is>(e))...);
			//return std::make_pair(std::get<0>(a).first, expr::add_all(std::get<Is>(a).second...));
		}


		template<typename V, typename G0, exp_key_t X0, typename... Gs, exp_key_t... Xs>
		auto _make_terms(V value, Term<G0, X0> const& term0, Term<Gs, Xs> const&... rest)
		{
			return OpTerms(value, term0, rest...);
		}

		template<typename V, typename G0, typename... Gs, exp_key_t... Xs>
		auto _make_terms(V value, Term<G0, 0> const& term0, Term<Gs, Xs> const&... rest)
		{
			return expr::transform::sift_term(value, rest...);
		}

		template<size_t N, typename V, typename... Gs, exp_key_t... Xs, size_t... Is>
		auto _select_terms(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<Is...>, std::index_sequence<>)
		{
			return std::pair(OpIdentity{}, e);
		}

		template<size_t N, typename C, typename V, typename... Gs, exp_key_t... Xs, size_t... Is, size_t I>
		auto _select_terms(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<Is...>, std::index_sequence<I>)
		{
			// If the type can factor another type more than once (e.g. |k|^2 factors |k|^4 twice), then if we want to
			// factor |k|^2 from |k|^4, we need to divide |k|^4 by |k|^2.
			constexpr size_t N0 = expr::factor_count<C, symphas::lib::type_at_index<I - 1, Gs...>>::value - 1;

			auto factor_data = expr::get<I>(e).data();
			auto factor_term = (Term(factor_data) * ~(Term<C, 1>(factor_data).template pow<N0>())).template pow<N>();
			auto nonfactor_term = expr::get<I>(e) * (~factor_term);

			return std::make_pair(_make_terms(OpIdentity{}, factor_term), _make_terms(e.term, nonfactor_term, expr::get<Is>(e)...));

		}



		template<size_t N, typename C, typename V, typename... Gs, exp_key_t... Xs, size_t... Is, bool... fs>
		auto _factor(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<Is...>, std::integer_sequence<bool, fs...>)
		{
			using symphas::lib::seq_join_t;
			using pick_factor_t = seq_join_t<
				std::index_sequence<>,
				std::conditional_t<
					fs,
					std::index_sequence<Is + 1>,
					std::index_sequence<>>...
				>;
			using pick_nonfactor_t = seq_join_t<
				std::index_sequence<>,
				std::conditional_t<
					!fs,
					std::index_sequence<Is + 1>,
					std::index_sequence<>>...
				>;
			return _select_terms<N, C>(e, pick_nonfactor_t{}, pick_factor_t{});
		}

		template<size_t N, typename C, typename V, typename... Gs, exp_key_t... Xs>
		auto _factor(OpTerms<V, Term<Gs, Xs>...> const& e)
		{
			using seq_t = std::make_index_sequence<sizeof...(Gs)>;
			using mask_t = std::integer_sequence<bool, (factor_count<C, Term<Gs, Xs>>::value >= N && expr::is_combinable<Gs>)...>;
			return _factor<N, C>(e, seq_t{}, mask_t{});
		}


		template<size_t N, typename C, typename... Es, typename Enable>
		auto _factor(OpAdd<Es...> const& e)
		{
			return _factor_adds<N, C>(e, std::make_index_sequence<sizeof...(Es)>{});
		}

		template<size_t N, typename C, typename E1, typename E2, 
			size_t AN = expr::factor_count<C, E1>::value, typename std::enable_if_t<(N > AN), int> = 0>
		auto _factor_sift(OpBinaryMul<E1, E2> const& e)
		{
			auto a = _factor<AN, C>(e.a);
			auto b = _factor<N - AN, C>(e.b);
			return std::make_pair(a.first * b.first, a.second * b.second);
		}

		template<size_t N, typename C, typename E1, typename E2,
			size_t AN = expr::factor_count<C, E1>::value, typename std::enable_if_t<(N <= AN), int> = 0>
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
		struct is_variable_data_factor<Variable<Z>>
		{
			static const bool value = true;
		};

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


		// automatically selects the largest possible number of datas to factor
		template<typename C0, typename E,
			typename std::enable_if<(expr::grid_can_combine<C0>::value || is_variable_data_factor<C0>::value), int>::type = 0>
		auto _factor(OpExpression<E> const& e)
		{
			return _factor<expr::factor_count<C0, E>::value, C0>(*static_cast<E const*>(&e));
		}

		template<typename C0, typename E,
			typename std::enable_if<(!expr::grid_can_combine<C0>::value && !is_variable_data_factor<C0>::value), int>::type = 0>
		auto _factor(OpExpression<E> const& e)
		{
			return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
		}

	}


	//! Factor an expression by the given terms.
	/*!
	 * Factor an expression by the given terms which are represented by
	 * the given types, not passed explicitly.
	 * Returns the result of factoring the given term in the first entry of the pair,
	 * and the product of all the factored terms in the second entry of the pair.
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
	 * Returns the result of factoring the given term in the first entry of the pair,
	 * and the product of all the factored terms in the second entry of the pair.
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
	 * Returns the result of factoring the given term in the first entry of the pair,
	 * and the product of all the factored terms in the second entry of the pair.
	 *
	 * The types to factor by must all be unique. 
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam N The number of times to factor each type out.
	 * \tparam C0 The first term to factor.
	 * \tparam Cs... The remaining terms to factor.
	 */
	template<size_t N, typename C0, typename E, typename std::enable_if_t<(!factor_pred<C0>::value), int> = 0>
	auto factor(OpExpression<E> const& e)
	{
		constexpr size_t min_order = fixed_min<expr::factor_count<C0, E>::value, N>;
		if constexpr (min_order > 0)
		{
			auto a = _factor<min_order, C0>(*static_cast<const E*>(&e));
			return std::make_pair(a.first, a.second);
		}
		else
		{
			return std::make_pair(OpIdentity{}, *static_cast<const E*>(&e));
		}
	}

	//! Factor an expression by the given terms the given number of times.
	/*!
	 * Factor an expression by the given terms which are represented by
	 * the given types, not passed explicitly. Each term is factored out
	 * the given number of times.
	 * Returns the result of factoring the given term in the first entry of the pair,
	 * and the product of all the factored terms in the second entry of the pair.
	 *
	 * The types to factor by must all be unique.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam N The number of times to factor each type out.
	 * \tparam C0 The first term to factor.
	 * \tparam Cs... The remaining terms to factor.
	 */
	template<size_t N, typename C0, typename C1, typename... Cs, typename E, typename std::enable_if_t<(!factor_pred<C0, C1, Cs...>::value), int> = 0>
	auto factor(OpExpression<E> const& e)
	{
		constexpr size_t min_order = fixed_min<expr::factor_count<C0, E>::value, N>;
		if constexpr (min_order > 0)
		{
			auto a = _factor<min_order, C0>(*static_cast<const E*>(&e));
			auto b = factor<N, C1, Cs...>(a.second);
			return std::make_pair(a.first * b.first, b.second);
		}
		else
		{
			return factor<N, C1, Cs...>(*static_cast<const E*>(&e));
		}
	}

	//! Factor an expression by the given terms the given number of times.
	/*!
	 * Factor an expression by the given terms which are represented by
	 * the given types, not passed explicitly. Each term is factored out
	 * the given number of times.
	 * Returns the result of factoring the given term in the first entry of the pair,
	 * and the product of all the factored terms in the second entry of the pair.
	 *
	 * The types to factor by must all be unique.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam N The number of times to factor each type out.
	 * \tparam C0 The first term to factor.
	 * \tparam Cs... The remaining terms to factor.
	 */
	template<size_t N, size_t Z0, size_t... Zs, typename E, typename std::enable_if_t<(!factor_pred<Variable<Z0>, Variable<Zs>...>::value), int> = 0>
	auto factor(OpExpression<E> const& e)
	{
		return factor<N, Variable<Z0>, Variable<Zs>...>(*static_cast<E const*>(&e));
	}


	//! Make a list of the result of each factorization.
	/*!
	 * For each given term, the original given expression is factored. The
	 * results are concatenated into a list.
	 * Returns the result of factoring the given term in the first entry of the pair,
	 * and the product of all the factored terms in the second entry of the pair.
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
	 * Returns the result of factoring the given term in the first entry of the pair,
	 * and the product of all the factored terms in the second entry of the pair.
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
	 * The expression is factored only once by each given term.
	 * 
	 * All the factors that are removed from the expression are in the first entry
	 * of the returned pair, and the factored expression is given in the
	 * second entry.
	 * 
	 * That is: Given a list of factors and an expression *E*, this function returns *A*, *B*:
	 * - *A* is a product of the factors that have been taken from *E*.
	 * - *B* is the factored expression, such that *E* = *A* * *B*.
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
	 * The expression is factored only once by each given term.
	 *
	 * All the factors that are removed from the expression are in the first entry
	 * of the returned pair, and the factored expression is given in the
	 * second entry.
	 * 
	 * That is: Given a list of factors and an expression *E*, this function returns *A*, *B*:
	 * - *A* is a product of the factors that have been taken from *E*.
	 * - *B* is the factored expression, such that *E* = *A* * *B*.
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
	 * The expression is factored only once by each given term.
	 *
	 * All the factors that are removed from the expression are in the first entry
	 * of the returned pair, and the factored expression is given in the
	 * second entry.
	 * 
	 * That is: Given a list of factors and an expression *E*, this function returns *A*, *B*:
	 * - *A* is a product of the factors that have been taken from *E*.
	 * - *B* is the factored expression, such that *E* = *A* * *B*.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam Z0 The first variable index to factor.
	 * \tparam Cs... The remaining variable indices to factor.
	 */
	template<size_t Z0, size_t... Zs, typename E, typename std::enable_if_t<(sizeof...(Zs) == 0), int> = 0>
	auto factor_list(OpExpression<E> const& e)
	{
		return factor_list<Variable<Z0>, Variable<Zs>...>(*static_cast<const E*>(&e));
	}


}

namespace symphas::internal
{
	template<typename E1, typename E2>
	auto remove_factors(E1 const& e1, E2 const& e2, symphas::lib::types_list<>)
	{
		return std::make_pair(e1, e2);
	}

	template<size_t N01, typename G01, size_t... N1s, typename... G1s, typename E1, typename E2>
	auto remove_factors(E1 const& e1, E2 const& e2, symphas::lib::types_list<std::pair<std::index_sequence<N01>, G01>, std::pair<std::index_sequence<N1s>, G1s>...>)
	{
		auto f = expr::split::factor<N01, G01>(e1);
		auto g = expr::split::factor<N01, G01>(e2);
		return remove_factors(f.second, g.second, symphas::lib::types_list<std::pair<std::index_sequence<N1s>, G1s>...>{});
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
	struct divide_with_factors;

	//! Specialization of expr::divide_with_factors with no factors given.
	template<>
	struct divide_with_factors<symphas::lib::types_list<>>
	{
		template<typename E1, typename E2>
		auto operator()(OpExpression<E1> const& a, OpExpression<E2> const& b)
		{
			return symphas::internal::terminate_div(*static_cast<E1 const*>(&a), expr::inverse(*static_cast<E2 const*>(&b)));
		}
	};

	//! Specialization of expr::divide_with_factors a list of factors given.
	template<size_t N0, typename G0, size_t... Ns, typename... Gs>
	struct divide_with_factors<symphas::lib::types_list<std::pair<std::index_sequence<N0>, G0>, std::pair<std::index_sequence<Ns>, Gs>...>>
	{
		template<typename E1, typename E2>
		auto operator()(OpExpression<E1> const& a, OpExpression<E2> const& b)
		{
			using factors_t = symphas::lib::types_list<std::pair<std::index_sequence<N0>, G0>, std::pair<std::index_sequence<Ns>, Gs>...>;
			auto [numerator, denominator] = 
				symphas::internal::remove_factors(
					*static_cast<E1 const*>(&a), 
					*static_cast<E2 const*>(&b),
					factors_t{});
			
			return symphas::internal::terminate_div(numerator, expr::inverse(denominator));
		}
		
	};


}

/*
 *
 *
 * Division rules
 *
 ******************************************************************************/

//! The division operator is overloaded to apply factoring in general.
template<typename E1, typename E2>
auto operator/(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::divide_with_factors<typename expr::factor_list_all<E1, E2>::type>{}(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
}

/*
 *
 *
 * Additional overloads to avoid ambiguous overloads
 *
 ******************************************************************************/

namespace expr
{

	namespace
	{
		template<typename G0, typename E>
		struct derivative_apply
		{
			auto operator()(symphas::internal::wrap_base, E const& e)
			{
				return OpVoid{};
			}

			auto operator()(symphas::internal::wrap_f<func_cos<E>>, E const& e)
			{
				return -OpFuncApply<func_sin<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e));
			}

			auto operator()(symphas::internal::wrap_f<func_sin<E>>, E const& e)
			{
				return OpFuncApply<func_cos<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e));
			}

			auto operator()(symphas::internal::wrap_f<func_tan<E>>, E const& e)
			{
				return OpFuncApply<func_sec<E>, OpIdentity, E>(e) * OpFuncApply<func_sec<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e));
			}

			auto operator()(symphas::internal::wrap_f<func_csc<E>>, E const& e)
			{
				return -OpFuncApply<func_cot<E>, OpIdentity, E>(e) * OpFuncApply<func_csc<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e));
			}

			auto operator()(symphas::internal::wrap_f<func_sec<E>>, E const& e)
			{
				return OpFuncApply<func_tan<E>, OpIdentity, E>(e) * OpFuncApply<func_sec<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e));
			}

			auto operator()(symphas::internal::wrap_f<func_cot<E>>, E const& e)
			{
				return -OpFuncApply<func_csc<E>, OpIdentity, E>(e) * OpFuncApply<func_csc<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e));
			}

			auto operator()(symphas::internal::wrap_f<func_sqrt<E>>, E const& e)
			{
				return expr::make_fraction<1, 2>() *
					expr::inverse(OpFuncApply<func_sqrt<E>, OpIdentity, E>(e)) *
					apply_operators(expr::make_derivative<1, G0>(e));
			}

			auto operator()(symphas::internal::wrap_f<func_log<E>>, E const& e)
			{
				return expr::inverse(e);
			}
		};
	}



	template<size_t O, typename V0, auto f, typename V1, typename E, typename G0>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V0, OpFuncApply<f, V1, E>, SymbolicDerivative<G0>> const& e)
	{
		auto&& function = expr::get_enclosed_expression(e);
		auto&& expr = expr::get_enclosed_expression(function);

		return apply_operators(
			expr::make_derivative<O - 1, G0>(
				expr::make_literal(e.value) * expr::make_literal(function.value),
				derivative_apply<G0, E>{}(symphas::internal::wrap_f<f>{}, expr)
				));
	}


	template<size_t O, typename V, typename V1, typename... Gs, expr::exp_key_t... Xs, typename G0, size_t N, typename>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpTerms<V1, Term<Gs, Xs>...>, SymbolicDerivative<OpTerm<OpIdentity, G0>>> const& e)
	{
		auto f = expr::split::factor<O, G0>(expr::get_enclosed_expression(e));
		return expr::factorial<N, N - O>() * expr::make_literal(e.value) * f.second;
	}

	template<size_t O, typename V, typename V1, typename... Gs, expr::exp_key_t... Xs, typename G0, size_t N, typename>
	auto apply_operators(OpFuncDerivative<std::index_sequence<O>, V, OpTerms<V1, Term<Gs, Xs>...>, SymbolicDerivative<G0>> const& e)
	{
		auto f = expr::split::factor<O, G0>(expr::get_enclosed_expression(e));
		return expr::factorial<N, N - O>() * expr::make_literal(e.value) * f.second;
	}

}




namespace symphas::internal
{
	template<size_t I, typename T>
	using itype = T;

	template<typename E1, typename E2, size_t... Rs, size_t R = sizeof...(Rs)>
	auto dot_tensor_components(OpExpression<E1> const& a, OpExpression<E2> const& b, std::index_sequence<Rs...>)
	{
		return (
			((expr::make_row_vector<Rs, R>() * (*static_cast<E1 const*>(&a))) 
				* (expr::make_row_vector<Rs, R>() * (*static_cast<E2 const*>(&b))))
			+ ...);

	}

	template<typename E1, typename E2, size_t... Rs, size_t R = sizeof...(Rs)>
	auto mul_tensor_components_rc(OpExpression<E1> const& a, OpExpression<E2> const& b, std::index_sequence<Rs...>)
	{
		return (
			(((*static_cast<E1 const*>(&a)) * expr::make_column_vector<Rs, R>()) 
				* (expr::make_row_vector<Rs, R>() * (*static_cast<E2 const*>(&b))))
			+ ...);
	}

	template<size_t R, size_t P, typename E1, typename E2, size_t... Rs, size_t... Ps>
	auto _mul_tensor_components_cr(OpExpression<E1> const& a, OpExpression<E2> const& b, std::index_sequence<Rs...>, std::index_sequence<Ps...>)
	{
		return ((expr::make_tensor<Rs, Ps, R, P>() *
			((expr::make_row_vector<Rs, R>() * (*static_cast<E1 const*>(&a)))
				* ((*static_cast<E2 const*>(&b)) * expr::make_column_vector<Ps, P>())))
			+ ...);
	}

	template<typename E1, typename E2, size_t... Rs, size_t R = sizeof...(Rs)>
	auto mul_tensor_components_cr(OpExpression<E1> const& a, OpExpression<E2> const& b, std::index_sequence<Rs...>)
	{
		return _mul_tensor_components_cr<R, R>(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b), 
			symphas::lib::seq_join_t<std::index_sequence<>, itype<Rs, std::index_sequence<Rs...>>...>{},
			symphas::lib::seq_join_t<std::index_sequence<>, symphas::lib::seq_repeating_value_t<R, size_t, Rs>...>{});
	}

	template<size_t R0, size_t R, size_t P0, size_t P, typename E1, typename E2, size_t... Qs, size_t Q = sizeof...(Qs)>
	auto get_dot_product_at(OpExpression<E1> const& a, OpExpression<E2> const& b, std::index_sequence<Qs...>)
	{
		return ((((expr::make_row_vector<R0, R>() * (*static_cast<E1 const*>(&a))) * expr::make_column_vector<Qs, Q>())
			* (expr::make_row_vector<Qs, Q>() * ((*static_cast<E2 const*>(&b)) * expr::make_column_vector<P0, P>()))) + ...);
	}

	template<size_t R, size_t P, typename E1, typename E2, size_t... Rs, size_t... Qs, size_t... Ps, size_t Q = sizeof...(Qs)>
	auto _mul_tensor_components(OpExpression<E1> const& a, OpExpression<E2> const& b, 
		std::index_sequence<Rs...>, std::index_sequence<Qs...>, std::index_sequence<Ps...>)
	{
		return ((expr::make_tensor<Rs, Ps, R, P>() *
			get_dot_product_at<Rs, R, Ps, P>(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b), std::index_sequence<Qs...>{})) + ...);
	}

	template<typename E1, typename E2, size_t... Rs, size_t... Qs, size_t... Ps,
		size_t R = sizeof...(Rs), size_t Q = sizeof...(Qs), size_t P = sizeof...(Ps)>
	auto mul_tensor_components(OpExpression<E1> const& a, OpExpression<E2> const& b, 
		std::index_sequence<Rs...>, std::index_sequence<Qs...>, std::index_sequence<Ps...>)
	{
		return _mul_tensor_components<R, P>(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b),
			symphas::lib::seq_join_t<std::index_sequence<>, itype<Ps, std::index_sequence<Rs...>>...>{},
			std::index_sequence<Qs...>{},
			symphas::lib::seq_join_t<std::index_sequence<>, symphas::lib::seq_repeating_value_t<R, size_t, Ps>...>{});
	}

	template<typename E1, typename E2, 
		size_t R1 = expr::eval_type<E1>::rank, size_t R2 = expr::eval_type<E2>::rank,
		size_t Q1 = expr::eval_type<E1>::template rank_<1>, size_t Q2 = expr::eval_type<E2>::template rank_<1>>
	auto dot(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		// multiply a row type by a column type
		if constexpr (R1 == 1 && Q2 == 1 && Q1 == R2 && R2 > 1)
		{
			return mul_tensor_components_rc(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b), std::make_index_sequence<R2>{});
		}
		// multiply a column type by a row type
		else if constexpr (R2 == 1 && Q1 == 1 && Q2 == R1 && R1 > 1)
		{
			return mul_tensor_components_cr(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b), std::make_index_sequence<R1>{});
		}
		// dot product of two column vector types
		else if constexpr (R1 == R2 && Q1 == 1 && Q2 == 1)
		{
			return dot_tensor_components(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b), std::make_index_sequence<R1>{});
		}
		// matrix multiplication
		else if constexpr (Q1 == R2)
		{
			return mul_tensor_components(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b), 
				std::make_index_sequence<R1>{}, std::make_index_sequence<Q1>{}, std::make_index_sequence<Q2>{});
		}
		// multiplying incompatible types just gives their multiplication
		else
		{
			return expr::make_mul(a, b);
		}
	}
}

namespace expr
{
	template<typename E1, typename E2>
	auto dot(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return symphas::internal::dot(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}
}

/*
 *
 * Multiplication between two expressions that are vector type
 * is assumed to be the dot product.
 *
 ******************************************************************************/

template<typename E1, typename E2,
	typename std::enable_if_t<(expr::eval_type<E1>::rank > 0 && expr::eval_type<E2>::rank > 0), int> = 0>
	auto operator*(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::dot(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}


namespace symphas::internal
{

	template<size_t Z, typename G, typename... Gs>
	auto filter_variables(std::tuple<Variable<Z, G>, Gs...> const& data_list);
	template<size_t Z0, typename G0>
	auto sort_variables(std::tuple<Variable<Z0, G0>> const& data_list);

	inline auto filter_variables(std::tuple<> const& data_list)
	{
		return std::tuple();
	}

	template<typename G, typename... Gs>
	auto filter_variables(std::tuple<G, Gs...> const& data_list)
	{
		return filter_variables(symphas::lib::get_tuple_ge<1>(data_list));
	}

	template<size_t Z, typename G, typename... Gs>
	auto filter_variables(std::tuple<Variable<Z, G>, Gs...> const& data_list)
	{
		return std::tuple_cat(std::make_tuple(std::get<0>(data_list)),
			filter_variables(symphas::lib::get_tuple_ge<1>(data_list)));
	}

	inline auto sort_variables(std::tuple<> const& data_list)
	{
		return std::tuple();
	}

	template<size_t Z0, typename G0>
	auto sort_variables(std::tuple<Variable<Z0, G0>> const& data_list)
	{
		return data_list;
	}

	template<size_t Z0, typename G0>
	auto sort_variables(Variable<Z0, G0> const& data0, std::tuple<> const& data_list)
	{
		return std::make_tuple(data0);
	}

	template<size_t Z0, size_t Z1, typename G0, typename G1, typename... Gs>
	auto sort_variables(Variable<Z0, G0> const& data0, std::tuple<Variable<Z1, G1>, Gs...> const& data_list)
	{
		if constexpr (Z0 > Z1)
		{
			return std::tuple_cat(std::make_tuple(std::get<0>(data_list), data0), symphas::lib::get_tuple_ge<1>(data_list));
		}
		else
		{
			return std::tuple_cat(std::make_tuple(data0, std::get<0>(data_list)), symphas::lib::get_tuple_ge<1>(data_list));
		}
	}

	template<size_t Z0, size_t Z1, typename G0, typename G1, typename... Gs>
	auto sort_variables(std::tuple<Variable<Z0, G0>, Variable<Z1, G1>, Gs...> const& data_list)
	{
		if constexpr (Z0 > Z1)
		{
			return sort_variables(
				std::get<1>(data_list),
				sort_variables(std::get<0>(data_list), sort_variables(symphas::lib::get_tuple_ge<2>(data_list))));
		}
		else
		{
			return sort_variables(
				std::get<0>(data_list),
				sort_variables(std::get<1>(data_list), sort_variables(symphas::lib::get_tuple_ge<2>(data_list))));
		}
	}



	template<typename... Gs>
	auto index_variables(std::tuple<Gs...> const& data_list, std::index_sequence<>)
	{
		return std::tuple();
	}

	template<size_t I0, size_t... Is>
	auto index_variables(std::tuple<> const& data_list, std::index_sequence<I0, Is...>)
	{
		return std::make_tuple(OpVoid{}, itype<Is, OpVoid>{}...);
	}

	template<size_t Z, Axis ax, typename G>
	auto remove_component(Variable<Z, VectorComponent<ax, G>> const& data0)
	{
		return Variable<Z, G>(*static_cast<G const*>(&data0));
	}

	template<size_t Z, typename G>
	auto remove_component(Variable<Z, G> const& data0)
	{
		return data0;
	}

	template<size_t Z0, typename G0, typename... Gs, size_t I0, size_t... Is>
	auto index_variables(std::tuple<Variable<Z0, G0>, Gs...> const& data_list, std::index_sequence<I0, Is...>)
	{
		if constexpr (Z0 == I0)
		{
			return std::tuple_cat(
				std::make_tuple(remove_component(std::get<0>(data_list))),
				index_variables(symphas::lib::get_tuple_ge<1>(data_list), std::index_sequence<Is...>{}));
		}
		else if (Z0 < I0)
		{
			return std::tuple_cat(
				std::make_tuple(OpVoid{}),
				index_variables(symphas::lib::get_tuple_ge<1>(data_list), std::index_sequence<Is...>{}));
		}
		else
		{
			return std::tuple_cat(
				std::make_tuple(OpVoid{}),
				index_variables(data_list, std::index_sequence<Is...>{}));
		}
	}

	template<size_t... Zs, typename... Gs>
	auto index_variables(std::tuple<Variable<Zs, Gs>...> const& data_list)
	{
		constexpr size_t Zm = symphas::lib::seq_index_value<sizeof...(Zs) - 1, std::index_sequence<Zs...>>::value;
		return index_variables(data_list, std::make_index_sequence<Zm + 1>{});

	}
}

namespace expr
{
	//! Returns a list of all the variables from the expression, sorted.
	/*!
	 * A list of all the variables from the OpTerm elements in the expression is aggregrated,
	 * duplicates are removed and then the list is sorted by variable ID, i.e. `Z` is the ID
	 * in `Variable<Z, G>`. The variables are placed in a tuple according to their ID, and OpVoid
	 * is placed where no ID exists.
	 * 
	 * \param e The expression in which to search variables.
	 */
	template<typename E>
	auto get_indexed_variable_list(OpExpression<E> const& e)
	{
		return symphas::internal::index_variables(
			symphas::internal::sort_variables(
				symphas::internal::filter_variables(
					expr::data_list(*static_cast<E const*>(&e)))));
	}
}








