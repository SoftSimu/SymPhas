
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

//#include "expressions.h"
#include "expressionderivatives.h"
#include "expressionintegrals.h"
#include "expressionexponentials.h"
#include "expressionfunctions.h"
#include "expressionconvolution.h"
#include "expressionrules.h"
#include "symbolicdata.h"
//#include "expressionproperties.h"

namespace expr
{
	//! Prune an expression to update the state of all nested terms.
	/*! 
	 * Specifies algorithms for different expression types that will update
	 * the underlying data before the expression can be evaluated.
	 */
	namespace prune {}
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
	template<typename A2>
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
	template<typename Dd1, typename Dd2, typename V1, typename V2, typename E, typename Sp1, typename Sp2>
	auto apply_operators(OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp1>, Sp2> const& e);

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
	auto apply_operators(OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp>, Sp> const& e);

	//! Distribute operators so they are applied to individual expressions.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<size_t O1, typename V, size_t O2, typename V1, typename E, typename G00>
	auto apply_operators(OpDerivative<std::index_sequence<O1>, V, OpDerivative<std::index_sequence<O2>, V1, E, SymbolicDerivative<G00>>, SymbolicDerivative<G00>> const& e);

	//! Distribute operators so they are applied to individual expressions.
	/*!
	 * For expressions that are derivatives of derivatives, the outermost
	 * derivatives might need to be distributed to the rest of the expression.
	 * Additionally, all derivatives which are applied to expressions that
	 * are linear combinations are distributed.
	 *
	 * \param e The expression which is distributed.
	 */
	template<size_t O1, typename V, size_t O2, typename V1, typename E, typename G00, typename G01,
		typename = std::enable_if_t<!std::is_same<G01, OpDerivative<std::index_sequence<O2>, V1, E, SymbolicDerivative<G00>>>::value, int>>
	auto apply_operators(OpDerivative<std::index_sequence<O1>, V, OpDerivative<std::index_sequence<O2>, V1, E, SymbolicDerivative<G00>>, SymbolicDerivative<G01>> const& e);

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



	//! Implementation of the division rule.
	template<typename Dd, typename Sp, typename V, typename E1, typename E2,
		typename std::enable_if_t<(OpDerivative<Dd, V, OpBinaryMul<E1, E2>, Sp>::order > 0), int> = 0>
		auto apply_operators(OpDerivative<Dd, V, OpBinaryMul<E1, E2>, Sp> const& e);

	 //! Implementation of the product rule.
	template<typename Dd, typename Sp, typename V, typename E1, typename E2,
		typename std::enable_if_t<(OpDerivative<Dd, V, OpBinaryDiv<E1, E2>, Sp>::order > 0), int> = 0>
	auto apply_operators(OpDerivative<Dd, V, OpBinaryDiv<E1, E2>, Sp> const& e);

	////! Implementation of the product rule for symbolic derivatives.
	//template<size_t O, typename V, typename A1, typename A2, typename G>
	//auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpOperatorChain<A1, A2>, SymbolicDerivative<G>> const& e);

	////! Implementation of the product rule for symbolic derivatives.
	//template<size_t O, typename V, typename A1, typename A2, typename E, typename G>
	//auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpChain<A1, A2, E>, SymbolicDerivative<G>> const& e);

	////! Implementation of the product rule for symbolic derivatives.
	//template<size_t O, typename V, typename A1, typename A2, typename G>
	//auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpOperatorCombination<A1, A2>, SymbolicDerivative<G>> const& e);

	////! Implementation of the product rule for symbolic derivatives.
	//template<size_t O, typename V, typename A1, typename A2, typename E, typename G>
	//auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpCombination<A1, A2, E>, SymbolicDerivative<G>> const& e);

	//! Specialization based on expr::grid_dim.
	template<typename V, typename E, typename F>
	auto apply_operators(OpFunction<V, E, F, void> const& e);

	//! Specialization based on expr::grid_dim.
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto apply_operators(OpFunction<V, E, F, Arg0, Args...> const& e);

	//! Specialization based on expr::grid_dim.
	template<auto f, typename V, typename E>
	auto apply_operators(OpFunctionApply<f, V, E> const& e);

	template<typename V, typename V0, typename E, typename T, typename G>
	auto apply_operators(OpFunctionalDerivative<V, OpDomainIntegral<V0, E, T>, G> const& e);

	
	template<size_t O, typename V, size_t Z, int I0, int P0, typename V0, typename E, typename... Ts,
		int... I0s, int... P0s, typename A, typename B, typename... Vs>
	auto apply_operators(
		OpDerivative<std::index_sequence<O>, V, OpSum<V0, E,
			Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
            A, B, symphas::lib::types_list<Vs...>>,
			SymbolicDerivative<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>>> const& e);

	template<size_t O, typename V, typename V0, typename E, typename... Ts,
		int... I0s, int... P0s, typename A, typename B, typename... Vs, typename GG>
	auto apply_operators(
		OpDerivative<std::index_sequence<O>, V, OpSum<V0, E,
			Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
			A, B, symphas::lib::types_list<Vs...>>, SymbolicDerivative<GG>> const& e);

	
	template<typename V, typename E, typename... Ts,
		int... I0s, int... P0s, typename A, typename B, typename... Vs>
	auto apply_operators(OpSum<V, E,
		Substitution<SymbolicDataArray<Ts>...>,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
        A, B, symphas::lib::types_list<Vs...>> const& e);


	//! Implementation of the product rule for terms.
	//template<typename Dd, typename Sp, typename V, typename V0, typename G0, expr::exp_key_t X0, typename G1, expr::exp_key_t X1, typename... Gs, expr::exp_key_t... Xs,
	//	typename std::enable_if_t<(Dd::order > 0), int> = 0>
	//auto apply_operators(OpDerivative<Dd, V, OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>, Sp> const& e);
	
	template<typename Dd, typename Sp, typename V, typename V0, typename E,
		typename std::enable_if_t<(OpDerivative<Dd, V, OpExponential<V0, E>, Sp>::order > 0), int> = 0>
	auto apply_operators(OpDerivative<Dd, V, OpExponential<V0, E>, Sp> const& e);

	template<size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<expr::symbols::i_<I0, P0>, X0>>, SymbolicDerivative<GG>> const& e);

	template<size_t O, typename V, typename V1, int I0, int P0, typename G1, typename... Gs,
		expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t... Xs, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V1, Term<expr::symbols::i_<I0, P0>, X0>, Term<G1, X1>, Term<Gs, Xs>...>, SymbolicDerivative<GG>> const& e);

	template<size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0, size_t Z, typename GG,
		size_t N = expr::_Xk_t<X0>::N, typename = std::enable_if_t<(O > 0 && N >= O), int>>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>, SymbolicDerivative<Variable<Z, GG>>> const& e);

	template<size_t O, typename V, typename V0, int I0, int P0, size_t D, expr::exp_key_t X0, size_t Z, typename GG,
		size_t N = expr::_Xk_t<X0>::N, typename = std::enable_if_t<(O > 0 && N >= O), int>>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, D>, X0>>, SymbolicDerivative<Variable<Z, GG>>> const& e);

	template<size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0, typename GG, 
		size_t N = expr::_Xk_t<X0>::N, typename = std::enable_if_t<(O > 0 && N >= O), int>>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>, SymbolicDerivative<DynamicVariable<GG>>> const& e);

	template<size_t O, typename V, typename V0, int I0, int P0, size_t D, expr::exp_key_t X0, typename GG,
		size_t N = expr::_Xk_t<X0>::N, typename = std::enable_if_t<(O > 0 && N >= O), int>>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, D>, X0>>, SymbolicDerivative<DynamicVariable<GG>>> const& e);

	template<size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0, typename GG, 
		size_t N = expr::_Xk_t<X0>::N, typename = std::enable_if_t<(O > 0 && N >= O && !expr::is_expression<GG>), int>>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>, SymbolicDerivative<GG>> const& e);

	template<size_t O, typename V, typename V0, int I0, int P0, size_t D, expr::exp_key_t X0, typename GG,
		size_t N = expr::_Xk_t<X0>::N, typename = std::enable_if_t<(O > 0 && N >= O && !expr::is_expression<GG>), int>>
		auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, D>, X0>>, SymbolicDerivative<GG>> const& e);

	template<size_t O, typename V, typename V0, typename G0, expr::exp_key_t X0, typename GG, 
		size_t N = expr::factor_count_list<GG, Term<G0, X0>>::value,
		typename = std::enable_if_t<(O > 0 && N >= O && !expr::is_expression<GG>), int>>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>>, SymbolicDerivative<GG>> const& e);

	template<size_t O, typename V, typename V0, typename G0, typename G1, typename... Gs,
		expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t... Xs, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>, SymbolicDerivative<GG>> const& e);

	template<size_t O, typename V, typename V0, typename G0, expr::exp_key_t X0, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>>, SymbolicDerivative<OpTerm<OpIdentity, GG>>> const& e);

	template<size_t O, typename V, typename V0, typename G0, typename G1, typename... Gs,
		expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t... Xs, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>, SymbolicDerivative<OpTerm<OpIdentity, GG>>> const& e);

	template<size_t O, typename V0, auto f, typename V1, typename E, typename G0>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V0, OpFunctionApply<f, V1, E>, SymbolicDerivative<G0>> const& e);

	template<size_t O, typename V, typename... Es, typename G0>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpAdd<Es...>, SymbolicDerivative<G0>> const& e);

	template<size_t O, typename G, typename V, typename E, typename std::enable_if_t<!expr_has_deriv<E>::value, int> = 0>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G>> const& e)
	{
		return OpVoid{};
	}

	template<typename V, expr::exp_key_t X0, typename V0, typename E0, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<1>, V, OpPow<X0, V0, E0>, SymbolicDerivative<GG>> const& e);

	template<size_t O, typename V, expr::exp_key_t X0, typename V0, typename E0, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpPow<X0, V0, E0>, SymbolicDerivative<GG>> const& e);

	template<expr::exp_key_t X, typename V, typename E>
	auto apply_operators(OpPow<X, V, E> const& e);

	template<size_t O, typename V, typename E1, typename E2, typename G>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpBinaryDiv<E1, E2>, SymbolicDerivative<G>> const& e);

	//! Implementation of the product rule for symbolic derivatives.
	template<size_t O, typename V, typename E1, typename E2, typename G>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpBinaryMul<E1, E2>, SymbolicDerivative<G>> const& e);

	//! Implementation of the product rule for symbolic derivatives.
	template<size_t O, typename V, typename A1, typename E, typename G, std::enable_if_t<expr::is_coeff<A1>, int> = 0>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpChain<A1, OpIdentity, E>, SymbolicDerivative<G>> const& e);


	template<typename V, typename V1, typename E, typename Dd, typename Sp>
	auto apply_operators(OpDerivative<std::index_sequence<1>, V, OpDerivative<Dd, V1, E, Sp>, SymbolicDerivative<OpDerivative<Dd, OpIdentity, E, Sp>>> const& e)
	{
		auto expr = expr::get_enclosed_expression(e);
		return expr::coeff(e) * expr::coeff(expr);
	}

	template<typename V, typename E, typename Dd, typename Sp>
	auto apply_operators(OpDerivative<std::index_sequence<1>, V, OpDerivative<Dd, OpIdentity, E, Sp>, SymbolicDerivative<OpDerivative<Dd, OpIdentity, E, Sp>>> const& e)
	{
		auto expr = expr::get_enclosed_expression(e);
		return expr::coeff(e);
	}

	template<size_t O, typename V, typename V1, typename E, typename Dd, typename Sp, typename G0, 
		typename std::enable_if_t<!std::is_same<G0, OpDerivative<Dd, V1, E, Sp>>::value, int> = 0>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpDerivative<Dd, V1, E, Sp>, SymbolicDerivative<G0>> const& e)
	{
		return OpVoid{};
	}

	template<size_t O1, typename V, size_t O2, typename V1, typename E, typename G00>
	auto apply_operators(OpDerivative<std::index_sequence<O1>, V, OpDerivative<std::index_sequence<O2>, V1, E, SymbolicDerivative<G00>>, SymbolicDerivative<G00>> const& e)
	{
		auto expr = expr::get_enclosed_expression(e);
		return expr::coeff(e) * expr::coeff(expr) * apply_operators(expr::make_derivative<O1 + O2, G00>(expr::get_enclosed_expression(expr), e.solver));
	}

	template<size_t O1, typename V, size_t O2, typename V1, typename E, typename G00, typename G01, typename>
	auto apply_operators(OpDerivative<std::index_sequence<O1>, V, OpDerivative<std::index_sequence<O2>, V1, E, SymbolicDerivative<G00>>, SymbolicDerivative<G01>> const& e)
	{
		auto expr = expr::get_enclosed_expression(e);
		return expr::coeff(e) * apply_operators(expr::make_derivative<O1, G01>(apply_operators(expr), e.solver));
	}

	template<typename G, typename V, typename E>
	auto apply_operators(OpDerivative<std::index_sequence<0>, V, E, SymbolicDerivative<G>> const& e)
	{
		return expr::coeff(e) * expr::get_enclosed_expression(e);
	}


	template<typename V, typename E>
	auto apply_operators(OpDerivative<std::index_sequence<1>, V, E, SymbolicDerivative<E>> const& e)
	{
		return expr::coeff(e);
	}

	namespace
	{

		template<typename Dd, typename Sp, typename E1, typename E2>
		auto handle_apply_mul(Sp const& solver, OpExpression<E1> const& lhs0, OpExpression<E2> const& rhs0)
		{
			auto lhs = apply_operators(expr::make_derivative<Dd>(*static_cast<E1 const*>(&lhs0), solver) * (*static_cast<E2 const*>(&rhs0)));
			auto rhs = apply_operators((*static_cast<E1 const*>(&lhs0)) * expr::make_derivative<Dd>(*static_cast<E2 const*>(&rhs0), solver));
			return lhs + rhs;
		}

		template<typename Dd, typename Sp, typename E1, typename E2>
		auto handle_apply_mul(Sp const& solver, OpOperator<E1> const& lhs0, OpExpression<E2> const& rhs0)
		{
			return apply_operators(expr::make_derivative<Dd>(apply_operators((*static_cast<E1 const*>(&lhs0)) * (*static_cast<E2 const*>(&rhs0))), solver));
		}

		template<typename Dd, typename Sp, typename E1, typename E2>
		auto handle_apply_mul(Sp const& solver, OpExpression<E1> const& lhs0, OpOperator<E2> const& rhs0)
		{
			auto lhs = apply_operators(expr::make_derivative<Dd>(*static_cast<E1 const*>(&lhs0), solver) * (*static_cast<E2 const*>(&rhs0)));
			auto rhs = apply_operators((*static_cast<E1 const*>(&lhs0)) * expr::make_derivative<Dd>(*static_cast<E2 const*>(&rhs0), solver));
			return lhs + rhs;
		}

		template<typename Dd, typename Sp, typename E1, typename E2>
		auto handle_apply_mul(Sp const& solver, OpOperator<E1> const& lhs0, OpOperator<E2> const& rhs0)
		{
			return apply_operators(expr::make_derivative<Dd>(apply_operators((*static_cast<E1 const*>(&lhs0)) * (*static_cast<E2 const*>(&rhs0))), solver));
		}

		template<typename Dd, typename Sp, typename E10, typename E11, typename E2>
		auto handle_apply_mul(Sp const& solver, OpOperatorCombination<E10, E11> const& lhs0, OpExpression<E2> const& rhs0)
		{
			return apply_operators(expr::make_derivative<Dd>(lhs0.f * (*static_cast<E2 const*>(&rhs0)) + lhs0.g * (*static_cast<E2 const*>(&rhs0)), solver));
		}

		template<typename Dd, typename Sp, typename E10, typename E11, typename E2>
		auto handle_apply_mul(Sp const& solver, OpOperatorChain<E10, E11> const& lhs0, OpExpression<E2> const& rhs0)
		{
			return apply_operators(expr::make_derivative<Dd>(lhs0.f(lhs0.g * (*static_cast<E2 const*>(&rhs0))), solver));
		}

		template<typename Dd, typename Sp, typename E10, typename E11, typename E2>
		auto handle_apply_mul(Sp const& solver, OpOperatorCombination<E10, E11> const& lhs0, OpOperator<E2> const& rhs0)
		{
			return apply_operators(expr::make_derivative<Dd>(lhs0.f * (*static_cast<E2 const*>(&rhs0)) + lhs0.g * (*static_cast<E2 const*>(&rhs0)), solver));
		}

		template<typename Dd, typename Sp, typename E10, typename E11, typename E2>
		auto handle_apply_mul(Sp const& solver, OpOperatorChain<E10, E11> const& lhs0, OpOperator<E2> const& rhs0)
		{
			return apply_operators(expr::make_derivative<Dd>(lhs0.f(lhs0.g * (*static_cast<E2 const*>(&rhs0))), solver));
		}
	}

	//! Implementation of the quotient rule for symbolic derivatives.
	template<size_t O, typename V, typename E1, typename E2, typename G>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpBinaryDiv<E1, E2>, SymbolicDerivative<G>> const& e)
	{
		//auto expr = expr::get_enclosed_expression(e);
		//return expr::coeff(e) * handle_apply_mul(expr::make_operator_derivative<O, G>(e.solver), apply_operators(expr.a), apply_operators(expr.b));
		
		auto expr = expr::get_enclosed_expression(e);
		auto a = apply_operators(expr.a);
		auto b = apply_operators(expr.b);

		auto lhs = apply_operators(expr::make_derivative<O, G>(a, e.solver) * b);
		auto rhs = apply_operators(a * expr::make_derivative<O, G>(b, e.solver));
		return expr::coeff(e) * (lhs - rhs) / (b * b);
	}

	//! Implementation of the quotient rule for symbolic derivatives.
	template<size_t O, typename V, typename E1, typename E2, typename G>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpBinaryMul<E1, E2>, SymbolicDerivative<G>> const& e)
	{
		auto expr = expr::get_enclosed_expression(e);
		return expr::coeff(e) * handle_apply_mul<std::index_sequence<O>>(e.solver, apply_operators(expr.a), apply_operators(expr.b));

	}


	//! Implementation of the product rule for symbolic derivatives.
	template<size_t O, typename V, typename A1, typename E, typename G, std::enable_if_t<expr::is_coeff<A1>, int>>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpChain<A1, OpIdentity, E>, SymbolicDerivative<G>> const& e)
	{
		auto expr = expr::get_enclosed_expression(e);
		auto applied = apply_operators(expr::make_derivative<O, G>(expr::get_enclosed_expression(expr), e.solver));
		return OpOperatorChain(expr::coeff(e) * expr.combination.f, OpIdentity{})(applied);
	}

	//! Implementation of the product rule.
	template<typename Dd, typename Sp, typename V, typename E1, typename E2,
		typename std::enable_if_t<(OpDerivative<Dd, V, OpBinaryMul<E1, E2>, Sp>::order > 0), int>>
	auto apply_operators(OpDerivative<Dd, V, OpBinaryMul<E1, E2>, Sp> const& e)
	{

		auto expr = expr::get_enclosed_expression(e);
		return expr::coeff(e) * handle_apply_mul<Dd>(e.solver, apply_operators(expr.a), apply_operators(expr.b));

		//auto expr = expr::get_enclosed_expression(e);
		//auto a = apply_operators(expr.a);
		//auto b = apply_operators(expr.b);
		//
		//auto lhs = apply_operators(expr::make_derivative<Dd>(a, e.solver) * b);
		//auto rhs = apply_operators(a * expr::make_derivative<Dd>(b, e.solver));
		//return expr::coeff(e) * (lhs + rhs);
	}

	//! Implementation of the quotient rule.
	template<typename Dd, typename Sp, typename V, typename E1, typename E2, 
		typename std::enable_if_t<(OpDerivative<Dd, V, OpBinaryDiv<E1, E2>, Sp>::order > 0), int>>
	auto apply_operators(OpDerivative<Dd, V, OpBinaryDiv<E1, E2>, Sp> const& e)
	{
		auto expr = expr::get_enclosed_expression(e);
		auto a = apply_operators(expr.a);
		auto b = apply_operators(expr.b);

		auto lhs = apply_operators(expr::make_derivative<Dd>(a, e.solver) * b);
		auto rhs = apply_operators(a * expr::make_derivative<Dd>(b, e.solver));
		return expr::coeff(e) * (lhs - rhs) / (b * b);
	}

	////! Implementation of the product rule for symbolic derivatives.
	//template<size_t O, typename V, typename A1, typename A2, typename G>
	//auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpOperatorChain<A1, A2>, SymbolicDerivative<G>> const& e)
	//{
	//	auto&& expr = expr::get_enclosed_expression(e);
	//	auto lhs = apply_operators(expr::make_derivative<O, G>(expr.f, e.solver))(expr.g);
	//	auto rhs = expr.f(apply_operators(expr::make_derivative<O, G>(expr.g, e.solver)));
	//	return expr::coeff(e) * (lhs + rhs);
	//}

	////! Implementation of the product rule for symbolic derivatives.
	//template<size_t O, typename V, typename A1, typename A2, typename E, typename G>
	//auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpChain<A1, A2, E>, SymbolicDerivative<G>> const& e)
	//{
	//	auto&& expr = apply_operators(expr::get_enclosed_expression(e));
	//	auto lhs = apply_operators(expr::make_derivative<O, G>(expr.combination, e.solver)) * expr::get_enclosed_expression(expr);
	//	auto rhs = expr.combination * apply_operators(expr::make_derivative<O, G>(expr::get_enclosed_expression(expr), e.solver));
	//	return expr::coeff(e) * (lhs + rhs);
	//}

	////! Implementation of the product rule for symbolic derivatives.
	//template<size_t O, typename V, typename A1, typename A2, typename G>
	//auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpOperatorCombination<A1, A2>, SymbolicDerivative<G>> const& e)
	//{
	//	auto&& expr = expr::get_enclosed_expression(e);
	//	return expr::coeff(e) * (apply_operators(expr::make_derivative<O, G>(expr.f, e.solver)) + apply_operators(expr::make_derivative<O, G>(expr.g, e.solver)));
	//}

	////! Implementation of the product rule for symbolic derivatives.
	//template<size_t O, typename V, typename A1, typename A2, typename E, typename G>
	//auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpCombination<A1, A2, E>, SymbolicDerivative<G>> const& e)
	//{
	//	auto&& expr = apply_operators(expr::get_enclosed_expression(e));
	//	auto lhs = apply_operators(expr::make_derivative<O, G>(expr.combination, e.solver)) * expr::get_enclosed_expression(expr);
	//	auto rhs = expr.combination * apply_operators(expr::make_derivative<O, G>(expr::get_enclosed_expression(expr), e.solver));
	//	return expr::coeff(e) * (lhs + rhs);
	//}

	//! Specialization based on expr::grid_dim.
	template<typename V, typename E, typename F>
	auto apply_operators(OpFunction<V, E, F, void> const& e)
	{
		return expr::coeff(e) * expr::make_function(apply_operators(expr::get_enclosed_expression(e)), e.name, e.f);
	};

	//! Specialization based on expr::grid_dim.
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto apply_operators(OpFunction<V, E, F, Arg0, Args...> const& e)
	{
		return expr::coeff(e) * expr::make_function(apply_operators(expr::get_enclosed_expression(e)), e.name, e.f, e.args);
	};

	//! Specialization based on expr::grid_dim.
	template<auto f, typename V, typename E>
	auto apply_operators(OpFunctionApply<f, V, E> const& e)
	{
		return expr::coeff(e) * expr::make_function<f>(apply_operators(expr::get_enclosed_expression(e)));
	};

	//! Specialization based on expr::grid_dim.
	template<typename V, typename E, typename T>
	auto apply_operators(OpIntegral<V, E, T> const& e)
	{
		return expr::coeff(e) * expr::make_integral(apply_operators(expr::get_enclosed_expression(e)), e.domain);
	};

	//! Implementation of the product rule for terms.
	//template<typename Dd, typename Sp, typename V, typename V0, typename G0, expr::exp_key_t X0, typename G1, expr::exp_key_t X1, typename... Gs, expr::exp_key_t... Xs,
	//	typename std::enable_if_t<(OpDerivative<Dd, V, OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>, Sp>::order > 0), int>>
	//auto apply_operators(OpDerivative<Dd, V, OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>, Sp> const& e)
	//{
	//	auto&& terms = expr::get_enclosed_expression(e);
	//	auto coeff = expr::coeff(terms);
	//	auto a = OpTerms(OpIdentity{}, expr::get<1>(terms));
	//	auto b = OpTerms(OpIdentity{}, *static_cast<OpTerms<Term<G1, X1>, Term<Gs, Xs>...> const*>(&terms));
	//
	//	auto lhs = apply_operators(expr::make_derivative<Dd>(a, e.solver)) * b;
	//	auto rhs = a * apply_operators(expr::make_derivative<Dd>(b, e.solver));
	//	return (coeff * e.value) * (lhs + rhs);
	//}

	/*
	//! Implementation of the product rule for terms.
	template<typename Dd, typename Sp, typename V, typename V0, typename G0, expr::exp_key_t X0,
		typename std::enable_if_t<(OpDerivative<Dd, V, OpTerms<V0, Term<G0, X0>>, Sp>::order > 0 
			&& (expr::Xk_t<X0>::N > expr::Xk_t<X0>::D && !expr::Xk_t<X0>::sign)), int>>
	auto apply_operators(OpDerivative<Dd, V, OpTerms<V0, Term<G0, X0>>, Sp> const& e)
	{
		auto&& term = expr::get_enclosed_expression(e);

		constexpr size_t N = expr::_Xk_t<X0>::N;
		constexpr size_t D = expr::_Xk_t<X0>::D;
		constexpr bool sign = expr::_Xk_t<X0>::sign;

		auto dvar = symphas::internal::to_term_element<expr::Xk<N - D, D, sign>>(expr::get<1>(term).data());
		auto var1 = symphas::internal::to_term_element(expr::get<1>(term).data());

		auto coeff = expr::coeff(term);
		auto dterm = OpTerms(coeff, dvar);

		auto lhs = apply_operators(expr::make_derivative<Dd>(dvar, e.solver)) * var1;
		auto rhs = dvar * apply_operators(expr::make_derivative<Dd>(var1, e.solver));
		return (lhs + rhs);
	}*/

	template<typename Dd, typename Sp, typename V, typename V0, typename E,
		typename std::enable_if_t<(OpDerivative<Dd, V, OpExponential<V0, E>, Sp>::order > 0), int>>
	auto apply_operators(OpDerivative<Dd, V, OpExponential<V0, E>, Sp> const& e)
	{
		auto&& expr = expr::get_enclosed_expression(e);
		auto&& pow = expr::get_enclosed_expression(expr);
		return expr::coeff(e) * expr * apply_operators(expr::make_derivative<Dd>(pow, e.solver));
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

		template<typename V, typename E>
		auto apply_operators_chain(OpChain<V, OpIdentity, E> const& e)
		{
			return OpChain(OpOperatorChain(apply_operators(e.combination.f), OpIdentity{}), 
				apply_operators(expr::get_enclosed_expression(e)));
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
			return (apply_operators(expr::make_derivative<Dd>(expr::get<Is>(e), solver)) + ...);
		}

		template<size_t O, typename G0, typename... Es, size_t... Is>
		auto apply_operators_adds(SymbolicDerivative<G0> const& solver, OpAdd<Es...> const& e, std::index_sequence<Is...>)
		{
			return (apply_operators(expr::make_derivative<O, G0>(expr::get<Is>(e), solver)) + ...);
		}


		template<typename A1, typename A2>
		auto apply_operators_chain(OpOperatorChain<A1, A2> const& e)
		{
			return e;
		}

		template<typename E>
		auto apply_operators_chain(OpExpression<E> const& e)
		{
			return apply_operators(*static_cast<E const*>(&e));
		}

		template<typename E>
		auto apply_operators_chain(OpOperator<E> const& e)
		{
			return apply_operators(*static_cast<E const*>(&e));
		}
	}


	template<typename A1, typename A2>
	auto apply_operators(OpOperatorChain<A1, A2> const& e)
	{
		return apply_operators_chain(apply_operators(e.f)(apply_operators(e.g)));
	}

	template<typename A2>
	auto apply_operators(OpOperatorChain<OpIdentity, A2> const& e)
	{
		return OpOperatorChain(OpIdentity{}, apply_operators(e.g));
	}

	//template<typename A1, typename B1, typename B2>
	//auto apply_operators(OpOperatorChain<A1, OpOperatorChain<B1, B2>> const& e)
	//{
	//	return apply_operators(((apply_operators(e.f)(apply_operators(e.g.f))) * apply_operators(e.g.g)))
	//		+ apply_operators((apply_operators(e.g.f) * (apply_operators(e.f)(apply_operators(e.g.g)))));
	//}

	//template<typename A1, typename B1, typename B2>
	//auto apply_operators(OpOperatorChain<A1, OpBinaryMul<B1, B2>> const& e)
	//{
	//	return apply_operators(((apply_operators(e.f)(apply_operators(e.g.a))) * apply_operators(e.g.b)))
	//		+ apply_operators((apply_operators(e.g.a) * (apply_operators(e.f)(apply_operators(e.g.b)))));
	//}

	template<typename A1, typename A2>
	auto apply_operators(OpOperatorCombination<A1, A2> const& e)
	{
		return apply_operators(e.f) + apply_operators(e.g);
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


		template<typename Sp, typename Seq, typename Axes>
		struct mixed_derivative_type {};

		template<typename Sp, size_t... Os, Axis... axs>
		struct mixed_derivative_type<Sp, std::index_sequence<Os...>, symphas::lib::axis_list<axs...>>
		{
			template<Axis ax0, Axis ax1, size_t O>
			static const size_t p1 = (ax0 == ax1) ? (O % 2) : 0;
			
			template<Axis ax0, Axis ax1, size_t O>
			static const size_t pO = (ax0 == ax1) ? (O - (O % 2)) : 0;


			template<size_t O, Axis axd>
			//using type = typename Solver<Sp>::template mixed_derivative<(Os + pO<axm, axs, O> + p1<axs, O>)...>;
			using type = typename Solver<Sp>::template mixed_derivative<(Os + pO<axd, axs, O> + p1<axd, axs, O>)...>;
		};

		template<typename Sp>
		struct combine_mixed_derivatives
		{
			using solver_type = Solver<Sp>;
			
			template<size_t... Os>
			using mixed_derivative = typename solver_type::template mixed_derivative<Os...>;
			template<Axis ax, size_t O>
			using directional_derivative = typename solver_type::template directional_derivative<ax, O>;
			template<Axis ax, size_t O>
			using derivative = typename solver_type::template derivative<ax, O>;

			template<typename E, Axis ax1, size_t O1, Axis ax2, size_t O2, Axis... axs>
			auto operator()(
				OpExpression<E> const& enclosed,
				Sp const& solver,
				directional_derivative<ax1, O1>,
				directional_derivative<ax2, O2>,
				symphas::lib::axis_list<axs...>) const
			{
				if constexpr (ax1 == ax2)
				{
					using Dd = typename Solver<Sp>::template directional_derivative<ax1, O1 + O2>;
					return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver);
				}
				else
				{
					using Dd = typename Solver<Sp>::template mixed_derivative<((ax1 == axs) ? O1 : (ax2 == axs) ? O2 : 0)...>;
					return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver);
				}
			}

			template<typename E, Axis ax1, Axis ax2, Axis... axs>
			auto operator()(
				OpExpression<E> const& enclosed,
				Sp const& solver,
				derivative<ax1, 1>,
				derivative<ax2, 1>,
				symphas::lib::axis_list<axs...>) const
			{
				if constexpr (ax1 == ax2)
				{
					using Dd = typename Solver<Sp>::template directional_derivative<ax1, 2>;
					return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver);
				}
				else
				{
					using Dd = typename Solver<Sp>::template mixed_derivative<((ax1 == axs || ax2 == axs) ? 1 : 0)...>;
					return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver);
				}
			}

			template<size_t... O1s>
			struct deduction_redirection
			{
				template<typename E, size_t... O2s, Axis... axs>
				auto operator()(
					OpExpression<E> const& enclosed,
					Sp const& solver,
					mixed_derivative<O2s...>,
					symphas::lib::axis_list<axs...>) const
				{
					using Dd = typename Solver<Sp>::template mixed_derivative<(O1s + O2s)...>;
					return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver);
				}
			};

			template<typename E, size_t... O1s, typename mixed_other, Axis... axs>
			auto operator()(
				OpExpression<E> const& enclosed,
				Sp const& solver,
				mixed_derivative<O1s...>,
				mixed_other,
				symphas::lib::axis_list<axs...>) const
			{
				return deduction_redirection<O1s...>{}(*static_cast<E const*>(&enclosed), solver, mixed_other{}, symphas::lib::axis_list<axs...>{});
			}

			template<typename E, Axis ax, size_t O, size_t... O2s, Axis... axs>
			auto operator()(
				OpExpression<E> const& enclosed,
				Sp const& solver,
				directional_derivative<ax, O>,
				mixed_derivative<O2s...>,
				symphas::lib::axis_list<axs...>) const
			{
				using Dd = typename Solver<Sp>::template mixed_derivative<((ax == axs) ? O + O2s : O2s)...>;
				return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver);
			}

			template<typename E, Axis ax, size_t O, size_t... O1s, Axis... axs>
			auto operator()(
				OpExpression<E> const& enclosed,
				Sp const& solver,
				mixed_derivative<O1s...>,
				directional_derivative<ax, O>,
				symphas::lib::axis_list<axs...>) const
			{
				using Dd = typename Solver<Sp>::template mixed_derivative<((ax == axs) ? O + O1s : O1s)...>;
				return expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver);
			}

			template<typename E, Axis ax, size_t O, size_t... O2s, Axis... axs>
			auto operator()(
				OpExpression<E> const& enclosed,
				Sp const& solver,
				derivative<ax, O>,
				mixed_derivative<O2s...>,
				symphas::lib::axis_list<axs...>) const
			{
				return (expr::make_derivative<
					typename mixed_derivative_type<Sp, std::index_sequence<O2s...>, symphas::lib::axis_list<axs...>>::template type<O, ax>
				>(*static_cast<E const*>(&enclosed), solver) );
			}

			template<typename E, Axis ax, size_t O, size_t... O1s, Axis... axs>
			auto operator()(
				OpExpression<E> const& enclosed,
				Sp const& solver,
				mixed_derivative<O1s...>,
				derivative<ax, O>,
				symphas::lib::axis_list<axs...>) const
			{
				return (expr::make_derivative<
					typename mixed_derivative_type<Sp, std::index_sequence<O1s...>, symphas::lib::axis_list<axs...>>::template type<O, ax>
				>(*static_cast<E const*>(&enclosed), solver) );
			}

			template<typename E, Axis ax1, size_t O1, Axis ax2, size_t O2, Axis... axs>
			auto operator()(
				OpExpression<E> const& enclosed,
				Sp const& solver,
				derivative<ax1, O1>,
				directional_derivative<ax2, O2>,
				symphas::lib::axis_list<axs...>) const
			{
				//using mixed_type = typename Solver<Sp>::template mixed_derivative<((axs == ax1) ? O1 : (axs == ax2) ? O2 : 0)...>;
				//return expr::make_derivative<mixed_type>(*static_cast<E const*>(&enclosed), solver);
				return (expr::make_derivative<
					typename mixed_derivative_type<Sp, std::index_sequence<((ax2 == axs) ? O2 : 0)...>, symphas::lib::axis_list<axs...>>::template type<O1, ax1>
				>(*static_cast<E const*>(&enclosed), solver) );
			}

			template<typename E, Axis ax1, size_t O1, Axis ax2, size_t O2, Axis... axs>
			auto operator()(
				OpExpression<E> const& enclosed,
				Sp const& solver,
				directional_derivative<ax1, O1>,
				derivative<ax2, O2>,
				symphas::lib::axis_list<axs...>) const
			{
				//using mixed_type = typename Solver<Sp>::template mixed_derivative<((axs == ax1) ? O1 : (axs == ax2) ? O2 : 0)...>;
				//return expr::make_derivative<mixed_type>(*static_cast<E const*>(&enclosed), solver);
				return (expr::make_derivative<
					typename mixed_derivative_type<Sp, std::index_sequence<((ax1 == axs) ? O1 : 0)...>, symphas::lib::axis_list<axs...>>::template type<O2, ax2>
				>(*static_cast<E const*>(&enclosed), solver) );
			}

		};




		template<typename Dd1, typename Dd2, typename V1, typename V2, typename E, typename Sp1, typename Sp2>
		auto apply_operator_derivative_nested(OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp1>, Sp2> const& e)
		{
			return e;
		}

		template<typename Dd1, typename V1, typename E, typename Sp2>
		auto apply_operator_derivative_nested(OpDerivative<Dd1, V1, E, Sp2> const& e)
		{
			return expr::apply_operators(e);
		}
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
	//template<typename Dd1, size_t N, typename V1, typename V2, typename E, typename Sp1, typename Sp2>
	//auto apply_operators(OpDerivative<Dd1, V1, OpDerivative<std::index_sequence<N>, V2, E, SymbolicDerivative<Sp1>>, Sp2> const& e);


	//template<typename Dd1, size_t N, typename V1, typename V2, typename E, typename Sp1, typename Sp2>
	//auto apply_operators(OpDerivative<Dd1, V1, OpDerivative<std::index_sequence<N>, V2, E, SymbolicDerivative<Sp1>>, Sp2> const& e)
	//{
	//	return e;
	//}




	template<typename Dd1, typename Dd2, typename V1, typename V2, typename E, typename Sp1, typename Sp2>
	auto apply_operators(OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp1>, Sp2> const& e)
	{
		return apply_operator_derivative_nested(expr::coeff(e) * expr::make_derivative<Dd1>(expr::apply_operators(expr::get_enclosed_expression(e)), e.solver));
	}

	template<bool parity1, bool parity2, size_t order1, size_t order2, Axis ax1, Axis ax2>
	struct combine_derivatives
	{
		template<typename Sp, typename E, Axis... axs>
		auto operator()(Sp const& solver, OpExpression<E> const& enclosed, symphas::lib::axis_list<axs...>)
		{
			auto dd1 = expr::make_operator_directional_derivative<ax1, parity1>(solver)
				* break_up_derivative<order1 - parity1>(solver, symphas::lib::axis_list<axs...>{});
			auto dd2 = expr::make_operator_directional_derivative<ax2, parity2>(solver)
				* break_up_derivative<order2 - parity2>(solver, symphas::lib::axis_list<axs...>{});
			
			return apply_operators((dd1 * dd2) * *static_cast<E const*>(&enclosed));
		}
	};

	template<size_t order1, size_t order2, Axis ax1, Axis ax2>
	struct combine_derivatives<0, 0, order1, order2, ax1, ax2>
	{
		template<typename Sp, typename E, Axis... axs>
		auto operator()(Sp const& solver, OpExpression<E> const& enclosed, symphas::lib::axis_list<axs...>)
		{
			using Dd = typename Solver<Sp>::template derivative<Axis::X, order1 + order2>;
			return apply_operators(expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver));
			//return apply_operators(expr::make_operator_derivative<order1 + order2>(solver)(OpVoid{}));
		}
	};

	template<size_t order1, size_t order2, Axis ax1, Axis ax2>
	struct combine_derivatives<1, 0, order1, order2, ax1, ax2>
	{
		template<typename Sp, typename E, Axis... axs>
		auto operator()(Sp const& solver, OpExpression<E> const& enclosed, symphas::lib::axis_list<axs...>)
		{
			using Dd = typename Solver<Sp>::template derivative<ax1, order1 + order2>;
			return apply_operators(expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver));
		}
	};

	template<size_t order1, size_t order2, Axis ax1, Axis ax2>
	struct combine_derivatives<0, 1, order1, order2, ax1, ax2>
	{
		template<typename Sp, typename E, Axis... axs>
		auto operator()(Sp const& solver, OpExpression<E> const& enclosed, symphas::lib::axis_list<axs...>)
		{
			using Dd = typename Solver<Sp>::template derivative<ax2, order1 + order2>;
			return apply_operators(expr::make_derivative<Dd>(*static_cast<E const*>(&enclosed), solver));
		}
	};

	template<typename Dd1, typename Dd2, typename V1, typename V2, typename E, typename Sp>
	auto apply_operators(OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp>, Sp> const& e)
	{
		using d1t = OpDerivative<Dd1, V1, OpDerivative<Dd2, V2, E, Sp>, Sp>;
		using d2t = OpDerivative<Dd2, V2, E, Sp>;

		constexpr size_t D = expr::grid_dim<d1t>::value;

		constexpr size_t order1 = d1t::order;
		constexpr size_t order2 = d2t::order;

		auto enclosed1 = expr::get_enclosed_expression(e);
		auto enclosed2 = expr::get_enclosed_expression(enclosed1);
		auto enclosed = expr::coeff(e) * expr::coeff(enclosed1) * apply_operators(enclosed2);

		auto axis_list = symphas::lib::make_axis_list<D>();

		if constexpr (!Dd1::is_mixed && !Dd2::is_mixed && !Dd1::is_directional && !Dd2::is_directional)
		{
			return combine_derivatives<order1 % 2, order2 % 2, order1, order2, d1t::axis, d2t::axis>{}(e.solver, enclosed, axis_list);
			//if constexpr (order1 % 2 == 0 && order2 % 2 == 0)
			//{
			//	//using Dd = typename Solver<Sp>::template derivative<Axis::X, order1 + order2>;
			//	//return apply_operators(expr::make_derivative<Dd>(enclosed, e.solver));
			//
			//	return combine_apply(e.solver, enclosed, axis_list);
			//}
			//else if constexpr (order1 % 2 == 1)
			//{
			//	//using Dd = typename Solver<Sp>::template derivative<d1t::axis, order1 + order2>;
			//	//return apply_operators(expr::make_derivative<Dd>(enclosed, e.solver));
			//}
			//else if constexpr (order2 % 2 == 1)
			//{
			//	//using Dd = typename Solver<Sp>::template derivative<d2t::axis, order1 + order2>;
			//	//return apply_operators(expr::make_derivative<Dd>(enclosed, e.solver));
			//}
			//else
			//{
			//	return combine_apply(e.solver, enclosed, axis_list);
			//	//auto dd1 = expr::make_operator_directional_derivative<Dd1::axis, order1 % 2>(e.solver)
			//	//	* break_up_derivative<order1 - (order1 % 2)>(e.solver, symphas::lib::make_axis_list<D>());
			//	//
			//	//auto dd2 = expr::make_operator_directional_derivative<Dd2::axis, order2 % 2>(e.solver)
			//	//	* break_up_derivative<order2 - (order2 % 2)>(e.solver, symphas::lib::make_axis_list<D>());
			//	//return apply_operators((dd1 * dd2) * enclosed);
			//}
		}
		else
		{
			return apply_operators(combine_mixed_derivatives<Sp>{}(enclosed, e.solver, Dd1{}, Dd2{}, axis_list));
		}
	}

	namespace
	{
		template<size_t O, typename V, typename E, typename GG, typename E0>
		auto apply_operators_deriv(
			OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& d,
			OpExpression<E0> const& e)
		{
			return apply_operators(expr::make_derivative<O, GG>(*static_cast<E0 const*>(&e), d.solver));
		}

		template<size_t O, typename V, typename E, typename GG, typename A1, typename A2, typename E0>
		auto apply_operators_deriv(
			OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& d,
			OpCombination<A1, A2, E0> const& e)
		{
			return apply_operators(expr::make_derivative<O, GG>(e.combination.f(e.e), d.solver))
				+ apply_operators(expr::make_derivative<O, GG>(e.combination.g(e.e), d.solver));
		}

		template<size_t O, typename V, typename E, typename GG, typename E0>
		auto apply_operators_deriv(
			OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& d,
			OpOperatorChain<OpIdentity, E0> const& e)
		{
			return apply_operators(expr::make_derivative<O, GG>(e.g));
		}
		
		template<size_t O, typename V, typename E, typename GG, typename E0>
		auto apply_operators_deriv(
			OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& d,
			OpOperator<E0> const& e)
		{
			return OpVoid{};
		}

		//template<typename Dd, typename V, typename E, typename Sp, typename E0>
		//auto apply_operators_deriv(OpDerivative<Dd, V, E, Sp> const& d, OpOperator<E0> const& e)
		//{
		//	return apply_operators(expr::make_derivative<Dd>(*static_cast<E0 const*>(&e), d.solver));
		//}

		template<typename Dd, typename V, typename E, typename Sp, typename E0>
		auto apply_operators_deriv(OpDerivative<Dd, V, E, Sp> const& d, OpOperator<E0> const& e);

		template<typename Dd, typename V, typename E, typename Sp, typename E0,
			typename std::enable_if_t<expr_has_deriv<E>::value, int> = 0>
		auto apply_operators_deriv(OpDerivative<Dd, V, E, Sp> const& d, OpExpression<E0> const& e)
		{
			return apply_operators(expr::make_derivative<Dd>(*static_cast<E0 const*>(&e), d.solver));
		}

		template<typename Dd, typename V, typename E, typename Sp, typename E0, 
			typename std::enable_if_t<!expr_has_deriv<E>::value, int> = 0>
		auto apply_operators_deriv(OpDerivative<Dd, V, E, Sp> const& d, OpExpression<E0> const& e)
		{
			return expr::make_derivative<Dd>(*static_cast<E0 const*>(&e), d.solver);
		}
	}

	template<typename Dd, typename V, typename E, typename Sp, typename>
	auto apply_operators(OpDerivative<Dd, V, E, Sp> const& e)
	{
		return expr::coeff(e) * apply_operators_deriv(e, apply_operators(expr::get_enclosed_expression(e)));
	}

	template<typename Dd, typename V, typename... Es, typename Sp>
	auto apply_operators(OpDerivative<Dd, V, OpAdd<Es...>, Sp> const& e)
	{
		return expr::coeff(e) * apply_operators_adds<Dd>(e.solver,
			expr::get_enclosed_expression(e), 
			std::make_index_sequence<sizeof...(Es)>{});
	}

	template<size_t O, typename V, typename... Es, typename G0>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpAdd<Es...>, SymbolicDerivative<G0>> const& e)
	{
		return expr::coeff(e) * apply_operators_adds<O>(e.solver,
			expr::get_enclosed_expression(e), 
			std::make_index_sequence<sizeof...(Es)>{});
	}

	template<typename... Es>
	auto apply_operators(OpAdd<Es...> const& e)
	{
		return apply_operators_adds(e, std::make_index_sequence<sizeof...(Es)>{});
	}

	namespace
	{

		template<typename A, typename B, typename E2>
		auto apply_operators_mul(OpOperatorCombination<A, B> const& combination, OpExpression<E2> const& b)
		{
			return apply_operators(combination.f * (*static_cast<E2 const*>(&b))) + apply_operators(combination.g * (*static_cast<E2 const*>(&b)));
		}

		//! Apply the chain operation to an expression.
		template<typename A1, typename A2, typename E>
		auto apply_operators_mul(OpOperatorChain<A1, A2> const& combination, OpExpression<E> const& b)
		{
			return apply_operators(combination.f(apply_operators(combination.g * *static_cast<E const*>(&b))));
		}

		template<typename A, typename B, typename E2>
		auto apply_operators_mul(OpOperatorCombination<A, B> const& combination, OpOperator<E2> const& b)
		{
			return apply_operators(combination.f * (*static_cast<E2 const*>(&b))) + apply_operators(combination.g * (*static_cast<E2 const*>(&b)));
		}

		//! Apply the chain operation to an expression.
		template<typename A1, typename A2, typename E>
		auto apply_operators_mul(OpOperatorChain<A1, A2> const& combination, OpOperator<E> const& b)
		{
			return apply_operators(combination.f(apply_operators(combination.g * *static_cast<E const*>(&b))));
		}

		template<typename E1, typename E2>
		auto apply_operators_mul(OpOperator<E1> const& a, OpExpression<E2> const& b)
		{
			return apply_operators((*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b)));
		}

		template<typename E1, typename E2>
		auto apply_operators_mul(OpExpression<E1> const& a, OpOperator<E2> const& b)
		{
			return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
		}

		template<typename E1, typename E2>
		auto apply_operators_mul(OpOperator<E1> const& a, OpOperator<E2> const& b)
		{
			return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
		}

		template<typename E1, typename E2>
		auto apply_operators_mul(OpExpression<E1> const& a, OpExpression<E2> const& b)
		{
			return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
		}

		template<typename... Es, typename E2>
		auto apply_operators_mul(OpAdd<Es...> const& a, OpExpression<E2> const& b)
		{
			return apply_operators(a * *static_cast<E2 const*>(&b));
		}

		template<typename... Es, typename E2>
		auto apply_operators_mul(OpAdd<Es...> const& a, OpOperator<E2> const& b)
		{
			return apply_operators(a * *static_cast<E2 const*>(&b));
		}

		template<typename... Es, typename E2>
		auto apply_operators_mul(OpAdd<Es...> const& a, OpOperatorChain<OpIdentity, E2> const& b)
		{
			return (a * b);
		}

		//! Apply the chain operation to an expression.
		template<typename A2, typename E>
		auto apply_operators_mul(OpOperatorChain<OpIdentity, A2> const& a, OpExpression<E> const& b)
		{
			return (a * *static_cast<E const*>(&b));
		}

		template<typename E1, typename E2>
		auto apply_operators_mul(OpOperatorChain<OpIdentity, E1> const& a, OpOperatorChain<OpIdentity, E2> const& b)
		{
			return OpOperatorChain(OpIdentity{}, apply_operators(a.g * b.g));
		}

		template<typename V, typename E1, typename E2>
		auto apply_operators_mul(OpOperator<E1> const& a, OpChain<V, OpIdentity, E2> const& b)
		{
			return b.combination.f * apply_operators((*static_cast<E1 const*>(&a)) * expr::get_enclosed_expression(b));
		}
	}


	template<typename E1, typename E2>
	auto apply_operators(OpBinaryMul<E1, E2> const& e)
	{
		return apply_operators_mul(apply_operators(e.a), apply_operators(e.b));
	}


	template<typename E1, typename E2>
	auto apply_operators(OpBinaryDiv<E1, E2> const& e)
	{
		return apply_operators(e.a) / apply_operators(e.b);
	}

	namespace
	{
		template<Axis ax, typename G, typename Sp>
		using grad_term_t = OpDerivative<typename Solver<Sp>::template derivative<ax, 1>, OpIdentity, OpTerm<OpIdentity, G>, Sp>;

		template<typename E, typename Sp, typename G, Axis... axs, size_t... Is>
		auto euler_lagrange_deriv(E const& e, Sp const& solver,
			symphas::lib::types_list<G, symphas::lib::axis_list<axs...>>,
			std::index_sequence<Is...>)
		{
			return expr::apply_operators(
				expr::make_operator_derivative<1>(solver) * (
					expr::apply_operators(expr::make_column_vector<Is, sizeof...(Is)>() * expr::make_derivative<1, grad_term_t<axs, G, Sp>>(e)) + ...));
		}

		template<typename G, size_t D, typename E, typename Sp>
		auto euler_lagrange_deriv(E const& e, Sp const& solver)
		{
			return euler_lagrange_deriv(e, solver,
				symphas::lib::make_axis_list<D, G>(), std::make_index_sequence<D>{});
		}

		template<typename G, size_t D, typename E>
		auto euler_lagrange_deriv(E const& e, int)
		{
			return OpVoid{};
		}
	}

	template<typename G, typename E, typename Sp>
	auto euler_lagrange_apply(SymbolicFunctionalDerivative<G> const& symbol, OpExpression<E> const& e, Sp const& solver)
	{
		constexpr size_t D = expr::grid_dim<G>::value;

		auto interface_term = euler_lagrange_deriv<G, D>(*static_cast<E const*>(&e), solver);
		auto bulk_term = expr::apply_operators(expr::make_derivative<1, G>(*static_cast<E const*>(&e)));
		return bulk_term - interface_term;

	}


	template<typename G, typename E, typename Sp>
	auto euler_lagrange_apply(SymbolicFunctionalDerivative<DynamicVariable<G>> const& symbol, OpExpression<E> const& e, Sp const& solver)
	{
		constexpr size_t D = expr::grid_dim<G>::value;

		auto interface_term = euler_lagrange_deriv<DynamicVariable<G>, D>(*static_cast<E const*>(&e), solver);
		auto bulk_term = expr::apply_operators(expr::make_derivative<1, DynamicVariable<G>>(*static_cast<E const*>(&e), symbol.index));
		return bulk_term - interface_term;
		//using swap_t = GridSymbol<expr::eval_type_t<G>, D>;
		//return expr::transform::swap_grid<swap_t, expr::symbols::placeholder_N_symbol>(result, expr::make_term_dynamic, symbol.index);

	}

	template<typename G, size_t D, typename E, typename Sp>
	auto euler_lagrange_apply(SymbolicFunctionalDerivative<GridSymbol<G, D>>, OpExpression<E> const& e, Sp const& solver)
	{
		auto interface_term = euler_lagrange_deriv<GridSymbol<G, D>, D>(*static_cast<E const*>(&e), solver);
		auto bulk_term = expr::apply_operators(expr::make_derivative<1, GridSymbol<G, D>>(*static_cast<E const*>(&e)));
		return bulk_term - interface_term;
	} 

	template<typename V, typename V0, typename E, typename T, typename G>
	auto apply_operators(OpFunctionalDerivative<V, OpDomainIntegral<V0, E, T>, G> const& e)
	{
		auto integral = expr::get_enclosed_expression(e);
		auto expr = expr::get_enclosed_expression(integral);
		return expr::coeff(e) * expr::coeff(integral) * euler_lagrange_apply(e.solver, expr, e.implicit_solver);
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
			Block<T> result(expr::data_length(*static_cast<E const*>(&e)));
			expr::result(*static_cast<E const*>(&e), result.values, result.len);
			return result;
		}
	};


	//! Used for the grid swapping routines; swaps all indicies with the given replacement.
	/*!
	 * Defined using the expr::split::factor routine, implemented below. All the matching indices
	 * are factored from the expression, and then the replacements are combined into a single
	 * term and substituted with the replacement.
	 */
	template<typename E, int N0, int... P0s, typename G_F>
	auto swap_matching_i(OpExpression<E> const& e, symphas::lib::types_list<expr::symbols::i_<N0, P0s>...>, G_F&& g);
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

	//! Evaluating a literal as a grid returns just the value of the literal.
	/*!
	 * Evaluating a literal as a grid returns just the value of the literal.
	 * See expr::transform::to_grid(OpExpressioN<E>&)
	 */
	template<typename T, typename G>
	auto to_grid(OpTerm<T, G> const e)
	{
		return expr::BaseData<G>::get(expr::get<1>(e).data());
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


	//! Remove the terms for which the base_data has a match in the given list.
	template<typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs, typename... Is>
	auto filter_from_term(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& e, symphas::lib::types_list<Is...>)
	{
		if constexpr ((std::is_same<expr::base_data_t<G0>, Is>::value || ... || false))
		{
			return filter_from_term(expr::terms_after_first(e), symphas::lib::types_list<Is...>{});
		}
		else
		{
			return expr::make_term(expr::get<1>(e))
				* filter_from_term(expr::terms_after_first(e), symphas::lib::types_list<Is...>{});
		}
	}

	//! Remove the terms for which the base_data has a match in the given list.
	template<typename V, typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs, typename... Is>
	auto filter_from_term(OpTerms<V, Term<G0, X0>, Term<Gs, Xs>...> const& e, symphas::lib::types_list<Is...>)
	{
		return expr::coeff(e) * filter_from_term(expr::terms_after_first(e), symphas::lib::types_list<Is...>{});
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
	template<size_t D, typename V, typename Sp, size_t... Os>
	auto to_ft(OpOperatorMixedDerivative<V, Sp, Os...> const& e, double const* h, const len_type* dims);

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
	auto to_ft(OpConvolution<V, E1, E2> const& e, double const* h, const len_type* dims);

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
	auto to_ft(OpConvolution<V, GaussianSmoothing<D>, E> const& e, double const* h, const len_type* dims);

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
	auto to_ft(OpDerivative<Dd, V, E, Sp> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename V, typename sub_t, typename E, typename... Ts>
	auto to_ft(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename V, typename sub_t, typename E, typename... Ts>
	auto to_ft(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e, double const* h, const len_type* dims);

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
	template<size_t D, typename V, typename E, typename... Ts, int... I0s, int... P0s,
		typename A, typename... As, typename B, typename C>
	auto to_ft(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
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
	template<size_t D, typename V, size_t O>
	auto to_ft(OpTerm<V, k_grid_type<O, D>> const& e, double const* h, const len_type*)
	{
		//return expr::make_operator_derivative<O>(Solver<void>{});
	}

	template<size_t D, typename V, Axis ax, size_t O>
	auto to_ft(OpTerm<V, k_grid_axis_type<ax, O, D>> const& e, double const* h, const len_type*)
	{
		//return symphas::internal::nth_symbolic_derivative_function<ax, O, Sp>::template get(Solver<void>{});
	}

	template<size_t D, typename V, Axis ax, size_t O>
	auto to_ft(OpTerm<V, k_grid_component_type<ax, O, D>> const& e, double const* h, const len_type*)
	{
		//return symphas::internal::make_directional_operator_derivative<O>(Solver<void>{});
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
	template<size_t D, typename V, typename T, typename E, size_t... Ns, typename... Ts>
	auto to_ft(OpSymbolicEval<V, NoiseData<NoiseType::WHITE, T, D>, SymbolicFunction<E, Variable<Ns, Ts>...>> const& e, double const* h, const len_type* dims);

	template<size_t D, typename E>
	auto to_ft(OpExpression<E> const& e, double const* h, const len_type*)
	{
		return expr::make_fourier_map(*static_cast<E const*>(&e));
	}

	template<size_t D, typename E>
	auto to_ft(OpOperator<E> const& e, double const* h, const len_type*)
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
		return OpFractionLiteral<N, D>{} *to_ft<D>(OpIdentity{}, h, dims);
	}

	template<size_t D, size_t N, size_t D0>
	auto to_ft(OpNegFractionLiteral<N, D>, double const* h, const len_type* dims)
	{
		return OpNegFractionLiteral<N, D>{} *to_ft<D>(OpIdentity{}, h, dims);
	}

	template<size_t D, typename T, size_t... Ns>
	auto to_ft(OpTensor<T, Ns...> const& tensor, double const* h, const len_type* dims)
	{
		return tensor * to_ft<D>(OpIdentity{}, h, dims);
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

	namespace
	{

		template<size_t D, typename V, typename Sp, size_t... Os, Axis... axs>
		auto to_ft_mixed(OpOperatorMixedDerivative<V, Sp, Os...> const& e, double const* h, const len_type* dims, symphas::internal::axis_list<axs...>)
		{
			return expr::coeff(e) * (expr::make_term(k_grid_axis_type<axs, Os, D>(dims, h)) * ...);
		}
	}

	template<size_t D, typename V, typename Sp, size_t... Os>
	auto to_ft(OpOperatorMixedDerivative<V, Sp, Os...> const& e, double const* h, const len_type* dims)
	{
		return to_ft_mixed(e, h, dims, symphas::lib::make_axis_list<sizeof...(Os)>());
	}

	template<size_t D, typename V, typename E1, typename E2>
	auto to_ft(OpConvolution<V, E1, E2> const& e, double const* h, const len_type* dims)
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

	template<size_t D, typename V, typename sub_t, typename E, typename... Ts>
	auto to_ft(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e, double const* h, const len_type* dims)
	{
		auto ft = to_ft<D>(e.f.e, h, dims);
		auto f = (expr::function_of(Ts{}...) = ft);
		f.set_data_tuple(e.data);
		return symphas::internal::make_symbolic_eval(expr::coeff(e), e.data, f);
	}

	template<size_t D, typename V, typename E, typename... Ts, int... I0s, int... P0s,
		typename A, typename... As, typename B, typename C>
	auto to_ft(OpSum<V, E, Substitution<SymbolicDataArray<Ts>...>, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		SymbolicFunction<A, As...>, B, C> const& series,
		double const* h, const len_type* dims)
	{
		auto ft = to_ft<D>(series.data.e, h, dims);
		return expr::coeff(series) * expr::recreate_series(ft, series.data);
	}

	template<size_t D, typename V, typename E>
	auto to_ft(OpConvolution<V, GaussianSmoothing<D>, E> const& e, double const* h, const len_type* dims)
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

	namespace
	{
		template<size_t D, typename Sp, size_t... Os, Axis... axs>
		auto to_ft_mixed(Sp const& solver, typename Sp::template mixed_derivative<Os...>, symphas::internal::axis_list<axs...>, double const* h, const len_type* dims)
		{
			return (OpIdentity{} * ... * to_ft<D>(expr::make_operator_directional_derivative<axs, Os>(solver), h, dims));
		}

		template<size_t D, typename Sp, size_t... Os>
		auto to_ft_mixed(Sp const& solver, typename Sp::template mixed_derivative<Os...>, double const* h, const len_type* dims)
		{
			return to_ft_mixed<D>(solver, typename Sp::template mixed_derivative<Os...>{}, symphas::lib::make_axis_list<sizeof...(Os)>(), h, dims);
		}
	}

	template<size_t D, typename Dd, typename V, typename E, typename Sp>
	auto to_ft(OpDerivative<Dd, V, E, Sp> const& e, double const* h, const len_type* dims)
	{
		constexpr Axis axis = OpDerivative<Dd, V, E, Sp>::axis;
		constexpr size_t order = OpDerivative<Dd, V, E, Sp>::order;


		if constexpr (Dd::is_directional)
		{
			if constexpr (Dd::is_mixed)
			{
				return expr::coeff(e) * to_ft_mixed<D, Sp>(e.solver, Dd{}, h, dims)* to_ft<D>(expr::get_enclosed_expression(e), h, dims);
			}
			else
			{
				return expr::coeff(e) * expr::make_term(k_grid_axis_type<axis, order, D>(dims, h))
					* to_ft<D>(expr::get_enclosed_expression(e), h, dims);
			}
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

}


namespace symphas::internal
{

	template<typename expression_condition_t>
	struct expression_condition_impl
	{

		using this_type = expression_condition_impl<expression_condition_t>;

		template<typename... Es>
		auto operator()(OpAddList<Es...> const& e) const
		{
			return this_type{};
		}

		template<typename... Gs>
		auto operator()(OpTermsList<Gs...> const& e) const
		{
			return this_type{};
		}

		template<typename... Es>
		auto operator()(OpAdd<Es...> const& e) const
		{
			return this_type{};
		}

		template<typename... Gs>
		auto operator()(OpTerms<Gs...> const& e) const
		{
			return this_type{};
		}

		template<typename E>
		auto operator()(OpEvaluable<E> const& e) const
		{
			return this_type{};
		}

		auto operator()(...) const
		{
			return this_type{};
		}

		template<typename E>
		auto get_value() const
		{
			return cast().operator()(E{});
		}

		expression_condition_t& cast()
		{
			return *static_cast<expression_condition_t*>(this);
		}

		const expression_condition_t& cast() const
		{
			return *static_cast<expression_condition_t const*>(this);
		}

		template<typename E>
		using type = std::invoke_result_t<decltype(&this_type::template get_value<E>), this_type>;

		template<typename E>
		static const bool value = !std::is_same<type<E>, this_type>::value;

	};

	template<typename... expression_condition_ts>
	struct not_expression_condition {};

	template<typename... expression_condition_ts>
	struct and_expression_condition {};

	template<typename... expression_condition_ts>
	struct or_expression_condition {};

	template<typename E, typename expression_condition_t>
	constexpr bool expression_satisfies_condition = expression_condition_t::template value<E>;

	template<typename E>
	constexpr bool expression_satisfies_condition<E, void> = false;

	template<typename E, typename... expression_condition_ts>
	constexpr bool expression_satisfies_condition<E, not_expression_condition<expression_condition_ts...>> =
		(!expression_satisfies_condition<E, expression_condition_ts> && ... && true);

	template<typename E, typename... expression_condition_ts>
	constexpr bool expression_satisfies_condition<E, and_expression_condition<expression_condition_ts...>> =
		(expression_satisfies_condition<E, expression_condition_ts> && ... && true);

	template<typename E, typename... expression_condition_ts>
	constexpr bool expression_satisfies_condition<E, or_expression_condition<expression_condition_ts...>> =
		(expression_satisfies_condition<E, expression_condition_ts> || ... || (sizeof...(expression_condition_ts) == 0));

}



// **************************************************************************************
// Split an expression into groups for evaluating by individual regions
// **************************************************************************************

namespace expr
{
	template<typename condition_t>
	using expr_cond = symphas::internal::expression_condition_impl<condition_t>;
	template<typename... condition_ts>
	using not_ = symphas::internal::not_expression_condition<condition_ts...>;
	template<typename... condition_ts>
	using and_ = symphas::internal::and_expression_condition<condition_ts...>;
	template<typename... condition_ts>
	using or_ = symphas::internal::or_expression_condition<condition_ts...>;

	template<typename E, typename expression_condition_t>
	constexpr bool satisfies = symphas::internal::expression_satisfies_condition<E, expression_condition_t>;



	struct matches_any : expr_cond<matches_any>
	{
		using expr_cond<matches_any>::operator();
		template<typename E>
		auto operator()(OpEvaluable<E> const& e) const
		{
			return matches_any{};
		}
	};

	struct matches_series : expr_cond<matches_series>
	{
		using expr_cond<matches_series>::operator();
		template<typename V0, typename E, typename... Ts, int... I0s, int... P0s, typename A, typename B, typename... Vs>
		auto operator()(OpSum<V0, E,
			Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
			A, B, symphas::lib::types_list<Vs...>> const& sum) const
		{
			return matches_series{};
		}
	};

	struct matches_mul : expr_cond<matches_mul>
	{
		using expr_cond<matches_mul>::operator();
		template<typename A, typename B>
		auto operator()(OpBinaryMul<A, B> const& e) const
		{
			return matches_mul{};
		}
	};

	struct matches_div : expr_cond<matches_div>
	{
		using expr_cond<matches_div>::operator();
		template<typename A, typename B>
		auto operator()(OpBinaryDiv<A, B> const& e) const
		{
			return matches_div{};
		}
	};

	struct matches_term : expr_cond<matches_term>
	{
		using expr_cond<matches_term>::operator();
		template<typename V, typename... Gs, expr::exp_key_t... Xs>
		auto operator()(OpTerms<V, Term<Gs, Xs>...> const& e) const
		{
			return matches_term{};
		}
	};

	struct matches_integral : expr_cond<matches_integral>
	{
		using expr_cond<matches_integral>::operator();
		template<typename V, typename E, typename T>
		auto operator()(OpIntegral<V, E, T> const& e) const
		{
			return matches_integral{};
		}
	};

	struct matches_derivative : expr_cond<matches_derivative>
	{
		using expr_cond<matches_derivative>::operator();
		template<typename Dd, typename V, typename E, typename Sp>
		auto operator()(OpDerivative<Dd, V, E, Sp> const& e) const
		{
			return matches_derivative{};
		}
	};

	struct matches_operator : expr_cond<matches_operator>
	{
		using expr_cond<matches_operator>::operator();
		template<typename E>
		auto operator()(OpOperator<E> const& e) const
		{
			return matches_operator{};
		}
	};

	template<typename matches_t, typename A, typename B>
	struct matching_in_mul_apply
	{
		static const bool value =
			symphas::internal::expression_satisfies_condition<A, matches_t>
			|| symphas::internal::expression_satisfies_condition<B, matches_t>;
	};

	template<typename matches_t, typename A, typename B, typename C>
	struct matching_in_mul_apply<matches_t, A, OpBinaryMul<B, C>>
	{
		static const bool value = symphas::internal::expression_satisfies_condition<A, matches_t> || matching_in_mul_apply<matches_t, B, C>::value;
	};

	template<typename matches_t, typename A, typename B, typename C>
	struct matching_in_mul_apply<matches_t, OpBinaryMul<A, B>, C>
	{
		static const bool value = matching_in_mul_apply<matches_t, A, B>::value || symphas::internal::expression_satisfies_condition<C, matches_t>;
	};

	template<typename matches_t, typename A, typename B, typename C, typename D>
	struct matching_in_mul_apply<matches_t, OpBinaryMul<A, B>, OpBinaryMul<C, D>>
	{
		static const bool value = matching_in_mul_apply<matches_t, A, B>::value || matching_in_mul_apply<matches_t, C, D>::value;
	};

	template<typename matches_t>
	struct matching_in_mul : expr_cond<matching_in_mul<matches_t>>
	{
		using expr_cond<matching_in_mul<matches_t>>::operator();
		template<typename A, typename B>
		auto operator()(OpBinaryMul<A, B> const& e) const
		{
			if constexpr (matching_in_mul_apply<matches_t, A, B>::value)
			{
				return matches_operator{};
			}
			else
			{
				return expr_cond<matching_in_mul<matches_t>>::operator()(OpVoid{});
			}
		}
	};


	template<typename E>
	struct matches_with : expr_cond<matches_with<E>>
	{
		using expr_cond<matches_with<E>>::operator();
		auto operator()(E const& e) const
		{
			return matches_with<E>{};
		}
	};

	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	struct matches_with<OpTerms<V, Term<Gs, Xs>...>> : expr_cond<matches_with<OpTerms<V, Term<Gs, Xs>...>>>
	{
		using expr_cond<matches_with<OpTerms<V, Term<Gs, Xs>...>>>::operator();
		template<typename V0, typename... G0s, expr::exp_key_t... X0s>
		auto operator()(OpTerms<V0, Term<G0s, X0s>...> const& e) const
		{
			using namespace symphas::lib;

			using G_list = types_list<Gs...>;
			using G0_list = types_list<G0s...>;
			using X_seq = std::integer_sequence<expr::exp_key_t, Xs...>;
			using X0_seq = std::integer_sequence<expr::exp_key_t, X0s...>;
			if constexpr (
				types_list_size<expand_types_list<filter_types<G_list, G0_list>, filter_types<G0_list, G_list>>>::value == 0
				&& seq_join_t<filter_seq_t<X_seq, X0_seq>, filter_seq_t<X0_seq, X_seq>>::size() == 0
				&& sizeof...(Gs) == sizeof...(G0s))
			{
				return matches_with<OpTerms<V, Term<Gs, Xs>...>>{};
			}
			else
			{
				return expr_cond<matches_with<OpTerms<V, Term<Gs, Xs>...>>>::operator()(OpVoid{});
			}
		}
	};


	template<typename M>
	struct contains_satisfying : expr_cond<contains_satisfying<M>>
	{
		using self_match_t = contains_satisfying<M>;

		template<typename E, std::enable_if_t<satisfies<E, M>, int> = 0>
		auto operator()(OpEvaluable<E> const& e) const
		{
			return contains_satisfying<M>{};
		}

		template<typename E, std::enable_if_t<!satisfies<E, M>, int> = 0>
		auto operator()(OpEvaluable<E> const& e) const
		{
			return expr_cond<contains_satisfying<M>>::operator()(OpVoid{});
		}

		template<typename... Es, std::enable_if_t<(satisfies<Es, M> || ...), int> = 0>
		auto operator()(OpAdd<Es...> const& e) const
		{
			return contains_satisfying<M>{};
		}
	};

	template<typename M>
	struct contains_satisfying_anywhere : expr_cond<contains_satisfying_anywhere<M>>
	{
		using self_match_t = contains_satisfying_anywhere<M>;

		template<typename E, std::enable_if_t<satisfies<E, M>, int> = 0>
		auto operator()(OpEvaluable<E> const& e) const
		{
			return contains_satisfying_anywhere<M>{};
		}

		template<typename E, std::enable_if_t<!satisfies<E, M>, int> = 0>
		auto operator()(OpEvaluable<E> const& e) const
		{
			return expr_cond<contains_satisfying_anywhere<M>>::operator()(OpVoid{});
		}

		template<typename A, typename B, typename E, 
			std::enable_if_t<(satisfies<OpOperatorCombination<A, B>, self_match_t> || satisfies<E, self_match_t>), int> = 0>
		auto operator()(OpCombination<A, B, E> const& e) const
		{
			return contains_satisfying_anywhere<M>{};
		}

		template<typename A, typename B, typename E, 
			std::enable_if_t<(satisfies<OpOperatorChain<A, B>, self_match_t> || satisfies<E, self_match_t>), int> = 0>
		auto operator()(OpChain<A, B, E> const& e) const
		{
			return contains_satisfying_anywhere<M>{};
		}

		template<typename A, typename B, std::enable_if_t<(satisfies<A, self_match_t> || satisfies<B, self_match_t>), int> = 0>
		auto operator()(OpOperatorCombination<A, B> const& e) const
		{
			return contains_satisfying_anywhere<M>{};
		}

		template<typename A, typename B, std::enable_if_t<(satisfies<A, self_match_t> || satisfies<B, self_match_t>), int> = 0>
		auto operator()(OpOperatorChain<A, B> const& e) const
		{
			return contains_satisfying_anywhere<M>{};
		}

		template<typename A, typename B, std::enable_if_t<(satisfies<A, self_match_t> || satisfies<B, self_match_t>), int> = 0>
		auto operator()(OpBinaryDiv<A, B> const& e) const
		{
			return contains_satisfying_anywhere<M>{};
		}

		template<typename A, typename B, std::enable_if_t<(satisfies<A, self_match_t> || satisfies<B, self_match_t>), int> = 0>
		auto operator()(OpBinaryMul<A, B> const& e) const
		{
			return contains_satisfying_anywhere<M>{};
		}

		template<typename... Es, std::enable_if_t<(satisfies<Es, self_match_t> || ...), int> = 0>
		auto operator()(OpAdd<Es...> const& e) const
		{
			return contains_satisfying_anywhere<M>{};
		}
	};

	template<typename E>
	using contains_matching = contains_satisfying<matches_with<E>>;
	template<typename E>
	using contains_matching_anywhere = contains_satisfying_anywhere<matches_with<E>>;

	// I can introduce other criteria for terms like, is linear, or is operator
}

namespace symphas::internal
{

	template<typename E, typename M>
	constexpr bool expression_satisfies_condition<E, expr::matches_with<expr::matches_with<M>>> =
		(expression_satisfies_condition<E, expr::matches_with<M>>);

}

namespace expr::transform
{

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
	template<size_t Z, typename E, typename G_F>
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
	template<size_t Z, typename V, typename... Gs, exp_key_t... Xs, typename G_F>
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
	template<size_t Z, size_t O, typename V, typename Sp, typename G_F>
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
	template<size_t Z, Axis ax, size_t O, typename V, typename Sp, typename G_F>
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
	template<size_t Z, size_t... Os, typename V, typename Sp, typename G_F>
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
	template<size_t Z, typename A1, typename E, typename G_F>
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
	template<size_t Z, typename V, size_t D, typename E, typename G_F>
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
	template<size_t Z, typename Dd, typename V, typename E, typename Sp, typename G_F>
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
	template<size_t Z, size_t O, typename V, typename E, typename GG, typename G_F>
	auto swap_grid(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& e, G_F&& g);

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
	template<size_t Z, typename V, typename E, typename T, typename G_F>
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
	template<size_t Z, typename G, typename V, typename E, typename G_F>
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
	template<size_t Z, expr::exp_key_t X, typename V, typename E, typename G_F>
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
	template<size_t Z, typename V, typename E, typename F, typename Arg0, typename... Args, typename G_F>
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
	template<size_t Z, auto f, typename V, typename E, typename G_F>
	auto swap_grid(OpFunctionApply<f, V, E> const& e, G_F&& g);

	template<size_t Z, typename V, typename sub_t, typename E, typename... Ts, typename G_F>
	auto swap_grid(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e, G_F&& g);

	template<size_t Z, int N, int P, typename G_F>
	auto swap_grid(expr::symbols::i_<N, P> const& e, G_F&& g)
	{
		return expr::symbols::i_<N, P>{};
	}

	template<size_t Z, typename T, typename G_F>
	auto swap_grid(SymbolicData<T> const& e, G_F&& g)
	{
		return e;
	}

	namespace
	{

		template<typename V, typename... Gs, exp_key_t... Xs, typename... Ts, size_t I0, size_t... I0s>
		auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<I0, I0s...>, std::tuple<Ts...> const& subbed, std::index_sequence<>);
		template<typename V, typename... Gs, exp_key_t... Xs, typename... Ts, size_t I1, size_t... I1s>
		auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<>, std::tuple<Ts...> const& subbed, std::index_sequence<I1, I1s...>);
		template<typename V, typename... Gs, exp_key_t... Xs, typename... Ts>
		auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<>, std::tuple<Ts...> const& subbed, std::index_sequence<>);
		template<typename V, typename... Gs, exp_key_t... Xs, typename... Ts, size_t I0, size_t... I0s, size_t I1, size_t... I1s>
		auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<I0, I0s...>, std::tuple<Ts...> const& subbed, std::index_sequence<I1, I1s...>);

		template<typename V, typename... Gs, exp_key_t... Xs, typename... Ts, size_t I0, size_t... I0s>
		auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<I0, I0s...>, std::tuple<Ts...> const& subbed, std::index_sequence<>)
		{
			return expr::make_term(expr::get<I0 + 1>(e)) * recombine_terms(e, std::index_sequence<I0s...>{}, subbed, std::index_sequence<>{});
		}

		template<typename V, typename... Gs, exp_key_t... Xs, typename... Ts, size_t I1, size_t... I1s>
		auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<>, std::tuple<Ts...> const& subbed, std::index_sequence<I1, I1s...>)
		{
			return std::get<I1>(subbed) * recombine_terms(e, std::index_sequence<>{}, subbed, std::index_sequence<I1s...>{});
		}

		template<typename V, typename... Gs, exp_key_t... Xs, typename... Ts>
		auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<>, std::tuple<Ts...> const& subbed, std::index_sequence<>)
		{
			return OpIdentity{};
		}

		template<typename V, typename... Gs, exp_key_t... Xs, typename... Ts, size_t I0, size_t... I0s, size_t I1, size_t... I1s>
		auto recombine_terms(OpTerms<V, Term<Gs, Xs>...> const& e, std::index_sequence<I0, I0s...>, std::tuple<Ts...> const& subbed, std::index_sequence<I1, I1s...>)
		{
			if constexpr (I0 < I1)
			{
				return expr::make_term(expr::get<I0 + 1>(e)) * recombine_terms(e, std::index_sequence<I0s...>{}, subbed, std::index_sequence<I1, I1s...>{});
			}
			else
			{
				return std::get<I1>(subbed) * recombine_terms(e, std::index_sequence<I0, I0s...>{}, subbed, std::index_sequence<I1s...>{});
			}
		}


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

		template<typename V, typename... G1s, exp_key_t... X1s, size_t... Ns, size_t... Ms, exp_key_t X, typename G_F>
		auto _pick_terms(OpTerms<V, Term<G1s, X1s>...> const& a, 
			std::index_sequence<Ns...>, std::index_sequence<Ms...>, 
			std::integer_sequence<exp_key_t, X>, G_F&& g)
		{
			if constexpr (!expr::is_expression<G_F>)
			{
				if constexpr (is_simple_data<G_F>)
				{
					return sift_term(sift_term(g), expr::get<Ns>(a)..., Term(expr::make_literal(std::forward<G_F>(g))).template pow<X>(), expr::get<Ms>(a)...);
				}
				else
				{
					return sift_term(OpIdentity{}, expr::get<Ns>(a)..., Term(std::forward<G_F>(g)).template pow<X>(), expr::get<Ms>(a)...);
				}
			}
			else
			{
				if constexpr (_Xk_t<X>::D == 1)
				{
					if constexpr (_Xk_t<X>::N == 1)
					{
						return sift_term(OpIdentity{}, expr::get<Ns>(a)...)
							* std::forward<G_F>(g) 
							* sift_term(OpIdentity{}, expr::get<Ms>(a)...);
					}
					else
					{
						return sift_term(OpIdentity{}, expr::get<Ns>(a)...)
							* expr::pow<_Xk_t<X>::N>(std::forward<G_F>(g)) 
							* sift_term(OpIdentity{}, expr::get<Ms>(a)...);
					}
				}
				else
				{
					return sift_term(OpIdentity{}, expr::get<Ns>(a)...)
						* expr::exp(
							((_Xk_t<X>::sign) ? OpNegIdentity{} : OpIdentity{}) * expr::make_fraction<_Xk_t<X>::N, _Xk_t<X>::D>()
							* expr::log(std::forward<G_F>(g)))
						* sift_term(OpIdentity{}, expr::get<Ms>(a)...);
				}

			}
		}

		template<typename V, typename... G1s, exp_key_t... X1s, typename G_F>
		auto pick_terms(OpTerms<V, Term<G1s, X1s>...> const& e,
			std::index_sequence<>, G_F&& g)
		{
			return expr::make_term(expr::terms_after_first(e));
		}

		template<typename V, typename... G1s, exp_key_t... X1s, size_t N, typename G_F>
		auto pick_terms(OpTerms<V, Term<G1s, X1s>...> const& e,
			std::index_sequence<N>, G_F&& g)
		{
			using namespace symphas::lib;

			using pick_terms1_t = seq_offset_t<1, std::make_index_sequence<N>>;
			using pick_terms2_t = seq_offset_t<N + 2, std::make_index_sequence<sizeof...(G1s) - N - 1>>;
			using pick_power_t = symphas::lib::type_at_index<N, std::integer_sequence<exp_key_t, X1s>...>;

			return _pick_terms(e, pick_terms1_t{}, pick_terms2_t{}, pick_power_t{}, std::forward<G_F>(g));
		}

		template<typename V, typename... G1s, exp_key_t... X1s, size_t N0, size_t N1, size_t... Ns, typename G_F>
		auto pick_terms(OpTerms<V, Term<G1s, X1s>...> const& e,
			std::index_sequence<N0, N1, Ns...>, G_F&& g)
		{
			auto list = std::make_tuple(std::forward<G_F>(g), std::forward<G_F>(g), ([&] (auto) { return std::forward<G_F>(g);  })(Ns)...);
			using seq_t = symphas::lib::filter_seq_t<std::make_index_sequence<sizeof...(G1s)>, std::index_sequence<Ns...>>;
			return recombine_terms(e, seq_t{}, list, std::make_index_sequence<sizeof...(Ns) + 2>{});
		}


		template<typename V, typename... Gs, exp_key_t... Xs, size_t... Ns>
		auto select_index(OpTerms<V, Term<Gs, Xs>...> const& e,
			std::index_sequence<Ns...>, DynamicIndexSet g)
		{
			using namespace symphas::lib;

			((expr::get<Ns + 1>(const_cast<OpTerms<V, Term<Gs, Xs>...>&>(e)).data().fix(g)), ...);
			return e;
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

			return pick_terms(e, swap_seq_t{}, std::forward<G_F>(g));
		}

		template<size_t Z, typename... Es, typename G_F, size_t... Is>
		auto swap_grid_adds(OpAdd<Es...> const& e, G_F&& g, std::index_sequence<Is...>)
		{
			return (swap_grid<Z>(expr::get<Is>(e), std::forward<G_F>(g)) + ... + OpVoid{});
		}
	}

	template<size_t Z, typename E, typename G_F>
	auto swap_grid(OpExpression<E> const& e, G_F&&)
	{
		return *static_cast<E const*>(&e);
	}

	template<size_t Z, typename E, typename G_F>
	auto swap_grid(OpOperator<E> const& e, G_F&&)
	{
		return *static_cast<E const*>(&e);
	}

	template<size_t Z, typename V, typename... Gs, exp_key_t... Xs, typename G_F>
	auto swap_grid(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g)
	{
		using mask_t = std::integer_sequence<bool, (expr::factor_count<Variable<Z>, Gs>::value > 0)...>;
		return expr::coeff(e) * swap_terms(e, std::make_index_sequence<sizeof...(Gs)>{}, mask_t{}, std::forward<G_F>(g));
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
	auto swap_grid(OpConvolution<V, E1, E2> const& e, G_F&& g)
	{
		return expr::make_convolution(
			e.value,
			swap_grid<Z>(e.a, std::forward<G_F>(g)), 
			swap_grid<Z>(e.b, std::forward<G_F>(g)));
	}

	template<size_t Z, typename V, size_t D, typename E, typename G_F>
	auto swap_grid(OpConvolution<V, GaussianSmoothing<D>, E> const& e, G_F&& g)
	{
		return expr::make_convolution(
			e.value, 
			swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), 
			e.smoother);
	}

	template<size_t Z, size_t O, typename V, typename Sp, typename G_F>
	auto swap_grid(OpOperatorDerivative<O, V, Sp> const& e, G_F&& g)
	{
		return e;
	}

	template<size_t Z, Axis ax, size_t O, typename V, typename Sp, typename G_F>
	auto swap_grid(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e, G_F&& g)
	{
		return e;
	}

	template<size_t Z, size_t... Os, typename V, typename Sp, typename G_F>
	auto swap_grid(OpOperatorMixedDerivative<V, Sp, Os...> const& e, G_F&& g)
	{
		return e;
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

	template<size_t Z, typename A1, typename E, typename G_F>
	auto swap_grid(OpChain<A1, OpIdentity, E> const& e, G_F&& g)
	{
		return OpChain(
			OpOperatorChain(swap_grid<Z>(e.combination.f, std::forward<G_F>(g)), OpIdentity{}),
			swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
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
	auto swap_grid(OpDerivative<Dd, V, E, Sp> const& e, G_F&& g)
	{
		constexpr size_t order = OpDerivative<Dd, V, E, Sp>::order;
		constexpr Axis axis = OpDerivative<Dd, V, E, Sp>::axis;
		return symphas::internal::nth_derivative_apply<axis, order, Sp>::template get(
			e.value, swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), e.solver);
	}

	template<size_t Z, size_t O, typename V, typename E, typename GG, typename G_F>
	auto swap_grid(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& e, G_F&& g)
	{
		return expr::make_derivative<O, GG>(e.value, swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), e.solver);
	}

	template<size_t Z, typename V, typename E, typename T, typename G_F>
	auto swap_grid(OpIntegral<V, E, T> const& e, G_F&& g)
	{
		return expr::make_integral(expr::coeff(e), swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), e.domain);
	}

	template<size_t Z, typename G, typename V, typename E, typename G_F>
	auto swap_grid(OpMap<G, V, E> const& e, G_F&& g)
	{
		auto eg = swap_grid<Z>(e.e, std::forward<G_F>(g));
		return OpMap<G, V, decltype(eg)>(e.value, eg);
	}

	template<size_t Z, expr::exp_key_t X, typename V, typename E, typename G_F>
	auto swap_grid(OpPow<X, V, E> const& e, G_F&& g)
	{
		return expr::make_pow<X>(expr::coeff(e), swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
	}

	template<size_t Z, typename V, typename E, typename F, typename Arg0, typename... Args, typename G_F>
	auto swap_grid(OpFunction<V, E, F, Arg0, Args...> const& e, G_F&& g)
	{
		return OpFunction(e.name, e.value, swap_grid<Z>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), e.f, e.tt);
	}

	template<size_t Z, auto f, typename V, typename E, typename G_F>
	auto swap_grid(OpFunctionApply<f, V, E> const& e, G_F&& g)
	{
        auto eg = swap_grid<Z>(e.e, std::forward<G_F>(g));
		return expr::coeff(e) * expr::make_function<f>(eg);
	}

	namespace
	{

		template<size_t Z, typename E, size_t... ArgNs, typename... Ts, size_t... ArgMs, size_t... Is, typename G_F>
		auto swap_grid_symbolic(
			SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f,
			SymbolicTemplate<E, ArgNs...> const& tmpl, std::index_sequence<ArgMs...>,
			std::index_sequence<Is...>, G_F&& g)
		{
			auto swapped = swap_grid<Z>(tmpl(swap_grid<Z>(std::get<Is>(f.data), std::forward<G_F>(g))...), std::forward<G_F>(g));
			return (expr::function_of(as_variable<ArgMs>(swap_grid<Z>(std::get<Is>(f.data), std::forward<G_F>(g)))...) = swapped);
		}

		template<size_t Z, typename E, size_t... ArgNs, typename... Ts, size_t... ArgMs, size_t... Is, typename G_F>
		auto swap_grid_symbolic(
			SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f,
			OpExpression<E> const& e, std::index_sequence<ArgMs...>,
			std::index_sequence<Is...>, G_F&& g)
		{
			auto swapped = swap_grid<Z>(*static_cast<E const*>(&e), std::forward<G_F>(g));
			return (expr::function_of(as_variable<ArgMs>(swap_grid<Z>(std::get<Is>(f.data), std::forward<G_F>(g)))...) = swapped);
		}

		template<size_t Z, typename E, size_t... ArgNs, typename... Ts, typename G_F>
		auto swap_grid_symbolic(SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f, G_F&& g)
		{
			using var_g = decltype(expr::get_independent_variables(std::forward<G_F>(g)));
			using seq_t = std::make_index_sequence<fixed_max<sizeof...(Ts), var_g::size()> +1>;
			using seq_filt_t = symphas::lib::filter_seq_t<seq_t, var_g>;
			using seq_cut_t = symphas::lib::seq_lt_t<sizeof...(Ts), seq_filt_t>;
			return swap_grid_symbolic<Z>(
				f, expr::template_of(Variable<ArgNs>{}...) = f.e,
				seq_cut_t{}, std::make_index_sequence<sizeof...(Ts)>{}, std::forward<G_F>(g));
		}

		template<size_t Z, typename T, typename G_F>
		auto swap_grid_symbolic(SymbolicData<T> const& e, G_F&& g)
		{
			return e;
		}

		template<size_t Z, typename T, typename G_F>
		auto swap_grid_symbolic(SymbolicDataArray<T> const& e, G_F&& g)
		{
			return e;
		}


		template<bool flag, typename A, typename B>
		decltype(auto) switch_tuple_element(A&& a, B&& b)
		{
			if constexpr (flag)
			{
				return std::forward<A>(a);
			}
			else
			{
				return std::forward<B>(b);
			}
		}

		template<typename... Ts, bool... Bs, size_t... Is, typename G_F>
		auto swap_symbolic_data_array(SymbolicDataArray<std::tuple<Term<Ts>...>> const& data,
			std::integer_sequence<bool, Bs...>, std::index_sequence<Is...>, G_F&& g)
		{
			return std::make_tuple(switch_tuple_element<Bs>(std::forward<G_F>(g), std::get<Is>(data.get_data_tuple()))...);
		}

		template<size_t Z, size_t... Zs, typename... Ts, typename G_F>
		auto swap_grid_symbolic(SymbolicDataArray<std::tuple<Term<Variable<Zs, Ts>>...>> const& data, G_F&& g)
		{
			using mask_seq = std::integer_sequence<bool, (Z == Zs)...>;
			return SymbolicDataArray(swap_symbolic_data_array(data, mask_seq{}, std::make_index_sequence<sizeof...(Ts)>{}, std::forward<G_F>(g)));
		}

		template<size_t Z, typename... Ts, size_t... Is, typename G_F>
		auto swap_grid_symbolic(Substitution<Ts...> const& data, std::index_sequence<Is...>, G_F&& g)
		{
			return Substitution(swap_grid_symbolic<Z>(std::get<Is>(data), std::forward<G_F>(g))...);
		}

		template<size_t Z, typename... Ts, typename G_F>
		auto swap_grid_symbolic(Substitution<Ts...> const& data, G_F&& g)
		{
			return swap_grid_symbolic<Z>(data, std::make_index_sequence<sizeof...(Ts)>{}, std::forward<G_F>(g));
		}

		template<size_t Z, expr::NoiseType nt, typename T, size_t D, typename G_F>
		auto swap_grid_symbolic(NoiseData<nt, T, D> const& data, G_F&& g)
		{
			return data;
		}

		template<size_t Z, typename G_F>
		auto swap_grid_symbolic(DynamicIndex const& data, G_F&& g)
		{
			return data;
		}

		template<size_t Z, int N, int P, typename G_F>
		auto swap_grid_symbolic(expr::symbols::i_<N, P> const& data, G_F&& g)
		{
			return data;
		}

		template<size_t Z, typename T, typename E0, typename... T0s, typename G_F>
		auto swap_grid_symbolic(T const& data, SymbolicFunction<E0, T0s...> const& f, G_F&& g)
		{
			auto e = swap_grid_symbolic<Z>(f, std::forward<G_F>(g));
			auto substitution = swap_grid_symbolic<Z>(data, std::forward<G_F>(g));

			return symphas::internal::make_symbolic_eval(OpIdentity{}, substitution, f);
		}

		template<size_t Z, typename Op, typename... Ts, typename E, typename E0, typename... T0s, typename G_F>
		auto swap_grid_symbolic(SymbolicSeries<Op, Substitution<Ts...>, E> const& series, SymbolicFunction<E0, T0s...> const& f, G_F&& g)
		{
			auto e = swap_grid_symbolic<Z>(series.e, std::forward<G_F>(g));
			auto substitution = swap_grid_symbolic<Z>(series.substitution, std::forward<G_F>(g));

			return expr::series<Op>(e)(series.limits, substitution);
		}
	}


	template<size_t Z, typename V, typename sub_t, typename E, typename... Ts, typename G_F>
	auto swap_grid(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e, G_F&& g)
	{
		return expr::coeff(e) * swap_grid_symbolic<Z>(e.data, e.f, std::forward<G_F>(g));
	}

	template<size_t Z0, size_t Z1, size_t... Zs, typename E, typename G_F>
	auto swap_grid(OpExpression<E> const& e, G_F&& g)
	{
		return swap_grid<Z1, Zs...>(swap_grid<Z0>(*static_cast<E const*>(&e), std::forward<G_F>(g)), std::forward<G_F>(g));
	}

	template<size_t Z0, size_t Z1, size_t... Zs, typename E, typename G_F>
	auto swap_grid(OpOperator<E> const& e, G_F&& g)
	{
		return swap_grid<Z1, Zs...>(swap_grid<Z0>(*static_cast<E const*>(&e), std::forward<G_F>(g)), std::forward<G_F>(g));
	}

	template<size_t Z0, size_t Z1, size_t... Zs, typename E, typename G_F0, typename G_F1, typename... G_Fs>
	auto swap_grid(OpExpression<E> const& e, G_F0&& g0, G_F1&& g1, G_Fs&& ...gs)
	{
		return swap_grid<Z1, Zs...>(swap_grid<Z0>(*static_cast<E const*>(&e), std::forward<G_F0>(g0)), std::forward<G_F1>(g1), std::forward<G_Fs>(gs)...);
	}

	template<size_t Z0, size_t Z1, size_t... Zs, typename E, typename G_F0, typename G_F1, typename... G_Fs>
	auto swap_grid(OpOperator<E> const& e, G_F0&& g0, G_F1&& g1, G_Fs&& ...gs)
	{
		return swap_grid<Z1, Zs...>(swap_grid<Z0>(*static_cast<E const*>(&e), std::forward<G_F0>(g0)), std::forward<G_F1>(g1), std::forward<G_Fs>(gs)...);
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
	template<typename Sg, typename E, typename G_F,
		std::enable_if_t<!expr::satisfies<E, expr::matches_with<Sg>>, int> = 0>
	decltype(auto) swap_grid(OpEvaluable<E> const& e, G_F&& g)
	{
		return *static_cast<E const*>(&e);
	}
	
	template<typename Sg, typename E, typename G_F,
		std::enable_if_t<expr::satisfies<E, expr::matches_with<Sg>>, int> = 0>
	decltype(auto) swap_grid(OpEvaluable<E> const& e, G_F&& g)
	{
		return expr::coeff(*static_cast<E const*>(&e)) * std::forward<G_F>(g);
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
	//template<typename Sg, typename T, typename G, typename G_F, 
	//	typename std::enable_if_t<(expr::is_same_base<Sg, G> && !expr::is_expression<G_F>), int> = 0>
	//decltype(auto) swap_grid(OpTerm<T, G> const& v, G_F&& g);

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
	template<typename Sg, typename V, typename... Gs, exp_key_t... Xs, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpTerms<V, Term<Gs, Xs>...>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename... Es, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpAdd<Es...>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename E1, typename E2, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpBinaryMul<E1, E2>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename E1, typename E2, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpBinaryDiv<E1, E2>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename A1, typename A2, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpOperatorCombination<A1, A2>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename A1, typename A2, typename E, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpCombination<A1, A2, E>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, size_t O, typename V, typename Sp, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpOperatorDerivative<O, V, Sp>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, Axis ax, size_t O, typename V, typename Sp, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpOperatorDirectionalDerivative<ax, O, V, Sp>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, size_t... Os, typename V, typename Sp, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpOperatorMixedDerivative<V, Sp, Os...>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename A1, typename A2, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpOperatorChain<A1, A2>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename A1, typename A2, typename E, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpChain<A1, A2, E>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename A1, typename E, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpChain<A1, OpIdentity, E>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename V, typename E1, typename E2, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpConvolution<V, E1, E2>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename V, size_t D, typename E, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpConvolution<V, GaussianSmoothing<D>, E>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename Dd, typename V, typename E, typename Sp, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpDerivative<Dd, V, E, Sp>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, size_t O, typename V, typename E, typename GG, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>>, expr::matches_with<Sg>>, int>>
	auto swap_grid(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& e, G_F&& g);

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
	template<typename Sg, typename V, typename E, typename T, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpIntegral<V, E, T>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename G, typename V, typename E, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpMap<G, V, E>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, expr::exp_key_t X, typename V, typename E, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpPow<X, V, E>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename V, typename E, typename F, typename Arg0, typename... Args, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpFunction<V, E, F, Arg0, Args...>, expr::matches_with<Sg>>, int>>
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
	template<typename Sg, typename V, typename sub_t, typename E, typename... Ts, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>>, expr::matches_with<Sg>>, int>>
	auto swap_grid(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e, G_F&& g);

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
	template<typename Sg, typename V, int N, int P, typename E, typename... Ts, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpSymbolicEval<V, expr::symbols::i_<N, P>, SymbolicFunction<E, Ts...>>, expr::matches_with<Sg>>, int>>
	auto swap_grid(OpSymbolicEval<V, expr::symbols::i_<N, P>, SymbolicFunction<E, Ts...>> const& e, G_F&& g);

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
	template<typename Sg, typename V, typename E0, typename K, typename E, typename... Ts, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpSymbolicEval<V, SymbolicListIndex<E0, K>, SymbolicFunction<E, Ts...>>, expr::matches_with<Sg>>, int>>
	auto swap_grid(OpSymbolicEval<V, SymbolicListIndex<E0, K>, SymbolicFunction<E, Ts...>> const& e, G_F&& g);

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
	template<typename Sg, auto f, typename V, typename E, typename G_F,
		typename = std::enable_if_t<!expr::satisfies<OpFunctionApply<f, V, E>, expr::matches_with<Sg>>, int>>
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
	template<typename... Sgs, typename C, typename E, typename... G_Fs>
	auto swap_grid(SymbolicCaseSwap<C>, OpExpression<E> const& e, G_Fs&&... gs);

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
	template<typename... Sgs, typename C, typename E, typename... G_Fs>
	auto swap_grid(SymbolicCaseSwap<C>, OpOperator<E> const& e, G_Fs&&... gs);

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
	template<typename Sg, typename T, typename I, typename G_F>
	auto swap_grid(OpCoeff<T, I> const&, G_F&& g);

	template<typename Sg, typename T, typename I, size_t... Ns, typename G_F>
	auto swap_grid(OpTensor<OpCoeff<T, I>, Ns...> const& coeff, G_F&& g);

	template<typename Sg, typename G_F>
	auto swap_grid(DynamicIndex const& data, G_F&& g);

	template<typename T>
	using strip_qualifiers_t = std::remove_const_t<std::remove_reference_t<T>>;


	namespace
	{

		template<typename... As, typename... Bs, exp_key_t X, typename C, typename... Sgs, size_t... Is, typename... G_Fs>
		auto handle_case(Term<SymbolicCase<expr::case_entry<As, Bs>...>, X> const& term, 
			SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>, std::tuple<G_Fs...> const& gs, std::index_sequence<Is...>)
		{
			constexpr int N = symphas::lib::index_of_type<C, As...>;
			if constexpr (N >= 0)
			{
				if constexpr (sizeof...(G_Fs) > 0)
				{
					auto subbed_case = swap_grid<Sgs...>(std::get<size_t(N)>(term.data().cases), std::get<Is>(gs)...);
					return expr::pow_x<X>(subbed_case);
				}
				else
				{
					return expr::pow_x<X>(std::get<size_t(N)>(term.data().cases));
				}
			}
			else
			{
				return expr::make_term(term);
			}
		}

		template<typename V, typename... Gs, exp_key_t... Xs, typename C, typename... Sgs, size_t... Is, typename... G_Fs>
		auto handle_cases(OpTerms<V, Term<Gs, Xs>...> const& e, SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>, 
			std::index_sequence<Is...>, std::tuple<G_Fs...> const& gs)
		{
			return std::make_tuple(
				handle_case(
					expr::get<Is + 1>(e), 
					SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>{}, 
					gs,
					std::make_index_sequence<sizeof...(G_Fs)>{})...);
		}

		template<typename V, typename... Gs, exp_key_t... Xs, typename C, typename... Sgs, size_t... Is, bool... fs, typename... G_Fs>
		decltype(auto) swap_terms_case(OpTerms<V, Term<Gs, Xs>...> const& e, SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>,
			std::index_sequence<Is...>, std::integer_sequence<bool, fs...>, std::tuple<G_Fs...> const& gs)
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
				auto subbed_cases = handle_cases(e, SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>{}, swap_seq_t{}, gs);
				return recombine_terms(e, symphas::lib::filter_seq_t<std::make_index_sequence<sizeof...(Gs)>, swap_seq_t>{}, subbed_cases, std::make_index_sequence<swap_seq_t::size()>{});
			}
			else
			{
				return expr::make_term(expr::terms_after_first(e));
			}
		}
		
		template<typename V, typename... Gs, exp_key_t... Xs, size_t... Is, bool... fs, typename G_F>
		decltype(auto) swap_terms_case(OpTerms<V, Term<Gs, Xs>...> const& e, SymbolicCaseSwap<>,
			std::index_sequence<Is...>, std::integer_sequence<bool, fs...>, G_F&& g)
		{
			using namespace symphas::lib;

			using swap_seq_t = seq_join_t<
				std::index_sequence<>,
				std::conditional_t<
					fs,
					std::index_sequence<Is>,
					std::index_sequence<>>...
				>;

			if constexpr (swap_seq_t::size() > 0)
			{
				return pick_terms(e, swap_seq_t{}, std::forward<G_F>(g));
			}
			else
			{
				return expr::make_term(expr::terms_after_first(e));
			}
		}

		template<typename V, typename... Gs, exp_key_t... Xs, size_t... Is, bool... fs>
		decltype(auto) swap_terms_case(OpTerms<V, Term<Gs, Xs>...> const& e, DynamicIndexSet,
			std::index_sequence<Is...>, std::integer_sequence<bool, fs...>, DynamicIndexSet g)
		{
			using namespace symphas::lib;

			using swap_seq_t = seq_join_t<
				std::index_sequence<>,
				std::conditional_t<
					fs,
					std::index_sequence<Is>,
					std::index_sequence<>>...
				>;

			if constexpr (swap_seq_t::size() > 0)
			{
				return select_index(e, swap_seq_t{}, g);
			}
			else
			{
				return expr::make_term(expr::terms_after_first(e));
			}
		}

		template<typename I>
		struct swap_index
		{
			template<int N, int P, typename G_F>
			auto operator()(expr::symbols::i_<N, P> const& e, G_F&& g) const
			{
				return expr::symbols::i_<N, P>{};
			}
		};


		template<int N0, int P0>
		struct swap_index<expr::symbols::i_<N0, P0>>
		{
			template<int P, int N00, int P00>
			auto operator()(expr::symbols::i_<N0, P> const& e, expr::symbols::i_<N00, P00> const&) const
			{
				return expr::symbols::i_<N00, P>{};
			}

			template<int P, typename G_F>
			auto operator()(expr::symbols::i_<N0, P> const& e, G_F const& g) const
			{
				return g;
			}

			template<int N, int P, typename G_F>
			auto operator()(expr::symbols::i_<N, P> const& e, G_F const& g) const
			{
				return expr::symbols::i_<N, P>{};
			}
		};

	}



	template<typename Sg, int N, int P, typename G_F>
	auto swap_grid(expr::symbols::i_<N, P> const& e, G_F&& g)
	{
		return swap_index<Sg>{}(e, std::forward<G_F>(g) + DynamicIndex(P));
	}

	template<typename Sg, typename G_F>
	auto swap_grid(int n, G_F&& g)
	{
		return n;
	}

	template<typename Sg, typename T, typename G_F>
	auto swap_grid(SymbolicData<T> const& e, G_F&& g)
	{
		return e;
	}

	template<bool symbolic_case_flag, bool index_flag, typename Sg>
	struct swap_grid_terms_redirect;

	template<typename Sg>
	struct swap_grid_terms_redirect<false, false, Sg>
	{
		template<typename V, typename... Gs, exp_key_t... Xs, typename G_F>
		auto operator()(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g)
		{
			using mask_t = std::integer_sequence<bool, (expr::factor_count<Sg, Gs>::value > 0)...>;
			auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
			return c * swap_terms(e, std::make_index_sequence<sizeof...(Gs)>{}, mask_t{}, std::forward<G_F>(g));
		}
	};

	template<typename Sg>
	struct swap_grid_terms_redirect<true, false, Sg>
	{
		template<typename V, typename... Gs, exp_key_t... Xs, typename G_F>
		auto operator()(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g)
		{
			using mask_t = std::integer_sequence<bool, (expr::factor_count<Sg, Gs>::value > 0)...>;
			auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
			return c * swap_terms_case(e, Sg{}, std::make_index_sequence<sizeof...(Gs)>{}, mask_t{}, std::forward<G_F>(g));
		}
	};

	template<typename Sg>
	struct swap_grid_terms_redirect<false, true, Sg>
	{
		template<typename V, typename... Gs, exp_key_t... Xs, typename G_F>
		auto operator()(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g)
		{
			using matching_ids_t = symphas::internal::select_all_i_<Sg, op_types_t<OpTerms<V, Term<Gs, Xs>...>>>;
			return symphas::internal::swap_matching_i(e, matching_ids_t{}, std::forward<G_F>(g));
		}
	};

	//template<typename Sg>
	//struct swap_grid_terms_redirect<false, false, true, Sg>
	//{
	//	template<typename V, typename... Gs, exp_key_t... Xs, typename G_F>
	//	auto operator()(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g)
	//	{
	//		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
	//		return c * symphas::internal::fix_dynamic_indices(e, std::make_index_sequence<sizeof...(Gs)>{}, std::forward<G_F>(g));
	//	}
	//};

	template<typename Sg, typename V, typename... Gs, exp_key_t... Xs, typename G_F, typename>
	auto swap_grid(OpTerms<V, Term<Gs, Xs>...> const& e, G_F&& g)
	{
		constexpr bool symbolic_case_flag = expr::factor_count<SymbolicCaseSwap<>, Sg>::value > 0 || std::is_same<Sg, DynamicIndexSet>::value;
		constexpr bool index_flag = expr::has_selected_index<Sg, OpTerms<V, Term<Gs, Xs>...>>;
		return swap_grid_terms_redirect<symbolic_case_flag, index_flag, Sg>{}(e, std::forward<G_F>(g));


		//if constexpr (expr::factor_count<SymbolicCaseSwap<>, Sg>::value > 0)
		//{
		//	using mask_t = std::integer_sequence<bool, (expr::factor_count<Sg, Gs>::value > 0)...>;
		//	auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		//	return c * swap_terms_case(e, Sg{}, std::make_index_sequence<sizeof...(Gs)>{}, mask_t{}, std::forward<G_F>(g));
		//}
		//else if constexpr (expr::has_selected_index<Sg, OpTerms<V, Term<Gs, Xs>...>>)
		//{
		//	using matching_ids_t = symphas::internal::select_all_i_<Sg, op_types_t<OpTerms<V, Term<Gs, Xs>...>>>;
		//	return symphas::internal::swap_matching_i(e, matching_ids_t{}, std::forward<G_F>(g));
		//}
		//else
		//{
		//	using mask_t = std::integer_sequence<bool, (expr::factor_count<Sg, Gs>::value > 0)...>;
		//	auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		//	return c * swap_terms(e, std::make_index_sequence<sizeof...(Gs)>{}, mask_t{}, std::forward<G_F>(g));
		//}
	}

	namespace
	{
		template<typename Sg, typename... Es, typename G_F, size_t... Is>
		auto swap_grid_adds(OpAdd<Es...> const& e, G_F&& g, std::index_sequence<Is...>)
		{
			return (swap_grid<Sg>(expr::get<Is>(e), std::forward<G_F>(g)) + ... + OpVoid{});
		}


		template<typename Sg, typename G, typename G_F>
		auto swap_grid_solver(SymbolicFunctionalDerivative<DynamicVariable<G>> const& solver, G_F&& g)
		{
			return SymbolicFunctionalDerivative<DynamicVariable<G>>(swap_grid<Sg>(solver.index, std::forward<G_F>(g)));
		}


		template<typename Sg, typename G, typename G_F>
		auto swap_grid_solver(SymbolicDerivative<DynamicVariable<G>> const& solver, G_F&& g)
		{
			return SymbolicDerivative<DynamicVariable<G>>(swap_grid<Sg>(solver.index, std::forward<G>(g)));
		}

		template<typename Sg, typename S, typename G_F>
		decltype(auto) swap_grid_solver(S const& solver, G_F&& g)
		{
			return solver;
		}
	}

	template<typename Sg, typename... Es, typename G_F, typename>
	auto swap_grid(OpAdd<Es...> const& e, G_F&& g)
	{
		return swap_grid_adds<Sg>(e, std::forward<G_F>(g), std::make_index_sequence<sizeof...(Es)>{});
	}

	template<typename Sg, typename E1, typename E2, typename G_F, typename>
	auto swap_grid(OpBinaryMul<E1, E2> const& e, G_F&& g)
	{
		return swap_grid<Sg>(e.a, std::forward<G_F>(g)) * swap_grid<Sg>(e.b, std::forward<G_F>(g));
	}

	template<typename Sg, typename E1, typename E2, typename G_F, typename>
	auto swap_grid(OpBinaryDiv<E1, E2> const& e, G_F&& g)
	{
		return swap_grid<Sg>(e.a, std::forward<G_F>(g)) / swap_grid<Sg>(e.b, std::forward<G_F>(g));
	}

	template<typename Sg, typename V, typename E1, typename E2, typename G_F, typename>
	auto swap_grid(OpConvolution<V, E1, E2> const& e, G_F&& g)
	{
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * expr::make_convolution(
			swap_grid<Sg>(e.a, std::forward<G_F>(g)),
			swap_grid<Sg>(e.b, std::forward<G_F>(g)));
	}

	template<typename Sg, typename V, size_t D, typename E, typename G_F, typename>
	auto swap_grid(OpConvolution<V, GaussianSmoothing<D>, E> const& e, G_F&& g)
	{
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * expr::make_convolution(
			swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)),
			e.smoother);
	}

	template<typename Sg, size_t O, typename V, typename Sp, typename G_F, typename>
	auto swap_grid(OpOperatorDerivative<O, V, Sp> const& e, G_F&& g)
	{
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * expr::make_operator_derivative<O>(swap_grid_solver<Sg>(e.solver, std::forward<G_F>(g)));
	}

	template<typename Sg, Axis ax, size_t O, typename V, typename Sp, typename G_F, typename>
	auto swap_grid(OpOperatorDirectionalDerivative<ax, O, V, Sp> const& e, G_F&& g)
	{
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * expr::make_operator_directional_derivative<ax, O>(swap_grid_solver<Sg>(e.solver, std::forward<G_F>(g)));
	}

	template<typename Sg, size_t... Os, typename V, typename Sp, typename G_F, typename>
	auto swap_grid(OpOperatorMixedDerivative<V, Sp, Os...> const& e, G_F&& g)
	{
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * expr::make_operator_mixed_derivative<Os...>(swap_grid_solver<Sg>(e.solver, std::forward<G_F>(g)));
	}

	template<typename Sg, typename A1, typename A2, typename G_F, typename>
	auto swap_grid(OpOperatorChain<A1, A2> const& e, G_F&& g)
	{
		return OpOperatorChain(swap_grid<Sg>(e.f, std::forward<G_F>(g)), swap_grid<Sg>(e.g, std::forward<G_F>(g)));
	}

	template<typename Sg, typename A1, typename A2, typename E, typename G_F, typename>
	auto swap_grid(OpChain<A1, A2, E> const& e, G_F&& g)
	{
		return swap_grid<Sg>(e.combination, std::forward<G_F>(g))(swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
	}

	template<typename Sg, typename A1, typename E, typename G_F, typename>
	auto swap_grid(OpChain<A1, OpIdentity, E> const& e, G_F&& g)
	{
		return OpChain(
			OpOperatorChain(swap_grid<Sg>(e.combination.f, std::forward<G_F>(g)), OpIdentity{}),
			swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
	}

	template<typename Sg, typename A1, typename A2, typename G_F, typename>
	auto swap_grid(OpOperatorCombination<A1, A2> const& e, G_F&& g)
	{
		return OpOperatorCombination(swap_grid<Sg>(e.f, std::forward<G_F>(g)), swap_grid<Sg>(e.g, std::forward<G_F>(g)));
	}

	template<typename Sg, typename A1, typename A2, typename E, typename G_F, typename>
	auto swap_grid(OpCombination<A1, A2, E> const& e, G_F&& g)
	{
		return swap_grid<Sg>(e.combination, std::forward<G_F>(g))(swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
	}

	template<typename Sg, typename Dd, typename V, typename E, typename Sp, typename G_F, typename>
	auto swap_grid(OpDerivative<Dd, V, E, Sp> const& e, G_F&& g)
	{
		constexpr size_t order = OpDerivative<Dd, V, E, Sp>::order;
		constexpr Axis axis = OpDerivative<Dd, V, E, Sp>::axis;
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * expr::make_derivative<Dd>(
			swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), swap_grid_solver<Sg>(e.solver, std::forward<G_F>(g)));
	}

	template<typename Sg, size_t O, typename V, typename E, typename GG, typename G_F, typename>
	auto swap_grid(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<GG>> const& e, G_F&& g)
	{
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * expr::make_derivative<O, GG>(swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), swap_grid_solver<Sg>(e.solver, std::forward<G_F>(g)));
	}

	template<typename Sg, typename V, typename E, typename T, typename G_F, typename>
	auto swap_grid(OpIntegral<V, E, T> const& e, G_F&& g)
	{
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return expr::make_integral(c, swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)), e.domain);
	}

	template<typename Sg, typename G, typename V, typename E, typename G_F, typename>
	auto swap_grid(OpMap<G, V, E> const& e, G_F&& g)
	{
		auto eg = swap_grid<Sg>(e.e, std::forward<G_F>(g));
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * OpMap<G, OpIdentity, decltype(eg)>(OpIdentity{}, eg);
	}

	template<typename Sg, expr::exp_key_t X, typename V, typename E, typename G_F, typename>
	auto swap_grid(OpPow<X, V, E> const& e, G_F&& g)
	{
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * expr::make_pow<X>(swap_grid<Sg>(expr::get_enclosed_expression(e), std::forward<G_F>(g)));
	}

	template<typename Sg, typename V, typename E, typename F, typename Arg0, typename... Args, typename G_F, typename>
	auto swap_grid(OpFunction<V, E, F, Arg0, Args...> const& e, G_F&& g)
	{
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * ::OpFunction(e.name, OpIdentity{}, swap_grid<Sg>(e.e, std::forward<G_F>(g)), e.f, e.tt);
	}

	template<typename Sg, auto f, typename V, typename E, typename G_F, typename>
	auto swap_grid(OpFunctionApply<f, V, E> const& e, G_F&& g)
	{
		auto&& eg = swap_grid<Sg>(e.e, std::forward<G_F>(g));
		auto&& c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * expr::make_function<f>(OpIdentity{}, eg);
	}


	namespace
	{

		template<typename Sg, typename E, size_t... ArgNs, typename... Ts, size_t... ArgMs, size_t... Is, typename G_F>
		auto swap_grid_symbolic(
			SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f, 
			SymbolicTemplate<E, ArgNs...> const& tmpl, std::index_sequence<ArgMs...>,
			std::index_sequence<Is...>, G_F&& g)
		{
			auto swapped = swap_grid<Sg>(tmpl(swap_grid<Sg>(std::get<Is>(f.data), std::forward<G_F>(g))...), std::forward<G_F>(g));
			return (expr::function_of(as_variable<ArgMs>(swap_grid<Sg>(std::get<Is>(f.data), std::forward<G_F>(g)))...) = swapped);
		}

		template<typename Sg, typename E, size_t... ArgNs, typename... Ts, size_t... ArgMs, size_t... Is, typename G_F>
		auto swap_grid_symbolic(
			SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f, 
			OpExpression<E> const& e, std::index_sequence<ArgMs...>,
			std::index_sequence<Is...>, G_F&& g)
		{
			auto swapped = swap_grid<Sg>(*static_cast<E const*>(&e), std::forward<G_F>(g));
			return (expr::function_of(as_variable<ArgMs>(swap_grid<Sg>(std::get<Is>(f.data), std::forward<G_F>(g)))...) = swapped);
		}

		template<typename Sg, typename E, size_t... ArgNs, typename... Ts, typename G_F>
		auto swap_grid_symbolic(SymbolicFunction<E, Variable<ArgNs, Ts>...> const& f, G_F&& g)
		{
			using var_g = decltype(expr::get_independent_variables(std::forward<G_F>(g)));
			using seq_t = std::make_index_sequence<sizeof...(Ts) + var_g::size()>;
			using seq_filt_t = symphas::lib::filter_seq_t<seq_t, var_g>;
			using seq_cut_t = symphas::lib::seq_lt_t<sizeof...(Ts), seq_filt_t>;
			return swap_grid_symbolic<Sg>(
				f, expr::template_of(Variable<ArgNs>{}...) = f.e, 
				seq_cut_t{}, std::make_index_sequence<sizeof...(Ts)>{}, std::forward<G_F>(g));
		}

		template<typename Sg, typename T, typename G_F>
		auto swap_grid_symbolic(SymbolicData<T> const& e, G_F&& g)
		{
			return e;
		}

		template<typename Sg, typename T, typename G_F>
		auto swap_grid_symbolic(SymbolicDataArray<T> const& e, G_F&& g)
		{
			return e;
		}

		template<typename Sg, typename... Ts, typename G_F>
		auto swap_grid_symbolic(SymbolicDataArray<std::tuple<Term<Ts>...>> const& e, G_F&& g)
		{
			using mask_seq = std::integer_sequence<bool, (expr::factor_count<OpTerm<OpIdentity, Ts>, Sg>::value > 0)...>;
			return SymbolicDataArray(swap_symbolic_data_array(e, mask_seq{}, std::make_index_sequence<sizeof...(Ts)>{}, std::forward<G_F>(g)));
		}

		template<typename Sg, typename... Ts, size_t... Is, typename G_F>
		auto swap_grid_symbolic(Substitution<Ts...> const& e, std::index_sequence<Is...>, G_F&& g)
		{
			return Substitution(swap_grid_symbolic<Sg>(std::get<Is>(e), std::forward<G_F>(g))...);
		}

		template<typename Sg, typename... Ts, typename G_F>
		auto swap_grid_symbolic(Substitution<Ts...> const& e, G_F&& g)
		{
			return swap_grid_symbolic<Sg>(e, std::make_index_sequence<sizeof...(Ts)>{}, std::forward<G_F>(g));
		}


		template<typename Sg, expr::NoiseType nt, typename T, size_t D, typename G_F>
		auto swap_grid_symbolic(NoiseData<nt, T, D> const& data, G_F&& g)
		{
			return data;
		}

		template<typename Sg, typename G_F>
		auto swap_grid_symbolic(DynamicIndex const& data, G_F&& g)
		{
			return swap_grid<Sg>(data, std::forward<G_F>(g));
		}

		template<typename Sg, typename E0, typename K, typename G_F>
		auto swap_grid_symbolic(SymbolicListIndex<E0, K> const& data, G_F&& g)
		{
			return SymbolicListIndex{ swap_grid<Sg>(data.e, std::forward<G_F>(g)), K{} };
		}

		template<typename Sg, int N, int P, typename G_F>
		auto swap_grid_symbolic(expr::symbols::i_<N, P> const& data, G_F&& g)
		{
			return SymbolicListIndex{ swap_grid<Sg>(data, std::forward<G_F>(g)) };
		}


		template<typename Sg, typename A1, typename B1, typename A2, typename B2, typename G_F>
		auto swap_grid_symbolic(expr::series_limits<std::pair<A1, B1>, std::pair<A2, B2>> const& limit, G_F&& g)
		{
			return expr::series_limits(
				std::make_pair(
					swap_grid<Sg>(expr::limit_0(limit).first, std::forward<G_F>(g)),
					swap_grid<Sg>(expr::limit_0(limit).second, std::forward<G_F>(g))),
				std::make_pair(
					swap_grid<Sg>(expr::limit_1(limit).first, std::forward<G_F>(g)),
					swap_grid<Sg>(expr::limit_1(limit).second, std::forward<G_F>(g))));
		}

		template<typename Sg, typename A, typename B, typename T1, typename G_F>
		auto swap_grid_symbolic(expr::series_limits<T1, std::pair<A, B>> const& limit, G_F&& g)
		{
			return expr::series_limits(
				swap_grid<Sg>(expr::limit_0(limit), std::forward<G_F>(g)),
				std::make_pair(
					swap_grid<Sg>(expr::limit_1(limit).first, std::forward<G_F>(g)),
					swap_grid<Sg>(expr::limit_1(limit).second, std::forward<G_F>(g))));
		}

		template<typename Sg, typename A, typename B, typename T2, typename G_F>
		auto swap_grid_symbolic(expr::series_limits<std::pair<A, B>, T2> const& limit, G_F&& g)
		{
			return expr::series_limits(
				std::make_pair(
					swap_grid<Sg>(expr::limit_0(limit).first, std::forward<G_F>(g)),
					swap_grid<Sg>(expr::limit_0(limit).second, std::forward<G_F>(g))),
				swap_grid<Sg>(expr::limit_1(limit), std::forward<G_F>(g)));
		}

		template<typename Sg, typename T1, typename T2, typename G_F>
		auto swap_grid_symbolic(expr::series_limits<T1, T2> const& limit, G_F&& g)
		{
			return expr::series_limits(swap_grid<Sg>(expr::limit_0(limit), std::forward<G_F>(g)), swap_grid<Sg>(expr::limit_1(limit), std::forward<G_F>(g)));
		}

		template<typename Sg, typename... T1s, typename... T2s, size_t... Is, typename G_F>
		auto swap_grid_symbolic(std::tuple<expr::series_limits<T1s, T2s>...> const& limits, std::index_sequence<Is...>, G_F&& g)
		{
			return std::make_tuple(swap_grid_symbolic<Sg>(std::get<Is>(limits), std::forward<G_F>(g))...);
		}

		template<typename Sg, typename... T1s, typename... T2s, typename G_F>
		auto swap_grid_symbolic(std::tuple<expr::series_limits<T1s, T2s>...> const& limits, G_F&& g)
		{
			return swap_grid_symbolic<Sg>(limits, std::make_index_sequence<sizeof...(T1s)>{}, std::forward<G_F>(g));
		}

		template<typename Sg, typename T, typename E0, typename... T0s, typename G_F>
		auto swap_grid_symbolic(T const& data, SymbolicFunction<E0, T0s...> const& f, G_F&& g)
		{
			auto e = swap_grid_symbolic<Sg>(f, std::forward<G_F>(g));
			auto substitution = swap_grid_symbolic<Sg>(data, std::forward<G_F>(g));

			return symphas::internal::make_symbolic_eval(OpIdentity{}, substitution, e);
		}

		template<typename Sg, typename Op, typename... Ts, typename E, typename E0, typename... T0s, typename G_F>
		auto swap_grid_symbolic(SymbolicSeries<Op, Substitution<Ts...>, E> const& series, SymbolicFunction<E0, T0s...> const& f, G_F&& g)
		{
			auto limits = swap_grid_symbolic<Sg>(series.limits, std::forward<G_F>(g));
			auto e = swap_grid<Sg>(series.e, std::forward<G_F>(g));
			auto substitution = swap_grid_symbolic<Sg>(series.substitution, std::forward<G_F>(g));
			
			return expr::recreate_series(e, limits, series, substitution);
		}
	}

	template<typename Sg, typename V, int N, int P, typename E, typename... Ts, typename G_F, typename>
	auto swap_grid(OpSymbolicEval<V, expr::symbols::i_<N, P>, SymbolicFunction<E, Ts...>> const& e, G_F&& g)
	{
		if constexpr (expr::factor_count<Sg, expr::symbols::i_<N, P>>::value > 0)
		{
			auto substitution = swap_grid_symbolic<Sg>(e.data, std::forward<G_F>(g));
			auto ev = symphas::internal::make_symbolic_eval(OpIdentity{}, SymbolicListIndex{ substitution, expr::symbols::i_<N, P>{} }, e.f);

			auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
			return c * ev;
		}
		else
		{
			auto swapped = swap_grid_symbolic<Sg>(e.data, e.f, std::forward<G_F>(g));
			auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
			return c * swapped;
		}
	}

	template<typename Sg, typename V, typename E0, typename K, typename E, typename... Ts, typename G_F, typename>
	auto swap_grid(OpSymbolicEval<V, SymbolicListIndex<E0, K>, SymbolicFunction<E, Ts...>> const& e, G_F&& g)
	{
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		if constexpr (expr::factor_count<Sg, K>::value > 0 || expr::factor_count<Sg, SymbolicListIndex<E0, K>>::value > 0)
		{
			auto substitution = swap_grid_symbolic<Sg>(e.data, std::forward<G_F>(g));
			auto ev = symphas::internal::make_symbolic_eval(OpIdentity{}, substitution, e.f);
			return c * ev;
		}
		else
		{
			auto swapped = swap_grid_symbolic<Sg>(e.data, e.f, std::forward<G_F>(g));
			return c * swapped;
		}
	}

	template<typename Sg, typename V, typename sub_t, typename E, typename... Ts, typename G_F, typename>
	auto swap_grid(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e, G_F&& g)
	{
		auto swapped = swap_grid_symbolic<Sg>(e.data, e.f, std::forward<G_F>(g));
		auto c = swap_grid<Sg>(expr::coeff(e), std::forward<G_F>(g));
		return c * swapped;
	}

	template<typename... Sgs, typename C, typename E, typename... G_Fs>
	auto swap_grid(SymbolicCaseSwap<C>, OpExpression<E> const& e, G_Fs&&... gs)
	{
		return swap_grid<SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>>
			(*static_cast<E const*>(&e), std::make_tuple(std::forward<G_Fs>(gs)...));
	}

	template<typename... Sgs, typename C, typename E, typename... G_Fs>
	auto swap_grid(SymbolicCaseSwap<C>, OpOperator<E> const& e, G_Fs&&... gs)
	{
		return swap_grid<SymbolicCaseSwap<C, symphas::lib::types_list<Sgs...>>>
			(*static_cast<E const*>(&e), std::make_tuple(std::forward<G_Fs>(gs)...));
	}

	//template<typename... Sgs, typename E, typename... G_Fs>
	//auto swap_grid(std::tuple<>, OpExpression<E> const& e, G_Fs&&... gs)
	//{
	//	return *static_cast<E const*>(&e);
	//}

	//template<typename Sg, typename... Sgs, typename C0, typename... Cs, typename E, typename G_F0, typename... G_Fs>
	//auto swap_grid(std::tuple<SymbolicCaseSwap<C0>, SymbolicCaseSwap<Cs>...>, OpExpression<E> const& e, G_F0&& g0, G_Fs&&... gs)
	//{
	//	auto c0 = swap_grid<Sg>(SymbolicCaseSwap<C0>{}, * static_cast<E const*>(&e), std::forward<G_F0>(g0));
	//	return swap_grid<Sgs...>(std::tuple<SymbolicCaseSwap<Cs>...>{}, * static_cast<E const*>(&e), std::forward<G_Fs>(gs)...);
	//}

	namespace
	{
		/*template<typename T, typename I, typename E>
		auto handle_coeff_swap(OpCoeff<T, I> const& coeff, E const& swap)
		{
			return coeff[expr::eval(swap)];
		}*/

		template<typename T, size_t N0>
		auto handle_coeff_swap(OpCoeff<T, expr::symbols::placeholder_N_symbol_<N0>> const& coeff, DynamicIndex const& index)
		{
			return coeff(index);
		}

		template<typename T, int N, int P>
		auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, P>> const& coeff, int index)
		{
			return OpCoeff(coeff.data.data + P, coeff.data.len - P)[index];
		}

		template<typename T, int N>
		auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, 0>> const& coeff, int index)
		{
			return coeff[index];
		}

		template<typename T, int N, int P, size_t N0>
		auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, P>> const& coeff, expr::symbols::placeholder_N_symbol_<N0>)
		{
			return OpCoeff(coeff.data.data + P, coeff.data.len - P)(expr::symbols::placeholder_N_symbol_<N0>{});
		}

		template<typename T, int N, size_t N0>
		auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, 0>> const& coeff, expr::symbols::placeholder_N_symbol_<N0>)
		{
			return coeff(expr::symbols::placeholder_N_symbol_<N0>{});
		}

		template<typename T, int N, int P, size_t N0>
		auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, P>> const& coeff, expr::symbols::placeholder_N_<N0>)
		{
			return handle_coeff_swap(coeff, expr::symbols::placeholder_N_symbol_<N0>{});
		}

		template<typename T, int N, int P>
		auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, P>> const& coeff, DynamicIndex const& index)
		{
			return OpCoeff(coeff.data.data + P, coeff.data.len - P)(index);
		}

		template<typename T, int N>
		auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, 0>> const& coeff, DynamicIndex const& index)
		{
			return coeff(index);
		}

		template<typename T>
		auto handle_coeff_swap(OpCoeff<T, DynamicIndex> const& coeff, OpCoeffSwap<DynamicIndex> const&)
		{
			return const_cast<OpCoeff<T, DynamicIndex>&>(coeff).fix();
		}

		template<typename T>
		auto handle_coeff_swap(OpCoeff<T, DynamicIndex> const& coeff, DynamicIndexSet const& set)
		{
			return const_cast<OpCoeff<T, DynamicIndex>&>(coeff).fix(set);
		}

		template<typename T, typename I>
		auto handle_coeff_swap(OpCoeff<T, I> const& coeff, DynamicIndexSet const& set)
		{
			return coeff;
		}

		//template<typename T, int N, size_t N0>
		//auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, 0>> const& coeff, expr::symbols::placeholder_N_symbol_<N0> const&)
		//{
		//	return coeff(expr::symbols::placeholder_N_symbol_<N0>{});
		//}

		//template<typename T, int N, int P, size_t N0>
		//auto handle_coeff_swap(OpCoeff<T, expr::symbols::i_<N, P>> const& coeff, expr::symbols::placeholder_N_symbol_<N0> const&)
		//{
		//	return OpCoeff(coeff.data.data + P, coeff.data.len - P)(expr::symbols::placeholder_N_symbol_<N0>{});
		//}
	}

	template<typename Sg, typename G_F>
	auto swap_grid(DynamicIndex const& data, G_F&& g)
	{
		if constexpr (std::is_same<Sg, DynamicIndexSet>::value)
		{
			const_cast<DynamicIndex&>(data).fix(std::forward<G_F>(g));
		}
		return data;
	}

	template<typename Sg, typename T, typename I, typename G_F>
	auto swap_grid(OpCoeff<T, I> const& coeff, G_F&& g)
	{
		if constexpr (expr::factor_count<Sg, OpCoeff<T, I>>::value > 0 || std::is_same<Sg, DynamicIndexSet>::value)
		{
			return handle_coeff_swap(coeff, std::forward<G_F>(g));
		}
		else
		{
			return coeff;
		}
	}

	template<typename Sg, typename T, typename I, size_t... Ns, typename G_F>
	auto swap_grid(OpTensor<OpCoeff<T, I>, Ns...> const& coeff, G_F&& g)
	{
		if constexpr (expr::factor_count<Sg, OpCoeff<T, I>>::value > 0 || std::is_same<Sg, DynamicIndexSet>::value)
		{
			return expr::make_tensor<Ns...>(handle_coeff_swap(static_cast<OpCoeff<T, I> const&>(coeff), std::forward<G_F>(g)));
		}
		else
		{
			return coeff;
		}
	}



	template<typename Sg0, typename Sg1, typename... Sgs, typename E, typename G_F>
	auto swap_grid(OpEvaluable<E> const& e, G_F&& g)
	{
		return swap_grid<Sg1, Sgs...>(swap_grid<Sg0>(*static_cast<E const*>(&e), std::forward<G_F>(g)), std::forward<G_F>(g));
	}

	template<typename Sg0, typename Sg1, typename... Sgs, typename E, typename G_F0, typename G_F1, typename... G_Fs>
	auto swap_grid(OpEvaluable<E> const& e, G_F0&& g0, G_F1&& g1, G_Fs&& ...gs)
	{
		return swap_grid<Sg1, Sgs...>(swap_grid<Sg0>(*static_cast<E const*>(&e), std::forward<G_F0>(g0)), std::forward<G_F1>(g1), std::forward<G_Fs>(gs)...);
	}

	template<Axis ax, size_t Z, typename E, typename G_F>
	auto swap_grid(OpEvaluable<E> const& e, G_F&& g)
	{
		return swap_grid<Variable<Z, VectorComponent<ax>>>(*static_cast<E const*>(&e), std::forward<G_F>(g));
	}

	template<Axis ax, typename Sg, typename E, typename G_F>
	auto swap_grid(OpEvaluable<E> const& e, G_F&& g)
	{
		return swap_grid<VectorComponent<ax, Sg>>(*static_cast<E const*>(&e), std::forward<G_F>(g));
	}



	template<typename... Sgs, typename E, typename G_F>
	auto swap_grid(symphas::lib::types_list<Sgs...>, E&& e, G_F&& g)
	{
		return swap_grid<Sgs...>(std::forward<E>(e), std::forward<G_F>(g));
	}



	template<typename Eg, typename E, typename E_F>
	auto swap_expression(OpEvaluable<E> const& e, E_F&& g)
	{
		return swap_grid<expr::matches_with<Eg>>(*static_cast<E const*>(&e), std::forward<E_F>(g));
	}

}

namespace expr
{

	template<typename E>
	auto fix_index(OpExpression<E>& e, DynamicIndexSet set)
	{
		expr::transform::swap_grid<DynamicIndexSet>(*static_cast<E const*>(&e), set);
	}

	template<typename E>
	auto fix_index(OpOperator<E>& e, DynamicIndexSet set)
	{
		expr::transform::swap_grid<DynamicIndexSet>(*static_cast<E const*>(&e), set);
	}

	template<typename E>
	auto fix_coeff(OpExpression<E>& e)
	{
		expr::transform::swap_grid<OpCoeffSwap<DynamicIndex>>(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto fix_coeff(OpOperator<E>& e)
	{
		expr::transform::swap_grid<OpCoeffSwap<DynamicIndex>>(*static_cast<E const*>(&e));
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

		template<typename E>
		auto pack_left(OpOperator<E> const& e)
		{
			return std::make_pair(*static_cast<const E*>(&e), OpVoid{});
		}

		template<typename E>
		auto pack_right(OpOperator<E> const& e)
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
			return std::make_pair(pair0.first, (pair0.second + ... + pairs.second));
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
	template<typename E>
	auto by_linear(OpExpression<E> const& e)
	{
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
	template<typename E>
	auto by_linear(OpOperator<E> const& e)
	{
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

	//! Split all linear variables from nonlinear variables.
	/*!
	 * Using the available predicates to determine linearity/nonlinearity of
	 * subexpressions, a primary expression `E` is expanded into `L` and `NL`, where
	 * `E = L + NL`. The `L` part is totally linear, consistent with the predicate
	 * for linearity. Similarly, the `NL` part is nonlinear.
	 *
	 * \param e The expression which is split.
	 */
	template<typename G, typename V, typename E>
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
	template<typename V, typename sub_t, typename E, typename... Ts>
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
	template<typename V, typename E, typename... Ts, int... I0s, int... P0s,
		typename A, typename... As, typename B, typename C>
	auto by_linear(OpSum<V, E,
		Substitution<SymbolicDataArray<Ts>...>,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		SymbolicFunction<A, As...>, B, C> const& series);



	// derivative
	namespace
	{

		template<typename Dd, typename V, typename E, typename T,
			typename std::enable_if_t<(expression_predicate<E>::linear), int> = 0>
		auto by_linear_derivative(OpDerivative<Dd, V, E, T> const& e)
		{
			return pack_left(e);
		}

		template<typename Dd, typename V, typename E, typename T,
			typename std::enable_if_t<(!expression_predicate<E>::linear && expression_predicate<E>::combination), int> = 0>
		auto by_linear_derivative(OpDerivative<Dd, V, E, T> const& e)
		{
			return expr::split::by_linear(expr::apply_operators(e));
		}

		template<typename Dd, typename V, typename E, typename T,
			typename std::enable_if_t<(!expression_predicate<E>::linear && !expression_predicate<E>::combination), int> = 0>
		auto by_linear_derivative(OpDerivative<Dd, V, E, T> const& e)
		{
			return pack_right(e);
		}
	}

	template<typename Dd, typename V, typename E, typename Sp>
	auto by_linear(OpDerivative<Dd, V, E, Sp> const& e)
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
			if constexpr (expr::is_linear<A1>::value && expr::is_linear<A2>::value)
			{
				auto [l, r] = expr::split::by_linear(expr::get_enclosed_expression(e));
				return std::make_pair(e.combination * l, e.combination * r);
			}
			else
			{
				return pack_right(e);
			}
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
			if constexpr (expr::is_linear<A1>::value && expr::is_linear<A2>::value)
			{
				auto [l, r] = expr::split::by_linear(expr::get_enclosed_expression(e));
				return std::make_pair(e.combination * l, e.combination * r);
			}
			else
			{
				return pack_right(e);
			}
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
			typename std::enable_if_t<expression_predicate<OpConvolution<V, E1, E2>>::linear, int> = 0>
		auto by_linear_convolution(OpConvolution<V, E1, E2> const& e)
		{
			return pack_left(e);
		}

		template<typename V, typename E1, typename E2,
			typename std::enable_if_t<!expression_predicate<OpConvolution<V, E1, E2>>::linear, int> = 0>
		auto by_linear_convolution(OpConvolution<V, E1, E2> const& e)
		{
			return pack_right(e);
		}

	}

	template<typename V, typename E1, typename E2>
	auto by_linear(OpConvolution<V, E1, E2> const& e)
	{
		return by_linear_convolution(e);
	}

	// handling general expressions (recursion termination) and binary operations (recursive)

	template<typename V, typename G>
	auto by_linear(OpTerm<V, G> const& e)
	{
		return pack_left(e);
	}


	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	auto by_linear(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		using index_list_t = symphas::internal::select_unique_i_<expr::op_types_t<OpTerms<V, Term<Gs, Xs>...>>>;
		using term_list_t = symphas::lib::filter_types<symphas::lib::types_list<Gs...>, index_list_t>;

		constexpr int N = sum_factors<OpTerms<V, Term<Gs, Xs>...>, term_list_t>::value;
		if constexpr (N == 1)
		{
			return pack_left(e);
		}
		else
		{
			return pack_right(e);
		}
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

		template<typename... E0s, typename... Ts>
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

		template<typename... E0s, typename... Ts>
		auto filter_linear(std::tuple<E0s...> const& terms0, Ts const& ...ts)
		{
			return filter_linear(terms0, std::make_index_sequence<sizeof...(E0s)>{}, ts...);
		}

		template<typename E>
		auto get_linear_term(OpExpression<E> const& e)
		{
			auto&& [l, _] = expr::split::by_linear(*static_cast<E const*>(&e));
			return std::make_tuple(l);
		}

		template<typename E>
		auto get_linear_term(OpOperator<E> const& e)
		{
			auto&& [l, _] = expr::split::by_linear(*static_cast<E const*>(&e));
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
		template<typename E, typename A>
		auto check_term_in_num(OpOperator<E> const& e, OpExpression<A> const& a, std::index_sequence<>)
		{
			return std::make_tuple(OpVoid{});
		}

		// \p a are the terms from the numerator, and \p e is the term from the denominator which is being
		// checked that it can divide at least one of them.
		template<typename E, typename A>
		auto check_term_in_num(OpExpression<E> const& e, OpOperator<A> const& a, std::index_sequence<>)
		{
			return std::make_tuple(OpVoid{});
		}

		// \p a are the terms from the numerator, and \p e is the term from the denominator which is being
		// checked that it can divide at least one of them.
		template<typename E, typename A>
		auto check_term_in_num(OpOperator<E> const& e, OpOperator<A> const& a, std::index_sequence<>)
		{
			return std::make_tuple(OpVoid{});
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
		template<typename E, typename... As, size_t... Is>
		auto check_term_in_num(OpOperator<E> const& e, OpAdd<As...> const& a, std::index_sequence<Is...>)
		{
			return std::make_tuple(OpVoid{});
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

		template<typename T0, typename... As, typename B>
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
			return (std::get<Is>(linear_terms) + ... + OpVoid{});
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



	template<typename G, typename V, typename E>
	auto by_linear(OpMap<G, V, E> const& e)
	{
		if constexpr (expr::is_combination<OpMap<G, V, E>>::value)
		{
			auto [l, r] = by_linear(expr::get_enclosed_expression(e));
			auto lm = expr::make_map<G>(expr::coeff(e) * l);
			auto rm = expr::make_map<G>(expr::coeff(e) * r);
			return std::make_pair(lm, rm);
		}
		else
		{
			return pack_right(e);
		}
	}


	template<typename V, typename sub_t, typename E, typename... Ts>
	auto by_linear(OpSymbolicEval<V, sub_t, SymbolicFunction<E, Ts...>> const& e)
	{
		auto&& [l, nl] = by_linear(e.f.e);
		auto fl = (expr::function_of(Ts{}...) = l);
		auto fnl = (expr::function_of(Ts{}...) = nl);
		return std::make_pair(
			symphas::internal::make_symbolic_eval(expr::coeff(e), e.data, fl), 
			symphas::internal::make_symbolic_eval(expr::coeff(e), e.data, fnl));
	}

	template<typename V, typename E, typename... Ts, int... I0s, int... P0s, 
		typename A, typename... As, typename B, typename C>
	auto by_linear(OpSum<V, E, 
		Substitution<SymbolicDataArray<Ts>...>, 
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
		SymbolicFunction<A, As...>, B, C> const& series)
	{
		auto&& [l, nl] = by_linear(series.data.e);
		return std::make_pair(
			expr::coeff(series) * expr::recreate_series(l, series.data),
			expr::coeff(series) * expr::recreate_series(nl, series.data));
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
	template<size_t Z, typename Dd, typename V, typename G, typename Sp>
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
	template<size_t Z, typename V, typename E, typename T>
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
	template<size_t Z, typename V, size_t D, typename E>
	auto separate_var(OpConvolution<V, GaussianSmoothing<D>, E> const& e);

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
	template<size_t Z, typename G, typename V, typename E>
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
	template<size_t Z, auto f, typename V, typename E>
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
	template<size_t Z, typename V, typename E, typename F, typename Arg0, typename... Args>
	auto separate_var(OpFunction<V, E, F, Arg0, Args...> const& e);


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
		constexpr bool svm_pred = svc_pred<Z, E1, E2>;
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
			auto separate_var_convolution(OpConvolution<V, E1, E2> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, typename V, typename E1, typename E2,
			typename std::enable_if_t<!svc_pred<Z, E1, E2>, int> = 0>
			auto separate_var_convolution(OpConvolution<V, E1, E2> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename V, typename E1, typename E2>
	auto separate_var(OpConvolution<V, E1, E2> const& e)
	{
		return separate_var_convolution<Z>(e);
	}


	namespace
	{
		template<size_t Z, typename V, size_t D, typename E,
			typename std::enable_if_t<svcg_pred<Z, E>, int> = 0>
			auto separate_var_convolution_g(OpConvolution<V, GaussianSmoothing<D>, E> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, typename V, size_t D, typename E,
			typename std::enable_if_t<!svcg_pred<Z, E>, int> = 0>
			auto separate_var_convolution_g(OpConvolution<V, GaussianSmoothing<D>, E> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename V, size_t D, typename E>
	auto separate_var(OpConvolution<V, GaussianSmoothing<D>, E> const& e)
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
		auto separate_var_derivative(OpDerivative<Dd, V, E, T> const& e)
		{
			return pack_left(e);
		}


		template<size_t Z, typename Dd, typename V, typename E, typename T,
			typename std::enable_if_t<(!svd_pred_1<Z, E>&& svd_pred_2<Z, E>), int> = 0>
		auto separate_var_derivative(OpDerivative<Dd, V, E, T> const& e)
		{
			return separate_var<Z>(expr::apply_operators(e));

		}

		template<size_t Z, typename Dd, typename V, typename E, typename T,
			typename std::enable_if_t<!(svd_pred_1<Z, E> || svd_pred_2<Z, E>), int> = 0>
		auto separate_var_derivative(OpDerivative<Dd, V, E, T> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename Dd, typename V, typename E, typename Sp>
	auto separate_var(OpDerivative<Dd, V, E, Sp> const& e)
	{
		return separate_var_derivative<Z>(e);
	}



	namespace
	{
		template<size_t Z, typename Dd, typename V, typename G, typename T,
			typename std::enable_if_t<svd_pred_1<Z, G>, int> = 0>
		auto separate_var_derivative_lop(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, T> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, typename Dd, typename V, typename G, typename T,
			typename std::enable_if_t<!svd_pred_1<Z, G>, int> = 0>
		auto separate_var_derivative_lop(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, T> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename Dd, typename V, typename G, typename Sp>
	auto separate_var(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const& e)
	{
		return separate_var_derivative_lop<Z>(e);
	}


	namespace
	{
		template<size_t Z, typename V, typename E, typename T,
			typename std::enable_if_t<svd_pred_1<Z, E>, int> = 0>
		auto separate_var_integral(OpIntegral<V, E, T> const& e)
		{
			return pack_left(e);
		}

		template<size_t Z, typename V, typename E, typename T,
			typename std::enable_if_t<!svd_pred_1<Z, E>, int> = 0>
		auto separate_var_integral(OpIntegral<V, E, T> const& e)
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename V, typename E, typename T>
	auto separate_var(OpIntegral<V, E, T> const& e)
	{
		return separate_var_integral<Z>(e);
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

	template<size_t Z, typename G, typename V, typename E>
	auto separate_var(OpMap<G, V, E> const& e)
	{
		if constexpr (expr::is_combination<OpMap<G, V, E>>::value)
		{
			auto [l, r] = separate_var<Z>(expr::get_enclosed_expression(e));
			auto lm = expr::make_map<G>(expr::coeff(e) * l);
			auto rm = expr::make_map<G>(expr::coeff(e) * r);
			return std::make_pair(lm, rm);
		}
		else
		{
			if constexpr (svcg_pred<Z, E>)
			{
				return pack_left(e);
			}
			else
			{
				return pack_right(e);
			}
		}
	}

	template<size_t Z, expr::exp_key_t X, typename V, typename E>
	auto separate_var(OpPow<X, V, E> const& e)
	{
		if constexpr (svcg_pred<Z, E>)
		{
			return pack_left(e);
		}
		else
		{
			return pack_right(e);
		}
	}

	template<size_t Z, auto f, typename V, typename E>
	auto separate_var(OpFunctionApply<f, V, E> const& e)
	{
		if constexpr (svcg_pred<Z, E>)
		{
			return pack_left(e);
		}
		else
		{
			return pack_right(e);
		}
	}

	template<size_t Z, typename V, typename E, typename F, typename Arg0, typename... Args>
	auto separate_var(OpFunction<V, E, F, Arg0, Args...> const& e)
	{
		if constexpr (svcg_pred<Z, E>)
		{
			return pack_left(e);
		}
		else
		{
			return pack_right(e);
		}
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
	//	size_t O = OpDerivative<Dd, V, E, Sp>::order,
	//	typename std::enable_if<(O == OD), int>::type = 0>
	//auto factor_deriv(OpDerivative<Dd, V, E, Sp> const& e)
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
	//	size_t O = OpDerivative<Dd, V, E, Sp>::order,
	//	typename std::enable_if_t<(O != OD), int>::type = 0>
	//auto factor_deriv(OpDerivative<Dd, V, E, Sp> const& e)
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

	template<typename E>
	auto separate_operator(OpOperator<E> const& e)
	{
		return pack_left(*static_cast<E const*>(&e));
	}

	namespace
	{
		template<typename Sp, size_t... Os, Axis... axs>
		auto separate_mixed(Sp const& solver, typename Sp::template mixed_derivative<Os...>, symphas::internal::axis_list<axs...>)
		{
			return (OpIdentity{} * ... * expr::make_operator_directional_derivative<axs, Os>(solver));
		}

		template<typename Sp, size_t... Os>
		auto separate_mixed(Sp const& solver, typename Sp::template mixed_derivative<Os...>)
		{
			return separate_mixed<Sp>(solver, typename Sp::template mixed_derivative<Os...>{}, symphas::lib::make_axis_list<sizeof...(Os)>());
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
	template<typename Dd, typename V, typename E, typename Sp, size_t R = expr::eval_type<E>::rank>
	auto separate_operator(OpDerivative<Dd, V, E, Sp> const& e)
	{
		constexpr size_t order = OpDerivative<Dd, V, E, Sp>::order;
		constexpr Axis axis = OpDerivative<Dd, V, E, Sp>::axis;
		
		
		if constexpr (Dd::is_directional)
		{
			if constexpr (Dd::is_mixed)
			{
				auto op = separate_mixed<Sp>(e.solver, Dd{});
				return std::make_pair(op, expr::coeff(e) * expr::get_enclosed_expression(e));
			}
			else
			{
				auto op = expr::make_operator_directional_derivative<Dd::axis, Dd::order>(e.solver);
				return std::make_pair(op, expr::coeff(e) * expr::get_enclosed_expression(e));
			}
		}
		else
		{
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
				return std::make_pair(OpTerm<OpIdentity, G>(e.data), expr::coeff(e));
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
				return std::make_pair(OpTerm<OpIdentity, Variable<Z, G>>(e.data), expr::coeff(e));
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
		
		template<size_t N, typename C, typename E>
		auto _factor(OpOperator<E> const& e)
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

		template<size_t N, typename C, typename V, typename... Gs, exp_key_t... Xs, size_t... Is>
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
			auto factor_term = (Term(factor_data) * ~(Term(factor_data).template pow<N0>())).template pow<N>();
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
			constexpr size_t N = expr::factor_count<C0, E>::value;
			if constexpr (N > 0)
			{
				return _factor<N, C0>(*static_cast<E const*>(&e));
			}
			else
			{
				return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
			}
		}

		template<typename C0, typename E,
			typename std::enable_if<(!expr::grid_can_combine<C0>::value && !is_variable_data_factor<C0>::value), int>::type = 0>
		auto _factor(OpExpression<E> const& e)
		{
			return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
		}

		template<typename C0, typename E>
		auto _factor(OpOperator<E> const& e)
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
	template<typename C0, typename... Cs, typename E>
	auto factor(OpOperator<E> const& e)
	{
		return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
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
			return _factor<min_order, C0>(*static_cast<const E*>(&e));
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
	template<size_t N, typename C0, typename C1, typename... Cs, typename E>
	auto factor(OpOperator<E> const& e)
	{
		return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
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
	template<size_t N, size_t Z0, size_t... Zs, typename E>
	auto factor(OpOperator<E> const& e)
	{
		return std::make_pair(OpIdentity{}, *static_cast<E const*>(&e));
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
		auto b = factor_tuple<Cs...>(*static_cast<const E*>(&e));
		return std::tuple_cat(a, b);
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
	template<typename C0, typename... Cs, typename E>
	auto factor_tuple(OpOperator<E> const& e)
	{
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
	template<typename C0, typename... Cs, typename E>
	auto factor_list(OpOperator<E> const& e)
	{
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
	 * That is: Given a list of factors and an expression *E*, this function returns *A*, *B*:
	 * - *A* is a product of the factors that have been taken from *E*.
	 * - *B* is the factored expression, such that *E* = *A* * *B*.
	 *
	 * \param e The expression to factor.
	 *
	 * \tparam Z0 The first variable index to factor.
	 * \tparam Cs... The remaining variable indices to factor.
	 */
	template<size_t Z0, size_t... Zs, typename E>
	auto factor_list(OpOperator<E> const& e)
	{
		return factor_list<Variable<Z0>, Variable<Zs>...>(*static_cast<const E*>(&e));
	}



}


namespace expr::split
{
	template<typename... condition_ts, typename E>
	auto separate_by(OpExpression<E> const& e)
	{
		using namespace symphas::internal;
		if constexpr ((expression_satisfies_condition<E, condition_ts> || ... || false))
		{
			return pack_left(*static_cast<E const*>(&e));
		}
		else
		{

			return pack_right(*static_cast<E const*>(&e));
		}
	}

	template<typename... condition_ts, typename... Es, size_t... Is>
	auto separate_by(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return adds_expand_pair(separate_by<condition_ts...>(expr::get<Is>(e))...);
	}

	template<typename... condition_ts, typename... Es>
	auto separate_by(OpAdd<Es...> const& e)
	{
		return separate_by<condition_ts...>(e, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<typename... condition_ts, typename E>
	auto filter(OpExpression<E> const& e)
	{
		using namespace symphas::internal;
		if constexpr ((expression_satisfies_condition<E, condition_ts> || ... || false))
		{
			return *static_cast<E const*>(&e);
		}
		else
		{

			return OpVoid{};
		}
	}

	template<typename... condition_ts, typename... Es, size_t... Is>
	auto filter(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return (expr::get<Is>(e) + ... + OpVoid{});
	}

	template<typename... condition_ts, typename... Es>
	auto filter(OpAdd<Es...> const& e)
	{
		using namespace symphas::internal;
		using mask = std::integer_sequence<bool, expression_satisfies_condition<Es, expr::or_<condition_ts...>>...>;
		using seq = std::make_index_sequence<sizeof...(Es)>;
		return filter<condition_ts...>(e, symphas::lib::filter_seq_t<seq, mask>{});
	}
}

namespace expr
{



	template<typename E>
	auto result_sum_by_term(OpEvaluable<E> const& e, grid::region_interval<0> const& interval)
	{
		return result_sum_by_term(*static_cast<E const*>(&e), 1);
	}

	template<typename E, size_t D>
	auto result_sum_by_term(OpEvaluable<E> const& e, grid::region_interval<D> const& interval)
	{
		auto group = std::reduce(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par_unseq,
#endif
			static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
			static_cast<const E*>(&e)->end(symphas::it_grp, interval));
		return group;
	}

	template<typename E, size_t D>
	auto result_sum_by_term(OpEvaluable<E> const& e, grid::region_interval_multiple<D> const& regions)
	{
		expr::eval_type_t<E> sum{};
		for (grid::region_interval<D> region : regions)
		{
			sum += result_sum_by_term(*static_cast<E const*>(&e), region);
		}
		return sum;
	}




	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename E, typename assign_type>
	void result_of_matching(OpEvaluable<E> const& e, assign_type&& data)
	{
		using namespace symphas::internal;
		if constexpr ((expression_satisfies_condition<E, condition_ts> || ... || expression_satisfies_condition<E, condition_t>))
		{
			result_accumulate(*static_cast<E const*>(&e), std::forward<assign_type>(data));
		}
	}

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename... Es, typename assign_type, size_t... Is>
	void result_of_matching(OpAdd<Es...> const& e, assign_type&& data, std::index_sequence<Is...>)
	{
		(result_of_matching(expr::get<Is>(e), std::forward<assign_type>(data)), ...);
	}

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename... Es, typename assign_type>
	void result_of_matching(OpAdd<Es...> const& e, assign_type&& data)
	{
		result_of_matching(e, std::forward<assign_type>(data), std::make_index_sequence<sizeof...(Es)>{});
	}

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename E, typename assign_type, typename region_type>
	void result_of_matching(OpEvaluable<E> const& e, assign_type&& data, region_type&& region)
	{
		using namespace symphas::internal;
		if constexpr ((expression_satisfies_condition<E, condition_ts> || ... || expression_satisfies_condition<E, condition_t>))
		{
			result_accumulate(*static_cast<E const*>(&e), std::forward<assign_type>(data), std::forward<region_type>(region));
		}
	}

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename... Es, typename assign_type, typename region_type, size_t... Is>
	void result_of_matching(OpAdd<Es...> const& e, assign_type&& data, region_type&& region, std::index_sequence<Is...>)
	{
		(result_of_matching<condition_t, condition_ts...>(expr::get<Is>(e), std::forward<assign_type>(data), std::forward<region_type>(region)), ...);
	}

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename... Es, typename assign_type, typename region_type>
	void result_of_matching(OpAdd<Es...> const& e, assign_type&& data, region_type&& region)
	{
		result_of_matching(e, std::forward<assign_type>(data), std::forward<region_type>(region), std::make_index_sequence<sizeof...(Es)>{});
	}


	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename E, typename region_type>
	auto result_sum_of_matching(OpEvaluable<E> const& e, region_type&& region)
	{
		using namespace symphas::internal;
		if constexpr ((expression_satisfies_condition<E, condition_ts> || ... || expression_satisfies_condition<E, condition_t>))
		{
			return result_sum(*static_cast<E const*>(&e), std::forward<region_type>(region));
		}
		else
		{
			return expr::eval_type_t<E>{};
		}
	}

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename... Es, typename region_type, size_t... Is>
	auto result_sum_of_matching(OpAdd<Es...> const& e, region_type&& region, std::index_sequence<Is...>)
	{
		return (result_sum_of_matching(expr::get<Is>(e), std::forward<region_type>(region)), ...);
	}

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename... Es, typename region_type>
	auto result_sum_of_matching(OpAdd<Es...> const& e, region_type&& region)
	{
		return result_sum_of_matching(e, std::forward<region_type>(region), std::make_index_sequence<sizeof...(Es)>{});
	}


	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename E>
	auto result_sum_of_matching(OpEvaluable<E> const& e)
	{
		using namespace symphas::internal;
		if constexpr ((expression_satisfies_condition<E, condition_ts> || ... || expression_satisfies_condition<E, condition_t>))
		{
			return result_sum(*static_cast<E const*>(&e));
		}
		else
		{
			return expr::eval_type_t<E>{};
		}
	}

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename... Es, size_t... Is>
	auto result_sum_of_matching(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return (result_sum_of_matching(expr::get<Is>(e)), ...);
	}

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename... Es>
	auto result_sum_of_matching(OpAdd<Es...> const& e)
	{
		return result_sum_of_matching(e, std::make_index_sequence<sizeof...(Es)>{});
	}


	template<typename condition_t = void>
	struct result_by_term_apply
	{
		template<typename E, typename assign_type, typename region_type>
		void operator()(OpEvaluable<E> const& e, assign_type&& data, region_type&& region) const
		{
			result_of_matching<condition_t>(*static_cast<E const*>(&e), std::forward<assign_type>(data), std::forward<region_type>(region));
		}

		template<typename E, typename assign_type>
		void operator()(OpEvaluable<E> const& e, assign_type&& data) const
		{
			result_of_matching<condition_t>(*static_cast<E const*>(&e), std::forward<assign_type>(data));
		}
	};


	template<>
	struct result_by_term_apply<expr::matches_series> : result_by_term_apply<void>
	{
		using result_by_term_apply<void>::operator();

		template<typename V0, typename E, typename... Ts, int... I0s, int... P0s, typename A, typename B, typename... Vs,
			typename assign_type, typename region_type>
		void operator()(OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, symphas::lib::types_list<Vs...>> const& sum,
			assign_type&& data, region_type&& region) const
		{
			result_accumulate(sum, std::forward<assign_type>(data), std::forward<region_type>(region));
		}

		template<typename V0, typename E, typename... Ts, int... I0s, int... P0s, typename A, typename B, typename... Vs, typename assign_type>
		void operator()(OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, symphas::lib::types_list<Vs...>> const& sum,
			assign_type&& data) const
		{
			TIME_THIS_CONTEXT_LIFETIME(apply_summation);
			auto coeff = expr::coeff(sum);
			for (iter_type i = 0; i < sum.data.persistent.len; ++i)
			{
				auto region = expr::iterable_domain(sum.data.persistent[i].e);
				result_accumulate(coeff * sum.data.persistent[i].e, std::forward<assign_type>(data), region);
			}
		}
	};

	template<typename condition_t = void>
	struct result_sum_by_term_apply
	{
		template<typename E, typename region_type>
		auto operator()(OpEvaluable<E> const& e, region_type&& region) const
		{
			return result_sum_of_matching<condition_t>(*static_cast<E const*>(&e), std::forward<region_type>(region));
		}

		template<typename E>
		auto operator()(OpEvaluable<E> const& e) const
		{
			return result_sum_of_matching<condition_t>(*static_cast<E const*>(&e));
		}
	};


	template<>
	struct result_sum_by_term_apply<expr::matches_series> : result_sum_by_term_apply<void>
	{
		using result_sum_by_term_apply<void>::operator();

		static const len_type group_size = 6;


		template<typename V0, typename E, typename... Ts, int... I0s, int... P0s, typename A, typename B, typename... Vs,
			typename region_type>
		auto operator()(OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, symphas::lib::types_list<Vs...>> const& sum,
			region_type&& region) const
		{
			return result_sum_by_term(sum, std::forward<region_type>(region));
		}

		template<typename V0, typename E, typename... Ts, int... I0s, int... P0s, typename A, typename B, typename... Vs>
		auto operator()(OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, symphas::lib::types_list<Vs...>> const& sum) const
		{
			if (sum.data.persistent.len > 0)
			{
				auto reduced = result_sum(expr::coeff(sum) * sum.data.persistent[0].e);
				auto coeff = expr::coeff(sum);
				for (iter_type i = 1; i < sum.data.persistent.len; ++i)
				{
					reduced += result_sum(coeff * sum.data.persistent[i].e);
				}
				return reduced;
			}
			else
			{
				using expr_type = OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
					symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, symphas::lib::types_list<Vs...>>;
				return expr::eval_type_t<expr_type>{};
			}
		}
	};

	template<>
	struct result_sum_by_term_apply<expr::matching_in_mul<expr::matches_series>> : result_sum_by_term_apply<void>
	{
		using result_sum_by_term_apply<void>::operator();

		template<typename A, typename B, typename region_type,
			std::enable_if_t<symphas::internal::expression_satisfies_condition<OpBinaryMul<A, B>, expr::matching_in_mul<expr::matches_series>>, int> = 0>
		auto operator()(OpBinaryMul<A, B> const& e, region_type&& region) const
		{
			return result_sum_by_term(sum, std::forward<region_type>(region));
		}

		template<typename A, typename B,
			std::enable_if_t<symphas::internal::expression_satisfies_condition<OpBinaryMul<A, B>, expr::matching_in_mul<expr::matches_series>>, int> = 0>
		auto operator()(OpBinaryMul<A, B> const& e) const
		{
			return apply_mul(e.a, e.b);
		}


		template<typename V0, typename E, typename... Ts, int... I0s, int... P0s, typename A, typename B, typename... Vs, typename... Es>
		auto apply_mul(OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, symphas::lib::types_list<Vs...>> const& sum,
			Es&&... terms) const
		{
			using expr_type = mul_result_t<mul_result_t<Es...>, OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, symphas::lib::types_list<Vs...>>>;

			if (sum.data.persistent.len > 0)
			{
				auto e = (expr::coeff(sum) * ... * std::forward<Es>(terms));
				auto region = expr::iterable_domain(e);
				auto reduced = expr::eval_type_t<expr_type>{};
				
				TIME_THIS_CONTEXT_LIFETIME(reduce_mul_summation);
				for (iter_type i = 0; i < sum.data.persistent.len; ++i)
				{
					auto region0 = region / expr::iterable_domain(sum.data.persistent[i].e);
					TIME_THIS_EXPRESSION_LIFETIME(iterable_domain_reduce, auto r = expr::iterable_domain(sum.data.persistent[i].e);)

					if (!grid::is_empty(region0))
					{
						reduced += result_sum(expr::make_mul(e, sum.data.persistent[i].e), region0);
					}
				}
				return reduced;
			}
			else
			{
				return expr::eval_type_t<expr_type>{};
			}
		}

		template<typename E0, typename V0, typename E, typename... Ts, int... I0s, int... P0s, typename A, typename B, typename... Vs, typename... Es>
		auto apply_mul(OpEvaluable<E0> const& e, OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, symphas::lib::types_list<Vs...>> const& sum,
			Es&&... terms) const
		{
			return apply_mul(sum, *static_cast<E0 const*>(&e), std::forward<Es>(terms)...);
		}

		template<typename E0, typename V0, typename E, typename... Ts, int... I0s, int... P0s, typename A, typename B, typename... Vs, typename... Es>
		auto apply_mul(OpOperator<E0> const& e, OpSum<V0, E, Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, symphas::lib::types_list<Vs...>> const& sum,
			Es&&... terms) const
		{
			return apply_mul(sum, *static_cast<E0 const*>(&e), std::forward<Es>(terms)...);
		}

		template<typename A, typename B, typename E, typename... Es>
		auto apply_mul(OpEvaluable<E> const& e0, OpBinaryMul<A, B> const& e1, Es&&... terms) const
		{
			return apply_mul(e1.a, e1.b, *static_cast<E const*>(&e0), std::forward<Es>(terms)...);
		}

		template<typename A, typename B, typename E, typename... Es>
		auto apply_mul(OpOperator<E> const& e0, OpBinaryMul<A, B> const& e1, Es&&... terms) const
		{
			return apply_mul(e1.a, e1.b, *static_cast<E const*>(&e0), std::forward<Es>(terms)...);
		}

		template<typename A, typename B, typename E, typename... Es>
		auto apply_mul(OpBinaryMul<A, B> const& e0, OpEvaluable<E> const& e1, Es&&... terms) const
		{
			return apply_mul(*static_cast<E const*>(&e1), e0, std::forward<Es>(terms)...);
		}

		template<typename A, typename B, typename C, typename D, typename... Es>
		auto apply_mul(OpBinaryMul<A, B> const& e0, OpBinaryMul<C, D> const& e1, Es&&... terms) const
		{
			return apply_mul(e0 * e1.a, e1.b, std::forward<Es>(terms)...);
		}
	};


	


	template<typename E, typename assign_type>
	void result_by_group(OpEvaluable<E> const& e, assign_type&& data, symphas::lib::types_list<>)
	{
		result_accumulate(*static_cast<E const*>(&e), std::forward<assign_type>(data));
	}

	template<typename E, typename assign_type, typename condition_t, typename... condition_ts>
	void result_by_group(OpEvaluable<E> const& e, assign_type&& data, symphas::lib::types_list<condition_t, condition_ts...>)
	{
		auto&& [eval, rest] = expr::split::separate_by<condition_t>(*static_cast<E const*>(&e));
		result_accumulate(eval, std::forward<assign_type>(data));
		result_by_group(rest, std::forward<assign_type>(data), symphas::lib::types_list<condition_ts...>{});
	}

	template<typename condition_t, typename... condition_ts, typename E, typename assign_type, typename E0>
	void result_by_group(OpEvaluable<E> const& e, assign_type&& data, E0 const& init)
	{
		result(init, std::forward<assign_type>(data));
		result_by_group(*static_cast<E const*>(&e), std::forward<assign_type>(data), symphas::lib::types_list<condition_t, condition_ts...>{});
	}



	template<typename E, typename assign_type, typename region_type>
	void result_by_group(OpEvaluable<E> const& e, assign_type&& data, region_type&& region, symphas::lib::types_list<>)
	{
		result_accumulate(*static_cast<E const*>(&e), std::forward<assign_type>(data), std::forward<region_type>(region));
	}

	template<typename E, typename assign_type, typename region_type, typename condition_t, typename... condition_ts>
	void result_by_group(OpEvaluable<E> const& e, assign_type&& data, region_type&& region, symphas::lib::types_list<condition_t, condition_ts...>)
	{
		auto&& [eval, rest] = expr::split::separate_by<condition_t>(*static_cast<E const*>(&e));
		result_accumulate(eval, std::forward<assign_type>(data));
		result_by_group(rest, std::forward<assign_type>(data), std::forward<region_type>(region), symphas::lib::types_list<condition_ts...>{});
	}

	template<typename condition_t, typename... condition_ts, typename E, typename assign_type, typename region_type, typename E0>
	void result_by_group(OpEvaluable<E> const& e, assign_type&& data, region_type&& region, E0 const& init)
	{
		result(init, std::forward<assign_type>(data), std::forward<region_type>(region));
		result_by_group(*static_cast<E const*>(&e), std::forward<assign_type>(data), std::forward<region_type>(region), symphas::lib::types_list<condition_t, condition_ts...>{});
	}

	//! Evaluate only the terms of the expression matching the condition.
	/*!
	 * The expression must be iterable over the entire given length.
	 *
	 * \param e Expression that is evaluated.
	 * \param data The array containing the result of the expression.
	 * \param len The number of elements in the array.
	 */
	template<typename condition_t, typename... condition_ts, typename E, typename assign_type>
	void result_only(OpEvaluable<E> const& e, assign_type&& data, len_type len)
	{
		if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>)
		{
			symphas::data_iterator it(std::forward<assign_type>(data));

			std::transform(
#ifdef EXECUTION_HEADER_AVAILABLE
				std::execution::par_unseq,
#endif
				static_cast<const E*>(&e)->begin(),
				static_cast<const E*>(&e)->end(len), it,
				forward_value{});
		}
	}

	template<typename condition_t, typename... condition_ts, typename E, typename assign_type, size_t D>
	void result_only(OpEvaluable<E> const& e, assign_type&& data, grid::region_interval<D> const& interval)
	{
		if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>)
		{
			symphas::data_iterator_group it(std::forward<assign_type>(data), interval);

			std::transform(
#ifdef EXECUTION_HEADER_AVAILABLE
				std::execution::par_unseq,
#endif
				static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
				static_cast<const E*>(&e)->end(symphas::it_grp, interval), it,
				forward_value{});
		}
	}

	template<typename condition_t, typename... condition_ts, typename E, typename assign_type, size_t D>
	void result_only(OpEvaluable<E> const& e, assign_type&& data, grid::region_interval_multiple<D> const& regions)
	{
		for (grid::region_interval<D> region : regions)
		{
			result_only<condition_t, condition_ts...>(*static_cast<E const*>(&e), std::forward<assign_type>(data), region);
		}
	}


	template<typename condition_t, typename... condition_ts, typename E, typename assign_type>
	void result_only(OpEvaluable<E> const& e, assign_type&& data, grid::region_interval<0> const& interval)
	{
		auto data_region = expr::iterable_domain(std::forward<assign_type>(data));
		if (grid::length(data_region) > 1)
		{
			result_only<condition_t, condition_ts...>(*static_cast<E const*>(&e), std::forward<assign_type>(data), data_region);
		}
		else
		{
			result_only<condition_t, condition_ts...>(*static_cast<E const*>(&e), std::forward<assign_type>(data), 1);
		}
	}

	template<typename condition_t, typename... condition_ts, typename E, typename assign_type>
	void result_only(OpEvaluable<E> const& e, assign_type&& data, grid::region_empty) {}



	template<typename... Es, typename assign_type, size_t... Is>
	void result_only(OpAdd<Es...> const& e, assign_type&& data, len_type len, std::index_sequence<Is...>)
	{
		symphas::data_iterator it(std::forward<assign_type>(data));
		auto start = symphas::reduce_seq_iterator(expr::get<Is>(e).begin()...);
		auto end = symphas::reduce_seq_iterator(expr::get<Is>(e).end(len)...);

		std::transform(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par_unseq,
#endif
			start,
			end, it,
			forward_value{});
	}

	template<typename... Es, typename assign_type, size_t... Is>
	void result_only(OpAdd<Es...> const& e, assign_type&& data, grid::region_interval<0> const& interval, std::index_sequence<Is...>)
	{
		auto data_region = expr::iterable_domain(std::forward<assign_type>(data));
		if (grid::length(data_region) > 1)
		{
			result_only(e, std::forward<assign_type>(data), data_region, std::index_sequence<Is...>{});
		}
		else
		{
			result_only(e, std::forward<assign_type>(data), 1, std::index_sequence<Is...>{});
		}
	}

	template<typename... Es, typename assign_type, size_t D, size_t... Is>
	void result_only(OpAdd<Es...> const& e, assign_type&& data, grid::region_interval<D> const& interval, std::index_sequence<Is...>)
	{
		symphas::data_iterator_group it(std::forward<assign_type>(data), interval);
		auto start = symphas::reduce_iterator(expr::get<Is>(e).begin(symphas::it_grp, interval)...);
		auto end = symphas::reduce_iterator(expr::get<Is>(e).end(symphas::it_grp, interval)...);

		std::transform(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par_unseq,
#endif
			start,
			end, it,
			forward_value{});
	}

	template<typename... Es, typename assign_type, size_t D, size_t... Is>
	void result_only(OpAdd<Es...> const& e, assign_type&& data, grid::region_interval_multiple<D> const& regions, std::index_sequence<Is...>)
	{
		for (grid::region_interval<D> region : regions)
		{
			result_only(e, std::forward<assign_type>(data), region, std::index_sequence<Is...>{});
		}
	}

	template<typename... Es, typename assign_type>
	void result_only(OpAdd<Es...> const& e, assign_type&& data, std::index_sequence<>) 
	{
		result(OpVoid{}, std::forward<assign_type>(data), expr::iterable_domain(std::forward<assign_type>(data)));
	}

	template<typename... Es, typename assign_type, size_t I0, size_t... Is>
	void result_only(OpAdd<Es...> const& e, assign_type&& data, std::index_sequence<I0, Is...>)
	{
		auto region = (expr::iterable_domain(expr::get<Is>(e)) + ... + expr::iterable_domain(expr::get<I0>(e)));
		result_only(e, std::forward<assign_type>(data), region, std::index_sequence<I0, Is...>{});
	}

	template<typename condition_t, typename... condition_ts, typename... Es, typename assign_type>
	void result_only(OpAdd<Es...> const& e, assign_type&& data)
	{
		using namespace symphas::internal;
		using mask = std::integer_sequence<bool, expression_satisfies_condition<Es, expr::or_<condition_t, condition_ts...>>...>;
		using seq = std::make_index_sequence<sizeof...(Es)>;
		result_only(e, std::forward<assign_type>(data), symphas::lib::filter_seq_t<seq, mask>{});
	}

	template<typename condition_t, typename... condition_ts, typename E, typename assign_type>
	void result_only(OpEvaluable<E> const& e, assign_type&& data)
	{
		if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>)
		{
			result(*static_cast<E const*>(&e), std::forward<assign_type>(data));
		}
	}


	template<typename condition_t>
	struct result_by_term_apply;

	template<typename... Es, typename assign_type, size_t... Is>
	void result_by_term(OpAdd<Es...> const& e, assign_type&& data, symphas::lib::types_list<>, std::index_sequence<Is...>) {}

	template<typename... Es, typename assign_type, typename condition_t, typename... condition_ts, size_t... Is>
	void result_by_term(OpAdd<Es...> const& e, assign_type&& data, symphas::lib::types_list<condition_t, condition_ts...>, std::index_sequence<Is...>)
	{
		(result_by_term_apply<condition_t>{}(expr::get<Is>(e), std::forward<assign_type>(data)), ...);
		result_by_term(e, std::forward<assign_type>(data), symphas::lib::types_list<condition_ts...>{}, std::index_sequence<Is...>{});
	}

	template<typename condition_t, typename... condition_ts, typename... Es, typename assign_type>
	void result_by_term(OpAdd<Es...> const& e, assign_type&& data)
	{
		if constexpr ((expr::satisfies<Es, expr::or_<condition_t, condition_ts...>> || ...))
		{
			result_only<expr::not_<expr::or_<condition_t, condition_ts...>>>(e, std::forward<assign_type>(data));
			result_by_term(e, std::forward<assign_type>(data), symphas::lib::types_list<condition_t, condition_ts...>{}, std::make_index_sequence<sizeof...(Es)>{});
		}
		else
		{
			result(e, std::forward<assign_type>(data));
		}
	}

	template<typename condition_t, typename... condition_ts, typename E, typename assign_type>
	void result_by_term(OpEvaluable<E> const& e, assign_type&& data)
	{
		result(*static_cast<E const*>(&e), std::forward<assign_type>(data));
	}



	//! Evaluate only the terms of the expression matching the condition.
	/*!
	 * The expression must be iterable over the entire given length.
	 *
	 * \param e Expression that is evaluated.
	 * \param data The array containing the result of the expression.
	 * \param len The number of elements in the array.
	 */
	template<typename condition_t, typename... condition_ts, typename E, typename assign_type>
	auto result_sum_only(OpEvaluable<E> const& e, len_type len)
	{
		if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>)
		{
			return std::reduce(
#ifdef EXECUTION_HEADER_AVAILABLE
				std::execution::par_unseq,
#endif
				static_cast<const E*>(&e)->begin(),
				static_cast<const E*>(&e)->end(len));
		}
		else
		{
			return expr::eval_type_t<E>{};
		}
	}

	template<typename condition_t, typename... condition_ts, typename E, size_t D>
	auto result_sum_only(OpEvaluable<E> const& e, grid::region_interval<D> const& interval)
	{
		if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>)
		{
			return std::reduce(
#ifdef EXECUTION_HEADER_AVAILABLE
				std::execution::par_unseq,
#endif
				static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
				static_cast<const E*>(&e)->end(symphas::it_grp, interval));
		}
		else
		{
			return expr::eval_type_t<E>{};
		}
	}

	template<typename condition_t, typename... condition_ts, typename E, size_t D>
	auto result_sum_only(OpEvaluable<E> const& e, grid::region_interval_multiple<D> const& regions)
	{
		expr::eval_type_t<E> result{};
		for (grid::region_interval<D> region : regions)
		{
			result += result_sum_only<condition_t, condition_ts...>(*static_cast<E const*>(&e), region);
		}
		return result;
	}


	template<typename condition_t, typename... condition_ts, typename E>
	auto result_sum_only(OpEvaluable<E> const& e, grid::region_interval<0> const& interval)
	{
		return result_sum_only(*static_cast<E const*>(&e), 1);
	}

	template<typename condition_t, typename... condition_ts, typename E>
	auto result_sum_only(OpEvaluable<E> const& e, grid::region_empty) 
	{
		return expr::eval_type_t<E>{};
	}


	template<typename... Es, size_t... Is>
	auto result_sum_only(OpAdd<Es...> const& e, len_type len, std::index_sequence<Is...>)
	{
		auto start = symphas::reduce_seq_iterator(expr::get<Is>(e).begin()...);
		auto end = symphas::reduce_seq_iterator(expr::get<Is>(e).end(len)...);

		return std::reduce(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par_unseq,
#endif
			start, end);
	}

	template<typename... Es, size_t... Is>
	auto result_sum_only(OpAdd<Es...> const& e, grid::region_interval<0> const& interval, std::index_sequence<Is...>)
	{
		return result_sum_only(e, 1, std::index_sequence<Is...>{});
	}

	template<typename... Es, size_t D, size_t... Is>
	auto result_sum_only(OpAdd<Es...> const& e, grid::region_interval<D> const& interval, std::index_sequence<Is...>)
	{
		auto start = symphas::reduce_iterator(expr::get<Is>(e).begin(symphas::it_grp, interval)...);
		auto end = symphas::reduce_iterator(expr::get<Is>(e).end(symphas::it_grp, interval)...);

		return std::reduce(
#ifdef EXECUTION_HEADER_AVAILABLE
			std::execution::par_unseq,
#endif
			start, end);
	}

	template<typename... Es, size_t D, size_t... Is>
	auto result_sum_only(OpAdd<Es...> const& e, grid::region_interval_multiple<D> const& regions, std::index_sequence<Is...>)
	{
		expr::eval_type_t<OpAdd<symphas::lib::type_at_index<Is, Es...>...>> result{};
		for (grid::region_interval<D> region : regions)
		{
			result += result_sum_only(e, region, std::index_sequence<Is...>{});
		}
		return result;
	}

	template<typename... Es>
	auto result_sum_only(OpAdd<Es...> const& e, std::index_sequence<>)
	{
		return expr::eval_type_t<OpAdd<Es...>>{};
	}

	template<typename... Es, size_t I0, size_t... Is>
	auto result_sum_only(OpAdd<Es...> const& e, std::index_sequence<I0, Is...>)
	{
		auto region = (expr::iterable_domain(expr::get<Is>(e)) + ... + expr::iterable_domain(expr::get<I0>(e)));
		return result_sum_only(e, region, std::index_sequence<I0, Is...>{});
	}

	template<typename condition_t, typename... condition_ts, typename... Es>
	auto result_sum_only(OpAdd<Es...> const& e)
	{
		using namespace symphas::internal;
		using mask = std::integer_sequence<bool, expression_satisfies_condition<Es, expr::or_<condition_t, condition_ts...>>...>;
		using seq = std::make_index_sequence<sizeof...(Es)>;
		return result_sum_only(e, symphas::lib::filter_seq_t<seq, mask>{});
	}

	template<typename condition_t, typename... condition_ts, typename E>
	auto result_sum_only(OpEvaluable<E> const& e)
	{
		if constexpr (expr::satisfies<E, expr::or_<condition_t, condition_ts...>>)
		{
			return result(*static_cast<E const*>(&e));
		}
		else
		{
			return OpVoid{};
		}
	}



	template<typename condition_t>
	struct result_sum_by_term_apply;

	template<typename... Es, size_t... Is>
	auto result_sum_by_term(OpAdd<Es...> const& e, symphas::lib::types_list<>, std::index_sequence<Is...>)
	{
		return eval_type_t<OpAdd<Es...>>{};
	}

	template<typename condition_t, typename... condition_ts, typename... Es, size_t... Is>
	auto result_sum_by_term(OpAdd<Es...> const& e, symphas::lib::types_list<condition_t, condition_ts...>, std::index_sequence<Is...>)
	{
		return (result_sum_by_term_apply<condition_t>{}(expr::get<Is>(e)) + ...)
			+ result_sum_by_term(e, symphas::lib::types_list<condition_ts...>{}, std::index_sequence<Is...>{});
	}

	template<typename condition_t, typename... condition_ts, typename... Es>
	auto result_sum_by_term(OpAdd<Es...> const& e)
	{
		return result_sum_only<expr::not_<expr::or_<condition_t, condition_ts...>>>(e)
			+ result_sum_by_term(e, symphas::lib::types_list<condition_t, condition_ts...>{}, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<typename condition_t, typename... condition_ts, typename E>
	auto result_sum_by_term(OpEvaluable<E> const& e)
	{
		return result_sum(*static_cast<E const*>(&e));
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


	template<typename G_F>
	auto swap_matching_i(OpIdentity, G_F&& g)
	{
		return OpIdentity{};
	}

	//template<int N0, int... P0s, expr::exp_key_t... Xs, typename G_F>
	//auto swap_matching_i(OpTermsList<> const& e, G_F&& g)
	//{
	//	return expr::pow<X0>(expr::transform::swap_grid<expr::symbols::i_<N0, P0>>(expr::symbols::i_<N0, P0>{}, std::forward<G_F>(g)))
	//		* swap_matching_i(expr::terms_after_first(e), std::forward<G_F>(g));
	//}

	//template<int N0, int P0, expr::exp_key_t X0, typename G_F>
	//auto swap_matching_i(OpTermsList<Term<expr::symbols::i_<N0, P0>, X0>> const& e, G_F&& g)
	//{
	//	return expr::pow_x<X0>(expr::transform::swap_grid<expr::symbols::i_<N0, P0>>(expr::symbols::i_<N0, P0>{}, std::forward<G_F>(g)));
	//}

	//template<int N0, int P0, int... P0s, expr::exp_key_t X0, expr::exp_key_t... Xs, typename G_F>
	//auto swap_matching_i(OpTermsList<Term<expr::symbols::i_<N0, P0s>, Xs>...> const& e, G_F&& g)
	//{
	//	return (expr::pow_x<Xs>(expr::transform::swap_grid<expr::symbols::i_<N0, P0s>>(expr::symbols::i_<N0, P0>{}, std::forward<G_F>(g))) * ...);
	//		//* swap_matching_i(expr::terms_after_first(e), std::forward<G_F>(g));
	//}

	template<int N0, int... P0s, expr::exp_key_t... Xs, typename G_F>
	auto swap_matching_i(OpTerms<OpIdentity, Term<expr::symbols::i_<N0, P0s>, Xs>...> const& e, G_F&& g)
	{
		return (expr::pow_x<Xs>(expr::transform::swap_grid<expr::symbols::i_<N0, P0s>>(expr::symbols::i_<N0, P0s>{}, std::forward<G_F>(g))) * ...);
		//return swap_matching_i(*static_cast<OpTermsList<Term<expr::symbols::i_<N0, P0s>, Xs>...> const*>(&e), std::forward<G_F>(g));
	}

	template<typename E, int N0, int... P0s, typename G_F>
	auto swap_matching_i(OpExpression<E> const& e, symphas::lib::types_list<expr::symbols::i_<N0, P0s>...>, G_F&& g)
	{
		auto [ind, r] = expr::split::factor<expr::symbols::i_<N0, P0s>...>(*static_cast<E const*>(&e));
		return expr::coeff(r) * swap_matching_i(ind, std::forward<G_F>(g)) * r;
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

		template<typename E1, typename E2>
		auto operator()(OpExpression<E1> const& a, OpOperator<E2> const& b)
		{
			return expr::make_div(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
		}

		template<typename E1, typename E2>
		auto operator()(OpOperator<E1> const& a, OpExpression<E2> const& b)
		{
			return expr::make_div(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
		}

		template<typename E1, typename E2>
		auto operator()(OpOperator<E1> const& a, OpOperator<E2> const& b)
		{
			return expr::make_div(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
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

		template<typename E1, typename E2>
		auto operator()(OpExpression<E1> const& a, OpOperator<E2> const& b)
		{
			return expr::make_div(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
		}

		template<typename E1, typename E2>
		auto operator()(OpOperator<E1> const& a, OpExpression<E2> const& b)
		{
			return expr::make_div(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
		}

		template<typename E1, typename E2>
		auto operator()(OpOperator<E1> const& a, OpOperator<E2> const& b)
		{
			return expr::make_div(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
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

		template<typename Dd, typename V, typename E, typename Sp, typename E0>
		auto apply_operators_deriv(OpDerivative<Dd, V, E, Sp> const& d, OpOperator<E0> const& e)
		{
			auto [dd, ee] = expr::split::separate_operator(d);
			return OpOperatorChain(dd, ee);
		}

		template<typename G0, typename E>
		struct symbolic_derivative_function
		{
			template<typename Sp>
			auto operator()(symphas::internal::wrap_base, E const& e)
			{
				return OpVoid{};
			}

			template<typename Sp>
			auto operator()(symphas::internal::wrap_f<func_cos<E>>, E const& e, Sp const& solver)
			{
				return -OpFunctionApply<func_sin<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e, solver));
			}

			template<typename Sp>
			auto operator()(symphas::internal::wrap_f<func_sin<E>>, E const& e, Sp const& solver)
			{
				return OpFunctionApply<func_cos<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e, solver));
			}

			template<typename Sp>
			auto operator()(symphas::internal::wrap_f<func_tan<E>>, E const& e, Sp const& solver)
			{
				return OpFunctionApply<func_sec<E>, OpIdentity, E>(e) * OpFunctionApply<func_sec<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e, solver));
			}

			template<typename Sp>
			auto operator()(symphas::internal::wrap_f<func_csc<E>>, E const& e, Sp const& solver)
			{
				return -OpFunctionApply<func_cot<E>, OpIdentity, E>(e) * OpFunctionApply<func_csc<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e, solver));
			}

			template<typename Sp>
			auto operator()(symphas::internal::wrap_f<func_sec<E>>, E const& e, Sp const& solver)
			{
				return OpFunctionApply<func_tan<E>, OpIdentity, E>(e) * OpFunctionApply<func_sec<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e, solver));
			}

			template<typename Sp>
			auto operator()(symphas::internal::wrap_f<func_cot<E>>, E const& e, Sp const& solver)
			{
				return -OpFunctionApply<func_csc<E>, OpIdentity, E>(e) * OpFunctionApply<func_csc<E>, OpIdentity, E>(e) *
					apply_operators(expr::make_derivative<1, G0>(e, solver));
			}

			template<typename Sp>
			auto operator()(symphas::internal::wrap_f<func_sqrt<E>>, E const& e, Sp const& solver)
			{
				return expr::make_fraction<1, 2>() *
					expr::inverse(OpFunctionApply<func_sqrt<E>, OpIdentity, E>(e)) *
					apply_operators(expr::make_derivative<1, G0>(e, solver));
			}

			template<typename Sp>
			auto operator()(symphas::internal::wrap_f<func_log<E>>, E const& e, Sp const& solver)
			{
				return expr::inverse(e);
			}
		};
	}


	template<size_t O, typename V0, auto f, typename V1, typename E, typename G0>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V0, OpFunctionApply<f, V1, E>, SymbolicDerivative<G0>> const& e)
	{
		auto&& function = expr::get_enclosed_expression(e);
		auto&& expr = expr::get_enclosed_expression(function);

		return apply_operators(
			expr::make_derivative<O - 1, G0>(
				expr::coeff(e) * expr::coeff(function),
				symbolic_derivative_function<G0, E>{}(symphas::internal::wrap_f<f>{}, expr, e.solver),
				e.solver));
	}

	template<size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<expr::symbols::i_<I0, P0>, X0>>, SymbolicDerivative<GG>> const& e)
	{
		return OpVoid{};
	}

	template<size_t O, typename V, typename V1, int I0, int P0, typename G1, typename... Gs,
		expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t... Xs, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V1, Term<expr::symbols::i_<I0, P0>, X0>, Term<G1, X1>, Term<Gs, Xs>...>, SymbolicDerivative<GG>> const& e)
	{
		auto&& terms = expr::get_enclosed_expression(e);
		auto coeff = expr::coeff(terms);
		auto a = OpTerms(OpIdentity{}, expr::get<1>(terms));
		auto b = expr::terms_after_n<1>(terms);// OpTerms(OpIdentity{}, *static_cast<OpTerms<Term<G1, X1>, Term<Gs, Xs>...> const*>(&terms));

		return (coeff * expr::coeff(e) * a) * expr::apply_operators(expr::make_derivative<O, GG>(expr::make_term(b), e.solver));
	}

	template<size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0, size_t Z, typename GG, size_t N, typename>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>, SymbolicDerivative<Variable<Z, GG>>> const& e)
	{
		using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
		auto f = expr::split::factor<O, factor_t>(expr::get_enclosed_expression(e));

		SymbolicCase c(expr::symbols::i_<I0, P0>{} = Variable<Z, GG>{}, f.second, OpVoid{});

		return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
	}

	template<size_t O, typename V, typename V0, int I0, int P0, size_t D, expr::exp_key_t X0, size_t Z, typename GG, size_t N, typename>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, D>, X0>>, SymbolicDerivative<Variable<Z, GG>>> const& e)
	{
		using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
		auto f = expr::split::factor<O, factor_t>(expr::get_enclosed_expression(e));

		SymbolicCase c(expr::symbols::i_<I0, P0>{} = Variable<Z, GG>{}, f.second, OpVoid{});

		return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
	}

	template<size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0, typename GG, size_t N, typename>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>, SymbolicDerivative<DynamicVariable<GG>>> const& e)
	{
		using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
		auto f = expr::split::factor<O, factor_t>(expr::get_enclosed_expression(e));

		SymbolicCase c(expr::symbols::i_<I0, P0>{} = DynamicVariable<GG>{}, f.second, OpVoid{});

		return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
	}

	template<size_t O, typename V, typename V0, int I0, int P0, size_t D, expr::exp_key_t X0, typename GG, size_t N, typename>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, D>, X0>>, SymbolicDerivative<DynamicVariable<GG>>> const& e)
	{
		using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
		auto f = expr::split::factor<O, factor_t>(expr::get_enclosed_expression(e));

		SymbolicCase c(expr::symbols::i_<I0, P0>{} = DynamicVariable<GG>{}, f.second, OpVoid{});

		return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
	}

	template<size_t O, typename V, typename V0, int I0, int P0, expr::exp_key_t X0, typename GG, size_t N, typename>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, X0>>, SymbolicDerivative<GG>> const& e)
	{
		using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
		auto f = expr::split::factor<O, factor_t>(expr::get_enclosed_expression(e));

		SymbolicCase c(expr::symbols::i_<I0, P0>{} = GG{}, f.second, OpVoid{});

		return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
	}

	template<size_t O, typename V, typename V0, int I0, int P0, size_t D, expr::exp_key_t X0, typename GG, size_t N, typename>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, D>, X0>>, SymbolicDerivative<GG>> const& e)
	{
		using factor_t = expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>;
		auto f = expr::split::factor<O, factor_t>(expr::get_enclosed_expression(e));

		SymbolicCase c(expr::symbols::i_<I0, P0>{} = GG{}, f.second, OpVoid{});

		return (expr::factorial<N, N - O>() * expr::coeff(e)) * expr::make_term(c);
	}

	template<size_t O, typename V, typename V0, typename G0, expr::exp_key_t X0, typename GG, size_t N, typename>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>>, SymbolicDerivative<GG>> const& e)
	{
		auto f = expr::split::factor<O, GG>(expr::get_enclosed_expression(e));
		return expr::factorial<N, N - O>() * expr::coeff(e) * f.second;
	}

	template<size_t O, typename V, typename V0, typename G0, typename G1, typename... Gs,
		expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t... Xs, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>, SymbolicDerivative<GG>> const& e)
	{
		auto&& terms = expr::get_enclosed_expression(e);
		auto coeff = expr::coeff(terms);
		auto a = OpTerms(OpIdentity{}, expr::get<1>(terms));
		auto b = expr::terms_after_n<1>(terms);

		auto lhs = apply_operators(expr::make_derivative<O, GG>(a, e.solver)) * expr::make_term(b);
		auto rhs = a * apply_operators(expr::make_derivative<O, GG>(expr::make_term(b), e.solver));
		return (coeff * expr::coeff(e)) * (lhs + rhs);
	}


	template<size_t O, typename V, typename V0, typename G0, expr::exp_key_t X0, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>>, SymbolicDerivative<OpTerm<OpIdentity, GG>>> const& e)
	{
		return apply_operators(expr::coeff(e) * expr::make_derivative<O, GG>(expr::get_enclosed_expression(e), e.solver));
	}

	template<size_t O, typename V, typename V0, typename G0, typename G1, typename... Gs,
		expr::exp_key_t X0, expr::exp_key_t X1, expr::exp_key_t... Xs, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpTerms<V0, Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>, SymbolicDerivative<OpTerm<OpIdentity, GG>>> const& e)
	{
		return apply_operators(expr::coeff(e) * expr::make_derivative<O, GG>(expr::get_enclosed_expression(e), e.solver));
	}
	

	template<typename V, expr::exp_key_t X0, typename V0, typename E0, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<1>, V, OpPow<X0, V0, E0>, SymbolicDerivative<GG>> const& e)
	{
		auto power = expr::get_enclosed_expression(e);
		auto p = apply_operators(expr::get_enclosed_expression(power));


		constexpr size_t N0 = expr::_Xk_t<X0>::N;
		constexpr size_t D0 = expr::_Xk_t<X0>::D;
		constexpr bool sign = expr::_Xk_t<X0>::sign;
		
		constexpr expr::exp_key_t _X = (sign) ? N0 + D0 : (N0 < D0) ? D0 - N0 : N0 - D0;
		constexpr expr::exp_key_t _sign = (sign) ? sign : (N0 < D0) ? true : false;

		auto result = expr::coeff(e) * expr::coeff(power) * expr::make_fraction<N0, D0>() * 
			expr::dot(expr::make_pow<_X>(p), apply_operators(expr::make_derivative<1, GG>(p, e.solver)));
		
		if constexpr (sign)
		{
			return -result;
		}
		else
		{
			return result;
		}
	}

	template<size_t O, typename V, expr::exp_key_t X0, typename V0, typename E0, typename GG>
	auto apply_operators(OpDerivative<std::index_sequence<O>, V, OpPow<X0, V0, E0>, SymbolicDerivative<GG>> const& e)
	{
		return expr::coeff(e) * apply_operators(expr::make_derivative<O - 1, GG>(
			apply_operators(expr::make_derivative<1, GG>(expr::get_enclosed_expression(e), e.solver))));
	}

	template<expr::exp_key_t X, typename V, typename E>
	auto apply_operators(OpPow<X, V, E> const& e)
	{
		auto p = apply_operators(expr::get_enclosed_expression(e));
		return expr::coeff(e) * expr::make_pow<X>(p);
		//if constexpr (expr::_Xk<X> > 1)
		//{
		//	constexpr size_t N0 = expr::_Xk_t<X>::N;
		//	constexpr size_t D0 = expr::_Xk_t<X>::D;
		//	constexpr bool sign = expr::_Xk_t<X>::sign;
		//
		//	return expr::coeff(e) * apply_operators(expr::make_pow<expr::Xk<N0 - D0, D0, sign>>(p)) * p;
		//}
		//else
		//{
		//	return expr::coeff(e) * expr::make_pow<X>(p);
		//}
	}


	namespace
	{

		template<typename... Ts>
		struct filter_for_case;

		template<typename... Seqs, typename... Ts>
		struct filter_for_case<symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>>
		{
			using type = symphas::lib::types_list<Ts...>;
		};

		template<typename... Seqs, typename... Ts, typename G0, typename... Gs>
		struct filter_for_case<symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>, G0, Gs...>
		{
			using type = typename filter_for_case<symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>, Gs...>::type;
		};

		template<typename... Seqs, typename... Ts, int I0, int P0, size_t Z0, typename G0, typename L, typename R, typename... Gs>
		struct filter_for_case<
			std::integer_sequence<int, -1>,
			symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
			SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>, Gs...>
		{
			using type = typename filter_for_case<
				symphas::lib::types_list<Seqs..., expr::symbols::i_<I0, P0>>,
				symphas::lib::types_list<Ts..., symphas::lib::types_list<SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>>>, Gs...>::type;
		};

		template<int I, typename... Seqs, typename... Ts, int I0, int P0, size_t Z0, typename G0, typename L, typename R, typename... Gs>
		struct filter_for_case<
			std::integer_sequence<int, I>,
			symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
			SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>, Gs...>
		{
			static const size_t N = size_t(I);
			using type_N = symphas::lib::type_at_index<N, Ts...>;
			
			using type_combined = symphas::lib::expand_types_list<type_N, SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>>;
			using seq_combined = symphas::lib::expand_types_list<symphas::lib::types_before_index<N, Seqs...>, expr::symbols::i_<I0, P0>, symphas::lib::types_after_at_index<N + 1, Seqs...>>;

			using type = typename filter_for_case<
				seq_combined,
				symphas::lib::expand_types_list<symphas::lib::types_before_index<N, Ts...>, symphas::lib::types_list<type_combined>, symphas::lib::types_after_at_index<N + 1, Ts...>>,
				Gs...>::type;
		};

		template<typename... Seqs, typename... Ts, int I0, int P0, size_t Z0, typename G0, typename L, typename R, typename... Gs>
		struct filter_for_case<
			symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>, 
			SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>, Gs...>
		{
			static const int ind = symphas::lib::index_of_type<expr::symbols::i_<I0, P0>, Seqs...>;

			using type = typename filter_for_case<
				std::integer_sequence<int, ind>,
				symphas::lib::types_list<Seqs...>,
				symphas::lib::types_list<Ts...>, SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, L, R>, Gs...>::type;
		};

		//template<typename... Seqs, typename... Ts, int I0, int P0, typename G0, typename L, typename R, typename... Gs>
		//struct filter_for_case<
		//	std::integer_sequence<int, -1>,
		//	symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
		//	SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, L, R>, Gs...>
		//{
		//	using type = typename filter_for_case<
		//		symphas::lib::types_list<Seqs..., expr::symbols::i_<I0, P0>>,
		//		symphas::lib::types_list<Ts..., symphas::lib::types_list<SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, L, R>>>, Gs...>::type;
		//};

		//template<int I, typename... Seqs, typename... Ts, int I0, int P0, typename G0, typename L, typename R, typename... Gs>
		//struct filter_for_case<
		//	std::integer_sequence<int, I>,
		//	symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
		//	SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, L, R>, Gs...>
		//{
		//	static const size_t N = size_t(I);
		//	using type_N = symphas::lib::type_at_index<N, Ts...>;

		//	using type_combined = symphas::lib::expand_types_list<type_N, SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, L, R>>;
		//	using seq_combined = symphas::lib::expand_types_list<symphas::lib::types_before_index<N, Seqs...>, expr::symbols::i_<I0, P0>, symphas::lib::types_after_at_index<N + 1, Seqs...>>;

		//	using type = typename filter_for_case<
		//		seq_combined,
		//		symphas::lib::expand_types_list<symphas::lib::types_before_index<N, Ts...>, symphas::lib::types_list<type_combined>, symphas::lib::types_after_at_index<N + 1, Ts...>>,
		//		Gs...>::type;
		//};

		//template<typename... Seqs, typename... Ts, int I0, int P0, typename G0, typename L, typename R, typename... Gs>
		//struct filter_for_case<
		//	symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
		//	SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, L, R>, Gs...>
		//{
		//	static const int ind = symphas::lib::index_of_type<expr::symbols::i_<I0, P0>, Seqs...>;

		//	using type = typename filter_for_case<
		//		std::integer_sequence<int, ind>,
		//		symphas::lib::types_list<Seqs...>,
		//		symphas::lib::types_list<Ts...>, SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, L, R>, Gs...>::type;
		//};


		template<typename... Seqs, typename... Ts, int I0, int P0, typename G0, typename L, typename R, typename... Gs>
		struct filter_for_case<
			std::integer_sequence<int, -1>,
			symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
			SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>, Gs...>
		{
			using type = typename filter_for_case<
				symphas::lib::types_list<Seqs..., expr::symbols::i_<I0, P0>>,
				symphas::lib::types_list<Ts..., symphas::lib::types_list<SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>>>, Gs...>::type;
		};

		template<int I, typename... Seqs, typename... Ts, int I0, int P0, typename G0, typename L, typename R, typename... Gs>
		struct filter_for_case<
			std::integer_sequence<int, I>,
			symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
			SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>, Gs...>
		{
			static const size_t N = size_t(I);
			using type_N = symphas::lib::type_at_index<N, Ts...>;

			using type_combined = symphas::lib::expand_types_list<type_N, SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>>;
			using seq_combined = symphas::lib::expand_types_list<symphas::lib::types_before_index<N, Seqs...>, expr::symbols::i_<I0, P0>, symphas::lib::types_after_at_index<N + 1, Seqs...>>;

			using type = typename filter_for_case<
				seq_combined,
				symphas::lib::expand_types_list<symphas::lib::types_before_index<N, Ts...>, symphas::lib::types_list<type_combined>, symphas::lib::types_after_at_index<N + 1, Ts...>>,
				Gs...>::type;
		};


		template<typename... Seqs, typename... Ts, int I0, int P0, typename G0, typename L, typename R, typename... Gs>
		struct filter_for_case<
			symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>,
			SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>, Gs...>
		{
			static const int ind = symphas::lib::index_of_type<expr::symbols::i_<I0, P0>, Seqs...>;

			using type = typename filter_for_case<
				std::integer_sequence<int, ind>,
				symphas::lib::types_list<Seqs...>,
				symphas::lib::types_list<Ts...>, SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, L, R>, Gs...>::type;
		};

		template<typename... Seqs, typename... Ts, typename... Gs>
		struct filter_for_case<symphas::lib::types_list<Seqs...>, symphas::lib::types_list<Ts...>, symphas::lib::types_list<Gs...>>
		{
			using type = typename filter_for_case<symphas::lib::types_list<Ts...>, Gs...>::type;
		};

		template<typename... Gs>
		struct filter_for_case<symphas::lib::types_list<Gs...>>
		{
			using type = typename filter_for_case<symphas::lib::types_list<>, symphas::lib::types_list<>, Gs...>::type;
		};

		template<size_t N0, typename T>
		struct pos_in_substitution;

		template<size_t N0, size_t... Ns, typename... Gs>
		struct pos_in_substitution<N0, SymbolicDataArray<std::tuple<Variable<Ns, Gs>...>>>
		{
			static const int value = symphas::lib::index_of_value<size_t, N0, Ns...>;
			using type = typename SymbolicDataArray<std::tuple<Variable<Ns, Gs>...>>::storage_type;
		};

		template<size_t N0, size_t... Ns, typename... Gs>
		struct pos_in_substitution<N0, SymbolicDataArray<std::tuple<NamedData<Variable<Ns, Gs>>...>>>
		{
			static const int value = symphas::lib::index_of_value<size_t, N0, Ns...>;
			using type = typename SymbolicDataArray<std::tuple<NamedData<Variable<Ns, Gs>>...>>::storage_type;
		};

		template<size_t N0, int N, int P>
		struct pos_in_substitution<N0, SymbolicDataArray<expr::symbols::v_id_type<expr::symbols::i_<N, P>>>>
		{
			static const int value = -1;
			using type = expr::symbols::Symbol;
		};

		template<size_t N0, typename T>
		struct pos_in_substitution<N0, SymbolicDataArray<T>>
		{
			static const int value = N0;
			using type = T;
		};

		template<size_t Z0, size_t... Ns, typename... Gs>
		auto get_arg(SymbolicDataArray<std::tuple<Variable<Ns, Gs>...>> const& substitution)
		{
			return expr::make_term<Z0>(*substitution.data[symphas::lib::index_of_value<size_t, Z0, Ns...>].data);
		}

		template<size_t Z0, size_t... Ns, typename... Gs>
		auto get_arg(SymbolicDataArray<std::tuple<NamedData<Variable<Ns, Gs>>...>> const& substitution)
		{
			auto el = substitution.data[symphas::lib::index_of_value<size_t, Z0, Ns...>];
			return expr::make_term<Z0>(NamedData(std::ref(*el.data), el.name));
		}

		template<size_t Z0, int N, int P>
		auto get_arg(SymbolicDataArray<expr::symbols::v_id_type<expr::symbols::i_<N, P>>> const& substitution)
		{
			return expr::symbols::v_<expr::symbols::i_<N, P>>{};
		}

		template<size_t Z0, typename T>
		auto get_arg(SymbolicDataArray<T> const& substitution)
		{
			//auto el = substitution.data[Z0];
			//return expr::make_term(NamedData(std::ref(*el.data), el.name));
			constexpr size_t D = expr::grid_dim<T>::value;
			return expr::make_term<Z0>(GridSymbol<expr::eval_type_t<OpTerm<OpIdentity, T>>, D>{});
		}

		template<size_t Z0, typename... Ts, size_t... Ns>
		auto get_args(Substitution<Ts...> const& substitution, std::index_sequence<Ns...>)
		{
			return std::make_tuple(get_arg<Z0>(std::get<Ns>(substitution))...);
		}

		template<typename T>
		auto get_arg(SymbolicDataArray<T> const& substitution)
		{
			constexpr size_t D = expr::grid_dim<T>::value;
			return expr::make_term(GridSymbol<expr::eval_type_t<OpTerm<OpIdentity, T>>, D>{});
		}

		template<typename T>
		auto get_arg(SymbolicDataArray<T> const& substitution, DynamicIndex const& index)
		{
			return expr::make_term_dynamic(index, substitution.data);
		}

		template<typename T>
		auto get_arg(SymbolicDataArray<NamedData<T*>> const& substitution, DynamicIndex const& index)
		{
			return expr::make_term_dynamic(index, NamedData(substitution.data, substitution.name));
		}

		template<typename... Ts, size_t... Ns>
		auto get_args(Substitution<Ts...> const& substitution, std::index_sequence<Ns...>)
		{
			return std::make_tuple(get_arg(std::get<Ns>(substitution))...);
		}

		template<typename... Ts, size_t... Ns>
		auto get_args(Substitution<Ts...> const& substitution, DynamicIndex const& index, std::index_sequence<Ns...>)
		{
			return std::make_tuple(get_arg(std::get<Ns>(substitution), index)...);
		}



		template<typename T, size_t D>
		auto get_arg(GridSymbol<T, D> const& substitution)
		{
			return substitution;
		}

		template<typename T, size_t D>
		auto get_arg(T const& substitution)
		{
			return substitution;
		}


		template<typename E, typename... Ts, int... I0s, int... P0s, typename B, typename C,
			typename E0, int I0, int P0, size_t Z0, typename G0, typename... Ls, typename... Rs>
		auto apply_operators_sum_case(
			SymbolicSum<E, Substitution<Ts...>, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B, C> const& series,
			OpExpression<E0> const& applied,
			SymbolicDerivative<Variable<Z0, G0>> const& symbol,
			symphas::lib::types_list<SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, Ls, Rs>...>)
		{
			constexpr int ind_N = symphas::lib::index_of_value<int, I0, I0s...>;
			using pos_t = pos_in_substitution<Z0, symphas::lib::type_at_index<size_t(ind_N), Ts...>>;
			constexpr int arg_N = pos_t::value;

			using v_type = expr::symbols::v_<expr::symbols::i_<I0, P0>>;

			auto arg = get_arg<Z0>(std::get<size_t(ind_N)>(series.substitution));

			auto left_case = expr::transform::swap_grid<v_type>
				(SymbolicCase(expr::symbols::i_<I0, P0>{} = Variable<Z0, G0>{}), *static_cast<E0 const*>(&applied), arg);
			auto right_case = expr::transform::swap_grid<v_type>
				(SymbolicCase(expr::symbols::i_<I0, P0>{} != Variable<Z0, G0>{}), *static_cast<E0 const*>(&applied), arg);

			auto left = expr::transform::swap_grid<SymbolicCaseSwap<>>(left_case, OpVoid{});
			auto right = expr::transform::swap_grid<SymbolicCaseSwap<>>(right_case, OpVoid{});
			
			auto series_left = expr::recreate_series(expr::symbols::i_<I0, P0>{} = expr::val<arg_N>, left, series);
			auto series_right = OpVoid{};// expr::recreate_series(expr::symbols::i_<I0, 0>{} != expr::val<arg_N>, expr_right, series);
			return series_left + series_right;
		}
		

		template<typename E, typename... Ts, int... I0s, int... P0s, typename B, typename C,
			typename E0, int I0, int P0, typename G0, typename... Ls, typename... Rs, std::enable_if_t<!expr::is_id_variable<G0>, int> = 0>
		auto apply_operators_sum_case(
			SymbolicSum<E, Substitution<Ts...>, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B, C> const& series,
			OpExpression<E0> const& applied,
			SymbolicDerivative<G0> const& symbol,
			symphas::lib::types_list<SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, Ls, Rs>...>)
		{
			constexpr int ind_N = symphas::lib::index_of_value<int, I0, I0s...>;

			using v_type = expr::symbols::v_<expr::symbols::i_<I0, P0>>;

			//auto arg = get_arg(std::get<size_t(ind_N)>(series.substitution));
			auto arg = get_arg(G0{});

			auto left_case = expr::transform::swap_grid<v_type>
				(SymbolicCase(expr::symbols::i_<I0, P0>{} = G0{}), *static_cast<E0 const*>(&applied), arg);
			auto right_case = expr::transform::swap_grid<v_type>
				(SymbolicCase(expr::symbols::i_<I0, P0>{} != G0{}), *static_cast<E0 const*>(&applied), arg);

			auto left = expr::transform::swap_grid<SymbolicCaseSwap<>>(left_case, OpVoid{});
			auto right = expr::transform::swap_grid<SymbolicCaseSwap<>>(right_case, OpVoid{});
			
			auto series_left = expr::recreate_series(expr::symbols::i_<I0, P0>{} = expr::symbols::placeholder_N{}, left, series);
			auto series_right = OpVoid{};// expr::recreate_series(expr::symbols::i_<I0, 0>{} != expr::val<arg_N>, expr_right, series);
			return series_left + series_right;
		}


		template<typename E, typename... Ts, int... I0s, int... P0s, typename B, typename C,
			typename E0, int I0, int P0, typename G0, typename... Ls, typename... Rs>
		auto apply_operators_sum_case(
			SymbolicSum<E, Substitution<Ts...>, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B, C> const& series,
			OpExpression<E0> const& applied,
			SymbolicDerivative<DynamicVariable<G0>> const& symbol,
			symphas::lib::types_list<SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, Ls, Rs>...>)
		{
			constexpr int ind_N = symphas::lib::index_of_value<int, I0, I0s...>;

			using v_type = expr::symbols::v_<expr::symbols::i_<I0, P0>>;

			auto arg = get_arg(std::get<size_t(ind_N)>(series.substitution), symbol.index);

			auto left_case = expr::transform::swap_grid<v_type>
				(SymbolicCase(expr::symbols::i_<I0, P0>{} = DynamicVariable<G0>{}), *static_cast<E0 const*>(&applied), arg);
			auto right_case = expr::transform::swap_grid<v_type>
				(SymbolicCase(expr::symbols::i_<I0, P0>{} != DynamicVariable<G0>{}), *static_cast<E0 const*>(&applied), arg);

			auto left = expr::transform::swap_grid<SymbolicCaseSwap<>>(left_case, OpVoid{});
			auto right = expr::transform::swap_grid<SymbolicCaseSwap<>>(right_case, OpVoid{});

			auto series_left = expr::recreate_series(expr::symbols::i_<I0, P0>{} = expr::symbols::placeholder_N{}, left, series);
			auto series_right = OpVoid{};// expr::recreate_series(expr::symbols::i_<I0, 0>{} != expr::val<arg_N>, expr_right, series);
			return expr::transform::swap_grid<expr::symbols::placeholder_N_symbol>(series_left + series_right, symbol.index);
		}

		template<typename E, typename... Ts, int... I0s, int... P0s, typename B, typename C,
			typename E0, int I0, int P0, typename G0, typename... Ls, typename... Rs>
		auto apply_operators_sum_case(
			SymbolicSum<E, Substitution<Ts...>, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B, C> const& series,
			OpExpression<E0> const& applied,
			SymbolicFunctionalDerivative<G0> const& symbol,
			symphas::lib::types_list<SymbolicTernaryCase<expr::symbols::i_<I0, P0>, G0, Ls, Rs>...> cases)
		{
			return apply_operators_sum_case(series, *static_cast<E0 const*>(&applied), SymbolicDerivative<G0>{}, cases);
		}

		template<typename E, typename... Ts, int... I0s, int... P0s, typename B, typename C,
			typename E0, int I0, int P0, typename G0, typename... Ls, typename... Rs>
		auto apply_operators_sum_case(
			SymbolicSum<E, Substitution<Ts...>, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B, C> const& series,
			OpExpression<E0> const& applied,
			SymbolicFunctionalDerivative<DynamicVariable<G0>> const& symbol,
			symphas::lib::types_list<SymbolicTernaryCase<expr::symbols::i_<I0, P0>, DynamicVariable<G0>, Ls, Rs>...> cases)
		{
			return apply_operators_sum_case(series, *static_cast<E0 const*>(&applied), SymbolicDerivative<DynamicVariable<G0>>(symbol.index), cases);
		}

		template<typename E, typename... Ts, int... I0s, int... P0s, typename B, typename C,
			typename E0, int I0, int P0, size_t Z0, typename G0, typename... Ls, typename... Rs>
		auto apply_operators_sum_case(
			SymbolicSum<E, Substitution<Ts...>, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, B, C> const& series,
			OpExpression<E0> const& applied,
			SymbolicFunctionalDerivative<Variable<Z0, G0>> const& symbol,
			symphas::lib::types_list<SymbolicTernaryCase<expr::symbols::i_<I0, P0>, Variable<Z0, G0>, Ls, Rs>...> cases)
		{
			return apply_operators_sum_case(series, *static_cast<E0 const*>(&applied), SymbolicDerivative<Variable<Z0, G0>>{}, cases);
		}

		template<typename E, typename T, typename Seq, typename B, typename C, typename E0, typename G>
		auto apply_operators_sum(
			SymbolicSum<E, T, Seq, B, C> const& series,
			OpExpression<E0> const& applied,
			SymbolicDerivative<G> const& symbol,
			symphas::lib::types_list<>)
		{
			return *static_cast<E0 const*>(&applied);
		}

		//template<typename V, typename E, typename... Ts,
		//	int... I0s, int... P0s, typename A, typename B, typename C, typename... Rest>
		//auto apply_operators_sum_recurse(
		//	OpSum<V, E,
		//		Substitution<SymbolicDataArray<Ts>...>,
		//		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, C> const& sum,
		//	symphas::lib::types_list<Rest...>);

		template<typename E, typename T, typename Seq, typename B, typename C,
			typename E0, typename G, typename... Cs, typename... Rest>
		auto apply_operators_sum(
			SymbolicSum<E, T, Seq, B, C> const& series,
			OpExpression<E0> const& applied,
			SymbolicDerivative<G> const& symbol,
			symphas::lib::types_list<symphas::lib::types_list<Cs...>, Rest...>)
		{
			return (
				apply_operators_sum_case(series, *static_cast<E0 const*>(&applied), symbol, symphas::lib::types_list<Cs...>{})
				+ ... + apply_operators_sum_case(series, *static_cast<E0 const*>(&applied), symbol, Rest{}));
		}

		template<typename E, typename T, typename Seq, typename B, typename C, typename E0, typename G>
		auto apply_operators_sum(
			SymbolicSum<E, T, Seq, B, C> const& series,
			OpExpression<E0> const& applied,
			SymbolicDerivative<G> const& symbol)
		{
			using ops_t = typename filter_for_case<expr::op_types_t<E0>>::type;
			return apply_operators_sum(series, *static_cast<E0 const*>(&applied), symbol, ops_t{});
		}


		//template<typename V, typename E, typename... Ts,
		//	int... I0s, int... P0s, typename A, typename B, typename C, typename... Rest>
		//auto apply_operators_sum_recurse(
		//	OpSum<V, E,
		//		Substitution<SymbolicDataArray<Ts>...>,
		//		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, C> const& sum,
		//	symphas::lib::types_list<Rest...>)
		//{
		//	auto expr = expr::coeff(sum) * sum.data.substitute_placeholders(sum.f);
		//	return apply_operators_sum(sum.data, expr, symphas::lib::types_list<Rest...>{});
		//}

		template<typename... Vs, size_t O, typename V, typename V0, typename E, typename... Ts,
			int... I0s, int... P0s, typename A, typename B, typename C, typename... GGs, typename... Gs, size_t... Ns>
		auto apply_operators_sum_multiple(
			std::index_sequence<O>, V const& coeff,
            OpSum<V0, E,
				Substitution<SymbolicDataArray<Ts>...>,
				symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, A, B, C> const& sum,
			std::tuple<SymbolicDerivative<GGs>...> const& symbols,
			std::tuple<Gs...> const& args, std::index_sequence<Ns...>)
		{
			auto expr = sum.data.substitute_placeholders(sum.f);
			return coeff * expr::coeff(sum) * 
				(expr::transform::swap_grid<Vs>(
					apply_operators_sum(sum.data, apply_operators(expr::make_derivative<O, GGs>(expr, std::get<Ns>(symbols))), std::get<Ns>(symbols)), 
					std::get<Ns>(args))
					+ ...);
		}

        template<size_t O, typename V, typename V0, typename E, typename... Ts,
            int... I0s, int... P0s, typename A, typename B, typename... Vs, size_t Z, typename GG>
        auto apply_operators_sum(
            std::index_sequence<O>, V const& coeff, 
            OpSum<V0, E,
                Substitution<SymbolicDataArray<Ts>...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
                A, B, symphas::lib::types_list<Vs...>> const& sum, 
            SymbolicDerivative<Variable<Z, GG>> const& symbol)
        {
			auto expr = sum.data.substitute_placeholders(sum.f);
			return coeff * expr::coeff(sum) * apply_operators_sum(sum.data, apply_operators(expr::make_derivative<O, Variable<Z, GG>>(expr, symbol)), symbol);
        }

        template<size_t O, typename V, typename V0, typename E, typename... Ts,
            int... I0s, int... P0s, typename A, typename B, typename... Vs, typename GG>
        auto apply_operators_sum(
            std::index_sequence<O>, V const& coeff, 
            OpSum<V0, E,
                Substitution<SymbolicDataArray<Ts>...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
                A, B, symphas::lib::types_list<Vs...>> const& sum, 
            SymbolicDerivative<DynamicVariable<GG>> const& symbol)
        {
			auto expr = sum.data.substitute_placeholders(sum.f);
			return coeff * expr::coeff(sum) * apply_operators_sum(sum.data, apply_operators(expr::make_derivative<O, DynamicVariable<GG>>(expr, symbol)), symbol);
        }

		template<typename T>
		struct is_variational_trait
		{
			static const bool value = false;
		};

		template<typename T>
		struct is_variational_trait<expr::variational_t<T>>
		{
			static const bool value = true;
		};

        template<size_t O, typename V, typename V0, typename E, typename... Ts,
            int... I0s, int... P0s, typename A, typename B, typename... Vs, typename GG,
            std::enable_if_t<(!expr::is_id_variable<GG> && !expr::is_functional_derivative<SymbolicDerivative<GG>>), int> = 0>
        auto apply_operators_sum(
            std::index_sequence<O>, V const& coeff, 
            OpSum<V0, E,
                Substitution<SymbolicDataArray<Ts>...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
                A, B, symphas::lib::types_list<Vs...>> const& sum, 
            SymbolicDerivative<GG> const& symbol)
        {
            if constexpr (expr::is_expression<GG>)
            {
                using vs_seq = vars_in_ops_t<op_types_t<GG>>;
                constexpr size_t L = vs_seq::size();

                // if multiple variables are in the expression with respect to differentiate
                if constexpr (L > 1)
                {
                    return OpVoid{};
                }
                else
                {
					constexpr size_t D = expr::grid_dim<GG>::value;
					if constexpr (vs_seq::size() > 0)
					{
						constexpr size_t Z0 = symphas::lib::seq_index_value<0, vs_seq>::value;

						auto dvs = std::make_tuple(
							SymbolicDerivative{ expr::transform::swap_grid<Z0>(GG{}, GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>, D>{}) }...);
						auto args = get_args<Z0>(sum.data.substitution, std::make_index_sequence<sizeof...(Ts)>{});
						return apply_operators_sum_multiple<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>...>(
							std::index_sequence<O>{}, coeff, sum, dvs, args, std::make_index_sequence<sizeof...(Ts)>{});
					}
					else
					{
						auto dvs = std::make_tuple(
							SymbolicDerivative{ expr::transform::swap_grid(op_types_t<GG>{}, GG{}, GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>, D>{}) }...);
						//auto dvs = symphas::lib::types_list<
						//	decltype(expr::transform::swap_grid(op_types_t<GG>{}, std::declval<GG>(), GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>, D>{}))... > {};
						auto args = get_args(sum.data.substitution, std::make_index_sequence<sizeof...(Ts)>{});
						return apply_operators_sum_multiple<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>...>(
							std::index_sequence<O>{}, coeff, sum, dvs, args, std::make_index_sequence<sizeof...(Ts)>{});
					}
                }
            }
            else
            {
				auto expr = sum.data.substitute_placeholders(sum.f);
				return coeff * expr::coeff(sum) * apply_operators_sum(sum.data, apply_operators(expr::make_derivative<O, GG>(expr, symbol)), symbol);
            }
        }
		
        template<size_t O, typename V, typename V0, typename E, typename... Ts,
            int... I0s, int... P0s, typename A, typename B, typename... Vs, size_t Z, typename GG>
        auto apply_operators_sum(
            std::index_sequence<O>, V const& coeff, 
            OpSum<V0, E,
                Substitution<SymbolicDataArray<Ts>...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
                A, B, symphas::lib::types_list<Vs...>> const& sum, 
            SymbolicFunctionalDerivative<Variable<Z, GG>> const& symbol)
        {
			constexpr size_t D = expr::grid_dim<GG>::value;
			auto symbols = std::make_tuple(SymbolicFunctionalDerivative<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>, D>>()...);
            auto args = get_args<Z>(sum.data.substitution, std::make_index_sequence<sizeof...(Ts)>{});

            return apply_operators_sum_multiple<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>...>(
                 std::index_sequence<O>{}, coeff, sum, symbols, args, std::make_index_sequence<sizeof...(Ts)>{});
        }
		
        template<size_t O, typename V, typename V0, typename E, typename... Ts,
            int... I0s, int... P0s, typename A, typename B, typename... Vs, typename GG>
        auto apply_operators_sum(
            std::index_sequence<O>, V const& coeff, 
            OpSum<V0, E,
                Substitution<SymbolicDataArray<Ts>...>,
                symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
                A, B, symphas::lib::types_list<Vs...>> const& sum, 
            SymbolicFunctionalDerivative<DynamicVariable<GG>> const& symbol)
        {
			constexpr size_t D = expr::grid_dim<GG>::value;
			auto symbols = std::make_tuple(SymbolicFunctionalDerivative<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>, D>>(symbol.index)...);
            auto args = get_args(sum.data.substitution, symbol.index, std::make_index_sequence<sizeof...(Ts)>{});

            auto result = apply_operators_sum_multiple<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>...>(
                std::index_sequence<O>{}, coeff, sum, symbols, args, std::make_index_sequence<sizeof...(Ts)>{});
			return expr::transform::swap_grid<expr::symbols::placeholder_N_symbol_<0>, OpCoeffSwap<expr::symbols::placeholder_N_symbol_<0>>>
				(result, symbol.index, symbol.index);
        }

	}

	template<size_t O, typename V, typename V0, typename E, typename... Ts,
		int... I0s, int... P0s, typename A, typename B, typename... Vs, typename GG>
	auto apply_operators(
		OpDerivative<std::index_sequence<O>, V, OpSum<V0, E,
			Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
            A, B, symphas::lib::types_list<Vs...>>, SymbolicDerivative<GG>> const& e)
	{
        return apply_operators_sum(std::index_sequence<O>{}, expr::coeff(e), expr::get_enclosed_expression(e), e.solver);
    }

	template<size_t O, typename V, size_t Z, int I0, int P0, typename V0, typename E, typename... Ts,
		int... I0s, int... P0s, typename A, typename B, typename... Vs>
	auto apply_operators(
		OpDerivative<std::index_sequence<O>, V, OpSum<V0, E,
			Substitution<SymbolicDataArray<Ts>...>,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
            A, B, symphas::lib::types_list<Vs...>>, 
		SymbolicDerivative<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>>> const& e)
	{
		auto sum = expr::get_enclosed_expression(e);
		return expr::coeff(e) * apply_operators_sum(sum.data, apply_operators(expr::make_derivative<O, expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>>(sum.data.e)), e.solver);
	}

	template<typename V, typename E, typename... Ts,
		int... I0s, int... P0s, typename A, typename B, typename... Vs>
	auto apply_operators(OpSum<V, E,
		Substitution<SymbolicDataArray<Ts>...>,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		A, B, symphas::lib::types_list<Vs...>> const& sum)
	{
		auto expr = apply_operators(sum.data.substitute_placeholders(sum.f));
		return expr::coeff(sum) * expr::recreate_series(expr, sum.data);
	}

}




namespace symphas::internal
{
	template<size_t I, typename T>
	using itype = T;
	template<typename I, typename T>
	using ttype = T;
	template<size_t I, typename T>
	constexpr size_t typei = I;

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
		// in the case when the expressions don't have tensors, but are simply scalar
		if constexpr (R1 == 0 && R2 == 0 && Q1 == 0 && Q2 == 0)
		{
			return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
		}
		// in the special case when 1D tensors are being multiplied
		else if constexpr (R1 == 1 && Q2 == 1 && Q1 == 1 && R1 == 1)
		{
			return (symphas::internal::tensor_cancel{} *(*static_cast<E1 const*>(&a)))* (symphas::internal::tensor_cancel{} *(*static_cast<E2 const*>(&b)));
		}
		// multiply a row type by a column type
		else if constexpr (R1 == 1 && Q2 == 1 && Q1 == R2 && R2 > 1)
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
			return expr::make_mul(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
		}
	}

	template<typename E1, typename E2>
	auto dot(OpExpression<E1> const& a, OpOperator<E2> const& b)
	{
		return expr::make_mul(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename E1, typename E2>
	auto dot(OpOperator<E1> const& a, OpExpression<E2> const& b)
	{
		return expr::make_mul(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename E1, typename E2>
	auto dot(OpOperator<E1> const& a, OpOperator<E2> const& b)
	{
		return expr::make_mul(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename E1, typename E2>
	auto dot(OpOperatorChain<OpIdentity, E1> const& a, OpOperatorChain<OpIdentity, E2> const& b)
	{
		return dot(a.g, b.g);
	}

	template<typename E1>
	auto dot(OpExpression<E1> const& a, OpVoid)
	{
		return OpVoid{};
	}

	template<typename E2>
	auto dot(OpVoid, OpExpression<E2> const& b)
	{
		return OpVoid{};
	}

	inline auto dot(OpVoid, OpVoid)
	{
		return OpVoid{};
	}
}

namespace expr
{
	template<typename E1, typename E2>
	auto dot(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return symphas::internal::dot(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename E1, typename E2>
	auto dot(OpExpression<E1> const& a, OpOperator<E2> const& b)
	{
		return symphas::internal::dot(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename E1, typename E2>
	auto dot(OpOperator<E1> const& a, OpExpression<E2> const& b)
	{
		return symphas::internal::dot(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename E1, typename E2>
	auto dot(OpOperator<E1> const& a, OpOperator<E2> const& b)
	{
		return symphas::internal::dot(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename T, size_t N, size_t D>
	auto transpose(OpTensor<T, N, D> const& tensor)
	{
		return expr::make_tensor<0, D, 1, N>(T(tensor));
	}	

	template<typename T, size_t N0, size_t N1, size_t D0, size_t D1>
	auto transpose(OpTensor<T, N0, N1, D0, D1> const& tensor)
	{
		return expr::make_tensor<N1, N0, D1, D0>(T(tensor));
	}

	template<size_t R1, size_t R2, size_t P, typename E, size_t... Qs>
	auto compute_transpose_row(OpExpression<E> const& e, std::index_sequence<Qs...>)
	{
		return ((expr::make_tensor<Qs, P, R2, R1>() * (expr::make_row_vector<P, R1>() * (*static_cast<E const*>(&e) * expr::make_column_vector<Qs, R2>()))) + ...);
	}

	template<size_t R1, size_t R2, typename E, size_t... Ps, size_t... Qs>
	auto compute_transpose(OpExpression<E> const& e, std::index_sequence<Ps...>, std::index_sequence<Qs...>)
	{
		return (compute_transpose_row<R1, R2, Ps>(*static_cast<E const*>(&e), std::index_sequence<Qs...>{}) + ...);
	}

	template<typename E, size_t R1 = expr::eval_type<E>::rank, size_t R2 = expr::eval_type<E>::template rank_<1>>
	auto transpose(OpExpression<E> const& e)
	{
		if constexpr ((R1 == 0 && R2 == 0) || (R1 == 1 && R2 == 1))
		{
			return *static_cast<E const*>(&e);
		}
		else
		{
			return compute_transpose<R1, R2>(*static_cast<E const*>(&e), std::make_index_sequence<R1>{}, std::make_index_sequence<R2>{});
		}
	}
}

namespace symphas::internal
{

	template<size_t Z, typename G, typename... Gs>
	auto filter_variables(std::tuple<Variable<Z, G>, Gs...> const& data_list);
	template<size_t Z0, typename G0>
	auto sort_variables(std::tuple<Variable<Z0, G0>> const& data_list);

	inline auto filter_variables(std::tuple<> const& data_list)
	{
		return std::make_tuple();
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
		return std::make_tuple();
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
		return std::make_tuple();
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
		else if constexpr (Z0 < I0)
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

	template<size_t Z0, typename G0, size_t... Zs, typename... Gs>
	auto index_variables(std::tuple<Variable<Z0, G0>, Variable<Zs, Gs>...> const& data_list)
	{
		constexpr size_t Zm = symphas::lib::seq_index_value<sizeof...(Zs), std::index_sequence<Z0, Zs...>>::value;
		return index_variables(data_list, std::make_index_sequence<Zm + 1>{});

	}

	inline auto index_variables(std::tuple<> const& data_list)
	{
		return std::make_tuple();

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

template<typename V0, typename... G0s, expr::exp_key_t... X0s, typename V1,
	typename W0, typename... G1s, expr::exp_key_t... X1s, typename W1,
	typename Dd, typename E, typename Sp,
	typename = std::enable_if_t<(
		expr::is_combinable<G0s...>
		&& symphas::lib::filter_seq_t<std::integer_sequence<expr::exp_key_t, X0s...>, std::integer_sequence<expr::exp_key_t, X1s...>>::size() == 0
		&& symphas::lib::filter_seq_t<std::integer_sequence<expr::exp_key_t, X1s...>, std::integer_sequence<expr::exp_key_t, X0s...>>::size() == 0
		&& ((symphas::lib::index_of_type<G0s, G1s...> >= 0) && ...)
		&& ((symphas::lib::index_of_type<G1s, G0s...> >= 0) && ...)
		), int>>
auto operator+(OpBinaryMul<OpTerms<V0, Term<G0s, X0s>...>, OpDerivative<Dd, V1, E, Sp>> const& a,
	OpBinaryMul<OpTerms<W0, Term<G1s, X1s>...>, OpDerivative<Dd, W1, E, Sp>> const& b)
{
	return (expr::coeff(b.a) * expr::coeff(b.b) + expr::coeff(a.a) * expr::coeff(a.b))
		* (OpTerms<OpIdentity, Term<G0s, X0s>...>(OpIdentity{}, expr::terms_after_first(a.a)) * expr::make_derivative<Dd>(expr::get_enclosed_expression(a.b), a.b.solver));
}

template<typename V0, typename... G0s, expr::exp_key_t... X0s, typename V1,
	typename W0, typename... G1s, expr::exp_key_t... X1s, typename W1,
	typename Dd, typename E, typename Sp,
	typename = std::enable_if_t<(
		expr::is_combinable<G0s...>
		&& symphas::lib::filter_seq_t<std::integer_sequence<expr::exp_key_t, X0s...>, std::integer_sequence<expr::exp_key_t, X1s...>>::size() == 0
		&& symphas::lib::filter_seq_t<std::integer_sequence<expr::exp_key_t, X1s...>, std::integer_sequence<expr::exp_key_t, X0s...>>::size() == 0
		&& ((symphas::lib::index_of_type<G0s, G1s...> >= 0) && ...)
		&& ((symphas::lib::index_of_type<G1s, G0s...> >= 0) && ...)
		), int>>
auto operator-(OpBinaryMul<OpTerms<V0, Term<G0s, X0s>...>, OpDerivative<Dd, V1, E, Sp>> const& a,
	OpBinaryMul<OpTerms<W0, Term<G1s, X1s>...>, OpDerivative<Dd, W1, E, Sp>> const& b)
{
	return (expr::coeff(b.a) * expr::coeff(b.b) - expr::coeff(a.a) * expr::coeff(a.b))
		* (OpTerms<OpIdentity, Term<G0s, X0s>...>(OpIdentity{}, expr::terms_after_first(a.a)) * expr::make_derivative<Dd>(expr::get_enclosed_expression(a.b), a.b.solver));
}


namespace expr::prune
{
	namespace
	{
		inline void _update(...) {}

		template<typename E, typename... condition_ts>
		inline void _update(OpExpression<E>& e, symphas::lib::types_list<condition_ts...>) {}

		template<typename T, typename... condition_ts>
		inline void _update(Block<T>&, symphas::lib::types_list<condition_ts...>) {}
		template<size_t N, typename T, typename... condition_ts>
		inline void _update(MultiBlock<N, T>&, symphas::lib::types_list<condition_ts...>) {}

		template<typename V, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpExponential<V, E>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpExponential<V, E>& e, symphas::lib::types_list<condition_ts...>);
		template<expr::exp_key_t X, typename V, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpPow<X, V, E>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpPow<X, V, E>& e, symphas::lib::types_list<condition_ts...>);
		template<typename V, typename E1, typename E2, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpConvolution<V, E1, E2>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpConvolution<V, E1, E2>& e, symphas::lib::types_list<condition_ts...>);
		template<typename V, typename E, typename F, typename... Args, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpFunction<V, E, F, Args...>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpFunction<V, E, F, Args...>& e, symphas::lib::types_list<condition_ts...>);
		template<auto f, typename V, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpFunctionApply<f, V, E>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpFunctionApply<f, V, E>& e, symphas::lib::types_list<condition_ts...>);
		template<typename V, size_t D, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpConvolution<V, GaussianSmoothing<D>, E>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpConvolution<V, GaussianSmoothing<D>, E>& e, symphas::lib::types_list<condition_ts...>);
		template<typename V, size_t D, typename G, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>& e, symphas::lib::types_list<condition_ts...>);
		template<size_t O, typename V, typename E, typename G, typename... condition_ts>
		inline void _update(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G>>& e, symphas::lib::types_list<condition_ts...>);
		template<typename Dd, typename V, typename G, typename Sp, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>& e, symphas::lib::types_list<condition_ts...>);
		template<typename Dd, typename V, typename E, typename Sp, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpDerivative<Dd, V, E, Sp>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpDerivative<Dd, V, E, Sp>& e, symphas::lib::types_list<condition_ts...>);
		template<typename V, typename E, typename T, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpIntegral<V, E, T>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpIntegral<V, E, T>& e, symphas::lib::types_list<condition_ts...>);
		template<typename... Es, typename... condition_ts, size_t... Is,
			std::enable_if_t<expr::satisfies<OpAddList<Es...>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpAddList<Es...>& e, symphas::lib::types_list<condition_ts...>, std::index_sequence<Is...>);
		template<typename... Es, typename... condition_ts,
			std::enable_if_t<(expr::satisfies<Es, expr::or_<condition_ts...>> || ...), int> = 0>
		inline void _update(OpAdd<Es...>& e, symphas::lib::types_list<condition_ts...>);
		template<typename A1, typename A2, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpBinaryMul<A1, A2>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpBinaryMul<A1, A2>& e, symphas::lib::types_list<condition_ts...>);
		template<typename A1, typename A2, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpBinaryDiv<A1, A2>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpBinaryDiv<A1, A2>& e, symphas::lib::types_list<condition_ts...>);
		template<typename A1, typename A2, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpOperatorCombination<A1, A2>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpOperatorCombination<A1, A2>& e, symphas::lib::types_list<condition_ts...>);
		template<typename A1, typename A2, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpCombination<A1, A2, E>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpCombination<A1, A2, E>& e, symphas::lib::types_list<condition_ts...>);
		template<typename A1, typename A2, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpOperatorChain<A1, A2>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpOperatorChain<A1, A2>& e, symphas::lib::types_list<condition_ts...>);
		template<typename A1, typename A2, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpChain<A1, A2, E>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpChain<A1, A2, E>& e, symphas::lib::types_list<condition_ts...>);
		template<typename G, typename V, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpMap<G, V, E>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpMap<G, V, E>& e, symphas::lib::types_list<condition_ts...>);
		template<typename... Ts, typename... condition_ts, size_t... Is,
			std::enable_if_t<expr::satisfies<OpTermsList<Ts...>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpTermsList<Ts...>& e, symphas::lib::types_list<condition_ts...>, std::index_sequence<Is...>);
		template<typename... Ts, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpTerms<Ts...>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpTerms<Ts...>& e, symphas::lib::types_list<condition_ts...>);
		template<typename G, expr::exp_key_t X, typename... condition_ts,
			std::enable_if_t<expr::satisfies<Term<G, X>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(Term<G, X>& e, symphas::lib::types_list<condition_ts...>);
		template<typename V, typename sub_t, typename eval_t, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpSymbolicEval<V, sub_t, eval_t>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpSymbolicEval<V, sub_t, eval_t>& e, symphas::lib::types_list<condition_ts...>);
		template<typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpOptimized<E>, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(OpOptimized<E> &e, symphas::lib::types_list<condition_ts...>);
		template<typename E, typename... Ts, typename... condition_ts,
			std::enable_if_t<expr::satisfies<E, expr::or_<condition_ts...>>, int> = 0>
		inline void _update(SymbolicFunction<E, Ts...>& e, symphas::lib::types_list<condition_ts...>);
		template<typename... Ts, typename... condition_ts>
		inline void _update(SymbolicCase<Ts...>& e, symphas::lib::types_list<condition_ts...>);

		//template<typename Op, typename E, typename Inds>
		//inline void _update(SymbolicSeries<Op, E, Inds>& e);
		//template<typename Op, typename... Ts, typename E, int... I0s, int... P0s,
		//	typename... T1s, typename... T2s, typename... Is>
		//inline void _update(SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
		//	symphas::lib::types_list<E,
		//		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		//		symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
		//		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>
		//	>>& e);

		//template<typename... Ts>
		//inline void _update(SymbolicCase<Ts...>& e);
		//template<size_t Z, typename G>
		//inline void _update(Variable<Z, G>& e);
		//template<typename G>
		//inline void _update(NamedData<G>& e);
		//template<Axis ax, typename G>
		//inline void _update(VectorComponent<ax, G>& e);
		//template<typename G>
		//inline void _update(symphas::ref<G>& e);


		template<expr::NoiseType nt, typename T, size_t D, typename... condition_ts>
		inline void _update(NoiseData<nt, T, D>& e, symphas::lib::types_list<condition_ts...>);
		template<typename T, typename... condition_ts>
		inline void _update(SymbolicData<T>&, symphas::lib::types_list<condition_ts...>);

		template<typename V, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpExponential<V, E>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpExponential<V, E>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(expr::get_enclosed_expression(e), symphas::lib::types_list<condition_ts...>{});
		}

		template<expr::exp_key_t X, typename V, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpPow<X, V, E>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpPow<X, V, E>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(expr::get_enclosed_expression(e), symphas::lib::types_list<condition_ts...>{});
		}

		template<typename V, typename E1, typename E2, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpConvolution<V, E1, E2>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpConvolution<V, E1, E2>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(e.a, symphas::lib::types_list<condition_ts...>{});
			_update(e.b, symphas::lib::types_list<condition_ts...>{});
			e.update(symphas::lib::types_list<condition_ts...>{});
		}

		template<typename V, typename E, typename F, typename... Args, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpFunction<V, E, F, Args...>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpFunction<V, E, F, Args...>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(expr::get_enclosed_expression(e), symphas::lib::types_list<condition_ts...>{});
		}

		template<auto f, typename V, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpFunctionApply<f, V, E>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpFunctionApply<f, V, E>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(expr::get_enclosed_expression(e), symphas::lib::types_list<condition_ts...>{});
		}

		template<typename V, size_t D, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpConvolution<V, GaussianSmoothing<D>, E>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpConvolution<V, GaussianSmoothing<D>, E>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(expr::get_enclosed_expression(e), symphas::lib::types_list<condition_ts...>{});
			e.update(symphas::lib::types_list<condition_ts...>{});
		}

		template<typename V, size_t D, typename G, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpConvolution<V, GaussianSmoothing<D>, OpTerm<OpIdentity, G>>& e, symphas::lib::types_list<condition_ts...>)
		{
			e.update(symphas::lib::types_list<condition_ts...>{});
		}


		/* derivative pruning
			*/

		template<size_t O, typename V, typename E, typename G, typename... condition_ts>
		inline void _update(OpDerivative<std::index_sequence<O>, V, E, SymbolicDerivative<G>>& e, symphas::lib::types_list<condition_ts...>) {}

		template<typename Dd, typename V, typename G, typename Sp, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp>& e, symphas::lib::types_list<condition_ts...>)
		{
			e.update(symphas::lib::types_list<condition_ts...>{});
		}

		template<typename Dd, typename V, typename E, typename Sp, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpDerivative<Dd, V, E, Sp>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpDerivative<Dd, V, E, Sp>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(expr::get_enclosed_expression(e), symphas::lib::types_list<condition_ts...>{});
			e.update(symphas::lib::types_list<condition_ts...>{});
		}


		template<typename V, typename E, typename T, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpIntegral<V, E, T>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpIntegral<V, E, T>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(expr::get_enclosed_expression(e), symphas::lib::types_list<condition_ts...>{});
			e.update(symphas::lib::types_list<condition_ts...>{});
		}

		/* binary op pruning
			*/

		template<typename... Es, typename... condition_ts, size_t... Is,
			std::enable_if_t<expr::satisfies<OpAddList<Es...>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpAddList<Es...>& e, symphas::lib::types_list<condition_ts...>, std::index_sequence<Is...>)
		{
			(_update(expr::get<Is>(e), symphas::lib::types_list<condition_ts...>{}), ...);
		}

		template<typename... Es, typename... condition_ts,
			std::enable_if_t<(expr::satisfies<Es, expr::or_<condition_ts...>> || ...), int>>
		inline void _update(OpAdd<Es...>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(*static_cast<OpAddList<Es...>*>(&e), symphas::lib::types_list<condition_ts...>{}, std::make_index_sequence<sizeof...(Es)>{});
		}

		template<typename A1, typename A2, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpBinaryMul<A1, A2>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpBinaryMul<A1, A2>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(e.a, symphas::lib::types_list<condition_ts...>{});
			_update(e.b, symphas::lib::types_list<condition_ts...>{});
		}

		template<typename A1, typename A2, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpBinaryDiv<A1, A2>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpBinaryDiv<A1, A2>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(e.a, symphas::lib::types_list<condition_ts...>{});
			_update(e.b, symphas::lib::types_list<condition_ts...>{});
		}

		/* operator pruning
			*/

		template<typename A1, typename A2, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpOperatorCombination<A1, A2>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpOperatorCombination<A1, A2>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(e.f, symphas::lib::types_list<condition_ts...>{});
			_update(e.g, symphas::lib::types_list<condition_ts...>{});
		}

		template<typename A1, typename A2, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpCombination<A1, A2, E>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpCombination<A1, A2, E>& e, symphas::lib::types_list<condition_ts...>)
		{
			//_update(e.combination);
			//_update(expr::get_enclosed_expression(e));
			//if constexpr (expr::has_state<E>::value)
			//{
			//}
			e.update(symphas::lib::types_list<condition_ts...>{});
		}

		template<typename A1, typename A2, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpOperatorChain<A1, A2>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpOperatorChain<A1, A2>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(e.f, symphas::lib::types_list<condition_ts...>{});
			_update(e.g, symphas::lib::types_list<condition_ts...>{});
		}

		template<typename A1, typename A2, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpChain<A1, A2, E>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpChain<A1, A2, E>& e, symphas::lib::types_list<condition_ts...>)
		{
			//_update(e.combination);
			//if constexpr (expr::has_state<E>::value)
			//{
			//	_update(expr::get_enclosed_expression(e));
			//}
			e.update(symphas::lib::types_list<condition_ts...>{});
		}

		template<typename G, typename V, typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpMap<G, V, E>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpMap<G, V, E>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(expr::get_enclosed_expression(e), symphas::lib::types_list<condition_ts...>{});
			e.update(symphas::lib::types_list<condition_ts...>{});
		}

		template<typename... Ts, typename... condition_ts, size_t... Is,
			std::enable_if_t<expr::satisfies<OpTermsList<Ts...>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpTermsList<Ts...>& e, symphas::lib::types_list<condition_ts...>, std::index_sequence<Is...>)
		{
			(_update(expr::get<Is>(e), symphas::lib::types_list<condition_ts...>{}), ...);
		}

		template<typename... Ts, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpTerms<Ts...>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpTerms<Ts...>& e, symphas::lib::types_list<condition_ts...>)
		{
			//_update(*static_cast<OpTermsList<Ts...>*>(&e), std::make_index_sequence<sizeof...(Ts)>{});
		}


		template<typename G, expr::exp_key_t X, typename... condition_ts,
			std::enable_if_t<expr::satisfies<Term<G, X>, expr::or_<condition_ts...>>, int>>
		inline void _update(Term<G, X>& e, symphas::lib::types_list<condition_ts...>)
		{
			//_update(e.data());
		}

		template<typename V, typename sub_t, typename eval_t, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpSymbolicEval<V, sub_t, eval_t>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpSymbolicEval<V, sub_t, eval_t>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(e.f, symphas::lib::types_list<condition_ts...>{});
			e.update(symphas::lib::types_list<condition_ts...>{});
		}

		template<typename E, typename... condition_ts,
			std::enable_if_t<expr::satisfies<OpOptimized<E>, expr::or_<condition_ts...>>, int>>
		inline void _update(OpOptimized<E>& e, symphas::lib::types_list<condition_ts...>)
		{
			e.update(symphas::lib::types_list<condition_ts...>{});
		}

		template<typename E, typename... Ts, typename... condition_ts,
			std::enable_if_t<expr::satisfies<E, expr::or_<condition_ts...>>, int>>
		inline void _update(SymbolicFunction<E, Ts...>& e, symphas::lib::types_list<condition_ts...>)
		{
			_update(e.e, symphas::lib::types_list<condition_ts...>{});
		}

		//template<typename Op, typename E, typename Inds>
		//inline void _update(SymbolicSeries<Op, E, Inds>& e) {}

		//template<typename Op, typename... Ts, typename E, int... I0s, int... P0s,
		//	typename... T1s, typename... T2s, typename... Is>
		//inline void _update(SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
		//	symphas::lib::types_list<E,
		//		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		//		symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
		//		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>
		//	>>& e)
		//{
		//	//e.update();
		//}

		template<typename... Ts, typename... condition_ts>
		inline void _update(SymbolicCase<Ts...>& e, symphas::lib::types_list<condition_ts...>) {}

		//template<size_t Z, typename G>
		//inline void _update(Variable<Z, G>& e)
		//{
		//    _update(*static_cast<G*>(&e));
		//}

		//template<typename G>
		//inline void _update(NamedData<G>& e)
		//{
		//    _update(*static_cast<G*>(&e));
		//}

		//template<Axis ax, typename G>
		//inline void _update(VectorComponent<ax, G>& e)
		//{
		//    _update(*static_cast<G*>(&e));
		//}

		//template<typename G>
		//inline void _update(symphas::ref<G>& e)
		//{
		//    _update(e.get());
		//}

		//template<expr::NoiseType nt, typename T, size_t D>
		//inline void _update(NoiseData<nt, T, D>& data)
		//{
		//    data.update();
		//}

	//     template<typename T>
	//     inline void _update(SymbolicData<T>& data)
	//     {
			//if (data.data != nullptr)
			//{
			//	_update(*data.data);
			//}
	//     }

	}

	//! Update underlying the given expression.
	/*!
		* For expressions which store intermediate data, such as derivatives, this
		* data must be updated before the expression can be evaluated. This will
		* traverse the expression tree and perform all necessary updating.
		*
		* \param e The expression that is updated.
		*/
	template<typename... condition_ts, typename E>
	inline void update(E& e)
	{
		TIME_THIS_CONTEXT_LIFETIME(expression_update)
		_update(e, symphas::lib::types_list<condition_ts...>{});
	}

}

namespace expr
{

	template<typename A, typename B>
	struct numeric_range_dist
	{
		using type = sub_result_t<B, A>;
	};

	template<size_t N, size_t D, size_t N0, size_t D0, bool first_neg, bool second_neg>
	struct numeric_range_count_impl;

	template<size_t N, size_t D, size_t N0, size_t D0>
	struct numeric_range_count_impl<N, D, N0, D0, false, false>
	{
		using type = sub_result_t<OpFractionLiteral<N0, D0>, OpFractionLiteral<N, D>>;
	};

	template<size_t N, size_t D, size_t N0, size_t D0>
	struct numeric_range_count_impl<N, D, N0, D0, false, true>
	{
		using type = sub_result_t<OpNegFractionLiteral<N0, D0>, OpFractionLiteral<N, D>>;
	};
	
	template<size_t N, size_t D, size_t N0, size_t D0>
	struct numeric_range_count_impl<N, D, N0, D0, true, false>
	{
		using type = sub_result_t<OpFractionLiteral<N0, D0>, OpNegFractionLiteral<N, D>>;
	};

	template<size_t N, size_t D, size_t N0, size_t D0>
	struct numeric_range_count_impl<N, D, N0, D0, true, true>
	{
		using type = sub_result_t<OpNegFractionLiteral<N0, D0>, OpNegFractionLiteral<N, D>>;
	};

	template<typename A>
	struct numeric_range_value
	{
		constexpr static int value = A{}.eval();
	};


	template<int N, size_t D, int N0, size_t D0>
	struct numeric_range_count
	{
		using type = typename numeric_range_count_impl<(N < 0) ? size_t(-N) : size_t(N), D, (N0 < 0) ? size_t(-N0) : size_t(N0), D0, N < 0, N0 < 0>::type;
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
	template<int N, size_t D, int N0, size_t D0, int C>
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
	template<int N, size_t D, int N0, size_t D0>
	struct numeric_interval_state
	{
		template<typename A, typename B>
		constexpr numeric_interval_state(A, B) {}


		constexpr auto operator&(OpIdentity) const
		{
			return numeric_range_state<N, D, N0, D0, 1>{};
		}

		constexpr auto operator&(OpNegIdentity) const
		{
			return numeric_range_state<N, D, N0, D0, -1>{};
		}

		template<size_t C>
		constexpr auto operator&(OpFractionLiteral<C, 1>) const
		{
			return numeric_range_state<N, D, N0, D0, int(C)>{};
		}

		template<size_t C>
		constexpr auto operator&(OpNegFractionLiteral<C, 1>) const
		{
			return numeric_range_state<N, D, N0, D0, -int(C)>{};
		}

		constexpr auto operator|(OpIdentity) const
		{
			constexpr int C = numeric_range_count<N, D, N0, D0>::value;
			return numeric_range_state<N, D, N0, D0, C>{};
		}

		constexpr auto operator|(OpNegIdentity) const
		{
			constexpr int C = numeric_range_count<N, D, N0, D0>::value;
			return numeric_range_state<N, D, N0, D0, C>{};
		}

		template<size_t dN, size_t dD>
		constexpr auto operator|(OpFractionLiteral<dN, dD>) const
		{
			constexpr int C = div_result_t<typename numeric_range_count<N, D, N0, D0>::type, OpFractionLiteral<dN, dD>>{}.eval();
			return numeric_range_state<N, D, N0, D0, C>{};
		}

		template<size_t dN, size_t dD>
		constexpr auto operator|(OpNegFractionLiteral<dN, dD>) const
		{
			constexpr int C = div_result_t<typename numeric_range_count<N, D, N0, D0>::type, OpFractionLiteral<dN, dD>>{}.eval();
			return numeric_range_state<N, D, N0, D0, C>{};
		}
	};


	numeric_interval_state(OpIdentity, OpIdentity) -> numeric_interval_state<1, 1, 1, 1>;
	numeric_interval_state(OpIdentity, OpNegIdentity) -> numeric_interval_state<1, 1, -1, 1>;
	numeric_interval_state(OpIdentity, OpVoid) -> numeric_interval_state<1, 1, 0, 1>;
	template<size_t N, size_t D>
	numeric_interval_state(OpIdentity, OpFractionLiteral<N, D>) -> numeric_interval_state<1, 1, int(N), D>;
	template<size_t N, size_t D>
	numeric_interval_state(OpIdentity, OpNegFractionLiteral<N, D>) -> numeric_interval_state<1, 1, -int(N), D>;

	numeric_interval_state(OpNegIdentity, OpIdentity)->numeric_interval_state<-1, 1, 1, 1>;
	numeric_interval_state(OpNegIdentity, OpNegIdentity)->numeric_interval_state<-1, 1, -1, 1>;
	numeric_interval_state(OpNegIdentity, OpVoid)->numeric_interval_state<-1, 1, 0, 1>;
	template<size_t N, size_t D>
	numeric_interval_state(OpNegIdentity, OpFractionLiteral<N, D>) -> numeric_interval_state<-1, 1, int(N), D>;
	template<size_t N, size_t D>
	numeric_interval_state(OpNegIdentity, OpNegFractionLiteral<N, D>) -> numeric_interval_state<-1, 1, -int(N), D>;

	numeric_interval_state(OpVoid, OpIdentity)->numeric_interval_state<0, 1, 1, 1>;
	numeric_interval_state(OpVoid, OpNegIdentity)->numeric_interval_state<0, 1, -1, 1>;
	numeric_interval_state(OpVoid, OpVoid)->numeric_interval_state<0, 1, 0, 1>;
	template<size_t N, size_t D>
	numeric_interval_state(OpVoid, OpFractionLiteral<N, D>) -> numeric_interval_state<0, 1, int(N), D>;
	template<size_t N, size_t D>
	numeric_interval_state(OpVoid, OpNegFractionLiteral<N, D>) -> numeric_interval_state<0, 1, -int(N), D>;

	template<size_t N, size_t D>
	numeric_interval_state(OpFractionLiteral<N, D>, OpIdentity) -> numeric_interval_state<int(N), D, 1, 1>;
	template<size_t N, size_t D>
	numeric_interval_state(OpFractionLiteral<N, D>, OpNegIdentity) -> numeric_interval_state<int(N), D, -1, 1>;
	template<size_t N, size_t D>
	numeric_interval_state(OpFractionLiteral<N, D>, OpVoid) -> numeric_interval_state<int(N), D, 0, 1>;
	template<size_t N, size_t D, size_t N0, size_t D0>
	numeric_interval_state(OpFractionLiteral<N, D>, OpFractionLiteral<N0, D0>) -> numeric_interval_state<int(N), D, int(N0), D0>;
	template<size_t N, size_t D, size_t N0, size_t D0>
	numeric_interval_state(OpFractionLiteral<N, D>, OpNegFractionLiteral<N0, D0>) -> numeric_interval_state<int(N), D, -int(N0), D0>;

	template<size_t N, size_t D>
	numeric_interval_state(OpNegFractionLiteral<N, D>, OpIdentity) -> numeric_interval_state<-int(N), D, 1, 1>;
	template<size_t N, size_t D>
	numeric_interval_state(OpNegFractionLiteral<N, D>, OpNegIdentity) -> numeric_interval_state<-int(N), D, -1, 1>;
	template<size_t N, size_t D>
	numeric_interval_state(OpNegFractionLiteral<N, D>, OpVoid) -> numeric_interval_state<-int(N), D, 0, 1>;
	template<size_t N, size_t D, size_t N0, size_t D0>
	numeric_interval_state(OpNegFractionLiteral<N, D>, OpFractionLiteral<N0, D0>) -> numeric_interval_state<-int(N), D, int(N0), D0>;
	template<size_t N, size_t D, size_t N0, size_t D0>
	numeric_interval_state(OpNegFractionLiteral<N, D>, OpNegFractionLiteral<N0, D0>) -> numeric_interval_state<-int(N), D, -int(N0), D0>;

	template<typename A, typename B>
	constexpr auto numeric_range(A, B)
	{
		return numeric_interval_state(A{}, B{});
	}

	//! Represents the beginning of a constant interval.
	/*!
	 * Represents the beginning of a constant interval.
	 *
	 * \tparam N The numerator of the starting value.
	 * \tparam D The denominator of the starting value.
	 */
	template<int N, size_t D>
	struct numeric_range_state_start
	{
		template<typename A>
		constexpr numeric_range_state_start(A) {}

		template<typename B>
		auto operator>(B) const
		{
			if constexpr (N < 0)
			{
				return numeric_range(-expr::make_fraction<size_t(N), D>(), B{});
			}
			else
			{
				return numeric_range(expr::make_fraction<size_t(N), D>(), B{});
			}
		}
	};

	numeric_range_state_start(OpIdentity)->numeric_range_state_start<1, 1>;
	numeric_range_state_start(OpVoid)->numeric_range_state_start<0, 1>;
	numeric_range_state_start(OpNegIdentity)->numeric_range_state_start<-1, 1>;
	template<size_t N, size_t D>
	numeric_range_state_start(OpFractionLiteral<N, D>) -> numeric_range_state_start<int(N), D>;
	template<size_t N, size_t D>
	numeric_range_state_start(OpNegFractionLiteral<N, D>) -> numeric_range_state_start<-int(N), D>;

	template<typename A>
	constexpr auto numeric_range(A)
	{
		return numeric_range_state_start(A{});
	}
}

constexpr auto OpVoid::operator--(int) const
{
	return expr::numeric_range(OpVoid{});
}

constexpr auto OpIdentity::operator--(int) const
{
	return expr::numeric_range(OpIdentity{});
}

constexpr auto OpNegIdentity::operator--(int) const
{
	return expr::numeric_range(OpNegIdentity{});
}

namespace symphas::internal
{
	template<typename G, typename V, typename E, typename T, typename working_grid>
	auto to_optimized(OpIntegral<V, E, T> const& e, OpTerm<OpIdentity, working_grid*> const& data)
	{
		return expr::transform::swap_expression<G>(e, data);
	}
}

template<typename V, typename E, typename T, typename G>
struct OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>> : OpExpression<OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>>>
{
	using this_type = OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>>;
	using working_grid = expr::storage_type_t<G>;

	using expr_t = std::invoke_result_t<decltype(&symphas::internal::to_optimized<G, V, E, T, working_grid>), OpIntegral<V, E, T>, OpTerm<OpIdentity, working_grid*>>;

	OpOptimized() : working{ 0 }, e{ }, term{ } {}

	OpOptimized(OpBinaryMul<OpIntegral<V, E, T>, G> const& e) : 
		working{ expr::data_dimensions(e.b) }, e{ symphas::internal::to_optimized<G>(e.a, expr::make_term(&working)) }, term{ e.b } 
	{
		symphas::internal::update_temporary_grid(working, term);
		expr::prune::update(this->e);
	}

	OpOptimized(OpOptimized const& other) : OpOptimized(other.get_expression()) {}
	OpOptimized(OpOptimized&& other) noexcept : OpOptimized() 
	{
		swap(*this, other);
	}
	OpOptimized& operator=(OpOptimized other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>>& first, OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>>& second)
	{
		using std::swap;
		swap(first.working, second.working);
		swap(first.e, second.e);
		swap(first.term, second.term);
	}

	auto eval(iter_type n) const
	{
		return e.eval(n) * expr::BaseData<working_grid>::get(working, n);
	}

	template<typename... condition_ts>
	void update(symphas::lib::types_list<condition_ts...>)
	{
		symphas::internal::update_temporary_grid(working, term);
		expr::result(term, working, grid::get_iterable_domain(working));
		expr::prune::update<condition_ts...>(e);
	}

	void update()
	{
		update(symphas::lib::types_list<>{});
	}

	auto operator-() const;

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return get_expression().print(out);
	}

	size_t print(char* out) const
	{
		return get_expression().print(out);
	}

	size_t print_length() const
	{
		return get_expression().print_length();
	}
#endif

	auto get_expression() const
	{
		return expr::transform::swap_grid<OpTerm<OpIdentity, working_grid*>>(e, term) * term;
	}

	working_grid working;
	expr_t e;
	G term;
};

template<typename E>
OpOptimized(E) -> OpOptimized<E>;

namespace expr::transform
{
	template<typename E>
	auto optimize(OpEvaluable<E> const& e)
	{
		return *static_cast<E const*>(&e);
	}

	template<typename V, typename E, typename T, typename G,
		typename = std::enable_if_t<expr::satisfies<E, expr::contains_matching_anywhere<G>>, int>>
	auto optimize(OpBinaryMul<OpIntegral<V, E, T>, G> const& e)
	{
		return OpOptimized(e);
	}

	template<typename... Es, size_t... Is>
	auto optimize(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return (optimize(expr::get<Is>(e)) + ...);
	}

	template<typename... Es>
	auto optimize(OpAdd<Es...> const& e)
	{
		return optimize(e, std::make_index_sequence<sizeof...(Es)>{});
	}
}



template<typename V, typename E, typename T, typename G>
auto OpOptimized<OpBinaryMul<OpIntegral<V, E, T>, G>>::operator-() const
{
	return expr::transform::optimize(-get_expression());
}


