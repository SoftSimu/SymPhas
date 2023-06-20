
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
 * PURPOSE: Defines expressions that are to be used as operators.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressions.h"
#include "gridpair.h"

//! \cond

#ifdef LATEX_PLOT
#define SYEX_COMBINATION_FMT_A "\\left("
#define SYEX_COMBINATION_FMT_B "\\right)"
#define SYEX_COMBINATION_FMT_SEP " + "
#define SYEX_CHAIN_FMT_SEP " "

#define SYEX_OP_MAP_FMT_A "\text{map}\\left("
#define SYEX_OP_MAP_FMT_B "\\right)"

#else
#define SYEX_COMBINATION_FMT_A "("
#define SYEX_COMBINATION_FMT_B ")"
#define SYEX_COMBINATION_FMT_SEP " + "
#define SYEX_CHAIN_FMT_SEP " * "

#define SYEX_OP_MAP_FMT_A "map("
#define SYEX_OP_MAP_FMT_B ")"
#endif

#define SYEX_CHAIN_FMT_A SYEX_COMBINATION_FMT_A
#define SYEX_CHAIN_FMT_B SYEX_COMBINATION_FMT_B


#define SYEX_COMBINATION_FMT SYEX_COMBINATION_FMT_A "%s" SYEX_COMBINATION_FMT_SEP "%s" SYEX_COMBINATION_FMT_B
#define SYEX_COMBINATION_FMT_LEN (STR_ARR_LEN(SYEX_COMBINATION_FMT_A SYEX_COMBINATION_FMT_SEP SYEX_COMBINATION_FMT_B) - 1)
#define SYEX_COMBINATION_APPLY_FMT_LEN (STR_ARR_LEN(SYEX_COMBINATION_FMT_A SYEX_COMBINATION_FMT_B) - 1)

#define SYEX_CHAIN_FMT SYEX_CHAIN_FMT_A "%s" SYEX_CHAIN_FMT_SEP "%s" SYEX_CHAIN_FMT_B
#define SYEX_CHAIN_FMT_LEN (STR_ARR_LEN(SYEX_CHAIN_FMT_A SYEX_CHAIN_FMT_SEP SYEX_CHAIN_FMT_B) - 1)
#define SYEX_CHAIN_APPLY_FMT_LEN (STR_ARR_LEN(SYEX_COMBINATION_FMT_A SYEX_COMBINATION_FMT_B) - 1)

#define SYEX_OP_MAP_FMT SYEX_OP_MAP_FMT_A "%s" SYEX_OP_MAP_FMT_B
#define SYEX_OP_MAP_FMT_LEN (STR_ARR_LEN(SYEX_OP_MAP_FMT_A SYEX_OP_MAP_FMT_B) - 1)


//! \endcond


/*
 */

 // ******************************************************************************************

//! An operator in the expression tree.
/*!
 * An operator is another type of expression, and is not a compound type
 * unless it is made concrete. Thus, there are always two types of any 
 * particular operator, the general type which must be applied to an expression
 * before it is evaluated, and then the concrete type that is 
 * applied to a specific expression.
 * 
 * The details associated with the actual evaluation of an operator (applying the operator
 * to an actual expression) is the more interesting part:
 * 
 * 1. It supports operations to "rechain" operators that are applied to it.
 * 2. It typically needs to be pruned in order to evaluate the system, before it is executed,
 * so that the base data is updated.
 * 3. **Most importantly** The operators that it applies need to have a specialized implementation
 * for only ::OpTerm, not just to avoid wasting time memory, but because
 * the combination implementation will apply operators to variables and
 * will NOT update each applied operator before evaluating the combination result. Therefore,
 * each operator needs to have an implementation where it computes the results
 * directly of a variable without any intermediate data that would need to be
 * updated.
 * 
 * \tparam E The specialized operator, used for the CRTP strategy.
 */
template<typename E>
struct OpOperator/* : OpExpression<E> */
{
	inline auto eval(iter_type) const
	{
		return expr::symbols::Symbol{};
	}

	template<typename E0>
	auto operator()(OpExpression<E0> const& e) const
	{
		return cast().apply_impl(*static_cast<E0 const*>(&e));
	}

	template<typename E0>
	auto operator()(OpOperator<E0> const& e) const
	{
		return apply(*static_cast<E0 const*>(&e));
	}

	template<typename E0>
	auto operator*(OpExpression<E0> const& e) const
	{
		return expr::make_mul(cast(), *static_cast<E0 const*>(&e));
		//return cast().apply_impl(*static_cast<E0 const*>(&e));
	}

	template<typename E0>
	auto operator*(OpOperator<E0> const& e) const
	{
		return expr::make_mul(cast(), *static_cast<E0 const*>(&e));
		//return OpOperatorChain(cast(), *static_cast<E0 const*>(&e));
	}

	//template<typename... Es>
	//auto operator*(OpAdd<Es...> const& e) const
	//{
	//	return cast().operator*(e);
	//}

	template<typename E0, 
		typename std::enable_if_t<(expr::is_symbol<expr::eval_type_t<E0>> && expr::grid_dim<E0>::value == 0), int> = 0>
	auto apply(OpExpression<E0> const& e) const
	{
		return OpOperatorChain(cast(), *static_cast<E0 const*>(&e));
	}

	template<typename E0, 
		typename std::enable_if_t<!(expr::is_symbol<expr::eval_type_t<E0>>&& expr::grid_dim<E0>::value == 0), int> = 0>
	auto apply(OpExpression<E0> const& e) const
	{
		return cast().apply_impl(*static_cast<E0 const*>(&e));
	}

	template<typename E0>
	auto apply(OpOperator<E0> const& e) const
	{
		return OpOperatorChain(cast(), *static_cast<E0 const*>(&e));
	}

	//! The addition of two operators creates a combination.
	template<typename F, std::enable_if_t<!std::is_same<F, E>::value, int> = 0>
	auto operator+(OpOperator<F> const& b) const
	{
		return OpOperatorCombination(cast(), *static_cast<F const*>(&b));
	}

	//! Specialized behaviour adding an operator to a combination operator.
	template<typename B1, typename B2, std::enable_if_t<!std::is_same<OpOperatorCombination<B1, B2>, E>::value, int> = 0>
	auto operator+(OpOperatorCombination<B1, B2> const& b) const
	{
		return OpOperatorCombination(b, cast());
	}

	//! Specialized behaviour adding an operator to a chain operator.
	template<typename B1, typename B2, std::enable_if_t<!std::is_same<OpOperatorChain<B1, B2>, E>::value, int> = 0>
	auto operator+(OpOperatorChain<B1, B2> const& b) const
	{
		return OpOperatorCombination(b, cast());
	}


	//! The subtraction of two operators creates a combination.
	template<typename F, std::enable_if_t<!std::is_same<F, E>::value, int> = 0>
	auto operator-(OpOperator<F> const& b) const
	{
		return OpOperatorCombination(cast(), -*static_cast<F const*>(&b));
	}

	//! Specialized behaviour subtracting operator to a combination operator.
	template<typename B1, typename B2, std::enable_if_t<!std::is_same<OpOperatorCombination<B1, B2>, E>::value, int> = 0>
	auto operator-(OpOperatorCombination<B1, B2> const& b) const
	{
		return OpOperatorCombination(-b, cast());
	}

	//! Specialized behaviour subtracting operator to a chain operator.
	template<typename B1, typename B2, std::enable_if_t<!std::is_same<OpOperatorChain<B1, B2>, E>::value, int> = 0>
	auto operator-(OpOperatorChain<B1, B2> const& b) const
	{
		return OpOperatorCombination(-b, cast());
	}

	//! An operator applied to a combination creates a chain operator.
	//template<typename B1, typename B2>
	//auto operator*(OpOperatorCombination<B1, B2> const& b) const
	//{
	//	//return OpOperatorChain(cast(), b);
	//	return expr::make_mul(cast(), b);
	//	//return cast() * b.f + cast() * b.g;
	//}

	 //! An operator applied to a combination creates a chain operator.
	template<typename B1, typename B2, typename F>
	auto operator*(OpCombination<B1, B2, F> const& b) const
	{
		//OpOperatorChain combination(cast(), b.combination);
		//return OpChain(combination, b.e);
		return expr::make_mul(cast(), b);
		//return (cast() * b.combination)(expr::get_enclosed_expression(b));
	}

	//! An operator applied to a chain builds the chain operator.
	template<typename B1, typename B2, typename F>
	auto operator*(OpChain<B1, B2, F> const& b) const
	{
		//OpOperatorChain combination(cast(), b.combination);
		//return OpChain(combination, b.e);
		return expr::make_mul(cast(), b);
		//return OpOperatorChain(cast(), b.combination)(expr::get_enclosed_expression(b));
	}

	E const& cast() const
	{
		return (*static_cast<E const*>(this));
	}

	E& cast()
	{
		return (*static_cast<E*>(this));
	}
};
//
//template<typename E, typename coeff_t, typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
//auto operator*(OpOperator<E> const& e, coeff_t const& v)
//{
//	return (*static_cast<E const*>(&e))(v);
//}


//! The addition of two operators creates a combination.
template<typename E>
auto operator+(OpOperator<E> const& a, OpOperator<E> const& b)
{
	if constexpr (std::is_same<expr::symbols::Symbol, expr::eval_type_t<E>>::value)
	{
		return OpOperatorChain(OpIdentity{} + OpIdentity{}, *static_cast<E const*>(&b));
	}
	else
	{
		return OpOperatorCombination(*static_cast<E const*>(&a), *static_cast<E const*>(&b));
	}
}

//! The addition of two operators creates a combination.
template<typename E>
auto operator-(OpOperator<E> const& a, OpOperator<E> const& b)
{
	if constexpr (std::is_same<expr::symbols::Symbol, expr::eval_type_t<E>>::value)
	{
		return OpVoid{};
	}
	else
	{
		return OpOperatorCombination(*static_cast<E const*>(&a), *static_cast<E const*>(&b));
	}
}

//! The addition of two same operators adds their coefficients.
template<typename coeff_t, typename E, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator+(OpOperatorChain<coeff_t, E> const& a, OpOperator<E> const& b)
{
	if constexpr (std::is_same<expr::symbols::Symbol, expr::eval_type_t<E>>::value)
	{
		using coeff2_t = decltype(coeff_t{} + OpIdentity{});
		if constexpr (std::is_same<coeff2_t, OpIdentity>::value)
		{
			return *static_cast<E const*>(&b);
		}
		else if constexpr (std::is_same<coeff2_t, OpVoid>::value)
		{
			return OpVoid{};
		}
		return OpOperatorChain(coeff2_t{}, *static_cast<E const*>(&b));
	}
	else
	{
		return OpOperatorCombination(a, *static_cast<E const*>(&b));
	}
}

//! The addition of two same operators adds their coefficients.
template<typename coeff_t, typename E, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator+(OpOperator<E> const& a, OpOperatorChain<coeff_t, E> const& b)
{
	if constexpr (std::is_same<expr::symbols::Symbol, expr::eval_type_t<E>>::value)
	{
		using coeff2_t = decltype(OpIdentity{} + coeff_t{});
		if constexpr (std::is_same<coeff2_t, OpIdentity>::value)
		{
			return *static_cast<E const*>(&b);
		}
		else if constexpr (std::is_same<coeff2_t, OpVoid>::value)
		{
			return OpVoid{};
		}
		return OpOperatorChain(coeff2_t{}, *static_cast<E const*>(&b));
	}
	else
	{
		return OpOperatorCombination(a, *static_cast<E const*>(&b));
	}
}

//! The addition of two same operators adds their coefficients.
template<typename coeff_t, typename E, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator-(OpOperatorChain<coeff_t, E> const& a, OpOperator<E> const& b)
{
	if constexpr (std::is_same<expr::symbols::Symbol, expr::eval_type_t<E>>::value)
	{
		using coeff2_t = decltype(coeff_t{} - OpIdentity{});
		if constexpr (std::is_same<coeff2_t, OpIdentity>::value)
		{
			return *static_cast<E const*>(&b);
		}
		else if constexpr (std::is_same<coeff2_t, OpVoid>::value)
		{
			return OpVoid{};
		}
		return OpOperatorChain(coeff2_t{}, *static_cast<E const*>(&b));
	}
	else
	{
		return OpOperatorCombination(a, *static_cast<E const*>(&b));
	}
}

//! The addition of two same operators adds their coefficients.
template<typename coeff_t, typename E, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator-(OpOperator<E> const& a, OpOperatorChain<coeff_t, E> const& b)
{
	if constexpr (std::is_same<expr::symbols::Symbol, expr::eval_type_t<E>>::value)
	{
		using coeff2_t = decltype(OpIdentity{} - coeff_t{});
		if constexpr (std::is_same<coeff2_t, OpIdentity>::value)
		{
			return *static_cast<E const*>(&b);
		}
		else if constexpr (std::is_same<coeff2_t, OpVoid>::value)
		{
			return OpVoid{};
		}
		return OpOperatorChain(coeff2_t{}, *static_cast<E const*>(&b));
	}
	else
	{
		return OpOperatorCombination(a, *static_cast<E const*>(&b));
	}
}



//template<typename E1, typename E2>
//auto operator+(OpOperator<E1> const& a, OpOperator<E2> const& b)
//{
//	return (*static_cast<E1 const*>(&a)).operator+(*static_cast<E2 const*>(&b));
//}
//
//template<typename E1, typename E2>
//auto operator-(OpOperator<E1> const& a, OpOperator<E2> const& b)
//{
//	return (*static_cast<E1 const*>(&a)).operator-(*static_cast<E2 const*>(&b));
//}

template<typename E1, typename E2>
auto operator/(OpOperator<E1> const& a, OpOperator<E2> const& b)
{
	return (*static_cast<E1 const*>(&a)).operator/(*static_cast<E2 const*>(&b));
}

/* rules when used with expressions
 */


 //! Adding an operator to a constant creates a combination.
template<typename A, typename coeff_t,
	typename = std::enable_if_t<(expr::is_fraction<coeff_t> || expr::is_identity<coeff_t>), int>>
auto operator+(OpOperator<A> const& a, coeff_t)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), coeff_t{});
}

//! Subtracting a constant from an operator creates a combination.
template<typename A, typename coeff_t,
	typename = std::enable_if_t<(expr::is_fraction<coeff_t> || expr::is_identity<coeff_t>), int>>
auto operator-(OpOperator<A> const& a, coeff_t)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), -coeff_t{});
}


//! Adding an operator to a constant creates a combination.
template<typename A, typename coeff_t,
	typename = std::enable_if_t<(expr::is_fraction<coeff_t> || expr::is_identity<coeff_t>), int>>
auto operator+(coeff_t, OpOperator<A> const& b)
{
	return OpOperatorCombination(coeff_t{}, *static_cast<A const*>(&b));
}

//! Subtracting an operator from a constant creates a combination.
template<typename A, typename coeff_t,
	typename = std::enable_if_t<(expr::is_fraction<coeff_t> || expr::is_identity<coeff_t>), int>>
auto operator-(coeff_t, OpOperator<A> const& b)
{
	return OpOperatorCombination(coeff_t{}, -*static_cast<A const*>(&b));
}




template<typename A, typename T>
auto operator+(OpOperator<A> const& a, OpLiteral<T> const& b)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), b);
}

//! Subtracting a constant from an operator creates a combination.
template<typename A, typename T>
auto operator-(OpOperator<A> const& a, OpLiteral<T> const& b)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), -b);
}


//! Adding an operator to a constant creates a combination.
template<typename T, typename A>
auto operator+(OpLiteral<T> const& a, OpOperator<A> const& b)
{
	return OpOperatorCombination(a, *static_cast<A const*>(&b));
}

//! Subtracting an operator from a constant creates a combination.
template<typename T, typename A>
auto operator-(OpLiteral<T> const& a, OpOperator<A> const& b)
{
	return OpOperatorCombination(a, -*static_cast<A const*>(&b));
}


template<typename A, typename T, typename I>
auto operator+(OpOperator<A> const& a, OpCoeff<T, I> const& b)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), b);
}

//! Subtracting a constant from an operator creates a combination.
template<typename A, typename T, typename I>
auto operator-(OpOperator<A> const& a, OpCoeff<T, I> const& b)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), -b);
}


//! Adding an operator to a constant creates a combination.
template<typename T, typename I, typename A>
auto operator+(OpCoeff<T, I> const& a, OpOperator<A> const& b)
{
	return OpOperatorCombination(a, *static_cast<A const*>(&b));
}

//! Subtracting an operator from a constant creates a combination.
template<typename T, typename I, typename A>
auto operator-(OpCoeff<T, I> const& a, OpOperator<A> const& b)
{
	return OpOperatorCombination(a, -*static_cast<A const*>(&b));
}

//! The operator as the second term of a multiplication is applied.
template<typename A1, typename A2, typename E>
auto operator*(OpBinaryMul<A1, OpOperator<E>> const& a, OpExpression<A2> const& b)
{
	return a.a * (a.b * (*static_cast<const A2*>(&b)));
}


template<typename E1, typename E2>
auto operator+(OpExpression<E1> const& a, OpOperator<E2> const& b)
{
	return expr::make_add(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
}

template<typename E1, typename E2>
auto operator+(OpOperator<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_add(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
}

template<typename E1, typename E2>
auto operator-(OpExpression<E1> const& a, OpOperator<E2> const& b)
{
	return expr::make_add(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
}

template<typename E1, typename E2>
auto operator-(OpOperator<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_add(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
}

template<typename E1, typename E2>
auto operator/(OpExpression<E1> const& a, OpOperator<E2> const& b)
{
	return expr::make_div(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
}

template<typename E1, typename E2>
auto operator/(OpOperator<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_div(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
}

// ******************************************************************************************
//
//
////! An operator is applied to an expression through multiplication.
//template<typename E, typename F>
//auto operator*(OpOperator<E> const& a, OpExpression<F> const& b)
//{
//	return (*static_cast<E const*>(&a)).apply(*static_cast<F const*>(&b));
//}
//
////! An operator is applied to an expression through multiplication.
//template<typename E, typename... Es>
//auto operator*(OpOperator<E> const& a, OpAdd<Es...> const& b)
//{
//	return (*static_cast<E const*>(&a)).apply(b);
//}
//
////! An operator is applied to an expression through multiplication.
//template<typename E, typename A, typename B>
//auto operator*(OpOperator<E> const& a, OpBinaryDiv<A, B> const& b)
//{
//	return (*static_cast<E const*>(&a)).apply(b);
//}
//
////! An operator is applied to an expression through multiplication.
//template<typename E>
//auto operator*(OpOperator<E> const& a, OpIdentity)
//{
//	return (*static_cast<E const*>(&a)).apply(OpIdentity{});
//}
//
////! An operator is applied to an expression through multiplication.
//template<typename E>
//auto operator*(OpOperator<E> const& a, OpNegIdentity)
//{
//	return (*static_cast<E const*>(&a)).apply(OpNegIdentity{});
//}
//
////! An operator is applied to an expression through multiplication.
//template<typename E, typename T>
//auto operator*(OpOperator<E> const& a, OpLiteral<T> const& v)
//{
//	return (*static_cast<E const*>(&a)).apply(v);
//}
//
////! An operator is applied to an expression through multiplication.
//template<typename E, size_t N, size_t D>
//auto operator*(OpOperator<E> const& a, OpFractionLiteral<N, D>)
//{
//	return (*static_cast<E const*>(&a)).apply(OpFractionLiteral<N, D>{});
//}
//
////! An operator is applied to an expression through multiplication.
//template<typename E, size_t N, size_t D>
//auto operator*(OpOperator<E> const& a, OpNegFractionLiteral<N, D>)
//{
//	return (*static_cast<E const*>(&a)).apply(OpNegFractionLiteral<N, D>{});
//}

// ******************************************************************************************

namespace expr
{


	//! The combination operator is applied to terms in a binary addition.
	/*!
	 * Distributing an operator will split apart the expression by add/subtraction to apply the
	 * operator to each term of the expression. The operator itself is not expanded.
	 */
	template<typename A1, typename A2, typename... Es>
	auto distribute_operator(OpOperatorCombination<A1, A2> const& a, OpAdd<Es...> const& b);

	//! The expression is applied to the combination operator.
	/*!
	 * Distributing an operator will split apart the expression by add/subtraction to apply the
	 * operator to each term of the expression. The operator itself is not expanded.
	 */
	template<typename A1, typename A2, typename E>
	auto distribute_operator(OpOperatorCombination<A1, A2> const& a, OpExpression<E> const& b);

	//! The expression is applied to the combination operator.
	/*!
	 * Distributing an operator will split apart the expression by add/subtraction to apply the
	 * operator to each term of the expression. The operator itself is not expanded.
	 */
	template<typename A1, typename A2, typename E>
	auto distribute_operator(OpOperatorCombination<A1, A2> const& a, OpOperator<E> const& b);

	//! The chain operator is applied to terms in a binary addition.
	/*!
	 * Distributing an operator will split apart the expression by add/subtraction to apply the
	 * operator to each term of the expression. The operator itself is not expanded.
	 */
	template<typename A1, typename A2, typename... Es>
	auto distribute_operator(OpOperatorChain<A1, A2> const& a, OpAdd<Es...> const& b);

	//! The expression is applied to the chain operator.
	/*!
	 * Distributing an operator will split apart the expression by add/subtraction to apply the
	 * operator to each term of the expression. The operator itself is not expanded.
	 */
	template<typename A1, typename A2, typename E>
	auto distribute_operator(OpOperatorChain<A1, A2> const& a, OpExpression<E> const& b);

	//! The expression is applied to the chain operator.
	/*!
	 * Distributing an operator will split apart the expression by add/subtraction to apply the
	 * operator to each term of the expression. The operator itself is not expanded.
	 */
	template<typename A1, typename A2, typename E>
	auto distribute_operator(OpOperatorChain<A1, A2> const& a, OpOperator<E> const& b);

	//! The expression is distributed among all terms in the expression.
	/*!
	 * Distributing an operator will split apart the expression by add/subtraction to apply the
	 * operator to each term of the expression. The operator itself is not expanded.
	 */
	template<typename E1, typename E2>
	auto distribute_operator(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
	}
	
	//! The expression is distributed among all terms in the expression.
	/*!
	 * Distributing an operator will split apart the expression by add/subtraction to apply the
	 * operator to each term of the expression. The operator itself is not expanded.
	 */
	template<typename E1, typename E2>
	auto distribute_operator(OpExpression<E1> const& a, OpOperator<E2> const& b)
	{
		return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
	}

	//! The expression is distributed among all terms in the expression.
	/*!
	 * Distributing an operator will split apart the expression by add/subtraction to apply the
	 * operator to each term of the expression. The operator itself is not expanded.
	 */
	template<typename E1, typename E2>
	auto distribute_operator(OpOperator<E1> const& a, OpExpression<E2> const& b)
	{
		return (*static_cast<E1 const*>(&a))(*static_cast<E2 const*>(&b));
	}

	//! The expression is distributed among all terms in the expression.
	/*!
	 * Distributing an operator will split apart the expression by add/subtraction to apply the
	 * operator to each term of the expression. The operator itself is not expanded.
	 */
	template<typename E1, typename E2>
	auto distribute_operator(OpOperator<E1> const& a, OpOperator<E2> const& b)
	{
		return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
	}

	template<typename A1, typename A2, typename... Es, size_t... Is>
	auto distribute_operator(OpOperatorCombination<A1, A2> const& a, OpAdd<Es...> const& b, std::index_sequence<Is...>)
	{
		return (expr::distribute_operator(a, expr::get<Is>(b)) + ...);
	}

	template<typename A1, typename A2, typename... Es>
	auto distribute_operator(OpOperatorCombination<A1, A2> const& a, OpAdd<Es...> const& b)
	{
		return expr::distribute_operator(a, b, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<typename A1, typename A2, typename E>
	auto distribute_operator(OpOperatorCombination<A1, A2> const& a, OpExpression<E> const& b)
	{
		return a(*static_cast<E const*>(&b));
	}

	template<typename A1, typename A2, typename E>
	auto distribute_operator(OpOperatorCombination<A1, A2> const& a, OpOperator<E> const& b)
	{
		return a * (*static_cast<E const*>(&b));
	}

	template<typename A1, typename A2, typename... Es, size_t... Is>
	auto distribute_operator(OpOperatorChain<A1, A2> const& a, OpAdd<Es...> const& b, std::index_sequence<Is...>)
	{
		return (expr::distribute_operator(a, expr::get<Is>(b)) + ...);
	}

	template<typename A1, typename A2, typename... Es>
	auto distribute_operator(OpOperatorChain<A1, A2> const& a, OpAdd<Es...> const& b)
	{
		return expr::distribute_operator(a, b, std::make_index_sequence<sizeof...(Es)>{});
	}

	template<typename A1, typename A2, typename E>
	auto distribute_operator(OpOperatorChain<A1, A2> const& a, OpExpression<E> const& b)
	{
		return a(*static_cast<E const*>(&b));
	}

	template<typename A1, typename A2, typename E>
	auto distribute_operator(OpOperatorChain<A1, A2> const& a, OpOperator<E> const& b)
	{
		return a * (*static_cast<E const*>(&b));
	}




	/*
	 * a distribution method for the operators to apply to combinations
	 */


	//! Expanding two expression is simply their multiplication.
	/*!
	 * Attempting to expand two expressions which are not operators
	 * will multiply the expressions together.
	 * 
	 * Expanding an operator will attempt to distribute each term of an operator
	 * to an expression. The expression is not changed (i.e. it is not split by add/subtraction).
	 * This is opposed to distributing an operator, which will attempt to
	 * apply each individual operator term of a chain or combination to an expression without
	 * breaking apart the expression.
	 *
	 * \param a The expression on the left hand side of the multiplication
	 * operator.
	 * \param b The expression on the right hand side of the multiplication
	 * operator.
	 */
	template<typename E1, typename E2>
	auto expand_operator(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
	}

	//! Expanding an operator into an expression means applying the operator.
	/*!
	 * Attempting to expand two expressions which are not operators
	 * will multiply the expressions together.
	 *
	 * Expanding an operator will attempt to distribute each term of an operator
	 * to an expression. The expression is not changed (i.e. it is not split by add/subtraction).
	 * This is opposed to distributing an operator, which will attempt to
	 * apply each individual operator term of a chain or combination to an expression without
	 * breaking apart the expression.
	 *
	 * \param a The expression on the left hand side of the multiplication
	 * operator.
	 * \param b The expression on the right hand side of the multiplication
	 * operator.
	 */
	template<typename E1, typename E2>
	auto expand_operator(OpOperator<E1> const& a, OpExpression<E2> const& b)
	{
		return (*static_cast<E1 const*>(&a))(*static_cast<E2 const*>(&b));
	}


	//! Expanding a combination into an expression.
	/*!
	 * Applies each of the constituent operators to the given expression.
	 * In this way, this distributes the operators to the expression.
	 * 
	 * Expanding an operator will attempt to distribute each term of an operator
	 * to an expression. The expression is not changed (i.e. it is not split by add/subtraction).
	 * This is opposed to distributing an operator, which will attempt to
	 * apply each individual operator term of a chain or combination to an expression without
	 * breaking apart the expression.
	 * 
	 * \param a The operator combination which is distributed.
	 * \param b The expression which forms the concrete operators once applied.
	 */
	template<typename A1, typename A2, typename E>
	auto expand_operator(OpOperatorCombination<A1, A2> const& a, OpExpression<E> const& b)
	{
		return expand_operator(a.f, *static_cast<E const*>(&b)) + expand_operator(a.g, *static_cast<E const*>(&b));
	}

	//! Expanding a chain operator into an expression.
	/*!
	 * Applies the chain operator to an expression, which will apply
	 * the expression into the nested operator and then apply the outer operator
	 * to that result.
	 *
	 * Expanding an operator will attempt to distribute each term of an operator
	 * to an expression. The expression is not changed (i.e. it is not split by add/subtraction).
	 * This is opposed to distributing an operator, which will attempt to
	 * apply each individual operator term of a chain or combination to an expression without
	 * breaking apart the expression.
	 * 
	 * \param a The chain operator which is applied to an expression.
	 * \param b The expression applied by the chain operator.
	 */
	template<typename A1, typename A2, typename E>
	auto expand_operator(OpOperatorChain<A1, A2> const& a, OpExpression<E> const& b)
	{
		return expand_operator(a.f, expand_operator(a.g, *static_cast<E const*>(&b)));
	}

}


// ******************************************************************************************


//! A linear combination of two general operators.
/*!
 * A linear combination of general operators (always addition, the coefficient 
 * of the second operator will indicate the correct sign).
 * There are two operators (additional operators are implemented recursively); 
 * this object is itself an operator.
 * 
 * This is the binary expression as related to operators. It is equivalent
 * to binary addition for usual expression objects.
 *
 * This operator is not necessarily linear, depending on its constituent 
 * operators.
 * 
 * \tparam A1 Type of the first operator.
 * \tparam A2 Type of the second operator.
 */
template<typename A1, typename A2>
struct OpOperatorCombination : OpOperator<OpOperatorCombination<A1, A2>>
{
	using parent_type = OpOperator<OpOperatorCombination<A1, A2>>;
	using parent_type::operator*;
	using parent_type::operator-;
	using parent_type::operator+;

	//! Create the combination of two operators.
	/*!
	 * Create the combination of two operators.
	 * 
	 * \param f The operator on the left hand side of the addition.
	 * \param g The operator on the right hand side of the addition.
	 */
	OpOperatorCombination(A1 const& f, A2 const& g) : f{ f }, g{ g } {}
	OpOperatorCombination() : f{}, g{} {}

	inline auto eval(iter_type n = 0) const
	{
		return f.eval(n) + g.eval(n);
	}

	auto operator-() const
	{
		return -f + (-g);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = 0;
		n += fprintf(out, SYEX_COMBINATION_FMT_A);
		n += f.print(out);
		n += fprintf(out, SYEX_COMBINATION_FMT_SEP);
		n += g.print(out);
		n += fprintf(out, SYEX_COMBINATION_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_COMBINATION_FMT_A);
		n += f.print(out + n);
		n += sprintf(out + n, SYEX_COMBINATION_FMT_SEP);
		n += g.print(out + n);
		n += sprintf(out + n, SYEX_COMBINATION_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return f.print_length() + g.print_length() + SYEX_COMBINATION_FMT_LEN;
	}


	//! Apply the chain operation to an expression.
	template<typename E, typename AA1 = A1, typename AA2 = A2, typename std::enable_if_t<
		!(std::is_same<mul_result_t<A1, E>, decltype(std::declval<AA1>().apply(std::declval<E>()))>::value/*std::invoke_result_t<decltype(&A1::template apply<E>), A1, E>>::value*/
			&& std::is_same<mul_result_t<A2, E>, decltype(std::declval<AA2>().apply(std::declval<E>()))>::value/*std::invoke_result_t<decltype(&A2::template apply<E>), A2, E>>::value*/), int> = 0>
	auto operator*(OpExpression<E> const& e)
	{
		return expr::make_mul(*this, *static_cast<E const*>(&e));
		//return a.f * (*static_cast<E const*>(&e)) + a.g * (*static_cast<E const*>(&e));
	}

	//! Apply the chain operation to an expression.
	template<typename E, typename AA1 = A1, typename AA2 = A2, typename std::enable_if_t<
		(std::is_same<mul_result_t<AA1, E>, decltype(std::declval<AA1>().apply(std::declval<E>()))>::value
			&& std::is_same<mul_result_t<AA2, E>, decltype(std::declval<AA2>().apply(std::declval<E>()))>::value), int> = 0>
	auto operator*(OpExpression<E> const& e)
	{
		return apply(*static_cast<E const*>(&e));
	}

	template<typename B1, typename B2>
	auto operator*(OpOperatorCombination<B1, B2> const& b)
	{
		return f * b + g * b;
	}


#endif

	//! Apply this operator an expression.
	/*!
	 * Apply this operator an expression and return the concrete form of the
	 * linear combination expression.
	 * 
	 * \param a The expression to which this operator is applied.
	 */
	template<typename E>
	auto apply_impl(OpExpression<E> const& e) const
	{
		//if constexpr (expr::is_symbol<expr::eval_type_t<E>>)
		//{
		//	return OpOperatorChain(*this, *static_cast<E const*>(&e));
		//}
		//else
		{
			return OpCombination(*this, *static_cast<E const*>(&e));
		}
	}

	
	A1 f; //!< Operator on the left of the plus sign.
	A2 g; //!< Operator on the right of the plus sign.
};
//
////! Apply the chain operation to an expression.
//template<typename E, typename A1, typename A2, typename std::enable_if_t<
//	!(std::is_same<mul_result_t<A1, E>, decltype(std::declval<A1>().apply(std::declval<E>()))>::value/*std::invoke_result_t<decltype(&A1::template apply<E>), A1, E>>::value*/
//		&& std::is_same<mul_result_t<A2, E>, decltype(std::declval<A2>().apply(std::declval<E>()))>::value/*std::invoke_result_t<decltype(&A2::template apply<E>), A2, E>>::value*/), int> = 0>
//auto operator*(OpOperatorCombination<A1, A2> const& a, OpExpression<E> const& e)
//{
//	return expr::make_mul(a, *static_cast<E const*>(&e));
//	//return a.f * (*static_cast<E const*>(&e)) + a.g * (*static_cast<E const*>(&e));
//}
//
////! Apply the chain operation to an expression.
//template<typename E, typename A1, typename A2, typename std::enable_if_t<
//	(std::is_same<mul_result_t<A1, E>, decltype(std::declval<A1>().apply(std::declval<E>()))>::value
//		&& std::is_same<mul_result_t<A2, E>, decltype(std::declval<A2>().apply(std::declval<E>()))>::value), int> = 0>
//auto operator*(OpOperatorCombination<A1, A2> const& a, OpExpression<E> const& e)
//{
//	return a.apply(*static_cast<E const*>(&e));
//}
//
//template<typename A1, typename A2, typename B1, typename B2>
//auto operator*(OpOperatorCombination<A1, A2> const& a, OpOperatorCombination<B1, B2> const& b)
//{
//	return a.f * b + a.g * b;
//}

//
//template<typename coeff_t, typename A1, typename A2, typename = std::enable_if_t<expr::is_coeff<coeff_t>, int>>
//auto operator*(coeff_t const& a, OpOperatorCombination<A1, A2> const& b)
//{
//	return a * b.f + a * b.g;
//}





//! A concrete linear combination of two general operators.
/*!
 * See ::OpOperatorCombination. Concrete implementation 
 * where the operator is applied to a specific expression.
 *
 * \tparam A1 Type of the first operator.
 * \tparam A2 Type of the second operator.
 * \tparam E The type of the expression the operator applies to.
 */
template<typename A1, typename A2, typename E>
struct OpCombination : OpExpression<OpCombination<A1, A2, E>>
{
	OpOperatorCombination<A1, A2> combination;	//!< The combination operator.

protected:

	template<typename A>
	static auto _make_eval_expr(OpOperator<A> const& a, OpExpression<E> const& e)
	{
		return (*static_cast<A const*>(&a))(*static_cast<E const*>(&e));
	}

	template<typename A>
	static auto _make_eval_expr(OpExpression<A> const& a, OpExpression<E> const& e)
	{
		return (*static_cast<A const*>(&a)) * (*static_cast<E const*>(&e));
	}

	template<typename A>
	static auto make_eval_expr(A const& a, E const& e)
	{
		return _make_eval_expr(a, e);
	}


	using expr_type_f = std::invoke_result_t<decltype(&OpCombination<A1, A2, E>::make_eval_expr<A1>), A1, E>;
	using expr_type_g = std::invoke_result_t<decltype(&OpCombination<A1, A2, E>::make_eval_expr<A2>), A2, E>;

	expr_type_f eval_expr_f;
	expr_type_g eval_expr_g;


public:

	OpCombination() : combination{}, eval_expr_f{}, eval_expr_g{}, e{} {}

	//! Create the combination of two operators applied to an expression.
	/*!
	 * Create the combination of two operators applied to an expression.
	 *
	 * \param combination The operator being applied.
	 * \param e The expression the operator is applied to.
	 */
	OpCombination(OpOperatorCombination<A1, A2> const& combination, E const& e) :
		combination{ combination },
		eval_expr_f{ make_eval_expr(combination.f, e) }, eval_expr_g{ make_eval_expr(combination.g, e) }, e{ e } {}

	inline auto update()
	{
		expr::prune::update(eval_expr_f);
		expr::prune::update(eval_expr_g);
	}

	inline auto eval(iter_type n) const
	{
		return eval_expr_f.eval(n) + eval_expr_g.eval(n);
	}

	auto operator-() const
	{
		return (-combination)(e);
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = combination.print(out);
		n += fprintf(out, SYEX_COMBINATION_FMT_A);
		n += e.print(out);
		n += fprintf(out, SYEX_COMBINATION_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = combination.print(out);
		n += sprintf(out + n, SYEX_COMBINATION_FMT_A);
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_COMBINATION_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return combination.print_length() + e.print_length() 
			+ SYEX_COMBINATION_APPLY_FMT_LEN;
	}

#endif

    template<typename A10, typename A20, typename E0>
	friend auto const& expr::get_enclosed_expression(OpCombination<A10, A20, E0> const&);
    template<typename A10, typename A20, typename E0>
	friend auto& expr::get_enclosed_expression(OpCombination<A10, A20, E0>&);

protected:

	E e;							//!< Expression to which this operator applies.
};


template<typename coeff_t, typename A1, typename A2, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t>), int> = 0>
auto operator*(coeff_t const& value, OpCombination<A1, A2, E> const& b)
{
	return (value * b.combination)(expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename A1, typename A2, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t>
		&& expr::is_tensor<coeff_t> && expr::eval_type<OpOperatorCombination<A1, A2>>::rank == 0), int> = 0>
auto operator*(coeff_t const& value, OpCombination<A1, A2, E> const& b)
{
	return (b.combination)(value * expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename A1, typename A2, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t>
		&& expr::is_tensor<coeff_t> && expr::eval_type<OpOperatorCombination<A1, A2>>::rank > 0), int> = 0>
	auto operator*(coeff_t const& value, OpCombination<A1, A2, E> const& b)
{
	return (value * b.combination)(expr::get_enclosed_expression(b));
}




// ******************************************************************************************


//! An expression representing an operator applied to another operator.
/*!
 * Represents the result of applying one operator to another.
 * There are two operators (additional operators are implemented recursively);
 * this object is itself an operator.
 *
 * This is the binary expression as related to operators. It is equivalent
 * to binary multiplication for usual expression objects.
 *
 * \tparam A1 Type of the first operator.
 * \tparam A2 Type of the second operator.
 */
template<typename A1, typename A2>
struct OpOperatorChain : OpOperator<OpOperatorChain<A1, A2>>
{
	using parent_type = OpOperator<OpOperatorChain<A1, A2>>;
	using parent_type::operator*;
	using parent_type::operator-;
	using parent_type::operator+;

	//! Create the chain of two operators.
	/*!
	 * Create the chain of two operators.
	 *
	 * \param f The operator on the left hand side of the addition.
	 * \param g The operator on the right hand side of the addition.
	 */
	OpOperatorChain(A1 const& f, A2 const& g) : f{ f }, g{ g } {}
	OpOperatorChain() : f{}, g{} {}

	inline auto eval(iter_type n = 0) const
	{
		return f.eval(n) * g.eval(n);
	}

	auto operator-() const
	{
		return ::OpOperatorChain(-f, g);
	}

	//template<typename E>
	//auto operator*(OpExpression<E> const& e) const
	//{
	//	return f(g * *static_cast<E const*>(&e));
	//}

	//! Apply the chain operation to an expression.
	template<typename E>
	auto apply_impl(OpExpression<E> const& e) const
	{
		//if constexpr (std::is_same<expr::symbols::Symbol, expr::eval_type_t<E>>::value)
		//{
		//	return ::OpOperatorChain(*this, *static_cast<E const*>(&e));
		//}
		//else
		{
			return OpChain(*this, *static_cast<E const*>(&e));
		}
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = 0;
		n += symphas::internal::mul_print(out, f, g);
		//n += fprintf(out, SYEX_CHAIN_FMT_A);
		//n += f.print(out);
		//n += fprintf(out, SYEX_CHAIN_FMT_SEP);
		//n += g.print(out);
		//n += fprintf(out, SYEX_CHAIN_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = 0;
		n += symphas::internal::mul_print(out, f, g);
		//n += sprintf(out + n, SYEX_CHAIN_FMT_A);
		//n += f.print(out + n);
		//n += sprintf(out + n, SYEX_CHAIN_FMT_SEP);
		//n += g.print(out + n);
		//n += sprintf(out + n, SYEX_CHAIN_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return f.print_length() + g.print_length() + /*SYEX_CHAIN_FMT_LEN*/SYEX_MUL_FMT_LEN;
	}

#endif



	A1 f; //!< Operator which is applied to the result of operator `g` on an expression.
	A2 g; //!< Operator which is applied first.
};
//
////! Overload with multiplication by zero.
//template<typename A>
//auto operator*(OpOperator<A> const& a, OpVoid)
//{
//	return OpVoid{};
//}
//
////! Overload with multiplication by zero.
//template<typename A>
//auto operator*(OpVoid, OpOperator<A> const& a)
//{
//	return OpVoid{};
//}
//
////! Apply the chain operation to zero.
//template<typename A1, typename A2>
//auto operator*(OpOperatorChain<A1, A2> const& chain, OpVoid)
//{
//	return OpVoid{};
//}
//
////! Apply the chain operation to zero.
//template<typename A1, typename A2>
//auto operator*(OpVoid, OpOperatorChain<A1, A2> const& chain)
//{
//	return OpVoid{};
//}
//
////! Apply the chain operation to an expression.
//template<typename A1, typename A2, typename E, typename std::enable_if_t<!expr::is_coeff<E>, int> = 0>
//auto operator*(OpOperatorChain<A1, A2> const& chain, E const& e)
//{
//	return OpOperatorChain(chain, e);//chain.f(chain.g * e);
//}
//
////! Apply the chain operation to an expression.
//template<typename A1, typename A2, typename coeff_t, 
//	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
//auto operator*(OpOperatorChain<A1, A2> const& chain, coeff_t const& e)
//{
//	return chain.f(chain.g * e);
//}
//
////! Apply the chain operation to an expression.
//template<typename A1, typename A2, typename B>
//auto operator*(OpOperatorChain<A1, A2> const& a, OpOperator<B> const& b)
//{
//	return OpOperatorChain(a, b);// a.f(a.g * (*static_cast<B const*>(&b)));
//}
//
////! Apply the chain operation to an expression.
//template<typename A, typename B1, typename B2>
//auto operator*(OpOperator<A> const& a, OpOperatorChain<B1, B2> const& b)
//{
//	return ((*static_cast<A const*>(&a)) * b.f)(b.g);
//}


template<typename coeff_t, typename A1, typename A2, 
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(coeff_t const& a, OpOperatorChain<A1, A2> const& b)
{
	return (a * b.f)(b.g);
}

template<typename coeff_t, typename A1, typename A2,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(coeff_t const& a, OpOperatorCombination<A1, A2> const& b)
{
	return a * b.f + a * b.g;
}



//! A concrete operator chain of two general operators.
/*!
 * See OpOperatorChain. Concrete implementation of a chain operator
 * applied to a specific expression.
 *
 * \tparam A1 Type of the first operator.
 * \tparam A2 Type of the second operator.
 * \tparam E The type of the expression the operator applies to.
 */
template<typename A1, typename A2, typename E>
struct OpChain : OpExpression<OpChain<A1, A2, E>>
{

protected:

	auto get_eval_expr(OpOperatorChain<A1, A2> const& combination, E e)
	{
		return combination.f(combination.g(e));
	}

	using expr_type = std::invoke_result_t<decltype(&OpChain<A1, A2, E>::get_eval_expr), OpChain<A1, A2, E>, OpOperatorChain<A1, A2>, E>;
	expr_type eval_expr;		//!< The result of applying the outer operator to the inner.

public:

	OpChain() : combination{}, eval_expr{}, e{} {}

	//! Create the combination of two operators applied to an expression.
	/*!
	 * Create the chain of two operators applied to an expression.
	 *
	 * \param combination The operator being applied.
	 * \param e The expression the operator is applied to.
	 */
	OpChain(OpOperatorChain<A1, A2> const& combination, E const& e) : 
		eval_expr{ get_eval_expr(combination, e) },
		combination{ combination }, e{ e } {}

	inline auto update()
	{
		expr::prune::update(eval_expr);
	}

	inline auto eval(iter_type n) const
	{
		return eval_expr.eval(n);
	}

	auto operator-() const
	{
		return combination(-e);
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = combination.print(out);
		n += fprintf(out, SYEX_CHAIN_FMT_A);
		n += e.print(out);
		n += fprintf(out, SYEX_CHAIN_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = combination.print(out);
		n += sprintf(out + n, SYEX_CHAIN_FMT_A);
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_CHAIN_FMT_B);
		return n;
}

	size_t print_length() const
	{
		return combination.print_length() + e.print_length()
			+ SYEX_CHAIN_APPLY_FMT_LEN;
	}


    template<typename A10, typename A20, typename E0>
	friend auto const& expr::get_enclosed_expression(OpChain<A10, A20, E0> const&);
    template<typename A10, typename A20, typename E0>
	friend auto& expr::get_enclosed_expression(OpChain<A10, A20, E0>&);

#endif

	OpOperatorChain<A1, A2> combination;	//!< The chain operator.
	E e;									//!< Expression to which this operator applies.

};


template<typename coeff_t, typename A1, typename A2, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t>), int> = 0>
auto operator*(coeff_t const& value, OpChain<A1, A2, E> const& b)
{
	return (value * b.combination)(b.e);
}

template<typename coeff_t, typename A1, typename A2, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> 
		&& expr::is_tensor<coeff_t> && expr::eval_type<OpOperatorChain<A1, A2>>::rank == 0), int> = 0>
auto operator*(coeff_t const& value, OpChain<A1, A2, E> const& b)
{
	return (b.combination)(value * b.e);
}

template<typename coeff_t, typename A1, typename A2, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t>
		&& expr::is_tensor<coeff_t> && expr::eval_type<OpOperatorChain<A1, A2>>::rank > 0), int> = 0>
auto operator*(coeff_t const& value, OpChain<A1, A2, E> const& b)
{
	return (value * b.combination)(b.e);
}


// *********************************************************************************************************************************

namespace symphas::internal
{
#ifdef PRINTABLE_EQUATIONS
	template<expr::exp_key_t X, typename E>
	size_t pow_print(FILE* out, OpExpression<E> const& e)
	{
		size_t n = 0;
		n += fprintf(out, SYEX_MUL_FMT_AB);
		n += static_cast<E const*>(&e)->print(out);
		n += fprintf(out, SYEX_MUL_FMT_BA);

		n += fprintf(out, SYEX_POW_SEP_A);
		if constexpr (expr::_Xk_t<X>::D == 1)
		{
			n += fprintf(out, "%u", expr::_Xk_t<X>::N);
		}
		else
		{
			n += fprintf(out, "%u" SYEX_POW_DIV_SEP "%u", expr::_Xk_t<X>::N, expr::_Xk_t<X>::D);
		}
		n += fprintf(out, SYEX_POW_SEP_B);
		return n;
	}

	template<expr::exp_key_t X, typename E>
	size_t pow_print(char* out, OpExpression<E> const& e)
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_MUL_FMT_AB);
		n += static_cast<E const*>(&e)->print(out);
		n += sprintf(out + n, SYEX_MUL_FMT_BA);

		n += sprintf(out + n, SYEX_POW_SEP_A);
		if constexpr (expr::_Xk_t<X>::D == 1)
		{
			n += sprintf(out + n, "%u", expr::_Xk_t<X>::N);
		}
		else
		{
			n += sprintf(out + n, "%u" SYEX_POW_DIV_SEP "%u", expr::_Xk_t<X>::N, expr::_Xk_t<X>::D);
		}
		n += sprintf(out + n, SYEX_POW_SEP_B);
		return n;
	}
#endif
}


//! Binary expression, the multiplication of two terms.
template<expr::exp_key_t X, typename V, typename E>
struct OpPow : OpExpression<OpPow<X, V, E>>
{
	OpPow(V const& value, E const& e) : value{ value }, e{ e } {}

	OpPow() : value{ V{} }, e{ E{} } {}

	inline auto eval(iter_type n) const
	{
		using std::pow;
		using symphas::math::pow;
		using expr::pow;

		if constexpr (expr::_Xk_t<X>::D == 1 && !expr::_Xk_t<X>::sign)
		{
			return expr::eval(value) * pow<expr::_Xk_t<X>::N>(e.eval(n));
		}
		else
		{
			return expr::eval(value) * pow(e.eval(n), expr::_Xk<X>);
		}
	}

	auto operator-() const;

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return symphas::internal::pow_print<X>(out, e);
	}

	size_t print(char* out) const
	{
		return symphas::internal::pow_print<X>(out, e);
	}

	size_t print_length() const
	{
		size_t n = STR_ARR_LEN(SYEX_POW_SEP_A SYEX_POW_SEP_B)
			+ symphas::lib::num_digits<expr::_Xk_t<X>::N>()
			+ e.print_length();

		if constexpr (expr::_Xk_t<X>::D > 1)
		{
			n += symphas::lib::num_digits<expr::_Xk_t<X>::D>() + 1;
		}

		return n;
	}

#endif

    template<expr::exp_key_t X0, typename V0, typename E0>
	friend auto const& expr::get_enclosed_expression(OpPow<X0, V0, E0> const&);
    template<expr::exp_key_t X0, typename V0, typename E0>
	friend auto& expr::get_enclosed_expression(OpPow<X0, V0, E0>&);

	V value;
	E e;		//!< Expression that is the base of the power.

};

namespace expr
{
	//! Constructs the expression representing an expression to a power.
	/*!
	 * Directly constructs the exponent expression of an
	 * expression without applying any rules.
	 *
	 * \param value The coefficient.
	 * \param e The expression.
	 */
	template<expr::exp_key_t X, typename V, typename E>
	auto make_pow(V const& value, E const& e)
	{
		if constexpr (X == expr::Xk<1>)
		{
			return e;
		}
		else
		{
			return OpPow<X, V, E>(value, e);
		}
	}

	//! Constructs the expression representing an expression to a power.
	/*!
	 * Directly constructs the exponent expression of an
	 * expression without applying any rules.
	 *
	 * \param value The coefficient.
	 * \param e The expression.
	 */
	template<expr::exp_key_t X, typename E>
	auto make_pow(E const& e)
	{
		return make_pow<X>(OpIdentity{}, e);
	}

	////! Specialization based on reevaluate(OpExpression<E> const&).
	//template<expr::exp_key_t X, typename V, typename E>
	//auto reevaluate(OpPow<X, V, E> const& e)
	//{
	//	constexpr auto N0 = expr::_Xk_t<X>::N;
	//	constexpr auto D0 = expr::_Xk_t<X>::D;
	//	auto p = expr::get_enclosed_expression(e);
	//
	//	constexpr expr::exp_key_t _X = (sign) ? N0 + D0 : (N0 < D0) ? D0 - N0 : N0 - D0;
	//	constexpr expr::exp_key_t _sign = (sign) ? sign : (N0 < D0) ? true : false;
	//	auto result = expr::coeff(e) * p * reevaluate(expr::make_pow<_X>(p));
	//	
	//	if constexpr (_sign)
	//	{
	//		return -result;
	//	}
	//	else
	//	{
	//		return result;
	//	}
	//}
}

template<typename coeff_t, expr::exp_key_t X, typename V, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V>), int> = 0>
auto operator*(coeff_t const& value, OpPow<X, V, E> const& e)
{
	return expr::make_pow<X>(value * e.value, expr::get_enclosed_expression(e));
}

template<typename coeff_t, expr::exp_key_t X, typename tensor_t, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t>&& expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpPow<X, tensor_t, E> const& e)
{
	return (value * e.value) * expr::make_pow<X>(expr::get_enclosed_expression(e));
}

template<expr::exp_key_t X, typename V, typename E>
auto OpPow<X, V, E>::operator-() const
{
	return expr::make_pow<X>(-value, e);
}



namespace symphas::internal
{

	//! Implementation of map expression construction.
	/*!
	 * Implementation of functions which generalize and simplify the way to
	 * construct map expressions. Wraps the template deduction necessary
	 * to initialize a map expression.
	 */
	template<typename G>
	struct make_map
	{
		//! Constructs the map with the identity coefficient.
		template<typename A>
		static auto get(A&& a);

		//! Constructs the map applied to an expression.
		template<typename V, typename E>
		static auto get(V v, OpExpression<E> const& e);
	};

#ifdef PRINTABLE_EQUATIONS

	template<typename G, typename V, typename E>
	size_t print_map(OpMap<G, V, E> const& map, FILE* out)
	{
		size_t n = expr::print_with_coeff(out, map.value);
		n += fprintf(out, SYEX_OP_MAP_FMT_A);
		n += map.e.print(out);
		n += fprintf(out, SYEX_OP_MAP_FMT_B);
		return n;
	}

	template<typename G, typename V, typename E>
	size_t print_map(OpMap<G, V, E> const& map, char* out)
	{
		size_t n = expr::print_with_coeff(out, map.value);
		n += sprintf(out + n, SYEX_OP_MAP_FMT_A);
		n += map.e.print(out + n);
		n += sprintf(out + n, SYEX_OP_MAP_FMT_B);
		return n;
	}

	template<typename G, typename V, typename E>
	size_t print_map_length(OpMap<G, V, E> const& map)
	{
		return expr::coeff_print_length(map.value) + map.e.print_length() + SYEX_OP_MAP_FMT_LEN;
	}


	template<typename S, typename T, size_t D, typename V, typename E>
	size_t print_map(OpMap<MapGridInverseFourier<S, T, D>, V, E> const& map, FILE* out)
	{
		size_t n = expr::print_with_coeff(out, map.value);
		n += fprintf(out, SYEX_IFT_OF_EXPR_FMT_A);
		n += map.e.print(out);
		n += fprintf(out, SYEX_IFT_OF_EXPR_FMT_B);
		return n;
	}

	template<typename S, typename T, size_t D, typename V0, typename V1, typename G>
	size_t print_map(OpMap<MapGridInverseFourier<S, T, D>, V0, OpTerm<V1, G>> const& map, FILE* out)
	{
		size_t n = expr::print_with_coeff(out, map.value);
		n += fprintf(out, SYEX_IFT_OF_OP_FMT_A);
		n += map.e.print(out);
		n += fprintf(out, SYEX_IFT_OF_OP_FMT_B);
		return n;
	}

	template<typename S, typename T, size_t D, typename V, typename E>
	size_t print_map(OpMap<MapGridFourier<S, T, D>, V, E> const& map, FILE* out)
	{
		size_t n = expr::print_with_coeff(out, map.value);
		n += fprintf(out, SYEX_FT_OF_EXPR_FMT_A);
		n += map.e.print(out);
		n += fprintf(out, SYEX_FT_OF_EXPR_FMT_B);
		return n;
	}

	template<typename S, typename T, size_t D, typename V0, typename V1, typename G>
	size_t print_map(OpMap<MapGridFourier<S, T, D>, V0, OpTerm<V1, G>> const& map, FILE* out)
	{
		size_t n = expr::print_with_coeff(out, map.value);
		n += fprintf(out, SYEX_FT_OF_OP_FMT_A);
		n += map.e.print(out);
		n += fprintf(out, SYEX_FT_OF_OP_FMT_B);
		return n;
	}

	template<typename S, typename T, size_t D, typename V, typename E>
	size_t print_map(OpMap<MapGridInverseFourier<S, T, D>, V, E> const& map, char* out)
	{
		size_t n = expr::print_with_coeff(out, map.value);
		n += sprintf(out + n, SYEX_IFT_OF_EXPR_FMT_A);
		n += map.e.print(out + n);
		n += sprintf(out + n, SYEX_IFT_OF_EXPR_FMT_B);
		return n;
	}

	template<typename S, typename T, size_t D, typename V0, typename V1, typename G>
	size_t print_map(OpMap<MapGridInverseFourier<S, T, D>, V0, OpTerm<V1, G>> const& map, char* out)
	{
		size_t n = expr::print_with_coeff(out, map.value);
		n += sprintf(out + n, SYEX_IFT_OF_OP_FMT_A);
		n += map.e.print(out + n);
		n += sprintf(out + n, SYEX_IFT_OF_OP_FMT_B);
		return n;
	}

	template<typename S, typename T, size_t D, typename V, typename E>
	size_t print_map(OpMap<MapGridFourier<S, T, D>, V, E> const& map, char* out)
	{
		size_t n = expr::print_with_coeff(out, map.value);
		n += sprintf(out + n, SYEX_FT_OF_EXPR_FMT_A);
		n += map.e.print(out + n);
		n += sprintf(out + n, SYEX_FT_OF_EXPR_FMT_B);
		return n;
	}

	template<typename S, typename T, size_t D, typename V0, typename V1, typename G>
	size_t print_map(OpMap<MapGridFourier<S, T, D>, V0, OpTerm<V1, G>> const& map, char* out)
	{
		size_t n = expr::print_with_coeff(out, map.value);
		n += sprintf(out + n, SYEX_FT_OF_OP_FMT_A);
		n += map.e.print(out + n);
		n += sprintf(out + n, SYEX_FT_OF_OP_FMT_B);
		return n;
	}


	template<typename S, typename T, size_t D, typename V, typename E>
	size_t print_map_length(OpMap<MapGridInverseFourier<S, T, D>, V, E> const& map)
	{
		return expr::coeff_print_length(map.value) + map.e.print_length() + SYEX_IFT_OF_EXPR_FMT_LEN;
	}

	template<typename S, typename T, size_t D, typename V0, typename V1, typename G>
	size_t print_map_length(OpMap<MapGridInverseFourier<S, T, D>, V0, OpTerm<V1, G>> const& map)
	{
		return expr::coeff_print_length(map.value) + map.e.print_length() + SYEX_IFT_OF_OP_FMT_LEN;
	}

	template<typename S, typename T, size_t D, typename V, typename E>
	size_t print_map_length(OpMap<MapGridFourier<S, T, D>, V, E> const& map)
	{
		return expr::coeff_print_length(map.value) + map.e.print_length() + SYEX_FT_OF_EXPR_FMT_LEN;
	}

	template<typename S, typename T, size_t D, typename V0, typename V1, typename G>
	size_t print_map_length(OpMap<MapGridFourier<S, T, D>, V0, OpTerm<V1, G>> const& map)
	{
		return expr::coeff_print_length(map.value) + map.e.print_length() + SYEX_FT_OF_OP_FMT_LEN;
	}

#endif
}

namespace expr
{

	template<typename G, typename V, typename E>
	auto make_map(V const& value, OpExpression<E> const& e)
	{
		return symphas::internal::make_map<G>::template get(value, *static_cast<E const*>(&e));
	}

	template<typename G, typename E>
	auto make_map(OpExpression<E> const& e)
	{
		return symphas::internal::make_map<G>::template get(*static_cast<E const*>(&e));
	}
}


//! An expression applying a GridPair type.
/*!
 * Uses a GridPair
 * transformation on the evaluated expression. Thus, also it needs to be 
 * updated before being evaluated.
 *
 * The primary difference between ::OpMap and ::OpFunctionApply is that OpMap is more
 * general and can deal with transformations over the whole grid, rather
 * than point-wise transformations which the function expression handles. 
 *
 * \tparam G The GridPair type which is used to transform the result
 * of the applied expression.
 * \tparam V The coefficient type.
 * \tparam E The expression type that is evaluated and then transformed.
 */
template<typename G, typename V, typename E>
struct OpMap : OpExpression<OpMap<G, V, E>>
{
	using result_type = expr::eval_type_t<E>;

	OpMap() : value{ V{} }, e{}, data{} {}

	//! Create a mapping expression.
	/*!
	 * Create an expression which maps the given expression through the
	 * the prescribed grid mapping type.
	 * 
	 * \param value The coefficient of the mapping expression.
	 * \param e The expression which is evaluated and mapped.
	 */
	OpMap(V value, E const& e) :
		value{ value }, e{ e }, data{ expr::data_dimensions(e) } {}

	inline auto eval(iter_type n) const
	{
		return expr::eval(value) * data[n];
	}

	//! Update the underlying data.
	/*!
	 * Evaluate the expression that the OpMap is applied to and store the
	 * result. Then update the GridPair object, which has been linked with the
	 * result in the constructor.
	 */
	void update()
	{
		expr::result(e, data.src, data.len);
		data.update();
	}

	auto operator-() const
	{
		return expr::make_map<G>(-value, e);
	}



#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return symphas::internal::print_map(*this, out);
	}

	size_t print(char* out) const
	{
		return symphas::internal::print_map(*this, out);
	}

	size_t print_length() const
	{
		return symphas::internal::print_map_length(*this);
	}

#endif

    template<typename G0, typename V0, typename E0>
	friend auto const& expr::get_enclosed_expression(OpMap<G0, V0, E0> const&);
    template<typename G0, typename V0, typename E0>
	friend auto& expr::get_enclosed_expression(OpMap<G0, V0, E0>&);

	V value;		//!< Coefficient of the map expression term.
	E e;			//!< Expression to which this operator applies.

protected:

	G data;

};


//! An expression applying a GridPair type.
/*!
 * Uses a GridPair
 * transformation on the evaluated expression. Thus, also it needs to be
 * updated before being evaluated.
 *
 * The primary difference between ::OpMap and ::OpFunctionApply is that OpMap is more
 * general and can deal with transformations over the whole grid, rather
 * than point-wise transformations which the function expression handles.
 *
 * \tparam G The GridPair type which is used to transform the result
 * of the applied expression.
 * \tparam V The coefficient type.
 * \tparam E The expression type that is evaluated and then transformed.
 */
template<typename V, typename E>
struct OpMap<void, V, E> : OpExpression<OpMap<void, V, E>>
{
	using result_type = expr::eval_type_t<E>;

	OpMap() : value{ V{} }, e{} {}

	//! Create a mapping expression.
	/*!
	 * Create an expression which maps the given expression through the
	 * the prescribed grid mapping type.
	 *
	 * \param value The coefficient of the mapping expression.
	 * \param e The expression which is evaluated and mapped.
	 */
	OpMap(V value, E const& e) :
		value{ value }, e{ e }
	{
		update();
	}

	inline auto eval(iter_type n) const
	{
		return expr::eval(value);
	}

	//! Update the underlying data.
	/*!
	 * Evaluate the expression that the OpMap is applied to and store the
	 * result. Then update the GridPair object, which has been linked with the
	 * result in the constructor.
	 */
	void update() {}

	auto operator-() const
	{
		return expr::make_map<void>(-value, e);
	}



#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return symphas::internal::print_map(*this, out);
	}

	size_t print(char* out) const
	{
		return symphas::internal::print_map(*this, out);
	}

	size_t print_length() const
	{
		return symphas::internal::print_map_length(*this);
	}

#endif

    template<typename V0, typename E0>
	friend auto const& expr::get_enclosed_expression(OpMap<void, V0, E0> const&);
    template<typename V0, typename E0>
	friend auto& expr::get_enclosed_expression(OpMap<void, V0, E0>&);

	V value;		//!< Coefficient of the map expression term.
	E e;			//!< Expression to which this operator applies.

};


template<>
struct symphas::internal::make_map<symphas::internal::HCTS>
{
	//! Constructs the map with the identity coefficient.
	template<typename E>
	static auto get(OpExpression<E> const& a);

	//! Constructs the map applied to an expression.
	template<typename = void>
	static auto get(OpVoid);

	//! Constructs the map applied to an expression.
	template<typename V, typename E>
	static auto get(V v, OpExpression<E> const& e);

	//! Constructs the map with the identity coefficient.
	template<typename E>
	static auto get(OpMap<symphas::internal::STHC, OpIdentity, E> const& a);
};


template<>
struct symphas::internal::make_map<symphas::internal::STHC>
{
	//! Constructs the map with the identity coefficient.
	template<typename E>
	static auto get(OpExpression<E> const& a);

	//! Constructs the map applied to an expression.
	template<typename = void>
	static auto get(OpVoid);

	//! Constructs the map applied to an expression.
	template<typename V, typename E>
	static auto get(V v, OpExpression<E> const& e);

	//! Constructs the map with the identity coefficient.
	template<typename E>
	static auto get(OpMap<symphas::internal::HCTS, OpIdentity, E> const& a);
};

namespace symphas::internal
{
	template<typename T, typename T0>
	T cast_to_result(T0&& result)
	{
		return std::forward<T0>(result);
	}

	template<>
	inline scalar_t cast_to_result<scalar_t, complex_t>(complex_t&& result)
	{
		return result.real();
	}
}

namespace expr
{

	//! Returns an index from the sequential grid using a half complex (FFTW format) grid.
	/*!
	 * \param n The index corresponding to the sequential grid.
	 * \param src The half complex grid the value is taken from.
	 * \param dims The underlying dimensions of the data.
	 */
	template<size_t D>
	struct eval_fftw_hcts 
	{
		template<typename E>
		auto operator()(iter_type n, OpExpression<E> const& e, const len_type* dims)
		{
			return static_cast<E const*>(&e)->eval(n);
		}
	};

	//! Specialization based on symphas::dft::get_fftw_hcts.
	template<>
	struct eval_fftw_hcts<1>
	{
		template<typename E>
		auto operator()(iter_type n, OpExpression<E> const& e, const len_type* dims) -> expr::eval_type_t<E>
		{
			using namespace expr;
			using symphas::math::conj;

			iter_type dn = dims[0] / 2 + 1;

			if (n < dn)
			{
				return static_cast<E const*>(&e)->eval(n);
			}
			else
			{
				return conj(static_cast<E const*>(&e)->eval(dn - (n - dn) - 2));
			}
		}
	};

	//! Specialization based on symphas::dft::get_fftw_hcts.
	template<>
	struct eval_fftw_hcts<2>
	{
		template<typename E>
		auto operator()(iter_type n, OpExpression<E> const& e, const len_type* dims) -> expr::eval_type_t<E>
		{
			using namespace expr;
			using symphas::math::conj;
			
			iter_type dn = dims[0] / 2 + 1;
			iter_type i0 = n % dims[0];
			iter_type j0 = n / dims[0];

			// Source half
			if (i0 < dn)
			{
				iter_type n0 = i0 + j0 * dn;
				return symphas::internal::cast_to_result<expr::eval_type_t<E>>(static_cast<E const*>(&e)->eval(n0));
			}
			// Hermitian half
			else
			{
				iter_type i00 = (dn - (i0 - dn) - 2);
				iter_type j00 = (j0 == 0) ? j0 : (dims[1] - j0);
				iter_type n00 = i00 + j00 * dn;
				return symphas::internal::cast_to_result<expr::eval_type_t<E>>(conj(static_cast<E const*>(&e)->eval(n00)));
			}
		}
	};

	//! Specialization based on symphas::dft::get_fftw_hcts.
	template<>
	struct eval_fftw_hcts<3>
	{
		template<typename E>
		auto operator()(iter_type n, OpExpression<E> const& e, const len_type* dims) -> expr::eval_type_t<E>
		{
			using namespace expr;
			using symphas::math::conj;

			iter_type dn = dims[0] / 2 + 1;
			iter_type i0 = n % dims[0];
			iter_type j0 = (n / dims[0]) % dims[1];
			iter_type k0 = n / (dims[0] * dims[1]);
			
			// Source half
			if (i0 < dn)
			{
				iter_type n0 = i0 + j0 * dn + k0 * dn * dims[1];
				return static_cast<E const*>(&e)->eval(n0);
			}
			// Hermitian half
			else
			{
				iter_type i00 = (dn - (i0 - dn) - 2);
				iter_type j00 = (j0 == 0) ? j0 : (dims[1] - j0);
				iter_type k00 = (k0 == 0) ? k0 : (dims[2] - k0);
				iter_type n00 = i00 + j00 * dn;
				return symphas::internal::cast_to_result<expr::eval_type_t<E>>(conj(static_cast<E const*>(&e)->eval(n00)));
			}
			
			if (i0 >= dn)
			{
				iter_type i00 = dims[0] - i0 - 1;
				iter_type j00 = dims[1] - j0 - 1;
				iter_type k00 = k0;
				return symphas::internal::cast_to_result<expr::eval_type_t<E>>
					(conj(static_cast<E const*>(&e)->eval(i00 + j00 * dn + k00 * dn * dims[1])));
			}
			else
			{
				return static_cast<E const*>(&e)->eval(i0 + j0 * dn + k0 * dn * dims[1]);
			}
		}
	};

	//! Returns an index from the half complex (FFTW format) grid using sequential grid.
	/*!
	 * \param n The index corresponding to the sequential grid.
	 * \param src The half complex grid the value is taken from.
	 * \param dims The underlying dimensions of the data.
	 */
	template<size_t D>
	struct eval_fftw_sthc 
	{
		template<typename E>
		auto operator()(iter_type n, OpExpression<E> const& e, const len_type* dims)
		{
			return static_cast<E const*>(&e)->eval(n);
		}
	};

	//! Specialization based on symphas::dft::get_fftw_sthc.
	template<>
	struct eval_fftw_sthc<1>
	{
		template<typename E>
		auto operator()(iter_type n, OpExpression<E> const& e, const len_type* dims)
		{
			return static_cast<E const*>(&e)->eval(n);
		}
	};

	//! Specialization based on symphas::dft::get_fftw_sthc.template<>
	template<>
	struct eval_fftw_sthc<2>
	{
		template<typename E>
		auto operator()(iter_type n, OpExpression<E> const& e, const len_type* dims)
		{
			iter_type dn = dims[0] / 2 + 1;
			iter_type i0 = n % dn;
			iter_type j0 = n / dn;

			return static_cast<E const*>(&e)->eval(i0 + j0 * dims[0]);
		}
	};

	//! Specialization based on symphas::dft::get_fftw_sthc.
	template<>
	struct eval_fftw_sthc<3>
	{
		template<typename E>
		auto operator()(iter_type n, OpExpression<E> const& e, const len_type* dims)
		{
			iter_type dn = dims[0] / 2 + 1;
			iter_type i0 = n % dn;
			iter_type j0 = (n / dn) % dims[1];
			iter_type k0 = n / (dn * dims[1]);

			return static_cast<E const*>(&e)->eval(i0 + j0 * dims[0] + k0 * dims[0] * dims[1]);
		}
	};
}
//
//template<typename E>
//struct OpExpression<OpMap<symphas::internal::HCTS, OpIdentity, E>> 
//{
//	using parent_type = OpMap<symphas::internal::HCTS, OpIdentity, E>;
//
//	explicit OpExpression(E const& rest) : e(rest) {}
//	explicit OpExpression(E&& rest) noexcept : e(std::move(rest)) {}
//
//	auto operator()(iter_type n) const
//	{
//		return cast().eval(n);
//	}
//
//	template<typename EE>
//	auto operator()(OpExpression<EE> const& e) const
//	{
//		return cast() * (*static_cast<EE const*>(&e));
//	}
//
//	symphas::internal::expression_iterator<parent_type> begin() const
//	{
//		return symphas::internal::expression_iterator<parent_type>(cast());
//	}
//
//	symphas::internal::expression_iterator<parent_type> end(len_type len) const
//	{
//		return symphas::internal::expression_iterator<parent_type>(cast(), len);
//	}
//
//	auto& cast() const
//	{
//		return *static_cast<parent_type const*>(this);
//	}
//
//	E e;
//};

//! Rearranges a complex-valued expression determined using FFTW algorithms.
/*!
 * Rearranges a complex-valued expression determined using FFTW algorithms.
 */
template<typename E>
struct OpMap<symphas::internal::HCTS, OpIdentity, E> : OpExpression<OpMap<symphas::internal::HCTS, OpIdentity, E>>
{
	//using parent_type = OpExpression<OpMap<symphas::internal::HCTS, OpIdentity, E>>;
	//using parent_type::e;

	using result_type = expr::eval_type_t<E>;
	static const size_t D = expr::grid_dim<E>::value;

	OpMap() : e{}, dims{} {}

	//! Create a mapping expression.
	/*!
	 * Create an expression which maps the given expression through the
	 * the prescribed grid mapping type.
	 *
	 * \param value The coefficient of the mapping expression.
	 * \param e The expression which is evaluated and mapped.
	 */
	OpMap(E const& e) : e{ e }, dims{} 
	{
		auto dims_copy = expr::data_dimensions(e);
		for (iter_type i = 0; i < D; ++i)
		{
			dims[i] = dims_copy[i];
		}
	}

	inline auto eval(iter_type n) const
	{
		return expr::eval_fftw_hcts<D>{}(n, e, dims);
	}

	void update() {}

	auto operator-() const
	{
		return expr::make_map<symphas::internal::HCTS>(-e);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return e.print(out);
	}

	size_t print(char* out) const
	{
		return e.print(out);
	}

	size_t print_length() const
	{
		return e.print_length();
	}

#endif

    template<typename E0>
	friend auto const& expr::get_enclosed_expression(OpMap<symphas::internal::HCTS, OpIdentity, E0> const&);
    template<typename E0>
	friend auto& expr::get_enclosed_expression(OpMap<symphas::internal::HCTS, OpIdentity, E0>&);

	E e;
	len_type dims[D];	//!< Dimensions of the expression.
};

//
//template<typename E>
//struct OpExpression<OpMap<symphas::internal::STHC, OpIdentity, E>> 
//{
//	using parent_type = OpMap<symphas::internal::STHC, OpIdentity, E>;
//
//	explicit OpExpression(E const& rest) : e(rest) {}
//	explicit OpExpression(E&& rest) noexcept : e(std::move(rest)) {}
//
//	auto operator()(iter_type n) const
//	{
//		return cast().eval(n);
//	}
//
//	template<typename EE>
//	auto operator()(OpExpression<EE> const& e) const
//	{
//		return cast() * (*static_cast<EE const*>(&e));
//	}
//
//	symphas::internal::expression_iterator<parent_type> begin() const
//	{
//		return symphas::internal::expression_iterator<parent_type>(cast());
//	}
//
//	symphas::internal::expression_iterator<parent_type> end(len_type len) const
//	{
//		return symphas::internal::expression_iterator<parent_type>(cast(), len);
//	}
//
//	auto& cast() const
//	{
//		return *static_cast<parent_type const*>(this);
//	}
//
//	E e;
//};


//! Rearranges a complex-valued expression determined using FFTW algorithms.
/*!
 * Rearranges a complex-valued expression determined using FFTW algorithms.
 */
template<typename E>
struct OpMap<symphas::internal::STHC, OpIdentity, E> : OpExpression<OpMap<symphas::internal::STHC, OpIdentity, E>>
{
	//using parent_type = OpExpression<OpMap<symphas::internal::STHC, OpIdentity, E>>;
	//using parent_type::e;

	using result_type = expr::eval_type_t<E>;
	static const size_t D = expr::grid_dim<E>::value;

	OpMap() : e{}, dims{} {}

	//! Create a mapping expression.
	/*!
	 * Create an expression which maps the given expression through the
	 * the prescribed grid mapping type.
	 *
	 * \param value The coefficient of the mapping expression.
	 * \param e The expression which is evaluated and mapped.
	 */
	OpMap(E const& e) : e{ e }, dims{}
	{
		auto dims_copy = expr::data_dimensions(e);
		for (iter_type i = 0; i < D; ++i)
		{
			dims[i] = dims_copy[i];
		}
	}

	inline auto eval(iter_type n) const
	{
		if (n >= symphas::dft::length<scalar_t, D>(dims))
		{
			n %= symphas::dft::length<scalar_t, D>(dims);
		}
		return expr::eval_fftw_sthc<D>{}(n, e, dims);
	}

	void update() {}

	auto operator-() const
	{
		return expr::make_map<symphas::internal::STHC>(-e);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return e.print(out);
	}

	size_t print(char* out) const
	{
		return e.print(out);
	}

	size_t print_length() const
	{
		return e.print_length();
	}

#endif

    template<typename E0>
	friend auto const& expr::get_enclosed_expression(OpMap<symphas::internal::HCTS, OpIdentity, E0> const&);
    template<typename E0>
	friend auto& expr::get_enclosed_expression(OpMap<symphas::internal::HCTS, OpIdentity, E0>&);

	E e;
	len_type dims[D];	//!< Dimensions of the expression.
};




template<typename coeff_t, typename G2, typename V2, typename E2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V2>), int> = 0>
auto operator*(coeff_t const& value, OpMap<G2, V2, E2> const& b)
{
	return expr::make_map<G2>(value * b.value, expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename tensor_t, typename G2, typename E2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpMap<G2, tensor_t, E2> const& b)
{
	return (value * b.value) * expr::make_map<G2>(expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(coeff_t const& value, OpMap<symphas::internal::HCTS, OpIdentity, E> const& b)
{
	return expr::make_map<symphas::internal::HCTS>(value * expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(coeff_t const& value, OpMap<symphas::internal::STHC, OpIdentity, E> const& b)
{
	return expr::make_map<symphas::internal::STHC>(value * expr::get_enclosed_expression(b));
}


//! Rearranges a complex-valued expression determined using FFTW algorithms.
/*!
 * Rearranges a complex-valued expression determined using FFTW algorithms.
 */
template<Axis ax, typename V, typename E>
struct OpMap<VectorComponent<ax>, V, E> : OpExpression<OpMap<VectorComponent<ax>, V, E>>
{
	OpMap() : e{} {}

	//! Create a mapping expression to get the component of an N-d expression.
	/*!
	 * Create an expression which maps the given expression to get the
	 * component of an N-d expression.
	 *
	 * \param value The coefficient of the mapping expression.
	 * \param e The expression which is evaluated and mapped.
	 */
	OpMap(V const& value, E const& e) : value{ value }, e { e } {}

	inline auto eval(iter_type n) const
	{
		return e.eval(n)[symphas::axis_to_index(ax)];
	}

	void update() {}

	auto operator-() const
	{
		return expr::to_axis<ax>(-e);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return e.print(out);
	}

	size_t print(char* out) const
	{
		return e.print(out);
	}

	size_t print_length() const
	{
		return e.print_length();
	}

#endif

	V value;
	E e;

	template<typename E0>
	friend auto const& expr::get_enclosed_expression(OpMap<symphas::internal::HCTS, OpIdentity, E0> const&);
	template<typename E0>
	friend auto& expr::get_enclosed_expression(OpMap<symphas::internal::HCTS, OpIdentity, E0>&);
};




namespace expr
{

	template<Axis ax, typename V, typename E>
	auto to_axis(V const& value, OpExpression<E> const& e)
	{
		return OpMap<VectorComponent<ax>, V, E>(value, *static_cast<E const*>(&e));
	}

	template<Axis ax, typename E>
	auto to_axis(OpExpression<E> const& e)
	{
		return to_axis<ax>(OpIdentity{}, *static_cast<E const*>(&e));;
	}

}



namespace symphas::internal
{
	template<typename T>
	struct real_space_type
	{
		using type = T;
	};

	template<typename T, size_t D>
	struct real_space_type<any_vector_t<T, D>>
	{
		using type = T;
	};
}


namespace expr
{
	namespace
	{
		template<typename... As, size_t... Is>
		auto hcts_adds(OpAdd<As...> const& e, std::index_sequence<Is...>)
		{
			return (hcts(expr::get<Is>(e)) + ...);
		}

		template<typename... As, size_t... Is>
		auto sthc_adds(OpAdd<As...> const& e, std::index_sequence<Is...>)
		{
			return (sthc(expr::get<Is>(e)) + ...);
		}
	}

	template<typename E>
	auto hcts(OpExpression<E> const& e)
	{
		if constexpr (expr::is_coeff<E>)
		{
			return *static_cast<E const*>(&e);
		}
		else
		{
			return expr::make_map<symphas::internal::HCTS>(*static_cast<E const*>(&e));
		}
	}

	template<typename E>
	auto sthc(OpExpression<E> const& e)
	{
		if constexpr (expr::is_coeff<E>)
		{
			return *static_cast<E const*>(&e);
		}
		else
		{
			return expr::make_map<symphas::internal::STHC>(*static_cast<E const*>(&e));
		}
	}

	//! Initializes a mapping expression based on the Fourier transformation.
	/*!
	 * Initializes an OpMap instance for the currently implemented grid pairs that
	 * can be used, which is currently limited to the MapGridFourier type. Will
	 * perform an arrangement based on FFTW scalar->complex array organization so
	 * that the result will be put in a sequential, or "full complex" format.
	 */
	template<typename V, typename E, size_t D = expr::grid_dim<E>::value, 
		size_t R = expr::eval_type<E>::rank, typename std::enable_if_t<(R == 0), int> = 0>
	auto make_fourier_map(V const& value, OpExpression<E> const& e)
	{
		static_assert(D > 0);

		using src_type = typename expr::eval_type<E>::type;
		using grid_pair_type = MapGridFourier<src_type, complex_t, D>;
		using rt = typename symphas::internal::real_space_type<src_type>::type;

		if constexpr (std::is_same<rt, scalar_t>::value)
		{
			return hcts(expr::make_map<grid_pair_type>(value, *static_cast<E const*>(&e)));
		}
		else
		{
			return expr::make_map<grid_pair_type>(value, *static_cast<E const*>(&e));
		}
	}

	template<typename V>
	auto make_fourier_map(V const& value, OpVoid)
	{
		return OpVoid{};
	}

	template<typename V, typename E, size_t D = expr::grid_dim<E>::value, size_t... Rs, size_t R = sizeof...(Rs)>
	auto make_fourier_map(V const& value, OpExpression<E> const& e, std::index_sequence<Rs...>)
	{
		return (make_fourier_map(expr::make_column_vector<Rs, R>() * value, expr::make_row_vector<Rs, R>() * (*static_cast<E const*>(&e))) + ...);
	}

	template<typename V, typename E, size_t D = expr::grid_dim<E>::value,
		size_t R = expr::eval_type<E>::rank, typename std::enable_if_t<(R > 0), int> = 0>
	auto make_fourier_map(V const& value, OpExpression<E> const& e)
	{
		return make_fourier_map(value, *static_cast<E const*>(&e), std::make_index_sequence<R>{});
	}

	template<typename E>
	auto make_fourier_map(OpExpression<E> const& e)
	{
		return expr::make_fourier_map(OpIdentity{}, *static_cast<E const*>(&e));
	}

	//! Initializes a mapping expression based on the Fourier transformation.
	/*!
	 * Initializes an OpMap instance for the currently implemented grid pairs that
	 * can be used, which is currently limited to the MapGridFourier type.
	 */
	template<typename target_type, typename V, typename E, size_t D = expr::grid_dim<E>::value,
		size_t R = expr::eval_type<E>::rank, typename std::enable_if_t<(R == 0), int> = 0>
	auto make_inv_fourier_map(V const& value, OpExpression<E> const& e)
	{
		static_assert(D > 0);
		using rt = typename symphas::internal::real_space_type<target_type>::type;

		using grid_pair_type = MapGridInverseFourier<complex_t, rt, D>;

		if constexpr (std::is_same<rt, scalar_t>::value)
		{
			return expr::make_map<grid_pair_type>(value, sthc(*static_cast<E const*>(&e)));
		}
		else
		{
			return expr::make_map<grid_pair_type>(value, *static_cast<E const*>(&e));
		}
	}

	template<typename target_type, typename V>
	auto make_inv_fourier_map(V const& value, OpVoid)
	{
		return OpVoid{};
	}


	template<typename target_type, typename V, typename E, size_t D = expr::grid_dim<E>::value, size_t... Rs, size_t R = sizeof...(Rs)>
	auto make_inv_fourier_map(V const& value, OpExpression<E> const& e, std::index_sequence<Rs...>)
	{
		return (make_inv_fourier_map<target_type>(expr::make_column_vector<Rs, R>() * value, 
			expr::make_row_vector<Rs, R>() * (*static_cast<E const*>(&e))) + ...);
	}

	template<typename target_type, typename V, typename E, size_t D = expr::grid_dim<E>::value,
		size_t R = expr::eval_type<E>::rank, typename std::enable_if_t<(R > 0), int> = 0>
	auto make_inv_fourier_map(V const& value, OpExpression<E> const& e)
	{
		return make_inv_fourier_map<target_type>(value, *static_cast<E const*>(&e), std::make_index_sequence<R>{});
	}


	//! Initializes a mapping expression based on the Fourier transformation.
	/*!
	 * Initializes an OpMap instance for the currently implemented grid pairs that
	 * can be used, which is currently limited to the MapGridFourier type.
	 */
	template<typename target_type, typename E>
	auto make_inv_fourier_map(OpExpression<E> const& e)
	{
		return expr::make_inv_fourier_map<target_type>(OpIdentity{}, *static_cast<E const*>(&e));
	}
}


template<typename G>
template<typename V, typename E>
auto symphas::internal::make_map<G>::get(V value, OpExpression<E> const& e)
{
	return OpMap<G, V, E>(value, *static_cast<E const*>(&e));
}

template<typename G>
template<typename A>
auto symphas::internal::make_map<G>::get(A&& a)
{
	return get(OpIdentity{}, std::forward<A>(a));
}

template<typename>
 auto symphas::internal::make_map<symphas::internal::HCTS>::get(OpVoid)
{
	return OpVoid{};
}

template<typename V, typename E>
auto symphas::internal::make_map<symphas::internal::HCTS>::get(V value, OpExpression<E> const& e)
{
	return get(value * *static_cast<E const*>(&e));
}

template<typename E>
auto symphas::internal::make_map<symphas::internal::HCTS>::get(OpExpression<E> const& e)
{
	return OpMap<symphas::internal::HCTS, OpIdentity, E>(*static_cast<E const*>(&e));
}

template<typename E>
auto symphas::internal::make_map<symphas::internal::HCTS>::get(OpMap<symphas::internal::STHC, OpIdentity, E> const& e)
{
	return expr::get_enclosed_expression(e);
}

template<typename>
inline auto symphas::internal::make_map<symphas::internal::STHC>::get(OpVoid)
{
	return OpVoid{};
}

template<typename V, typename E>
auto symphas::internal::make_map<symphas::internal::STHC>::get(V value, OpExpression<E> const& e)
{
	return get(value * *static_cast<E const*>(&e));
}

template<typename E>
auto symphas::internal::make_map<symphas::internal::STHC>::get(OpExpression<E> const& e)
{
	return OpMap<symphas::internal::STHC, OpIdentity, E>(*static_cast<E const*>(&e));
}

template<typename E>
auto symphas::internal::make_map<symphas::internal::STHC>::get(OpMap<symphas::internal::HCTS, OpIdentity, E> const& e)
{
	return expr::get_enclosed_expression(e);
}




template<typename source_type, typename V, typename E>
using OpFourierTransform = OpMap<MapGridFourier<source_type, complex_t, expr::grid_dim<E>::value>, V, E>;

template<typename target_type, typename V, typename E>
using OpInverseFourierTransform = OpMap<MapGridInverseFourier<complex_t, target_type, expr::grid_dim<E>::value>, V, E>;



// Overloads to allow symbolic rules to work betwewen STHC expressions.

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator+(OpMap<symphas::internal::STHC, OpIdentity, E> const& a, coeff_t const& b)
{
	return expr::sthc(expr::get_enclosed_expression(a) + b);
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator+(coeff_t const& a, OpMap<symphas::internal::STHC, OpIdentity, E> const& b)
{
	return expr::sthc(a + expr::get_enclosed_expression(b));
}

template<typename E1, typename E2>
auto operator+(OpMap<symphas::internal::STHC, OpIdentity, E1> const& a, OpMap<symphas::internal::STHC, OpIdentity, E2> const& b)
{
	return expr::sthc(expr::get_enclosed_expression(a) + expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator-(OpMap<symphas::internal::STHC, OpIdentity, E> const& a, coeff_t const& b)
{
	return expr::sthc(expr::get_enclosed_expression(a) - b);
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator-(coeff_t const& a, OpMap<symphas::internal::STHC, OpIdentity, E> const& b)
{
	return expr::sthc(a - expr::get_enclosed_expression(b));
}

template<typename E1, typename E2>
auto operator-(OpMap<symphas::internal::STHC, OpIdentity, E1> const& a, OpMap<symphas::internal::STHC, OpIdentity, E2> const& b)
{
	return expr::sthc(expr::get_enclosed_expression(a) - expr::get_enclosed_expression(b));
}

template<typename E1, typename E2>
auto operator*(OpMap<symphas::internal::STHC, OpIdentity, E1> const& a, OpMap<symphas::internal::STHC, OpIdentity, E2> const& b)
{
	return expr::sthc(expr::get_enclosed_expression(a) * expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator/(OpMap<symphas::internal::STHC, OpIdentity, E> const& a, coeff_t const& b)
{
	return expr::sthc(expr::get_enclosed_expression(a) / b);
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator/(coeff_t const& a, OpMap<symphas::internal::STHC, OpIdentity, E> const& b)
{
	return expr::sthc(a / expr::get_enclosed_expression(b));
}

template<typename E1, typename E2>
auto operator/(OpMap<symphas::internal::STHC, OpIdentity, E1> const& a, OpMap<symphas::internal::STHC, OpIdentity, E2> const& b)
{
	return expr::sthc(expr::get_enclosed_expression(a) / expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator+(OpMap<symphas::internal::HCTS, OpIdentity, E> const& a, coeff_t const& b)
{
	return expr::hcts(expr::get_enclosed_expression(a) + b);
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator+(coeff_t const& a, OpMap<symphas::internal::HCTS, OpIdentity, E> const& b)
{
	return expr::hcts(a + expr::get_enclosed_expression(b));
}

template<typename E1, typename E2>
auto operator+(OpMap<symphas::internal::HCTS, OpIdentity, E1> const& a, OpMap<symphas::internal::HCTS, OpIdentity, E2> const& b)
{
	return expr::hcts(expr::get_enclosed_expression(a) + expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator-(OpMap<symphas::internal::HCTS, OpIdentity, E> const& a, coeff_t const& b)
{
	return expr::hcts(expr::get_enclosed_expression(a) - b);
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator-(coeff_t const& a, OpMap<symphas::internal::HCTS, OpIdentity, E> const& b)
{
	return expr::hcts(a - expr::get_enclosed_expression(b));
}

template<typename E1, typename E2>
auto operator-(OpMap<symphas::internal::HCTS, OpIdentity, E1> const& a, OpMap<symphas::internal::HCTS, OpIdentity, E2> const& b)
{
	return expr::hcts(expr::get_enclosed_expression(a) - expr::get_enclosed_expression(b));
}

template<typename E1, typename E2>
auto operator*(OpMap<symphas::internal::HCTS, OpIdentity, E1> const& a, OpMap<symphas::internal::HCTS, OpIdentity, E2> const& b)
{
	return expr::hcts(expr::get_enclosed_expression(a) * expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator/(OpMap<symphas::internal::HCTS, OpIdentity, E> const& a, coeff_t const& b)
{
	return expr::hcts(expr::get_enclosed_expression(a) / b);
}

template<typename coeff_t, typename E,
	typename std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator/(coeff_t const& a, OpMap<symphas::internal::HCTS, OpIdentity, E> const& b)
{
	return expr::hcts(a / expr::get_enclosed_expression(b));
}

template<typename E1, typename E2>
auto operator/(OpMap<symphas::internal::HCTS, OpIdentity, E1> const& a, OpMap<symphas::internal::HCTS, OpIdentity, E2> const& b)
{
	return expr::hcts(expr::get_enclosed_expression(a) / expr::get_enclosed_expression(b));
}

#include "expressionexponentials.h"

namespace expr
{
	template<typename E>
	auto exp(OpMap<symphas::internal::STHC, OpIdentity, E> const& e)
	{
		return sthc(exp(get_enclosed_expression(e)));
	}

	template<typename E>
	auto exp(OpMap<symphas::internal::HCTS, OpIdentity, E> const& e)
	{
		return hcts(exp(get_enclosed_expression(e)));
	}

	template<typename E>
	constexpr auto coeff(OpMap<symphas::internal::STHC, OpIdentity, E> const& e)
	{
		return coeff(e.e);
	}

	template<typename E>
	constexpr auto coeff(OpMap<symphas::internal::HCTS, OpIdentity, E> const& e)
	{
		return coeff(e.e);
	}
}


#undef SYEX_COMBINATION_FMT
#undef SYEX_COMBINATION_FMT_LEN
#undef SYEX_COMBINATION_FMT_LATEX
#undef SYEX_COMBINATION_FMT_LATEX_LEN
#undef SYEX_CHAIN_FMT
#undef SYEX_CHAIN_FMT_LEN
#undef SYEX_CHAIN_FMT_LATEX
#undef SYEX_CHAIN_FMT_LATEX_LEN
#undef SYEX_OP_APPLY_FMT
#undef SYEX_OP_APPLY_FMT_LEN
#undef SYEX_OP_APPLY_FMT_LATEX
#undef SYEX_OP_APPLY_FMT_LATEX_LEN
#undef SYEX_OP_MAP_FMT
#undef SYEX_OP_MAP_FMT_LEN
#undef SYEX_OP_MAP_FMT_LATEX
#undef SYEX_OP_MAP_FMT_LATEX_LEN



