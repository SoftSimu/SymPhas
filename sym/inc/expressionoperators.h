
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

#include "expressionprototypes.h"
#include "gridpair.h"

//! \cond

#ifdef LATEX_PLOT
#define SYEX_COMBINATION_FMT_A "\\left("
#define SYEX_COMBINATION_FMT_B "\\right)"
#define SYEX_COMBINATION_FMT_SEP " + "
#define SYEX_CHAIN_FMT_SEP " \\cdot "

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

#define SYEX_MAP_FMT SYEX_OP_MAP_FMT_A "%s" SYEX_OP_MAP_FMT_B
#define SYEX_MAP_FMT_LEN (STR_ARR_LEN(SYEX_OP_MAP_FMT_A SYEX_OP_MAP_FMT_B) - 1)


#ifdef LATEX_PLOT
#define SYEX_FT_OF_EXPR_FMT_A "\\hat{\\mathcal{F}}\\left\\{"
#define SYEX_FT_OF_EXPR_FMT_B "\\right\\}"

#define SYEX_FT_OF_OP_FMT_A "\\hat{"
#define SYEX_FT_OF_OP_FMT_B "}"
#else
#define SYEX_FT_OF_EXPR_FMT_A "F^{"
#define SYEX_FT_OF_EXPR_FMT_B "}"
#define SYEX_FT_OF_OP_FMT_A SYEX_FT_OF_EXPR_FMT_A
#define SYEX_FT_OF_OP_FMT_B SYEX_FT_OF_EXPR_FMT_B
#endif

#define SYEX_FT_OF_EXPR_FMT SYEX_FT_OF_EXPR_FMT_A "%s" SYEX_FT_OF_EXPR_FMT_B
#define SYEX_FT_OF_EXPR_FMT_LEN (STR_ARR_LEN(SYEX_FT_OF_EXPR_FMT_A SYEX_FT_OF_EXPR_FMT_B) - 1)
#define SYEX_FT_OF_OP_FMT SYEX_FT_OF_OP_FMT_A "%s" SYEX_FT_OF_OP_FMT_B
#define SYEX_FT_OF_OP_FMT_LEN (STR_ARR_LEN(SYEX_FT_OF_OP_FMT_A SYEX_FT_OF_OP_FMT_B) - 1)

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
 * for only ::OpLVariable, not just to avoid wasting time memory, but because
 * the combination implementation will apply operators to variables and
 * will NOT update each applied operator before evaluating the combination result. Therefore,
 * each operator needs to have an implementation where it computes the results
 * directly of a variable without any intermediate data that would need to be
 * updated.
 * 
 * \tparam E The specialized operator, used for the CRTP strategy.
 */
template<typename E>
struct OpOperator : OpExpression<E> {};



//! Specialized behaviour applying an operator to a combination operator.
template<typename E, typename B1, typename B2>
auto operator*(OpOperator<E> const& a, OpOperatorCombination<B1, B2> const& b)
{
	return OpOperatorChain(*static_cast<E const*>(&a), b);
}

//! Specialized behaviour adding an operator to a combination operator.
template<typename E, typename B1, typename B2>
auto operator+(OpOperator<E> const& a, OpOperatorCombination<B1, B2> const& b)
{
	return OpOperatorCombination(*static_cast<E const*>(&a), b);
}

//! Specialized behaviour subtracting operator to a combination operator.
template<typename E, typename B1, typename B2>
auto operator-(OpOperator<E> const& a, OpOperatorCombination<B1, B2> const& b)
{
	return OpOperatorCombination(*static_cast<E const*>(&a), -b);
}

//! Specialized behaviour applying an operator to a chain operator.
template<typename E, typename B1, typename B2>
auto operator*(OpOperator<E> const& a, OpOperatorChain<B1, B2> const& b)
{
	return OpOperatorChain(*static_cast<E const*>(&a), b);
}

//! Specialized behaviour adding an operator to a chain operator.
template<typename E, typename B1, typename B2>
auto operator+(OpOperator<E> const& a, OpOperatorChain<B1, B2> const& b)
{
	return OpOperatorCombination(*static_cast<E const*>(&a), b);
}

//! Specialized behaviour subtracting operator to a chain operator.
template<typename E, typename B1, typename B2>
auto operator-(OpOperator<E> const& a, OpOperatorChain<B1, B2> const& b)
{
	return OpOperatorCombination(*static_cast<E const*>(&a), -b);
}

/*
 * applied to another combination or chain object
 */

//! An operator applied to a combination creates a chain operator.
template<typename E, typename B1, typename B2, typename F>
auto operator*(OpOperator<E> const& a, OpCombination<B1, B2, F> const& b)
{
	return OpOperatorChain(*static_cast<E const*>(&a), b.combination)* expr::compound_get::template expr(b);
}

//! An operator applied to a chain builds the chain operator.
template<typename E, typename B1, typename B2, typename F>
auto operator*(OpOperator<E> const& a, OpChain<B1, B2, F> const& b)
{
	return OpOperatorChain(*static_cast<E const*>(&a), b.combination)* expr::compound_get::template expr(b);
}


//! An operator applied to another one creates a chain operator.
template<typename E, typename F>
auto operator*(OpOperator<E> const& a, OpOperator<F> const& b)
{
	return OpOperatorChain(*static_cast<E const*>(&a), *static_cast<F const*>(&b));
}

//! The addition of two operators creates a combination.
template<typename E, typename F>
auto operator+(OpOperator<E> const& a, OpOperator<F> const& b)
{
	return OpOperatorCombination(*static_cast<E const*>(&a), *static_cast<F const*>(&b));
}

//! The subtraction of two operators creates a combination.
template<typename E, typename F>
auto operator-(OpOperator<E> const& a, OpOperator<F> const& b)
{
	return OpOperatorCombination(*static_cast<E const*>(&a), -*static_cast<F const*>(&b));
}


/* rules when used with expressions
 */


 //! Adding an operator to a constant creates a combination.
template<typename A>
auto operator+(OpOperator<A> const& a, OpIdentity const)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), OpIdentity{});
}

//! Subtracting an operator from a constant creates a combination.
template<typename A>
auto operator-(OpOperator<A> const& a, OpIdentity const)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), OpIdentity{});
}

//! Adding an operator to a constant creates a combination.
template<typename A>
auto operator+(OpOperator<A> const& a, OpNegIdentity const)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), OpNegIdentity{});
}

//! Subtracting an operator from a constant creates a combination.
template<typename A>
auto operator-(OpOperator<A> const& a, OpNegIdentity const)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), OpNegIdentity{});
}


//! Adding a constant to an operator creates a combination.
template<typename A>
auto operator+(OpIdentity const, OpOperator<A> const& b)
{
	return OpOperatorCombination(OpIdentity{}, *static_cast<A const*>(&b));
}

//! Subtracting a constant from an operator creates a combination.
template<typename A>
auto operator-(OpIdentity const, OpOperator<A> const& b)
{
	return OpOperatorCombination(OpIdentity{}, *static_cast<A const*>(&b));
}

//! Adding a constant to an operator creates a combination.
template<typename A>
auto operator+(OpNegIdentity const, OpOperator<A> const& b)
{
	return OpOperatorCombination(OpNegIdentity{}, *static_cast<A const*>(&b));
}

//! Subtracting a constant from an operator creates a combination.
template<typename A>
auto operator-(OpNegIdentity const, OpOperator<A> const& b)
{
	return OpOperatorCombination(OpNegIdentity{}, *static_cast<A const*>(&b));
}


//! Adding a constant to an operator creates a combination.
template<typename A, typename T>
auto operator+(OpOperator<A> const& a, OpLiteral<T> const& b)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), b);
}

//! Subtracting a constant from an operator creates a combination.
template<typename A, typename T>
auto operator-(OpOperator<A> const& a, OpLiteral<T> const& b)
{
	return OpOperatorCombination(*static_cast<A const*>(&a), b);
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

//! The operator as the second term of a multiplication is applied.
template<typename A1, typename A2, typename E>
auto operator*(OpBinaryMul<A1, OpOperator<E>> const& a, OpExpression<A2> const& b)
{
	return a.a * (a.b * (*static_cast<const A2*>(&b)));
}



// ******************************************************************************************


//! An operator is applied to an expression through multiplication.
template<typename E, typename F>
auto operator*(OpOperator<E> const& a, OpExpression<F> const& b)
{
	return (*static_cast<E const*>(&a)).apply(*static_cast<F const*>(&b));
}

//! An operator is applied to an expression through multiplication.
template<typename E>
auto operator*(OpOperator<E> const& a, OpIdentity const)
{
	return (*static_cast<E const*>(&a)).apply(OpIdentity{});
}

//! An operator is applied to an expression through multiplication.
template<typename E>
auto operator*(OpOperator<E> const& a, OpNegIdentity const)
{
	return -((*static_cast<E const*>(&a)) * OpIdentity{});
}

//! An operator is applied to an expression through multiplication.
template<typename E, typename T>
auto operator*(OpOperator<E> const& a, OpLiteral<T> const& v)
{
	return v * ((*static_cast<E const*>(&a)) * OpIdentity{});
}


// ******************************************************************************************

namespace expr
{

	//! The combination operator is applied to terms in a binary addition.
	template<typename A1, typename A2, typename E1, typename E2>
	auto distribute_operator(OpOperatorCombination<A1, A2> const& a, OpBinaryAdd<E1, E2> const& b)
	{
		return expr::distribute_operator(a, b.a) + expr::distribute_operator(a, b.b);
	}

	//! The combination operator is applied to terms in a binary subtraction.
	template<typename A1, typename A2, typename E1, typename E2>
	auto distribute_operator(OpOperatorCombination<A1, A2> const& a, OpBinarySub<E1, E2> const& b)
	{
		return expr::distribute_operator(a, b.a) - expr::distribute_operator(a, b.b);
	}

	//! The expression is applied to the combination operator.
	template<typename A1, typename A2, typename E>
	auto distribute_operator(OpOperatorCombination<A1, A2> const& a, OpExpression<E> const& b)
	{
		return a * (*static_cast<E const*>(&b));
	}


	//! The chain operator is applied to terms in a binary addition.
	template<typename A1, typename A2, typename E1, typename E2>
	auto distribute_operator(OpOperatorChain<A1, A2> const& a, OpBinaryAdd<E1, E2> const& b)
	{
		return expr::distribute_operator(a, b.a) + expr::distribute_operator(a, b.b);
	}

	//! The chain operator is applied to terms in a binary subtraction.
	template<typename A1, typename A2, typename E1, typename E2>
	auto distribute_operator(OpOperatorChain<A1, A2> const& a, OpBinarySub<E1, E2> const& b)
	{
		return expr::distribute_operator(a, b.a) - expr::distribute_operator(a, b.b);
	}

	//! The expression is applied to the chain operator.
	template<typename A1, typename A2, typename E>
	auto distribute_operator(OpOperatorChain<A1, A2> const& a, OpExpression<E> const& b)
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

	//! Expanding a combination into an expression.
	/*!
	 * Applies each of the constituent operators to the given expression.
	 * In this way, this distributes the operators to the expression.
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
	 * \param a The chain operator which is applied to an expression.
	 * \param b The expression applied by the chain operator.
	 */
	template<typename A1, typename A2, typename E>
	auto expand_operator(OpOperatorChain<A1, A2> const& a, OpExpression<E> const& b)
	{
		return expand_operator(a.g, expand_operator(a.f, *static_cast<E const*>(&b)));
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
	//! Create the combination of two operators.
	/*!
	 * Create the combination of two operators.
	 * 
	 * \param f The operator on the left hand side of the addition.
	 * \param g The operator on the right hand side of the addition.
	 */
	OpOperatorCombination(A1 const& f, A2 const& g) : f{ f }, g{ g } {}

	inline auto eval(iter_type) const {}

	auto operator-() const
	{
		return ::OpOperatorCombination(-f, -g);
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

#endif

	//! Apply this operator an expression.
	/*!
	 * Apply this operator an expression and return the concrete form of the
	 * linear combination expression.
	 * 
	 * \param a The expression to which this operator is applied.
	 */
	template<typename E>
	auto apply(OpExpression<E> const& a) const
	{
		return OpCombination(*this, *static_cast<E const*>(&a));
	}

	A1 f; //!< Operator on the left of the plus sign.
	A2 g; //!< Operator on the right of the plus sign.
};



template<typename T, typename A1, typename A2>
auto operator*(OpLiteral<T> const& a, OpOperatorCombination<A1, A2> const& b)
{
	return OpOperatorCombination(a * b.f, a * b.g);
}



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
	using G = typename expr::grid_type<E>::type;

	OpOperatorCombination<A1, A2> combination;	//!< The combination operator.
	G data;										//!< Data which stores the result of the expression.

protected:


	auto make_eval_expr()
	{
		return combination.f * expr::make_op(data) + combination.g * expr::make_op(data);
	}

	using eval_expr_type = std::invoke_result_t<decltype(&OpCombination<A1, A2, E>::make_eval_expr), OpCombination<A1, A2, E>>;

public:

	//! Create the combination of two operators applied to an expression.
	/*!
	 * Create the combination of two operators applied to an expression.
	 *
	 * \param combination The operator being applied.
	 * \param e The expression the operator is applied to.
	 */
	OpCombination(OpOperatorCombination<A1, A2> const& combination, E const& e) :
		combination{ combination }, data{ expr::property::data_dimensions(e) }, e{ e }, eval_expr{ make_eval_expr() } { update(); }

	inline auto eval(iter_type n) const
	{
		return eval_expr.eval(n);
	}

	//! Save the result of the expression in memory to be evaluated.
	void update()
	{
		expr::result(e, data.values, data.len);
	}

	auto operator-() const
	{
		return combination * -e;
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

	friend struct expr::compound_get;



protected:

	E e;							//!< Expression to which this operator applies.
	eval_expr_type eval_expr;		//!< The variable to apply the system.
};


//! See OpCombination.
template<typename A1, typename A2, typename T, typename G>
struct OpCombination<A1, A2, OpLVariable<T, G>> : OpExpression<OpCombination<A1, A2, OpLVariable<T, G>>>
{
	using E = OpLVariable<T, G>;

	OpOperatorCombination<A1, A2> combination;	//!< The combination operator.
	OpLVariable<T, G> e;						//!< Variable to which this operator applies.

protected:


	auto make_eval_expr()
	{
		return combination.f * e + combination.g * e;
	}

	using eval_expr_type = std::invoke_result_t<decltype(&OpCombination<A1, A2, E>::make_eval_expr), OpCombination<A1, A2, E>>;

public:

	//! Create the combination of two operators applied to a variable.
	/*!
	 * Create the combination of two operators applied to a variable.
	 *
	 * \param combination The operator being applied.
	 * \param e The variable the operator is applied to.
	 */
	OpCombination(OpOperatorCombination<A1, A2> const& combination, E const& e) : combination{ combination }, e{ e }, eval_expr{ make_eval_expr() } {}

	inline auto eval(iter_type n) const
	{
		return eval_expr.eval(n);
	}

	void update() {}

	auto operator-() const
	{
		return combination * -e;
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

	friend struct expr::compound_get;


	eval_expr_type eval_expr;					//!< The variable to apply the system.


};



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
	//! Create the chain of two operators.
	/*!
	 * Create the chain of two operators.
	 *
	 * \param f The operator on the left hand side of the addition.
	 * \param g The operator on the right hand side of the addition.
	 */
	OpOperatorChain(A1 const& f, A2 const& g) : f{ f }, g{ g } {}

	inline auto eval(iter_type) const {}

	auto operator-() const
	{
		return OpOperatorChain(-f, g);
	}


	//! Apply the chain operation to an expression.
	template<typename E>
	auto apply(OpExpression<E> const& a) const
	{
		return OpChain(*this, *static_cast<E const*>(&a));
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = 0;
		n += fprintf(out, SYEX_CHAIN_FMT_A);
		n += f.print(out);
		n += fprintf(out, SYEX_CHAIN_FMT_SEP);
		n += g.print(out);
		n += fprintf(out, SYEX_CHAIN_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = 0;
		n += sprintf(out + n, SYEX_CHAIN_FMT_A);
		n += f.print(out + n);
		n += sprintf(out + n, SYEX_CHAIN_FMT_SEP);
		n += g.print(out + n);
		n += sprintf(out + n, SYEX_CHAIN_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return f.print_length() + g.print_length() + SYEX_CHAIN_FMT_LEN;
	}

#endif


	A1 f; //!< Operator which is applied to the result of operator `g` on an expression.
	A2 g; //!< Operator which is applied first.
};


/* when multiplied by a literal
 */

template<typename T, typename A1, typename A2>
auto operator*(OpLiteral<T> const& a, OpOperatorChain<A1, A2> const& b)
{
	return OpOperatorChain(a * b.f, b.g);
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
	using inner_type = mul_result_t<A2, E>;

public:
	inner_type inner; //!< The result of applying the inner operator to the expression.

	auto inner_grid()
	{
		return expr::make_op(expr::compound_get::template grid(inner));
	}

	using inner_var = typename std::invoke_result_t<decltype(&OpChain<A1, A2, E>::inner_grid), OpChain<A1, A2, E>>;
	using outer_type = mul_result_t<A1, inner_var>;

public:

	outer_type outer; //!< The result of applying the outer operator to the inner.

	//! Create the combination of two operators applied to an expression.
	/*!
	 * Create the chain of two operators applied to an expression.
	 *
	 * \param combination The operator being applied.
	 * \param e The expression the operator is applied to.
	 */
	OpChain(OpOperatorChain<A1, A2> const& combination, E const& e) : 
		combination{ combination }, e{ e }, inner{ combination.g * e }, outer{ combination.f * expr::make_op(expr::compound_get::template grid(inner)) } {}

	inline auto eval(iter_type n) const
	{
		return outer.eval(n);
	}

	auto operator-() const
	{
		return combination * -e;
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

#endif


	friend struct expr::compound_get;


	OpOperatorChain<A1, A2> combination;	//!< The chain operator.
	E e;									//!< Expression to which this operator applies.

};


//! A concrete operator chain of two general operators.
/*!
 * See OpChain. Specializes applying a chain operator when the chain operator
 * is applied to a variable.
 *
 * \tparam A1 Type of the first operator.
 * \tparam A2 Type of the second operator.
 * \tparam T The type of the variable .
 */
template<typename A1, typename A2, typename T, typename G>
struct OpChain<A1, A2, OpLVariable<T, G>> : OpExpression<OpChain<A1, A2, OpLVariable<T, G>>>
{
	using E = OpLVariable<T, G>;

protected:
	using inner_type = mul_result_t<A2, E>;

public:
	inner_type inner; //!< The result of applying the inner operator to the variable.

protected:
	auto inner_grid()
	{
		return expr::make_op(expr::compound_get::template grid(inner));
	}

	using inner_var = typename std::invoke_result_t<decltype(&OpChain<A1, A2, OpLVariable<T, G>>::inner_grid), OpChain<A1, A2, OpLVariable<T, G>>>;
	using outer_type = mul_result_t<A1, inner_var>;

public:

	outer_type outer; //!< The result of applying the outer operator to the result of the inner.


	//! Create the combination of two operators applied to an expression.
	/*!
	 * Create the chain of two operators applied to an expression.
	 *
	 * \param combination The operator being applied.
	 * \param e The expression the operator is applied to.
	 */
	OpChain(OpOperatorChain<A1, A2> const& combination, E const& e) : combination{ combination }, e{ e }, inner{ combination.g * e }, outer{ combination.f * inner_grid() } {}
	inline auto eval(iter_type n) const
	{
		return outer.eval(n);
	}

	auto operator-() const
	{
		return combination * -e;
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

#endif

	friend struct expr::compound_get;



	OpOperatorChain<A1, A2> combination;	//!< The chain operator.
	OpLVariable<T, G> e;					//!< Variable to which this operator applies.


};


// *********************************************************************************************************************************

namespace symphas::internal
{

	//! Implementation of map expression construction.
	/*!
	 * Implementation of functions which generalize and simplify the way to
	 * construct map expressions. Wraps the template deduction necessary
	 * to initialize a map expression.
	 */
	struct make_map
	{
		//! Constructs the map with the identity coefficient.
		template<typename G, typename A>
		static auto get(A&& a);

		//! Constructs the map applied to an expression.
		template<typename G, typename V, typename E>
		static auto get(V v, OpExpression<E> const& e);
	};

}




//! An expression applying a GridPair type.
/*!
 * Uses a GridPair
 * transformation on the evaluated expression. Thus, also it needs to be 
 * updated before being evaluated.
 *
 * The primary difference between ::OpMap and ::OpFuncApply is that OpMap is more
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
	using result_grid = typename expr::grid_type<E>::type;

	//! Create a mapping expression.
	/*!
	 * Create an expression which maps the given expression through the
	 * the prescribed grid mapping type.
	 * 
	 * \param value The coefficient of the mapping expression.
	 * \param e The expression which is evaluated and mapped.
	 */
	OpMap(V value, E const& e) : value{ value }, e{ e }, result{ expr::property::data_dimensions(e) }, data{ result } { update(); }

	inline auto eval(iter_type n) const
	{
		return value * data[n];
	}

	//! Update the underlying data.
	/*!
	 * Evaluate the expression that the OpMap is applied to and store the
	 * result. Then update the GridPair object, which has been linked with the
	 * result in the constructor.
	 */
	void update()
	{
		expr::result(e, result, result.len);
		data.update();
	}


	auto operator-() const
	{
		return symphas::internal::make_map::get(-value, e);
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += fprintf(out, SYEX_OP_MAP_FMT_A);
		n += e.print(out);
		n += fprintf(out, SYEX_OP_MAP_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += sprintf(out + n, SYEX_OP_MAP_FMT_A);
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_OP_MAP_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + e.print_length()
			+ SYEX_CHAIN_APPLY_FMT_LEN;
	}

#endif


	friend struct expr::compound_get;

	V value;		//!< Coefficient of the map expression term.
	E e;			//!< Expression to which this operator applies.

protected:

	result_grid result;
	G data;

};

template<typename S1, typename G, typename V, typename E>
auto operator*(OpLiteral<S1> const& a, OpMap<G, V, E> const& b)
{
	return symphas::internal::make_map::template get<G>(a.value * b.value, expr::compound_get::template expr(b));
}




template<typename G, typename A>
inline auto symphas::internal::make_map::get(A&& a)
{
	return get<G>(OpIdentity{}, std::forward<A>(a));
}

template<typename G, typename V, typename E>
inline auto symphas::internal::make_map::get(V v, OpExpression<E> const& e)
{
	return OpMap<G, V, E>(v, *static_cast<const E*>(&e));
}







namespace expr
{
	//! Initializes a mapping expression based on the Fourier transformation.
	/*!
	 * Initializes an OpMap instance for the currently implemented grid pairs that
	 * can be used, which is currently limited to the MapGridFourier type.
	 */
	template<typename V, typename E, size_t D = expr::grid_dim<E>::dimension>
	auto make_fourier_map(V value, OpExpression<E> const& e)
	{
		using src_type = typename expr::eval_type<E>::type;
		using grid_pair_type = MapGridFourier<src_type, complex_t, D>;
		return OpMap<grid_pair_type, V, E>(value, *static_cast<E const*>(&e));

	}

	//! Initializes a mapping expression based on the Fourier transformation.
	/*!
	 * Initializes an OpMap instance for the currently implemented grid pairs that
	 * can be used, which is currently limited to the MapGridFourier type.
	 */
	template<typename E, size_t D = expr::grid_dim<E>::dimension>
	auto make_fourier_map(OpExpression<E> const& e)
	{
		return expr::make_fourier_map(OpIdentity{}, e);
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



