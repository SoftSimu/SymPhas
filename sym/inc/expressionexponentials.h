
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
 * PURPOSE: Implements derivatives of the symbolic algebra library.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressionlib.h"

//! \cond

#ifdef LATEX_PLOT
#define SYEX_EXP_FMT_A "e^{"
#define SYEX_EXP_FMT_B "}"
#else
#define SYEX_EXP_FMT_A "e^{"
#define SYEX_EXP_FMT_B "}"
#endif

#define SYEX_EXP_FMT SYEX_EXP_FMT_A "%s" SYEX_EXP_FMT_B
#define SYEX_EXP_FMT_LEN (STR_ARR_LEN(SYEX_EXP_FMT_A SYEX_EXP_FMT_B) - 1)
//! \endcond

// *************************************************************************************

namespace symphas::internal
{
	//! Implementation of exponential expression construction.
	/*!
	 * Implementation of functions which generalize and simplify the way to
	 * construct exponential expressions. Wraps the template deduction necessary
	 * to initialize a exponential expression.
	 */
	struct make_exponential
	{
		//! Constructs the exponential with the identity coefficient.
		template<typename A>
		static auto get(A&&);

		template<typename V, typename E>
		static auto get(V const&, OpExpression<E> const&);

		template<typename V, typename S>
		static auto get(V const& v, OpLiteral<S> const& e);

		template<typename V>
		static auto get(V const& v, OpVoid const e);

		template<typename V>
		static auto get(V const& v, OpIdentity const e);

		template<typename V>
		static auto get(V const& v, OpNegIdentity const e);

		template<typename V, size_t N, size_t D>
		static auto get(V const& v, OpFractionLiteral<N, D> const e);

		template<typename V, size_t N, size_t D>
		static auto get(V const& v, OpNegFractionLiteral<N, D> const e);

		template<typename V, typename E>
		static auto get(OpLiteral<V> const& v, E&& e)
		{
			return get(v.value, std::forward<E>(e));
		}
	};
}


//! Exponential function of an expression.
/*!
 * The expression is natural number \f$e\f$ to the power of the given
 * expression.
 * 
 * \tparam V The type of the coefficient.
 * \tparam E The type of the expression which is exponentiated.
 */
template<typename V, typename E>
struct OpExponential : OpExpression<OpExponential<V, E>>
{
	OpExponential() : value{ V{} }, e{} {}

	//! Generate an exponential term with the given expression.
	/*!
	 * Generate an exponential term with the given expression.
	 * 
	 * \param value The coefficient of the exponential.
	 * \param e The expression that is exponentiated.
	 */
	OpExponential(V value, E const& e) : value{ value }, e{ e } {}

	inline auto eval(iter_type n) const
	{
		using namespace std;
		using namespace symphas::math;
		return expr::eval(value) * exp(e.eval(n));
	}

	auto operator-() const
	{
		return symphas::internal::make_exponential::get(-value, e);
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += fprintf(out, SYEX_EXP_FMT_A);
		n += e.print(out);
		n += fprintf(out, SYEX_EXP_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += sprintf(out + n, SYEX_EXP_FMT_A);
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_EXP_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + SYEX_EXP_FMT_LEN + e.print_length();
	}

#endif

	V value;						//!< Value multiplying the result of this convolution.
	E e;							//!< The expression which is exponentiated.


};


namespace symphas::internal
{

	template<typename A>
	inline auto make_exponential::get(A&& a)
	{
		return get(OpIdentity{}, std::forward<A>(a));
	}

	template<typename V, typename E>
	inline auto make_exponential::get(V const& v, OpExpression<E> const& a)
	{
		static_assert(!std::is_same<V, OpLiteral<double>>::value);
		return OpExponential<V, E>(v, *static_cast<const E*>(&a));
	}

	template<typename V, typename S>
	inline auto make_exponential::get(V const& v, OpLiteral<S> const& e)
	{
		using namespace std;
		using namespace symphas::math;
		return expr::make_literal(v * exp(e.value));
	}

	template<typename V>
	inline auto make_exponential::get(V const& v, OpVoid const)
	{
		return v;
	}

	template<typename V>
	inline auto make_exponential::get(V const& v, OpIdentity const e)
	{
		return get(v, expr::make_literal(e.eval()));
	}

	template<typename V>
	inline auto make_exponential::get(V const& v, OpNegIdentity const e)
	{
		return get(v, expr::make_literal(e.eval()));
	}

	template<typename V, size_t N, size_t D>
	inline auto make_exponential::get(V const& v, OpFractionLiteral<N, D> const e)
	{
		return get(v, expr::make_literal(e.eval()));
	}

	template<typename V, size_t N, size_t D>
	inline auto make_exponential::get(V const& v, OpNegFractionLiteral<N, D> const e)
	{
		return get(v, expr::make_literal(e.eval()));
	}
}


/* multiplication between exponential expressions
 */


template<typename V1, typename E1, typename V2, typename E2>
auto operator*(OpExponential<V1, E1> const& a, OpExponential<V2, E2> const& b)
{
	return symphas::internal::make_exponential::template get(a.value * b.value, a.e + b.e);
}

template<typename coeff_t, typename V, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V>), int> = 0>
auto operator*(coeff_t const& value, OpExponential<V, E> const& b)
{
	return symphas::internal::make_exponential::template get(value * b.value, expr::get_enclosed_expression(b));
}

template<typename coeff_t, typename tensor_t, typename E,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpExponential<tensor_t, E> const& b)
{
	return (value * b.value) * symphas::internal::make_exponential::template get(OpIdentity{}, expr::get_enclosed_expression(b));
}



namespace expr
{
	//! Create the expression for an exponential of an expression.
	/*!
	 * Create the expression for an exponential of an expression.
	 * 
	 * \param a The expression which is exponentiated.
	 */
	template<typename E>
	auto exp(OpExpression<E> const& a)
	{
		return symphas::internal::make_exponential::template get(*static_cast<E const*>(&a));
	}
}




