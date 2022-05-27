
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
 * PURPOSE: Sets up a framework to obtain the indefinite integral of
 * an expression.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressionsprint.h"



namespace symphas::internal
{

	//! Implementation of convolution expression construction.
	/*!
	 * Implementation of functions which generalize and simplify the way to
	 * construct convolution expressions. Wraps the template deduction necessary
	 * to initialize a convolution expression.
	 */
	template<size_t Z>
	struct make_integral
	{
		template<typename A, typename B>
		static auto get(A&& a);

		template<typename V, typename E>
		static auto get(V v, OpExpression<E> const& e);

	};

}



namespace expr
{
	template<size_t Z, typename A>
	auto make_integral(A&& a)
	{
		return symphas::internal::make_integral<Z>::template get(std::forward<A>(a));
	}

	template<size_t Z, typename V, typename A>
	auto make_integral(V&& v, A&& a)
	{
		return symphas::internal::make_integral<Z>::template get(std::forward<V>(v), std::forward<A>(a));
	}
}



//! Represents the integration of an expression.
/*!
 * The variable of integration is chosen by the variable index specified
 * by the first template parameter.
 * 
 * \tparam Z The index of the variable of integration.
 * \tparam V The type of the coefficient.
 * \tparam E The type of the expression which is being integrated.
 */
template<size_t Z, typename V, typename E>
struct OpFuncIntegral : OpExpression<OpFuncIntegral<Z, V, E>>
{
	OpFuncIntegral(V value, E const& e) : e{ e }, value{ value } {}

	inline auto eval(iter_type n) const
	{
		return OpVoid{};
	}

	auto operator-() const
	{
		return make_integral::get(-value, e);
	}

	auto apply()
	{

	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::print_with_coeff(out, "", value);
	}

	size_t print(char* out) const
	{
		return expr::print_with_coeff(out, "", value);
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value);
	}

#endif

	V value;							// value multiplying the result of this derivative

protected:
	E e;								// expression object specifying grid values

};


template<typename S1, typename V2, typename E2>
auto operator*(OpLiteral<S1> const& a, OpFuncIntegral<V2, E2> const& b)
{
	return make_integral::template get(a.value * b.value, b.e);
}

template<typename S2, typename V1, typename E1>
auto operator*(OpFuncIntegral<V1, E1> const& a, OpLiteral<S2> const& b)
{
	return make_integral::template get(b.value * a.value, a.e);
}



namespace symphas::internal
{

	template<size_t Z>
	template<typename V, typename E>
	inline auto make_integral<Z>::get(V v, OpExpression<E> const& e)
	{
		return OpFuncIntegral<Z, V, E>(v, e);
	}

	template<size_t Z>
	template<typename A, typename B>
	inline auto make_integral<Z>::get(A&& a)
	{
		return get(OpIdentity{}, std::forward<A>(a));
	}

}



