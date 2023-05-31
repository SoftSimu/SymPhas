
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
 * PURPOSE: Defines almost all the rules present in the symbolic algebra
 * functionality.
 * 
 * ***************************************************************************
 */

#pragma once

#include <cassert>

#include "expressionaggregates.h"

#include "symbolicdata.h"

// ******************************************************************************************


/*
 *
 * Multiplication between basic values, particularly identities.
 *
 ******************************************************************************/



//! Multiplication between integer and complex types.
inline auto operator*(int const a, complex_t const& b)
{
	return static_cast<scalar_t>(a) * b;
}

//! Multiplication between complex and integer types.
inline auto operator*(complex_t const& a, int const b)
{
	return a * static_cast<scalar_t>(b);
}

//! Multiplication between vector and complex types.
template<size_t D>
auto operator*(vector_t<D> const &a, complex_t const& b)
{
	any_vector_t<mul_result_t<scalar_t, complex_t>, D> result;
	for (iter_type i = 0; i < D; ++i)
	{
		result[i] = a[i] * b;
	}
	return result;
}

//! Multiplication between complex and vector types.
template<size_t D>
auto operator*(complex_t const& a, vector_t<D> const &b)
{
	any_vector_t<mul_result_t<scalar_t, complex_t>, D> result;
	for (iter_type i = 0; i < D; ++i)
	{
		result[i] = a * b[i];
	}
	return result;
}


//! Dividing with a vector in the denominator.
template<typename E, typename T, size_t... Ns>
auto operator/(E const& a, OpTensor<T, Ns...> const& b) = delete;


//
////! Multiplication between the identity expression and an integer value.
//template<typename E>
//inline auto operator*(OpIdentity, E const& e)
//{
//	return e;
//}
//
////! Multiplication between the identity expression and a vector.
//template<typename E>
//auto operator*(E const& e, OpIdentity)
//{
//	return e;
//}


//! Multiplication between the identity expression and a vector.
template<typename T, size_t D>
auto operator*(expr::symbols::Symbol, any_vector_t<T, D> const& b)
{
	return any_vector_t<expr::symbols::Symbol, D>{};
}

//! Multiplication between the identity expression and a vector.
template<typename T, size_t D>
auto operator*(any_vector_t<T, D> const& a, expr::symbols::Symbol)
{
	return any_vector_t<expr::symbols::Symbol, D>{};
}

//! Multiplication between the identity expression and a vector.
template<typename T, size_t D>
auto operator*(OpIdentity, any_vector_t<T, D> const& b)
{
	return b;
}

//! Multiplication between the identity expression and a scalar.
inline auto operator*(OpIdentity, scalar_t const& b)
{
	return b;
}

//! Multiplication between the identity expression and a complex value.
inline auto operator*(OpIdentity, complex_t const& b)
{
	return b;
}

//! Multiplication between the identity expression and an integer value.
inline auto operator*(OpIdentity, int b)
{
	return b;
}

//! Multiplication between the identity expression and a vector.
template<typename T, size_t D>
auto operator*(any_vector_t<T, D> const& a, OpIdentity)
{
	return a;
}

//! Multiplication between the identity expression and a scalar.
inline auto operator*(scalar_t const& a, OpIdentity)
{
	return a;
}

//! Multiplication between the identity expression and a complex value.
inline auto operator*(complex_t const& a, OpIdentity)
{
	return a;
}

//! Multiplication between the identity expression and an integer value.
inline auto operator*(int a, OpIdentity)
{
	return a;
}

//! Multiplication between the identity expression and a vector.
template<typename T, size_t D>
auto operator*(OpNegIdentity, any_vector_t<T, D> const& b)
{
	return -b;
}

//! Multiplication between the identity expression and a scalar.
inline auto operator*(OpNegIdentity, scalar_t const& b)
{
	return -b;
}

//! Multiplication between the identity expression and a complex value.
inline auto operator*(OpNegIdentity, complex_t const& b)
{
	return -b;
}

//! Multiplication between the identity expression and an integer value.
inline auto operator*(OpNegIdentity, int b)
{
	return -b;
}

//! Multiplication between the identity expression and a vector.
template<typename T, size_t D>
auto operator*(any_vector_t<T, D> const& a, OpNegIdentity)
{
	return -a;
}

//! Multiplication between the identity expression and a scalar.
inline auto operator*(scalar_t const& a, OpNegIdentity)
{
	return -a;
}

//! Multiplication between the identity expression and a complex value.
inline auto operator*(complex_t const& a, OpNegIdentity)
{
	return -a;
}

//! Multiplication between the identity expression and an integer value.
inline auto operator*(int a, OpNegIdentity)
{
	return -a;
}

//! Multiplication between the identity expression and a scalar.
template<typename T, size_t D, typename S>
auto operator*(any_vector_t<T, D> const& a, OpLiteral<S> const& b)
{
	return expr::make_literal(a * b.value);
}

//! Multiplication between the identity expression and a scalar.
template<typename T>
auto operator*(scalar_t a, OpLiteral<T> const& b)
{
	return expr::make_literal(a * b.value);
}

//! Multiplication between the identity expression and a complex value.
template<typename T>
auto operator*(complex_t a, OpLiteral<T> const& b)
{
	return expr::make_literal(a * b.value);
}

//! Multiplication between the identity expression and an integer value.
template<typename T>
auto operator*(int a, OpLiteral<T> const& b)
{
	return expr::make_literal(a * b.value);
}


//! Multiplication between the identity expression and a scalar.
template<typename T, size_t D, typename S>
auto operator*(OpLiteral<S> const& a, any_vector_t<T, D> const& b)
{
	return expr::make_literal(a.value * b);
}

//! Multiplication between the identity expression and a scalar.
template<typename T>
auto operator*(OpLiteral<T> const& a, scalar_t b)
{
	return expr::make_literal(a.value * b);
}

//! Multiplication between the identity expression and a complex value.
template<typename T>
auto operator*(OpLiteral<T> const& a, complex_t b)
{
	return expr::make_literal(a.value * b);
}

//! Multiplication between the identity expression and an integer value.
template<typename T>
auto operator*(OpLiteral<T> const& a, int b)
{
	return expr::make_literal(a.value * b);
}



//! Multiplication between the identity expression and anything.
template<typename T>
auto operator*(OpIdentity, OpLiteral<T> const& b)
{
	return b;
}

//! Multiplication between anything and the identity expression.
template<typename T>
auto operator*(OpLiteral<T> const& a, OpIdentity)
{
	return a;
}

//! Multiplication between the identity expression and anything.
template<typename T>
auto operator*(OpNegIdentity, OpLiteral<T> const& b)
{
	return -b;
}

//! Multiplication between anything and the identity expression.
template<typename T>
auto operator*(OpLiteral<T> const& a, OpNegIdentity)
{
	return -a;
}

template<typename T1, typename I, typename T2>
auto operator+(OpCoeff<T1, I> const& a, OpLiteral<T2> const& b)
{
	auto result = expr::init_coeff_from<add_result_t<T1, T2>, I>(a);
	for (iter_type i = 0; i < a.data.len; ++i)
	{
		result.data[i] = a.data[i] + b;
	}
	return result;
}

template<typename T1, typename I, typename T2>
auto operator-(OpCoeff<T1, I> const& a, OpLiteral<T2> const& b)
{
	auto result = expr::init_coeff_from<sub_result_t<T1, T2>, I>(a);
	for (iter_type i = 0; i < a.data.len; ++i)
	{
		result.data[i] = a.data[i] - b;
	}
	return result;
}

template<typename T1, typename I, typename T2>
auto operator*(OpCoeff<T1, I> const& a, OpLiteral<T2> const& b)
{
	auto result = expr::init_coeff_from<mul_result_t<T1, T2>, I>(a);
	for (iter_type i = 0; i < a.data.len; ++i)
	{
		result.data[i] = a.data[i] * b;
	}
	return result;
}

template<typename T1, typename I, typename T2>
auto operator/(OpCoeff<T1, I> const& a, OpLiteral<T2> const& b)
{
	auto result = expr::init_coeff_from<div_result_t<T1, T2>, I>(a);
	for (iter_type i = 0; i < a.data.len; ++i)
	{
		result.data[i] = a.data[i] / b;
	}
	return result;
}

template<typename T, typename I, typename coeff_t, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator+(OpCoeff<T, I> const& a, coeff_t const& b)
{
	auto result = expr::init_coeff_from<T, I>(a);
	for (iter_type i = 0; i < a.data.len; ++i)
	{
		result.data[i] = a.data[i] + b;
	}
	return result;
}

template<typename T, typename I, typename coeff_t, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator-(OpCoeff<T, I> const& a, coeff_t const& b)
{
	auto result = expr::init_coeff_from<T, I>(a);
	for (iter_type i = 0; i < a.data.len; ++i)
	{
		result.data[i] = a.data[i] - b;
	}
	return result;
}

template<typename T, typename I, typename coeff_t, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(OpCoeff<T, I> const& a, coeff_t const& b)
{
	auto result = expr::init_coeff_from<T, I>(a);
	for (iter_type i = 0; i < a.data.len; ++i)
	{
		result.data[i] = a.data[i] * b;
	}
	return result;
}

template<typename T, typename I, typename coeff_t, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator/(OpCoeff<T, I> const& a, coeff_t const& b)
{
	auto result = expr::init_coeff_from<T, I>(a);
	for (iter_type i = 0; i < a.data.len; ++i)
	{
		result.data[i] = a.data[i] / b;
	}
	return result;
}

template<typename T, typename I, typename coeff_t, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator+(coeff_t const& a, OpCoeff<T, I> const& b)
{
	auto result = expr::init_coeff_from<T, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = a + b.data[i];
	}
	return result;
}

template<typename T, typename I, typename coeff_t, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator-(coeff_t const& a, OpCoeff<T, I> const& b)
{
	auto result = expr::init_coeff_from<T, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = a - b.data[i];
	}
	return result;
}

template<typename T, typename I, typename coeff_t, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator*(coeff_t const& a, OpCoeff<T, I> const& b)
{
	auto result = expr::init_coeff_from<T, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = a * b.data[i];
	}
	return result;
}

template<typename T, typename I, typename coeff_t, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
auto operator/(coeff_t const& a, OpCoeff<T, I> const& b)
{
	auto result = expr::init_coeff_from<T, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = a / b.data[i];
	}
	return result;
}

template<typename T, typename I, typename V, size_t... Ns>
auto operator*(OpTensor<V, Ns...> const& a, OpCoeff<T, I> const& b)
{
	auto result = expr::init_coeff_from<mul_result_t<V, T>, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = V(a) * b.data[i];
	}
	return expr::make_tensor<Ns...>(result);
}


template<typename T, typename I, typename V>
auto operator/(OpLiteral<V> const& a, OpCoeff<T, I> const& b)
{
	auto result = expr::init_coeff_from<div_result_t<V, T>, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = a / b.data[i];
	}
	return result;
}

template<typename T, typename I>
auto operator/(OpIdentity, OpCoeff<T, I> const& b)
{
	auto result = expr::init_coeff_from<T, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = 1. / b.data[i];
	}
	return result;
}

template<typename T, typename I>
auto operator/(OpNegIdentity, OpCoeff<T, I> const& b)
{
	auto result = expr::init_coeff_from<T, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = -1. / b.data[i];
	}
	return result;
}

template<size_t N, size_t D, typename T, typename I>
auto operator/(OpFractionLiteral<N, D> const& a, OpCoeff<T, I> const& b)
{
	auto result = expr::init_coeff_from<T, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = a / b.data[i];
	}
	return result;
}

template<size_t N, size_t D, typename T, typename I>
auto operator/(OpNegFractionLiteral<N, D> const& a, OpCoeff<T, I> const& b)
{
	auto result = expr::init_coeff_from<T, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = a / b.data[i];
	}
	return result;
}

template<typename T, typename I>
auto operator/(OpVoid, OpCoeff<T, I> const& b)
{
	return OpVoid{};
}

template<typename T, typename I>
auto operator/(OpCoeff<T, I> const& b, OpVoid) = delete;




template<typename T1, typename I, typename T2>
auto operator+(OpCoeff<T1, I> const& a, OpCoeff<T2, I> const& b)
{
	auto result = expr::init_coeff_from<add_result_t<T1, T2>, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = a.data[i] + b.data[i];
	}
	return result;
}

template<typename T1, typename I, typename T2>
auto operator-(OpCoeff<T1, I> const& a, OpCoeff<T2, I> const& b)
{
	auto result = expr::init_coeff_from<sub_result_t<T1, T2>, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = a.data[i] - b.data[i];
	}
	return result;
}

template<typename T1, typename I, typename T2>
auto operator*(OpCoeff<T1, I> const& a, OpCoeff<T2, I> const& b)
{
	auto result = expr::init_coeff_from<mul_result_t<T1, T2>, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = a.data[i] * b.data[i];
	}
	return result;
}

template<typename T1, typename I, typename T2>
auto operator/(OpCoeff<T1, I> const& a, OpCoeff<T2, I> const& b)
{
	auto result = expr::init_coeff_from<div_result_t<T1, T2>, I>(b);
	for (iter_type i = 0; i < b.data.len; ++i)
	{
		result.data[i] = a.data[i] / b.data[i];
	}
	return result;
}



//! Multiplication between the identity expression and a scalar.
template<typename T, typename I, typename coeff_t, std::enable_if_t<!expr::is_expression<coeff_t>, int> = 0>
auto operator+(OpCoeff<T, I> const& a, coeff_t const& b)
{
	return a + expr::make_literal(b);
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I, typename coeff_t, std::enable_if_t<!expr::is_expression<coeff_t>, int> = 0>
auto operator-(OpCoeff<T, I> const& a, coeff_t const& b)
{
	return a - expr::make_literal(b);
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I, typename coeff_t, std::enable_if_t<!expr::is_expression<coeff_t>, int> = 0>
auto operator*(OpCoeff<T, I> const& a, coeff_t const& b)
{
	return a * expr::make_literal(b);
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I, typename coeff_t, std::enable_if_t<!expr::is_expression<coeff_t>, int> = 0>
auto operator/(OpCoeff<T, I> const& a, coeff_t const& b)
{
	return a / expr::make_literal(b);
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I, typename coeff_t, std::enable_if_t<!expr::is_expression<coeff_t>, int> = 0>
auto operator+(coeff_t const& a, OpCoeff<T, I> const& b)
{
	return expr::make_literal(a) + b;
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I, typename coeff_t, std::enable_if_t<!expr::is_expression<coeff_t>, int> = 0>
auto operator-(coeff_t const& a, OpCoeff<T, I> const& b)
{
	return expr::make_literal(a) - b;
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I, typename coeff_t, std::enable_if_t<!expr::is_expression<coeff_t>, int> = 0>
auto operator*(coeff_t const& a, OpCoeff<T, I> const& b)
{
	return expr::make_literal(a) * b;
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I, typename coeff_t, std::enable_if_t<!expr::is_expression<coeff_t>, int> = 0>
auto operator/(coeff_t const& a, OpCoeff<T, I> const& b)
{
	return expr::make_literal(a) / b;
}



//! Multiplication between the identity expression and a scalar.
template<typename T, typename I>
auto operator+(OpCoeff<T, I> const& a, expr::symbols::Symbol)
{
	return expr::symbols::Symbol{};
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I>
auto operator-(OpCoeff<T, I> const& a, expr::symbols::Symbol)
{
	return expr::symbols::Symbol{};
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I>
auto operator*(OpCoeff<T, I> const& a, expr::symbols::Symbol)
{
	return expr::symbols::Symbol{};
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I>
auto operator/(OpCoeff<T, I> const& a, expr::symbols::Symbol)
{
	return expr::symbols::Symbol{};
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I>
auto operator+(expr::symbols::Symbol, OpCoeff<T, I> const& b)
{
	return expr::symbols::Symbol{};
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I>
auto operator-(expr::symbols::Symbol, OpCoeff<T, I> const& b)
{
	return expr::symbols::Symbol{};
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I>
auto operator*(expr::symbols::Symbol, OpCoeff<T, I> const& b)
{
	return expr::symbols::Symbol{};
}

//! Multiplication between the identity expression and a scalar.
template<typename T, typename I>
auto operator/(expr::symbols::Symbol, OpCoeff<T, I> const& b)
{
	return expr::symbols::Symbol{};
}



//! Multiplication between the identity expression and anything.
template<typename E>
auto operator*(OpIdentity, OpExpression<E> const& e)
{
	return (*static_cast<E const*>(&e));
}

//! Multiplication between anything and the identity expression.
template<typename E>
auto operator*(OpExpression<E> const& e, OpIdentity)
{
	return (*static_cast<E const*>(&e));
}

//! Multiplication between the negative identity expression and anything.
template<typename E>
auto operator*(OpNegIdentity, OpExpression<E> const& e)
{
	return -(*static_cast<E const*>(&e));
}

//! Multiplication between the negative identity expression and anything.
template<typename E>
auto operator*(OpExpression<E> const& e, OpNegIdentity)
{
	return -(*static_cast<E const*>(&e));
}

//! Multiplication between two identity expressions.
inline auto operator*(OpIdentity, OpIdentity)
{
	return OpIdentity{};
}

//! Multiplication between two negative identity expressions.
inline auto operator*(OpNegIdentity, OpNegIdentity)
{
	return OpIdentity{};
}

//! Multiplication between the negative identity and identity expressions.
inline auto operator*(OpNegIdentity, OpIdentity)
{
	return OpNegIdentity{};
}

//! Multiplication between the identity and negative identity expressions.
inline auto operator*(OpIdentity, OpNegIdentity)
{
	return OpNegIdentity{};
}

//! Multiplication between anything and the 0 identity.
template<typename E>
auto operator*(E&&, OpVoid)
{
	return OpVoid{};
}

//! Multiplication between the 0 identity and anything.
template<typename E>
auto operator*(OpVoid, E&&)
{
	return OpVoid{};
}

//! Multiplication between the identity expression and 0 identity.
inline auto operator*(OpIdentity, OpVoid)
{
	return OpVoid{};
}

//! Multiplication between the 0 identity and identity expression.
inline auto operator*(OpVoid, OpIdentity)
{
	return OpVoid{};
}

//! Multiplication between the negative identity and 0 identity.
inline auto operator*(OpNegIdentity, OpVoid)
{
	return OpVoid{};
}

//! Multiplication between the 0 identity and negative identity.
inline auto operator*(OpVoid, OpNegIdentity)
{
	return OpVoid{};
}

//! Multiplication between two 0 identities.
inline auto operator*(OpVoid, OpVoid)
{
	return OpVoid{};
}



inline int& operator*=(int& lhs, OpIdentity)
{
	return lhs;
}

inline int& operator*=(int& lhs, OpNegIdentity)
{
	lhs = -lhs;
	return lhs;
}

inline scalar_t& operator*=(scalar_t& lhs, OpIdentity) 
{
	return lhs;
}

inline scalar_t& operator*=(scalar_t& lhs, OpNegIdentity)
{
	lhs = -lhs;
	return lhs;
}

inline complex_t& operator*=(complex_t& lhs, OpIdentity)
{
	return lhs;
}

inline complex_t& operator*=(complex_t& lhs, OpNegIdentity)
{
	lhs = -lhs;
	return lhs;
}

template<typename T, size_t D>
any_vector_t<T, D>& operator*=(any_vector_t<T, D>& lhs, OpIdentity)
{
	return lhs;
}

template<typename T, size_t D>
any_vector_t<T, D>& operator*=(any_vector_t<T, D>& lhs, OpNegIdentity)
{
	lhs = -lhs;
	return lhs;
}

template<typename T>
int& operator*=(int& lhs, OpLiteral<T> const& v)
{
	lhs = lhs * v;
	return lhs;
}

template<typename T>
scalar_t& operator*=(scalar_t& lhs, OpLiteral<T> const& v)
{
	lhs = lhs * v;
	return lhs;
}

template<typename T>
complex_t& operator*=(complex_t& lhs, OpLiteral<T> const& v)
{
	lhs = lhs * v;
	return lhs;
}

template<typename T1, size_t D, typename T2>
any_vector_t<T1, D>& operator*=(any_vector_t<T1, D>& lhs, OpLiteral<T2> const& v)
{
	lhs = lhs * v;
	return lhs;
}




// ******************************************************************************************



namespace expr
{
	template<typename T, size_t... Ns>
	auto inverse(OpTensor<T, Ns...> const& tensor) = delete;


	//! Apply an inverse to a scalar value.
	inline auto inverse(scalar_t e)
	{
		return symphas::lib::get_identity<scalar_t>() / e;
	}

	//! Apply an inverse to an integer value.
	inline auto inverse(int e)
	{
		return symphas::lib::get_identity<scalar_t>() / e;
	}

	//! Apply an inverse to a complex value.
	inline auto inverse(complex_t const& e)
	{
		return e / std::norm(e);
	}

	//! Apply an inverse to an expression.
	template<typename E>
	auto inverse(OpExpression<E> const& e)
	{
		return expr::make_div(OpIdentity{}, * static_cast<const E*>(&e));
	}

	//! Apply an inverse to an expression.
	template<typename A, typename B>
	auto inverse(OpBinaryDiv<A, B> const& e)
	{
		return e.b / e.a;
	}

	//! Apply an inverse to an expression literal.
	template<typename T>
	auto inverse(OpLiteral<T> const& e)
	{
		return expr::make_literal(inverse(e.value));
	}

	//! Apply an inverse to an expression literal.
	template<size_t N, size_t D>
	constexpr auto inverse(OpFractionLiteral<N, D>)
	{
		return OpFractionLiteral<D, N>{};
	}

	//! Apply an inverse to an expression literal.
	template<size_t N, size_t D>
	constexpr auto inverse(OpNegFractionLiteral<N, D>)
	{
		return OpNegFractionLiteral<D, N>{};
	}

	//! Apply an inverse to an expression literal.
	constexpr inline auto inverse(OpIdentity)
	{
		return OpIdentity{};
	}

	//! Apply an inverse to an expression literal.
	constexpr inline auto inverse(OpNegIdentity)
	{
		return OpNegIdentity{};
	}

	//! Apply an inverse to an expression.
	template<typename V, typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
	auto inverse(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& e)
	{
		return expr::make_div(OpIdentity{}, OpTerms(OpIdentity{}, e));
	}

	//! Apply an inverse to an expression.
	template<typename V, typename... Gs, expr::exp_key_t... Xs>
	auto inverse(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		return expr::make_div(expr::inverse(expr::coeff(e)), OpTerms(OpIdentity{}, expr::terms_after_first(e)));
	}

	//! Apply an inverse to an expression.
	template<typename V, typename E>
	auto inverse(OpExponential<V, E> const& e)
	{
		return exp(e.value, -expr::get_enclosed_expression(e));
	}

	//! Apply an inverse to an expression.
	template<expr::exp_key_t X, typename V, typename E>
	auto inverse(OpPow<X, V, E> const& e)
	{
		constexpr expr::exp_key_t Xx = Xk<_Xk_t<X>::N, _Xk_t<X>::D, !_Xk_t<X>::sign>;
		return make_pow<Xx>(e.value, expr::get_enclosed_expression(e));
	}

	//template<typename T, size_t... Ns>
	//auto inverse(OpTensor<T, Ns...> const& tensor)
	//{
	//	return symphas::internal::tensor_cancel{};
	//}
}


/*
 *
 * Division between basic values, particularly identities.
 *
 ******************************************************************************/

 //! Division between anything and the identity expression.
template<typename E>
decltype(auto) operator/(E&& a, OpIdentity)
{
	return std::forward<E>(a);
}

//! Division between the identity expression and anything.
template<typename E>
decltype(auto) operator/(OpIdentity, E&& b)
{
	return expr::inverse(std::forward<E>(b));
}

//! Division between the identity expression and anything.
inline auto operator/(OpIdentity, OpVoid) = delete;

//! Division between anything and the identity expression.
template<typename E>
decltype(auto) operator/(E&& a, OpNegIdentity)
{
	return -std::forward<E>(a);
}

//! Division between the negative identity expression and anything.
template<typename E>
auto operator/(OpNegIdentity, E&& b)
{
	return -expr::inverse(std::forward<E>(b));
}

//! Division between two identities.
constexpr inline auto operator/(OpIdentity, OpIdentity)
{
	return OpIdentity{};
}

//! Division between two negative identities.
constexpr inline auto operator/(OpNegIdentity, OpNegIdentity)
{
	return OpIdentity{};
}

//! Division between the identity and the negative identity.
constexpr inline auto operator/(OpIdentity, OpNegIdentity)
{
	return OpNegIdentity{};
}

//! Division between the negative identity and the identity.
constexpr inline auto operator/(OpNegIdentity, OpIdentity)
{
	return OpNegIdentity{};
}

//! Division between anything and the 0 identity.
template<typename E>
auto operator/(E&&, OpVoid const) = delete;

//! Division between 0 and 0.
inline auto operator/(OpVoid const, OpVoid const) = delete;

//! Division between 0 identity and anything.
template<typename E>
constexpr auto operator/(OpVoid const, E const&)
{
	return OpVoid{};
}

//! Division between 0 identity and identity expression.
constexpr inline auto operator/(OpVoid const, OpIdentity)
{
	return OpVoid{};
}

//! Division between 0 identity and negative identity expression.
constexpr inline auto operator/(OpVoid const, OpNegIdentity)
{
	return OpVoid{};
}

//! Division between a value constant and anything.
template<typename T, typename E>
auto operator/(OpLiteral<T> const& a, E const& b)
{
	return expr::inverse(expr::inverse(a) * b);
}

//! Division between two value constants.
template<typename T>
auto operator/(OpLiteral<T> const& a, OpLiteral<T> const& b)
{
	return a * expr::inverse(b);
}

//! Division between two value constants.
template<typename T, typename S>
auto operator/(OpLiteral<T> const& a, OpLiteral<S> const& b)
{
	return a * expr::inverse(b);
}

//! Division between anything and a value constant.
template<typename T, typename E>
auto operator/(E const& a, OpLiteral<T> const& b)
{
	return expr::inverse(b) * a;
}

//! Division between anything and a tensor is deleted.
template<typename E, typename T, size_t... Ns>
auto operator/(E&&, OpTensor<T, Ns...> const& tensor) = delete;





/*
 *
 * Addition between basic values, particularly identities.
 *
 ******************************************************************************/


 //! Addition between integer and complex types.
inline auto operator+(int const a, complex_t const& b)
{
	return static_cast<scalar_t>(a) + b;
}

//! Addition between complex and integer types.
inline auto operator+(complex_t const& a, int const b)
{
	return a + static_cast<scalar_t>(b);
}



//! Addition between anything and the 0 identity.
template<typename E>
decltype(auto) operator+(E&& a, OpVoid)
{
	return std::forward<E>(a);
}

//! Addition between the 0 identity and anything.
template<typename E>
decltype(auto) operator+(OpVoid, E&& b)
{
	return std::forward<E>(b);
}

//! Addition between two 0 identities.
inline auto operator+(OpVoid, OpVoid)
{
	return OpVoid{};
}

//! Addition between two identities.
inline auto operator+(OpIdentity, OpIdentity)
{
	return expr::make_fraction<2, 1>();
}

//! Addition between two negative identities.
inline auto operator+(OpNegIdentity, OpNegIdentity)
{
	return -expr::make_fraction<2, 1>();
}

//! Addition between the identity expression and negative identity.
inline OpVoid operator+(OpIdentity, OpNegIdentity)
{
	return OpVoid{};
}

//! Addition between the negative identity expression and identity.
inline OpVoid operator+(OpNegIdentity, OpIdentity)
{
	return OpVoid{};
}

inline auto operator+(OpVoid, OpNegIdentity)
{
	return OpNegIdentity{};
}

inline auto operator+(OpVoid, OpIdentity)
{
	return OpIdentity{};
}




inline auto operator+(OpVoid, expr::symbols::Symbol)
{
	return expr::symbols::Symbol{};
}

inline auto operator+(expr::symbols::Symbol, OpVoid)
{
	return expr::symbols::Symbol{};
}

inline auto operator-(OpVoid, expr::symbols::Symbol)
{
	return expr::symbols::Symbol{};
}

inline auto operator-(expr::symbols::Symbol, OpVoid)
{
	return expr::symbols::Symbol{};
}

inline auto operator*(OpVoid, expr::symbols::Symbol)
{
	return OpVoid{};
}

inline auto operator*(expr::symbols::Symbol, OpVoid)
{
	return OpVoid{};
}

inline auto operator/(OpVoid, expr::symbols::Symbol)
{
	return OpVoid{};
}

inline auto operator/(expr::symbols::Symbol, OpVoid) = delete;





//! Addition between a primitive value and the identity expression.
inline auto operator+(scalar_t const a, OpIdentity)
{
	return a + OpIdentity{}.eval();
}

//! Addition between the identity expression and a primitive value.
inline auto operator+(OpIdentity, scalar_t const b)
{
	return OpIdentity{}.eval() + b;
}

//! Addition between a primitive value and the identity expression.
inline auto operator+(int const a, OpIdentity)
{
	return a + OpIdentity{}.eval();
}

//! Addition between the identity expression and a primitive value.
inline auto operator+(OpIdentity, int const b)
{
	return OpIdentity{}.eval() + b;
}

//! Addition between a primitive complex type and the identity expression.
inline auto operator+(complex_t const& a, OpIdentity)
{
	return a + static_cast<complex_t>(OpIdentity{}.eval());
}

//! Addition between the identity expression and a primitive complex type.
inline auto operator+(OpIdentity, complex_t const b)
{
	return static_cast<complex_t>(OpIdentity{}.eval()) + b;
}




//! Addition between a primitive value and the negative identity expression.
inline auto operator+(scalar_t const a, OpNegIdentity)
{
	return a + OpNegIdentity{}.eval();
}

//! Addition between the negative identity expression and a primitive value.
inline auto operator+(OpNegIdentity, scalar_t const b)
{
	return OpNegIdentity{}.eval() + b;
}

//! Addition between a primitive value and the negative identity expression.
inline auto operator+(int const a, OpNegIdentity)
{
	return a + OpNegIdentity{}.eval();
}

//! Addition between the negative identity expression and a primitive value.
inline auto operator+(OpNegIdentity, int const b)
{
	return OpNegIdentity{}.eval() + b;
}

//! Addition between a primitive complex type and the negative identity.
inline auto operator+(complex_t const& a, OpNegIdentity)
{
	return a + static_cast<complex_t>(OpNegIdentity{}.eval());
}

//! Addition between the negative identity and a primitive complex type.
inline auto operator+(OpNegIdentity, complex_t const b)
{
	return static_cast<complex_t>(OpNegIdentity{}.eval()) + b;
}



/*
 *
 * Subtraction between basic values, particularly identities.
 *
 ******************************************************************************/

 //! Subtraction between anything and the 0 identity.
inline auto operator-(int const a, complex_t const& b)
{
	return static_cast<scalar_t>(a) - b;
}

//! Subtraction between complex and integer types.
inline auto operator-(complex_t const& a, int const b)
{
	return a - static_cast<scalar_t>(b);
}




//! Subtraction between anything and the 0 identity.
template<typename E>
decltype(auto) operator-(E&& a, OpVoid const)
{
	return std::forward<E>(a);
}

//! Subtraction between the 0 identity and anything.
template<typename E>
decltype(auto) operator-(OpVoid const, E&& b)
{
	return -std::forward<E>(b);
}

//! Subtraction between two 0 identities.
inline auto operator-(OpVoid const, OpVoid const)
{
	return OpVoid{};
}

//! Subtraction between two identities.
inline auto operator-(OpIdentity, OpIdentity)
{
	return OpVoid{};
}

//! Subtraction between two negative identities.
inline auto operator-(OpNegIdentity, OpNegIdentity)
{
	return OpVoid{};
}

//! Subtraction between two identities.
inline auto operator-(OpIdentity, OpNegIdentity)
{
	return OpIdentity{} + OpIdentity{};
}

//! Subtraction between two negative identities.
inline auto operator-(OpNegIdentity, OpIdentity)
{
	return OpNegIdentity{} + OpNegIdentity{};
}

inline auto operator-(OpVoid, OpNegIdentity)
{
	return OpNegIdentity{};
}

inline auto operator-(OpVoid, OpIdentity)
{
	return OpIdentity{};
}

//! Subtraction between a primitive value and the identity expression.
inline auto operator-(scalar_t const a, OpIdentity)
{
	return a - OpIdentity{}.eval();
}

//! Subtraction between the identity expression and a primitive value.
inline auto operator-(OpIdentity, scalar_t const b)
{
	return OpIdentity{}.eval() - b;
}

//! Subtraction between a primitive value and the identity expression.
inline auto operator-(int const a, OpIdentity)
{
	return a - OpIdentity{}.eval();
}

//! Subtraction between the identity expression and a primitive value.
inline auto operator-(OpIdentity, int const b)
{
	return OpIdentity{}.eval() - b;
}

//! Subtraction between a primitive complex type and the identity expression.
inline auto operator-(complex_t const a, OpIdentity)
{
	return a - static_cast<complex_t>(OpIdentity{}.eval());
}

//! Subtraction between the identity expression and a primitive complex type.
inline auto operator-(OpIdentity, complex_t const b)
{
	return static_cast<complex_t>(OpIdentity{}.eval()) - b;
}





//! Subtraction between a primitive value and the negative identity expression.
inline auto operator-(scalar_t const a, OpNegIdentity)
{
	return a - OpNegIdentity{}.eval();
}

//! Subtraction between the negative identity expression and a primitive value.
inline auto operator-(OpNegIdentity, scalar_t const b)
{
	return OpNegIdentity{}.eval() - b;
}

//! Subtraction between a primitive value and the negative identity expression.
inline auto operator-(int const a, OpNegIdentity)
{
	return a - OpNegIdentity{}.eval();
}

//! Subtraction between the negative identity expression and a primitive value.
inline auto operator-(OpNegIdentity, int const b)
{
	return OpNegIdentity{}.eval() - b;
}

//! Subtraction between a primitive complex type and the negative identity.
inline auto operator-(complex_t const& a, OpNegIdentity)
{
	return a - static_cast<complex_t>(OpNegIdentity{}.eval());
}

//! Subtraction between the negative identity and a primitive complex type.
inline auto operator-(OpNegIdentity, complex_t const& b)
{
	return static_cast<complex_t>(OpNegIdentity{}.eval()) - b;
}

//template<typename E>
//inline auto operator-(E&& e, OpNegIdentity)
//{
//	return std::forward<E>(e) + OpIdentity{};
//}
//
//template<typename E>
//inline auto operator+(E&& e, OpNegIdentity)
//{
//	return std::forward<E>(e) - OpIdentity{};
//}
//
//template<typename E, size_t... Ns>
//inline auto operator-(E&& e, OpTensor<OpNegIdentity, Ns...>)
//{
//	return std::forward<E>(e) + OpTensor<OpIdentity, Ns...>{};
//}
//
//template<typename E, size_t... Ns>
//inline auto operator+(E&& e, OpTensor<OpNegIdentity, Ns...>)
//{
//	return std::forward<E>(e) - OpTensor<OpIdentity, Ns...>{};
//}
//
//template<size_t N0, size_t N1, size_t N>
//inline auto operator+(OpTensor<OpNegIdentity, N0, N>, OpTensor<OpNegIdentity, N1, N>)
//{
//	return expr::make_add(OpTensor<OpNegIdentity, N0, N>{}, OpTensor<OpNegIdentity, N1, N>{});
//}
//
//template<size_t N0, size_t N1, size_t N>
//inline auto operator-(OpTensor<OpNegIdentity, N0, N>, OpTensor<OpNegIdentity, N1, N>)
//{
//	return expr::make_add(OpTensor<OpNegIdentity, N0, N>{}, OpTensor<OpIdentity, N1, N>{});
//}



//! Subtraction between a primitive value and the identity expression.
template<typename T>
auto operator-(OpLiteral<T> const a, OpIdentity)
{
	return expr::make_literal(a.value - OpIdentity{}.eval());
}

//! Subtraction between the identity expression and a primitive value.
template<typename T>
auto operator-(OpIdentity, OpLiteral<T> const b)
{
	return expr::make_literal(OpIdentity{}.eval() - b.value);
}

//! Subtraction between a primitive value and the identity expression.
template<typename T>
auto operator-(OpLiteral<T> const a, OpNegIdentity)
{
	return expr::make_literal(a.value - OpNegIdentity{}.eval());
}

//! Subtraction between the identity expression and a primitive value.
template<typename T>
auto operator-(OpNegIdentity, OpLiteral<T> const b)
{
	return expr::make_literal(OpNegIdentity{}.eval() - b.value);
}
	
//! Addition between a primitive value and the identity expression.
template<typename T>
auto operator+(OpLiteral<T> const a, OpIdentity)
{
	return expr::make_literal(a.value + OpIdentity{}.eval());
}

//! Addition between the identity expression and a primitive value.
template<typename T>
auto operator+(OpIdentity, OpLiteral<T> const b)
{
	return expr::make_literal(OpIdentity{}.eval() + b.value);
}

//! Addition between a primitive value and the identity expression.
template<typename T>
auto operator+(OpLiteral<T> const a, OpNegIdentity)
{
	return expr::make_literal(a.value + OpNegIdentity{}.eval());
}

//! Addition between the identity expression and a primitive value.
template<typename T>
auto operator+(OpNegIdentity, OpLiteral<T> const b)
{
	return expr::make_literal(OpNegIdentity{}.eval() + b.value);
}

/*
 *
 * Operations between fraction terms.
 *
 ******************************************************************************/




template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator-(OpFractionLiteral<N1, D1>, OpFractionLiteral<N2, D2>)
{
	constexpr size_t Na = N1 * D2;
	constexpr size_t Nb = N2 * D1;
	constexpr size_t D = D1 * D2;

	if constexpr (Na > Nb)
	{
		return expr::make_fraction<Na - Nb, D>();
	}
	else
	{
		return -expr::make_fraction<Nb - Na, D>();
	}

}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator+(OpFractionLiteral<N1, D1>, OpFractionLiteral<N2, D2>)
{
	constexpr size_t N = N1 * D2 + N2 * D1;
	constexpr size_t D = D1 * D2;

	return expr::make_fraction<N, D>();
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator-(OpFractionLiteral<N1, D1>, OpNegFractionLiteral<N2, D2>)
{
	return OpFractionLiteral<N1, D1>{} + OpFractionLiteral<N2, D2>{};
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator+(OpFractionLiteral<N1, D1>, OpNegFractionLiteral<N2, D2>)
{
	return OpFractionLiteral<N1, D1>{} - OpFractionLiteral<N2, D2>{};
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator-(OpNegFractionLiteral<N1, D1>, OpFractionLiteral<N2, D2>)
{
	return -(OpFractionLiteral<N1, D1>{} + OpFractionLiteral<N2, D2>{});
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator+(OpNegFractionLiteral<N1, D1>, OpFractionLiteral<N2, D2>)
{
	return OpFractionLiteral<N2, D2>{} - OpFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator-(OpNegFractionLiteral<N1, D1>, OpNegFractionLiteral<N2, D2>)
{
	return OpFractionLiteral<N2, D2>{} - OpFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator+(OpNegFractionLiteral<N1, D1>, OpNegFractionLiteral<N2, D2>)
{
	return -(OpFractionLiteral<N1, D1>{} + OpFractionLiteral<N2, D2>{});
}


// Add and subtract with identities

template<size_t N1, size_t D1>
constexpr auto operator-(OpNegFractionLiteral<N1, D1>, OpIdentity)
{
	return OpNegFractionLiteral<N1, D1>{} - OpFractionLiteral<1, 1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator+(OpNegFractionLiteral<N1, D1>, OpIdentity)
{
	return OpNegFractionLiteral<N1, D1>{} + OpFractionLiteral<1, 1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator-(OpFractionLiteral<N1, D1>, OpIdentity)
{
	return OpFractionLiteral<N1, D1>{} - OpFractionLiteral<1, 1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator+(OpFractionLiteral<N1, D1>, OpIdentity)
{
	return OpFractionLiteral<N1, D1>{} + OpFractionLiteral<1, 1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator-(OpNegFractionLiteral<N1, D1>, OpNegIdentity)
{
	return OpNegFractionLiteral<N1, D1>{} + OpFractionLiteral<1, 1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator+(OpNegFractionLiteral<N1, D1>, OpNegIdentity)
{
	return OpNegFractionLiteral<N1, D1>{} - OpFractionLiteral<1, 1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator-(OpFractionLiteral<N1, D1>, OpNegIdentity)
{
	return OpFractionLiteral<N1, D1>{} + OpFractionLiteral<1, 1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator+(OpFractionLiteral<N1, D1>, OpNegIdentity)
{
	return OpFractionLiteral<N1, D1>{} - OpFractionLiteral<1, 1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator-(OpIdentity, OpNegFractionLiteral<N1, D1>)
{
	return OpFractionLiteral<1, 1>{} + OpFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator+(OpIdentity, OpNegFractionLiteral<N1, D1>)
{
	return OpFractionLiteral<1, 1>{} + OpNegFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator-(OpIdentity, OpFractionLiteral<N1, D1>)
{
	return OpFractionLiteral<1, 1>{} + OpNegFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator+(OpIdentity, OpFractionLiteral<N1, D1>)
{
	return OpFractionLiteral<1, 1>{} + OpFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator-(OpNegIdentity, OpNegFractionLiteral<N1, D1>)
{
	return OpNegFractionLiteral<1, 1>{} + OpFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator+(OpNegIdentity, OpFractionLiteral<N1, D1>)
{
	return OpNegFractionLiteral<1, 1>{} + OpFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator+(OpNegIdentity, OpNegFractionLiteral<N1, D1>)
{
	return OpNegFractionLiteral<1, 1>{} - OpFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator-(OpNegIdentity, OpFractionLiteral<N1, D1>)
{
	return OpNegFractionLiteral<1, 1>{} - OpFractionLiteral<N1, D1>{};
}



template<size_t N1, size_t D1>
constexpr auto operator*(OpIdentity, OpNegFractionLiteral<N1, D1>)
{
	return OpNegFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator*(OpIdentity, OpFractionLiteral<N1, D1>)
{
	return OpFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator*(OpNegIdentity, OpNegFractionLiteral<N1, D1>)
{
	return OpFractionLiteral<N1, D1>{};
}

template<size_t N1, size_t D1>
constexpr auto operator*(OpNegIdentity, OpFractionLiteral<N1, D1>)
{
	return OpNegFractionLiteral<N1, D1>{};
}



// add and subtract with normal scalar values


template<size_t N, size_t D, typename T, 
	typename = std::enable_if_t<is_simple_data<T>, int>>
auto operator-(OpFractionLiteral<N, D> const&, OpLiteral<T> const& value)
{
	return static_cast<scalar_t>(OpFractionLiteral<N, D>{}) - value;
}

template<size_t N, size_t D, typename T,
	typename = std::enable_if_t<is_simple_data<T>, int>>
auto operator+(OpFractionLiteral<N, D> const&, OpLiteral<T> const& value)
{
	return static_cast<scalar_t>(OpFractionLiteral<N, D>{}) + value;
}

template<size_t N, size_t D, typename T, 
	typename = std::enable_if_t<is_simple_data<T>, int>>
auto operator-(OpLiteral<T> const& value, OpFractionLiteral<N, D> const&)
{
	return value - static_cast<scalar_t>(OpFractionLiteral<N, D>{});
}

template<size_t N, size_t D, typename T,
	typename = std::enable_if_t<is_simple_data<T>, int>>
auto operator+(OpLiteral<T> const& value, OpFractionLiteral<N, D> const&)
{
	return value + static_cast<scalar_t>(OpFractionLiteral<N, D>{});
}

template<size_t N, size_t D, typename T,
	typename = std::enable_if_t<is_simple_data<T>, int>>
	auto operator-(OpNegFractionLiteral<N, D> const&, OpLiteral<T> const& value)
{
	return static_cast<scalar_t>(OpNegFractionLiteral<N, D>{}) - value;
}

template<size_t N, size_t D, typename T,
	typename = std::enable_if_t<is_simple_data<T>, int>>
	auto operator+(OpNegFractionLiteral<N, D> const&, OpLiteral<T> const& value)
{
	return static_cast<scalar_t>(OpNegFractionLiteral<N, D>{}) + value;
}

template<size_t N, size_t D, typename T,
	typename = std::enable_if_t<is_simple_data<T>, int>>
	auto operator-(OpLiteral<T> const& value, OpNegFractionLiteral<N, D> const&)
{
	return value - static_cast<scalar_t>(OpNegFractionLiteral<N, D>{});
}

template<size_t N, size_t D, typename T,
	typename = std::enable_if_t<is_simple_data<T>, int>>
	auto operator+(OpLiteral<T> const& value, OpNegFractionLiteral<N, D> const&)
{
	return value + static_cast<scalar_t>(OpNegFractionLiteral<N, D>{});
}

// divide with normal scalar values

template<size_t N, size_t D, typename T,
	typename = std::enable_if_t<is_simple_data<T>, int>>
auto operator/(OpFractionLiteral<N, D> const&, OpLiteral<T> const& value)
{
	return static_cast<scalar_t>(OpFractionLiteral<N, D>{}) * expr::inverse(value);
}

template<size_t N, size_t D, typename T,
	typename = std::enable_if_t<is_simple_data<T>, int>>
auto operator/(OpLiteral<T> const& value, OpFractionLiteral<N, D> const&)
{
	return value * static_cast<scalar_t>(expr::inverse(OpFractionLiteral<N, D>{}));
}

template<size_t N, size_t D, typename T,
	typename = std::enable_if_t<is_simple_data<T>, int>>
auto operator/(OpNegFractionLiteral<N, D> const&, OpLiteral<T> const& value)
{
	return static_cast<scalar_t>(OpNegFractionLiteral<N, D>{}) * expr::inverse(value);
}

template<size_t N, size_t D, typename T,
	typename = std::enable_if_t<is_simple_data<T>, int>>
auto operator/(OpLiteral<T> const& value, OpNegFractionLiteral<N, D> const&)
{
	return value * static_cast<scalar_t>(expr::inverse(OpNegFractionLiteral<N, D>{}));
}


// multiplication and division with expressions

template<typename E, size_t N, size_t D>
auto operator*(OpExpression<E> const& e, OpFractionLiteral<N, D>)
{
	return OpFractionLiteral<N, D>{} * *static_cast<E const*>(&e);
}

template<typename E, size_t N, size_t D>
auto operator*(OpExpression<E> const& e, OpNegFractionLiteral<N, D>)
{
	return OpNegFractionLiteral<N, D>{} * *static_cast<E const*>(&e);
}

template<typename E, size_t N, size_t D>
auto operator/(OpExpression<E> const& e, OpFractionLiteral<N, D>)
{
	return OpFractionLiteral<D, N>{} * *static_cast<E const*>(&e);
}

template<typename E, size_t N, size_t D>
auto operator/(OpExpression<E> const& e, OpNegFractionLiteral<N, D>)
{
	return OpNegFractionLiteral<D, N>{} * *static_cast<E const*>(&e);
}

template<typename E, size_t N, size_t D>
auto operator/(OpFractionLiteral<N, D>, OpExpression<E> const& e)
{
	return OpFractionLiteral<N, D>{} * expr::inverse(*static_cast<E const*>(&e));
}

template<typename E, size_t N, size_t D>
auto operator/(OpNegFractionLiteral<N, D>, OpExpression<E> const& e)
{
	return OpNegFractionLiteral<N, D>{} * expr::inverse(*static_cast<E const*>(&e));
}


// multiplication and division between fractions

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator*(OpFractionLiteral<N1, D1>, OpFractionLiteral<N2, D2>)
{
	return expr::make_fraction<N1* N2, D1* D2>();
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator*(OpNegFractionLiteral<N1, D1>, OpFractionLiteral<N2, D2>)
{
	return -expr::make_fraction<N1* N2, D1* D2>();
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator*(OpFractionLiteral<N1, D1>, OpNegFractionLiteral<N2, D2>)
{
	return -expr::make_fraction<N1* N2, D1* D2>();
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator*(OpNegFractionLiteral<N1, D1>, OpNegFractionLiteral<N2, D2>)
{
	return expr::make_fraction<N1* N2, D1* D2>();
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator/(OpFractionLiteral<N1, D1>, OpFractionLiteral<N2, D2>)
{
	return expr::make_fraction<N1* D2, D1* N2>();
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator/(OpNegFractionLiteral<N1, D1>, OpFractionLiteral<N2, D2>)
{
	return -expr::make_fraction<N1* D2, D1* N2>();
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator/(OpFractionLiteral<N1, D1>, OpNegFractionLiteral<N2, D2>)
{
	return -expr::make_fraction<N1* D2, D1* N2>();
}

template<size_t N1, size_t D1, size_t N2, size_t D2>
constexpr auto operator/(OpNegFractionLiteral<N1, D1>, OpNegFractionLiteral<N2, D2>)
{
	return expr::make_fraction<N1 * D2, D1 * N2>();
}

template<size_t N1, size_t D1>
constexpr auto operator/(OpFractionLiteral<N1, D1>, OpIdentity)
{
	return OpFractionLiteral<N1, D1>{};
}

template<size_t N2, size_t D2>
constexpr auto operator/(OpIdentity, OpFractionLiteral<N2, D2>)
{
	return OpFractionLiteral<D2, N2>{};
}

template<size_t N1, size_t D1>
constexpr auto operator/(OpFractionLiteral<N1, D1>, OpNegIdentity)
{
	return OpNegFractionLiteral<N1, D1>{};
}

template<size_t N2, size_t D2>
constexpr auto operator/(OpNegIdentity, OpFractionLiteral<N2, D2>)
{
	return OpNegFractionLiteral<D2, N2>{};
}

template<size_t N1, size_t D1>
constexpr auto operator/(OpNegFractionLiteral<N1, D1>, OpIdentity)
{
	return OpNegFractionLiteral<N1, D1>{};
}

template<size_t N2, size_t D2>
constexpr auto operator/(OpIdentity, OpNegFractionLiteral<N2, D2>)
{
	return OpNegFractionLiteral<D2, N2>{};
}

template<size_t N1, size_t D1>
constexpr auto operator/(OpNegFractionLiteral<N1, D1>, OpNegIdentity)
{
	return OpFractionLiteral<N1, D1>{};
}

template<size_t N2, size_t D2>
constexpr auto operator/(OpNegIdentity, OpNegFractionLiteral<N2, D2>)
{
	return OpFractionLiteral<D2, N2>{};
}

// ******************************************************************************************



/*
 *
 * Collection of like terms.
 *
 ******************************************************************************/

namespace expr
{
	//! Returns true if the base type matches between two types.
	/*!
	 * The base type of a type is evaluated by expr::base_data_type, and 
	 * this will return true if the base type of two types is the same.
	 */
	template<typename G1, typename G2>
	constexpr bool is_same_base = std::is_same<
		typename expr::base_data_type<G1>::type,
		typename expr::base_data_type<G2>::type>::value;
}

//! Addition of two variables with data that can be combined.
template<typename A, typename B, typename Dd, typename G, typename Sp, typename std::enable_if<expr::is_combinable<G>, int>::type = 0>
auto operator+(OpDerivative<Dd, A, OpTerm<OpIdentity, G>, Sp> const& a, OpDerivative<Dd, B, OpTerm<OpIdentity, G>, Sp> const& b)
{
	return (expr::coeff(a) + expr::coeff(b)) * OpDerivative<Dd, OpIdentity, OpTerm<OpIdentity, G>, Sp>(OpIdentity{}, expr::get_enclosed_expression(a), a.solver);
}

//! Subtraction of two variables with data that can be combined.
template<typename A, typename B, typename Dd, typename G, typename Sp, typename std::enable_if<expr::is_combinable<G>, int>::type = 0>
auto operator-(OpDerivative<Dd, A, OpTerm<OpIdentity, G>, Sp> const& a, OpDerivative<Dd, B, OpTerm<OpIdentity, G>, Sp> const& b)
{
	return (expr::coeff(a) - expr::coeff(b)) * OpDerivative<Dd, OpIdentity, OpTerm<OpIdentity, G>, Sp>(OpIdentity{}, expr::get_enclosed_expression(a), a.solver);
}

//! Addition of two variables with data that can be combined.
template<typename A, typename B, typename E, typename T, typename std::enable_if<expr::is_combinable<E>, int>::type = 0>
auto operator+(OpIntegral<A, E, T> const& a, OpIntegral<B, E, T> const& b)
{
	return (expr::coeff(a) + expr::coeff(b)) * OpIntegral<OpIdentity, E, T>(OpIdentity{}, expr::get_enclosed_expression(a));
}

//! Subtraction of two variables with data that can be combined.
template<typename A, typename B, typename E, typename T, typename std::enable_if<expr::is_combinable<E>, int>::type = 0>
auto operator-(OpIntegral<A, E, T> const& a, OpIntegral<B, E, T> const& b)
{
	return (expr::coeff(a) - expr::coeff(b)) * OpIntegral<OpIdentity, E, T>(OpIdentity{}, expr::get_enclosed_expression(a));
}


/*
 *
 * Overloads to simplify expressions for division.
 *
 ******************************************************************************/









/* Division rules for some other objects but will certainly not be used in
 * practice.
 */

template<typename V1, typename V2, typename Dd, typename E, typename Sp>
auto operator/(OpDerivative<Dd, V1, E, Sp> const& a, OpDerivative<Dd, V2, E, Sp> const& b)
{
	return expr::make_literal(a.value * expr::inverse(b.value));
}




/*
 *
 * Distributivity rules.
 *
 ******************************************************************************/



namespace expr
{
	template<typename E0, typename... Es, size_t... Is>
	auto pick_like_terms(symphas::lib::types_list<Es...>, std::index_sequence<Is...>)
	{
		// True on the index if there is a like term. 
		using mask_t = std::integer_sequence<bool, !std::is_same<OpAdd<E0, Es>, add_result_t<E0, Es>>::value...>;

		//return std::make_pair(
		//	symphas::lib::seq_join_t<std::index_sequence<>, std::conditional_t<symphas::lib::seq_index_value<Is, mask_t>::value, std::index_sequence<Is>, std::index_sequence<>>...>{},
		//	symphas::lib::seq_join_t<std::index_sequence<>, std::conditional_t<!symphas::lib::seq_index_value<Is, mask_t>::value, std::index_sequence<Is>, std::index_sequence<>>...>{}
		//);
		return symphas::lib::seq_join_t<std::index_sequence<>, std::conditional_t<symphas::lib::seq_index_value<Is, mask_t>::value, std::index_sequence<Is>, std::index_sequence<>>...>{};
	}

	template<size_t N0, typename... Es, size_t... Is>
	auto pick_like_terms(symphas::lib::types_list<> const&, symphas::lib::types_list<Es...>, std::index_sequence<Is...>)
	{
		return symphas::lib::types_list<symphas::lib::types_list<>, symphas::lib::types_list<>>{};
	}

	// Returns a types list of sequences, such that for each sequence, its first entry represents the index from E0, E0s which
	// matches with the entries from Es using the other entries of the sequence. Is... indexes Es... .
	template<size_t N0, typename E0, typename... E0s, typename... Es, size_t... Is>
	auto pick_like_terms(symphas::lib::types_list<E0, E0s...> const&, symphas::lib::types_list<Es...>, std::index_sequence<Is...>)
	{
		using namespace symphas::lib;
		if constexpr (sizeof...(Es) == 0)
		{
			return types_list<types_list<>, types_list<>>{};
		}
		else
		{
			using match_seq = decltype(pick_like_terms<E0>(
				types_list<Es...>{},
				std::make_index_sequence<sizeof...(Es)>{}));

			using seq_choose = seq_join_t<select_types<match_seq, std::index_sequence<Is>...>>;
			using seq_next = filter_seq_t<std::index_sequence<Is...>, seq_choose>;

			using rest_seq = decltype(pick_like_terms<N0 + 1>(
				types_list<E0s...>{},
				filter_types_on_index<match_seq, Es...>{},
				seq_next{}));
			
			return types_list<
				expand_types_list<std::index_sequence<N0>, type_at_index<0, unroll_types_list<rest_seq>>>,
				expand_types_list<seq_choose, type_at_index<1, unroll_types_list<rest_seq>>>>{};
		}
	}

	template<typename E0, typename... Es, size_t... Is>
	auto apply_add(E0 const& e0, OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return (e0 + ... + expr::get<Is>(e));
	}

	template<typename... Es, size_t... Is>
	auto neg_terms(OpAdd<Es...> const& e, std::index_sequence<Is...>)
	{
		return expr::make_add(-expr::get<Is>(e)...);
	}

	template<typename... Es, size_t... Is, typename... E0s>
	auto compile_add(OpAdd<Es...> const& pick, std::index_sequence<Is...>, E0s&&... es)
	{
		if constexpr (sizeof...(Is) == sizeof...(Es))
		{
			return make_add(pick, std::forward<E0s>(es)...);
		}
		else
		{
			return make_add(expr::get<Is>(pick)..., std::forward<E0s>(es)...);
		}
	}

	template<typename... Es, size_t N0, size_t N1, size_t... Ns, typename... E0s>
	auto compile_add_remove(OpAdd<Es...> const& pick, std::index_sequence<N0, N1, Ns...>, E0s&&... es)
	{
		return compile_add(pick, symphas::lib::filter_seq_t<std::make_index_sequence<sizeof...(Es)>, std::index_sequence<N0, N1, Ns...>>{}, std::forward<E0s>(es)...);
	}

	template<typename... Es, size_t N0, typename... E0s>
	auto compile_add_remove(OpAdd<Es...> const& pick, std::index_sequence<N0>, E0s&&... es)
	{
		return make_add(expr::terms_after_n<N0>(pick), expr::terms_before_n<N0>(pick), std::forward<E0s>(es)...);
	}

	template<typename... Es, typename... E0s>
	auto compile_add_remove(OpAdd<Es...> const& pick, std::index_sequence<>, E0s&&... es)
	{
		return make_add(pick, std::forward<E0s>(es)...);
	}


	template<typename A, typename... Es>
	auto add_one(A const& add, OpAdd<Es...> const& e)
	{
		auto like_seq = pick_like_terms<A>(symphas::lib::types_list<Es...>{}, std::make_index_sequence<sizeof...(Es)>{});
		return compile_add_remove(e, like_seq, apply_add(add, e, like_seq));
	}

	template<typename B, typename... Es>
	auto add_one(OpAdd<Es...> const& e, B const& add)
	{
		auto like_seq = pick_like_terms<B>(symphas::lib::types_list<Es...>{}, std::make_index_sequence<sizeof...(Es)>{});
		return compile_add_remove(e, like_seq, apply_add(add, e, like_seq));
	}


	template<typename E0>
	auto add_all(E0 const& e0)
	{
		return e0;
	}

	template<typename E0, typename E1>
	auto add_all(E0 const& e0, E1 const& e1)
	{
		return e0 + e1;
	}

	template<typename E0, typename E1, typename E2, typename... Es>
	auto add_all(E0 const& e0, E1 const& e1, E2 const& e2, Es const& ...es)
	{
		return (e0 + e1 + (e2 + ... + es));
	}

	template<typename E0, typename E1, typename... Es, size_t... Is>
	auto add_all(OpAdd<E0> const& e, std::index_sequence<Is...>)
	{
		return expr::get<0>(e);
	}

	template<typename E0, typename E1, typename... Es, size_t... Is>
	auto add_all(OpAdd<E0, E1, Es...> const& e, std::index_sequence<Is...>)
	{
		return (expr::get<Is>(e) + ...);
	}

	template<typename... As, typename... E0s, size_t... Is, size_t... Js>
	auto compile_add_group(
		OpAdd<As...> const& add, std::index_sequence<Is...>,
		OpAdd<E0s...> const& add0, std::index_sequence<Js...>)
	{
		if constexpr (sizeof...(Is) == sizeof...(As))
		{
			return expr::make_add(add, expr::get<Js>(add0)...);
		}
		else
		{
			return expr::make_add(expr::get<Is>(add)..., expr::get<Js>(add0)...);
		}
	}

	template<typename... As, typename... E0s, size_t... Is, size_t... Js>
	auto add_group_apply(
		OpAdd<As...> const& add, std::index_sequence<Is...>,
		OpAdd<E0s...> const& add0, std::index_sequence<Js...>,
		std::index_sequence<>, symphas::lib::types_list<>)
	{
		using namespace symphas::lib;
		return expr::make_add(add, add0);
	}
	
	template<typename... As, typename... E0s, size_t... Is, size_t... Js, size_t P0, size_t... Ps, typename Seq0, typename... Seqs>
	auto add_group_apply(
		OpAdd<As...> const& add, std::index_sequence<Is...>, 
		OpAdd<E0s...> const& add0, std::index_sequence<Js...>,
		std::index_sequence<P0, Ps...>, symphas::lib::types_list<Seq0, Seqs...>)
	{
		using namespace symphas::lib;

		using unmatched_seq_a = seq_skip_indices<sizeof...(Is), std::index_sequence<P0, Ps...>>;
		using unmatched_seq_b = seq_skip_indices<sizeof...(Js), sorted_seq<seq_join_t<Seq0, Seqs...>>>;

		return expr::make_add(
			compile_add_group(add, unmatched_seq_a{}, add0, unmatched_seq_b{}),
			(expr::get<P0>(add) + compile_add(add0, Seq0{})), (expr::get<Ps>(add) + compile_add(add0, Seqs{}))...);
	}


	template<typename... As, typename... E0s, size_t... Is>
	auto add_group(OpAdd<As...> const& add, std::index_sequence<Is...>, OpAdd<E0s...> const& add0)
	{
		using namespace symphas::lib;

		using match_terms_t = decltype(pick_like_terms<0>(types_list<As...>{}, types_list<E0s...>{}, std::make_index_sequence<sizeof...(E0s)>{}));
		using pick_a = seq_join_t<type_at_index<0, unroll_types_list<match_terms_t>>>;
		using match_b = type_at_index<1, unroll_types_list<match_terms_t>>;

		return add_group_apply(
			add, std::index_sequence<Is...>{},
			add0, std::make_index_sequence<sizeof...(E0s)>{},
			pick_a{}, match_b{});
	}

	template<typename... As, typename... E0s, typename T0, typename... Ts, size_t... Is>
	auto add_group(OpAdd<As...> const& add, std::index_sequence<Is...>, OpAdd<E0s...> const& e0, T0&& e1, Ts&&... rest);

	template<typename E>
	auto add_group(E const& e)
	{
		return e;
	}

	inline auto add_group()
	{
		return OpVoid{};
	}

	template<typename... As, typename... E0s, typename... Ts>
	auto add_group(OpAdd<As...> const& add, OpAdd<E0s...> const& e, Ts&&... rest);

	template<typename E1, typename E2, typename... Ts>
	auto add_group(OpExpression<E1> const& e1, OpExpression<E2> const& e2, Ts&&... rest);

	template<typename E, typename... E0s, typename... Ts>
	auto add_group(OpExpression<E> const& add, OpAdd<E0s...> const& e, Ts&&... rest);

	template<typename E, typename... E0s, typename... Ts>
	auto add_group(OpAdd<E0s...> const& e, OpExpression<E> const& add, Ts&&... rest);

	template<typename... As, typename... E0s, typename T0, typename... Ts, size_t... Is>
	auto add_group(OpAdd<As...> const& add, std::index_sequence<Is...>, OpAdd<E0s...> const& e0, T0&& e1, Ts&&... rest);

	template<typename E, typename... E0s, typename... Ts>
	auto add_group(OpExpression<E> const& add, OpAdd<E0s...> const& e, Ts&&... rest)
	{
		return add_group(add_one(*static_cast<E const*>(&add), e), std::forward<Ts>(rest)...);
	}

	template<typename E, typename... E0s, typename... Ts>
	auto add_group(OpAdd<E0s...> const& e, OpExpression<E> const& add, Ts&&... rest)
	{
		return add_group(add_one(e, *static_cast<E const*>(&add)), std::forward<Ts>(rest)...);
	}

	template<typename... As, typename... E0s, typename... Ts>
	auto add_group(OpAdd<As...> const& add, OpAdd<E0s...> const& e, Ts&&... rest)
	{
		return add_group(add, std::make_index_sequence<sizeof...(As)>{}, e, std::forward<Ts>(rest)...);
	}

	template<typename E1, typename E2, typename... Ts>
	auto add_group(OpExpression<E1> const& e1, OpExpression<E2> const& e2, Ts&&... rest)
	{
		return add_group(*static_cast<E1 const*>(&e1) + *static_cast<E2 const*>(&e2), std::forward<Ts>(rest)...);
	}

	template<typename... As, typename... E0s, typename T0, typename... Ts, size_t... Is>
	auto add_group(OpAdd<As...> const& add, std::index_sequence<Is...>, OpAdd<E0s...> const& add0, T0&& e1, Ts&&... rest)
	{

		return add_group(add_group(add, std::index_sequence<Is...>{}, add0), std::forward<T0>(e1), std::forward<Ts>(rest)...);
	}


	//! Distributing addition expression into a addition expression.
	template<typename A0, typename... Bs, size_t... Js,
		typename std::enable_if_t<(std::is_same<mul_result_t<A0, Bs>, OpBinaryMul<A0, Bs>>::value && ...), int> = 0>
	auto distribute_adds(A0 const& a, OpAdd<Bs...> const& b, std::index_sequence<Js...>)
	{
		//return expr::make_add((a * expr::get<Js>(b))...);
		return make_add((a * expr::get<Js>(b))...);
	}

	//! Distributing addition expression into a addition expression.
	template<typename A0, typename... Bs, size_t... Js,
		typename std::enable_if_t<!(std::is_same<mul_result_t<A0, Bs>, OpBinaryMul<A0, Bs>>::value && ...), int> = 0>
	auto distribute_adds(A0 const& a, OpAdd<Bs...> const& b, std::index_sequence<Js...>)
	{
		return ((a * expr::get<Js>(b)) + ...);
	}

	//! Distributing addition expression into a addition expression.
	template<typename... As, typename B0, size_t... Is,
		typename std::enable_if_t<(std::is_same<mul_result_t<As, B0>, OpBinaryMul<As, B0>>::value && ...), int> = 0>
	auto distribute_adds(OpAdd<As...> const& a, B0 const& b, std::index_sequence<Is...>)
	{
		//return expr::make_add((expr::get<Is>(a) * b)...);
		return make_add((expr::get<Is>(a) * b)...);
	}

	//! Distributing addition expression into a addition expression.
	template<typename... As, typename B0, size_t... Is,
		typename std::enable_if_t<!(std::is_same<mul_result_t<As, B0>, OpBinaryMul<As, B0>>::value && ...), int> = 0>
	auto distribute_adds(OpAdd<As...> const& a, B0 const& b, std::index_sequence<Is...>)
	{
		return ((expr::get<Is>(a) * b) + ...);
	}

	//! Distributing addition expression into a addition expression.
	template<typename... As, typename... Bs, size_t... Is, size_t... Js>
	auto distribute_adds(OpAdd<As...> const& a, OpAdd<Bs...> const& b, std::index_sequence<Is...>, std::index_sequence<Js...>)
	{
		return add_group(distribute_adds(expr::get<Is>(a), b, std::index_sequence<Js...>{})...);
	}
}


//! Distributing addition expression into a addition expression.
template<typename... As, typename... Bs>
auto operator*(OpAdd<As...> const& a, OpAdd<Bs...> const& b)
{
	return expr::distribute_adds(a, b, std::make_index_sequence<sizeof...(As)>{}, std::make_index_sequence<sizeof...(Bs)>{});
}

//! Distributing an RHS expression between operands in addition expression.
template<typename... As, typename B2, 
	typename std::enable_if_t<(!expr::is_identity<B2> && !expr::is_fraction<B2> && !expr::is_coeff<OpAdd<As...>>), int> = 0>
auto operator*(OpAdd<As...> const& a, OpExpression<B2> const& b)
{
	return expr::distribute_adds(a, *static_cast<const B2*>(&b),
		std::make_index_sequence<sizeof...(As)>{});
}

//! Distributing an LHS expression between operands in addition expression.
template<typename A1, typename... Bs,
	typename std::enable_if_t<(!expr::is_identity<A1> && !expr::is_coeff<OpAdd<Bs...>>), int> = 0>
auto operator*(OpExpression<A1> const& a, OpAdd<Bs...> const& b)
{
	return expr::distribute_adds(*static_cast<const A1*>(&a), b,
		std::make_index_sequence<sizeof...(Bs)>{});
}

//! Distributing an LHS expression between operands in addition expression.
template<typename A1, typename A2, typename... Bs,
	typename std::enable_if_t<(!expr::is_identity<A1> && !expr::is_coeff<OpAdd<Bs...>>), int> = 0>
auto operator*(OpBinaryMul<A1, A2> const& a, OpAdd<Bs...> const& b)
{
	return expr::distribute_adds(a, b,
		std::make_index_sequence<sizeof...(Bs)>{});
}


//! Distributing an LHS expression between operands in addition expression.
template<typename A1, typename... Bs,
	typename std::enable_if_t<(expr::is_coeff<OpAdd<Bs...>>), int> = 0>
auto operator*(OpExpression<A1> const& a, OpAdd<Bs...> const& b)
{
	return expr::distribute_adds(*static_cast<const A1*>(&a), b,
		std::make_index_sequence<sizeof...(Bs)>{});
}

//! Distributing an LHS expression between operands in addition expression.
template<typename... Bs, typename std::enable_if_t<(expr::is_coeff<OpAdd<Bs...>>), int> = 0>
auto operator*(OpIdentity const& a, OpAdd<Bs...> const& b)
{
	return b;
}


//! Distributing a literal into the first operand in multiplication expression.
template<typename T, typename B1, typename B2>
auto operator*(OpLiteral<T> const& a, OpBinaryMul<B1, B2> const& b)
{
	return expr::make_mul(a * b.a, b.b);
}

//! Distributing a literal into the first operand in division expression.
template<typename T, typename B1, typename B2>
auto operator*(OpLiteral<T> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return expr::make_div(a * b.a, b.b);
}

//! Distributing a literal into the first operand in multiplication expression.
template<typename coeff_t, typename B1, typename B2, typename std::enable_if_t<(expr::is_fraction<coeff_t>), int> = 0>
auto operator*(coeff_t, OpBinaryMul<B1, B2> const& b)
{
	return expr::make_mul(coeff_t{} * b.a, b.b);
}

//! Distributing a literal into the first operand in division expression.
template<typename coeff_t, typename B1, typename B2, typename std::enable_if_t<(expr::is_fraction<coeff_t>), int> = 0>
auto operator*(coeff_t, OpBinaryDiv<B1, B2> const& b)
{
	return expr::make_div(coeff_t{} * b.a, b.b);
}

namespace symphas::internal
{

	 //! Terminates the multiplication recursion of expressions.
	 /*!
	  * These multiplication rules take into account multiplication between variables
	  * and existing multiplication objects, so that variables which don't multiply
	  * with each other can be multiplied into a new variable which allows putting
	  * them into a multi variable.
	  *
	  * Consequently this assumes all multiplication is associative.
	  *
	  * Another assumption is that multiplication between objects that don't
	  * multiply with each other is commutative. i.e.:
	  * consider objects \f$A\f$, \f$B\f$ and \f$C\f$. \f$A\f$ and \f$B\f$ can be
	  * multiplied with each other to produce a multi variable variable, but \f$C\f$
	  * can only be multiplied with itself to produce a multi variable. Then

template<typename E1, typename E2,
	typename std::enable_if_t<(expr::eval_type<E1>::rank == 0 || expr::eval_type<E2>::rank == 0), int> = 0>
auto operator*(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_mul(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}

template<typename E1, typename E2,
	typename std::enable_if_t<(expr::eval_type<E1>::rank == 0 || expr::eval_type<E2>::rank == 0), int> = 0>
auto operator*(OpExpression<E1> const& a, OpOperator<E2> const& b)
{
	return OpOperatorChain(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}

template<typename E1, typename E2,
	typename std::enable_if_t<(expr::eval_type<E1>::rank > 0 && expr::eval_type<E2>::rank > 0), int> = 0>
	auto operator*(OpExpression<E1> const& a, OpOperator<E2> const& b)
{
	return expr::dot(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}

template<typename A, typename B, typename E>
auto operator*(OpBinaryMul<A, B> const& a, OpOperator<E> const& b)
{
	return OpOperatorChain(a, *static_cast<const E*>(&b));
}

template<typename E1, typename E2>
auto operator*(OpOperator<E1> const& a, OpOperator<E2> const& b)
{
	return (*static_cast<const E1*>(&a)).operator*(*static_cast<const E2*>(&b));
}

	  * \f$A \cdot B = AB\f$, but \f$A \cdot C = A \cdot C\f$ and
	  * \f$B \cdot C = B \cdot C\f$ however, the last two expressions are assumed to
	  * commute, e.g. \f$A \cdot C = C \cdot A\f$.
	  */
	template<typename E1, typename E2>
	auto terminate_mul(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return expr::make_mul(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename A1, typename A2, typename E>
	auto terminate_mul(OpBinaryMul<A1, A2> const& a, OpExpression<E> const& b)
	{
		return expr::make_mul(a.a, a.b * (*static_cast<E const*>(&b)));
	}

}

/*
 *
 *
 *
 * Recursive behaviour of addition; adding many terms together requires
 * simplifying if possible and collection of like terms if possible. This is
 * performed by adding together newly added terms with all other ones in the
 * sum, recursively having all terms be added or subtracted with each other.
 *
 ******************************************************************************/

template<typename... As, typename E>
auto operator+(OpExpression<E> const& a, OpAdd<As...> const& b)
{
	return expr::add_one(*static_cast<E const*>(&a), b);
}

template<typename... As, typename E>
auto operator+(OpOperator<E> const& a, OpAdd<As...> const& b)
{
	return expr::add_one(*static_cast<E const*>(&a), b);
}

template<typename... As, typename E>
auto operator+(OpAdd<As...> const& a, OpExpression<E> const& b)
{
	return expr::add_one(a, *static_cast<E const*>(&b));
}

template<typename... As, typename E>
auto operator+(OpAdd<As...> const& a, OpOperator<E> const& b)
{
	return expr::add_one(a, *static_cast<E const*>(&b));
}

template<typename... Bs, typename E>
auto operator-(OpExpression<E> const& a, OpAdd<Bs...> const& b)
{
	return expr::add_one(*static_cast<E const*>(&a), expr::neg_terms(b, std::make_index_sequence<sizeof...(Bs)>{}));
}

template<typename... Bs, typename E>
auto operator-(OpOperator<E> const& a, OpAdd<Bs...> const& b)
{
	return expr::add_one(*static_cast<E const*>(&a), expr::neg_terms(b, std::make_index_sequence<sizeof...(Bs)>{}));
}

template<typename... As, typename E>
auto operator-(OpAdd<As...> const& a, OpExpression<E> const& b)
{
	return expr::add_one(a, -*static_cast<E const*>(&b));
}

template<typename... As, typename E>
auto operator-(OpAdd<As...> const& a, OpOperator<E> const& b)
{
	return expr::add_one(a, -*static_cast<E const*>(&b));
}

template<typename... As, typename... Bs>
auto operator+(OpAdd<As...> const& a, OpAdd<Bs...> const& b)
{
	return expr::add_group(a, std::make_index_sequence<sizeof...(As)>{}, b);
}

template<typename... As, typename... Bs>
auto operator-(OpAdd<As...> const& a, OpAdd<Bs...> const& b)
{
	return expr::add_group(a, std::make_index_sequence<sizeof...(As)>{}, expr::neg_terms(b, std::make_index_sequence<sizeof...(Bs)>{}));
}


/*
 *
 * Overloads to simplify expressions for multiplication.
 *
 ******************************************************************************/

//
//template<typename A1, typename A2, typename V, typename... Gs, expr::exp_key_t... Xs>
//auto operator*(OpBinaryMul<A1, A2> const& a, OpTerms<V, Term<Gs, Xs>...> const& b)
//{
//	return symphas::internal::terminate_mul(a.a, (a.b * expr::coeff(b)) * OpTerms(OpIdentity{}, expr::terms_after_first(b)));
//}
//
//template<typename V, typename... Gs, expr::exp_key_t... Xs, typename B1, typename B2>
//auto operator*(OpTerms<V, Term<Gs, Xs>...> const& a, OpBinaryMul<B1, B2> const& b)
//{
//	return symphas::internal::terminate_mul(a * b.a, b.b);
//}

template<typename A1, typename A2, typename E/*, typename std::enable_if_t<!expr::is_coeff<E>, int> = 0*/>
auto operator*(OpBinaryMul<A1, A2> const& a, OpExpression<E> const& b)
{
	return expr::make_mul(a.a, a.b * (*static_cast<E const*>(&b)));
}

template<typename A1, typename A2>
auto operator*(OpBinaryMul<A1, A2> const& a, OpIdentity)
{
	return a;
}

template<typename A1, typename A2>
auto operator*(OpBinaryMul<A1, A2> const& a, OpNegIdentity)
{
	return -a;
}

template<typename A1, typename A2, size_t N, size_t D>
auto operator*(OpBinaryMul<A1, A2> const& a, OpFractionLiteral<N, D>)
{
	return expr::make_mul(a.a, OpFractionLiteral<N, D>{} *a.b);
}

template<typename A1, typename A2, size_t N, size_t D>
auto operator*(OpBinaryMul<A1, A2> const& a, OpNegFractionLiteral<N, D>)
{
	return expr::make_mul(a.a, OpNegFractionLiteral<N, D>{} * a.b);
}

//
//template<typename A1, typename A2, typename E, typename std::enable_if_t<!expr::is_coeff<E>, int> = 0>
//auto operator*(OpExpression<E> const& a, OpBinaryMul<A1, A2> const& b)
//{
//	return expr::make_mul((*static_cast<E const*>(&a)) * b.a,  b.b);
//}

template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryMul<A1, A2> const& a, OpBinaryMul<B1, B2> const& b)
{
	return expr::make_mul(a.a, a.b * b);
}






namespace symphas::internal
{

	//! Terminates the recursion of the division operation.
	/*!
	 * Endpoint of the recursion for applying the division rules, which **multiplies**
	 * two reduced terms, with a key emphasis on reduced terms. There are several
	 * rules for handling various reduced forms. There is also the capability to
	 * automatically reduce terms. Makes the assumption that a division is always
	 * the first term given to the function.
	 */
	template<typename A1, typename A2, typename E,
		typename std::enable_if_t<(expr::factor_list_all<E, A2>::value == 0), int> = 0>
	auto terminate_div(OpBinaryDiv<A1, A2> const& a, OpExpression<E> const& b)
	{
		return expr::make_div(a.a * (*static_cast<const E*>(&b)), a.b);
	}

	template<typename A1, typename A2, typename E,
		typename std::enable_if_t<(expr::factor_list_all<E, A2>::value > 0), int> = 0>
	auto terminate_div(OpBinaryDiv<A1, A2> const& a, OpExpression<E> const& b)
	{
		return a.a * ((*static_cast<const E*>(&b)) / a.b);
	}

	template<typename B1, typename B2, typename E,
		typename std::enable_if_t<(expr::factor_list_all<E, B2>::value == 0), int> = 0>
	auto terminate_div(OpExpression<E> const& a, OpBinaryDiv<B1, B2> const& b)
	{
		return expr::make_div((*static_cast<const E*>(&a)) * b.a, b.b);
	}

	template<typename B1, typename B2, typename E,
		typename std::enable_if_t<(expr::factor_list_all<E, B2>::value > 0), int> = 0>
	auto terminate_div(OpExpression<E> const& a, OpBinaryDiv<B1, B2> const& b)
	{
		return ((*static_cast<const E*>(&a)) / b.b) * b.a;
	}

	template<typename A1, typename A2, typename B1, typename B2,
		typename std::enable_if_t<(expr::factor_list_all<A1, B2>::value == 0 && expr::factor_list_all<A2, B1>::value == 0), int> = 0>
	auto terminate_div(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
	{
		return expr::make_div(
			a.a * b.a,
			a.b * b.b);
	}

	template<typename A1, typename A2, typename B1, typename B2,
		typename std::enable_if_t<(expr::factor_list_all<A1, B2>::value > 0 || expr::factor_list_all<A2, B1>::value > 0), int> = 0>
	auto terminate_div(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
	{
		return (a.a * b.a) / (a.b * b.b);
	}
	
	template<typename A1, typename A2, typename B1, typename B2,
		typename std::enable_if_t<(expr::factor_list_all<A1, B2>::value > 0), int> = 0>
	auto terminate_div(OpBinaryMul<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
	{
		return (a.a / b.b) * a.b * b.a;
	}

	template<typename A1, typename A2, typename B1, typename B2,
		typename std::enable_if_t<(expr::factor_list_all<A2, B2>::value > 0), int> = 0>
		auto terminate_div(OpBinaryMul<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
	{
		return a.a * (a.b / b.b) * b.a;
	}

	template<typename A1, typename A2, typename B1, typename B2,
		typename std::enable_if_t<(expr::factor_list_all<OpBinaryMul<A1, A2>, B2>::value == 0), int> = 0>
	auto terminate_div(OpBinaryMul<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
	{
		return expr::make_div(a * b.a, b.b);
	}

	template<typename E, typename T>
	auto terminate_div(OpLiteral<T> const& a, OpExpression<E> const& e)
	{
		return a * (*static_cast<const E*>(&e));
	}

	template<typename E, typename T>
	auto terminate_div(OpExpression<E> const& e, OpLiteral<T> const& a)
	{
		return a * (*static_cast<const E*>(&e));
	}

	template<typename coeff_t, typename A, typename B,
		typename std::enable_if_t<(expr::is_identity<coeff_t> || expr::is_fraction<coeff_t>), int> = 0>
	auto terminate_div(coeff_t, OpBinaryDiv<A, B> const& b)
	{
		return coeff_t{} * expr::make_div(b.a, b.b);
	}

	template<typename coeff_t, typename A, typename B,
		typename std::enable_if_t<(expr::is_identity<coeff_t> || expr::is_fraction<coeff_t>), int> = 0>
	auto terminate_div(OpBinaryDiv<A, B> const& a, coeff_t)
	{
		return coeff_t{} * expr::make_div(a.a, a.b);
	}


	template<typename A1, typename A2, typename E>
	auto terminate_div(OpBinaryMul<A1, A2> const& a, OpExpression<E> const& b)
	{
		return a.a * (a.b * *static_cast<E const*>(&b));
	}

	template<typename B1, typename B2, typename E>
	auto terminate_div(OpExpression<E> const& a, OpBinaryMul<B1, B2> const& b)
	{
		return ((*static_cast<E const*>(&a)) * b.a) * b.b;
	}


	template<typename E1, typename E2>
	auto terminate_div(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return (*static_cast<E1 const*>(&a)) * (*static_cast<E2 const*>(&b));
	}


	template<typename E, typename B1, typename T, typename... Gs, expr::exp_key_t... Xs,
		typename std::enable_if_t<((expr::factor_count<Gs, E>::value == 0) && ...), int> = 0>
	auto terminate_div(OpExpression<E> const& a, OpBinaryDiv<B1, OpTerms<T, Term<Gs, Xs>...>> const& b)
	{
		return expr::make_div(
			expr::inverse(expr::coeff(b.b)) * (*static_cast<E const*>(&a)) * b.a,
			OpTerms(OpIdentity{}, expr::terms_after_first(b.b)));
	}

	template<typename E, typename A1, typename T, typename... Gs, expr::exp_key_t... Xs,
		typename std::enable_if_t<((expr::factor_count<Gs, E>::value == 0) && ...), int> = 0>
	auto terminate_div(OpBinaryDiv<A1, OpTerms<T, Term<Gs, Xs>...>> const& a, OpExpression<E> const& b)
	{
		return expr::make_div(
			expr::inverse(expr::coeff(a.b)) * a.a * (*static_cast<E const*>(&b)),
			OpTerms(OpIdentity{}, expr::terms_after_first(a.b)));
	}

}


/*
 * Recursively handle division between multiplications and divisions.
 */



template<typename A1, typename A2, typename E,
	typename std::enable_if<(expr::factor_list_all<OpBinaryMul<A1, A2>, E>::value == 0 && !expr::is_fraction<E> && !expr::is_identity<E>), int>::type = 0>
auto operator/(OpBinaryMul<A1, A2> const& a, OpExpression<E> const& b)
{
	return expr::make_div(a, *static_cast<E const*>(&b));
	//return a.a * (a.b / (*static_cast<const E*>(&b)));
}

template<typename B1, typename B2, typename E,
	typename std::enable_if<(expr::factor_list_all<E, OpBinaryMul<B1, B2>>::value == 0 && !expr::is_fraction<E> && !expr::is_identity<E>), int>::type = 0>
auto operator/(OpExpression<E> const& a, OpBinaryMul<B1, B2> const& b)
{
	return expr::make_div(*static_cast<E const*>(&a), b);
	//return ((*static_cast<const E*>(&a)) / b.a) * b.b;
}

template<typename A1, typename A2, typename E,
	typename std::enable_if<(expr::factor_list_all<A1, E>::value == 0 && !expr::is_fraction<E> && !expr::is_identity<E>), int>::type = 0>
auto operator/(OpBinaryDiv<A1, A2> const& a, OpExpression<E> const& b)
{
	return a.a / (a.b * (*static_cast<const E*>(&b)));
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator/(OpBinaryDiv<A1, A2> const& a, OpBinaryMul<B1, B2> const& b)
{
	return a.a / (a.b * b);
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator/(OpBinaryMul<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return (a * b.b) / b.a;
}

template<typename B1, typename B2, typename E>
auto operator/(OpExpression<E> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return ((*static_cast<const E*>(&a)) * b.b) / b.a;
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator/(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return (a.a * b.b) / (a.b * b.a);
}

template<typename... As, typename B1, typename B2>
auto operator/(OpAdd<As...> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return (a * b.b) / b.a;
}

template<typename A1, typename A2, typename... Bs>
auto operator/(OpBinaryDiv<A1, A2> const& a, OpAdd<Bs...> const& b)
{
	return a.a / (a.b * b);
}



template<typename E, typename B1, typename V, typename... Gs, expr::exp_key_t... Xs,
	typename std::enable_if_t<(((expr::factor_count<Gs, E>::value == 0) && ...) && !expr::is_identity<E>), int> = 0>
auto operator/(OpExpression<E> const& a, OpTerms<V, Term<Gs, Xs>...> const& b)
{
	return expr::make_div(
		expr::inverse(expr::coeff(b)) * (*static_cast<E const*>(&a)),
		OpTerms(OpIdentity{}, expr::terms_after_first(b)));
}

template<typename E, typename B1, typename V, typename... Gs, expr::exp_key_t... Xs,
	typename std::enable_if_t<(((expr::factor_count<Gs, E>::value == 0) && ...) && !expr::is_identity<E>), int> = 0>
auto operator*(OpExpression<E> const& a, OpBinaryDiv<B1, OpTerms<V, Term<Gs, Xs>...>> const& b)
{
	return expr::make_div(
		expr::inverse(expr::coeff(b.b)) * ((*static_cast<E const*>(&a)) * b.a),
		OpTerms(OpIdentity{}, expr::terms_after_first(b.b)));
}

template<typename E, typename B1, typename V, typename... Gs, expr::exp_key_t... Xs,
	typename std::enable_if_t<(((expr::factor_count<Gs, E>::value == 0) && ...) && !expr::is_identity<E>), int> = 0>
auto operator*(OpOperator<E> const& a, OpBinaryDiv<B1, OpTerms<V, Term<Gs, Xs>...>> const& b)
{
	return (*static_cast<E const*>(&a)).operator()(b);
}


/*
 *
 *
 * Multiplication rules in the context of division, used for canceling like 
 * terms. All the divisions that are parameters to these functions are assumed 
 * to be in reduced form. 
 *
 ******************************************************************************/

template<typename A1, typename A2, typename E, typename std::enable_if_t<!(expr::is_identity<E> || expr::is_fraction<E>), int> = 0>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpExpression<E> const& b)
{
	return (a.a * (*static_cast<const E*>(&b))) / a.b;
}

template<typename E, typename B1, typename B2, typename std::enable_if_t<!expr::is_identity<E>, int> = 0>
auto operator*(OpExpression<E> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return ((*static_cast<const E*>(&a)) * b.a) / b.b;
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return (a.a * b.a) / (a.b * b.b);
}

template<typename A1, typename A2, typename... Bs>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpAdd<Bs...> const& b)
{
	return expr::distribute_adds(a.a, b, std::make_index_sequence<sizeof...(Bs)>{}) / a.b;
}

template<typename... As, typename A2, typename... Bs>
auto operator*(OpBinaryDiv<OpAdd<As...>, A2> const& a, OpAdd<Bs...> const& b)
{
	return expr::distribute_adds(a.a, b, std::make_index_sequence<sizeof...(As)>{}, std::make_index_sequence<sizeof...(Bs)>{}) / a.b;
}

template<typename... As, typename B1, typename B2>
auto operator*(OpAdd<As...> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return expr::distribute_adds(a, b.a, std::make_index_sequence<sizeof...(As)>{}) / b.b;
}

template<typename... As, typename... Bs, typename B2>
auto operator*(OpAdd<As...> const& a, OpBinaryDiv<OpAdd<Bs...>, B2> const& b)
{
	return expr::distribute_adds(a, b.a, std::make_index_sequence<sizeof...(As)>{}, std::make_index_sequence<sizeof...(Bs)>{}) / b.b;
}


// ******************************************************************************************
//
//template<typename A1, typename T, typename G, typename B1, typename B2, 
//	typename std::enable_if_t<(expr::factor_count<G, B2>::value > 0), int> = 0>
//auto operator+(OpBinaryDiv<A1, OpTerm<T, G>> const& a, OpBinaryDiv<B1, B2> const& b)
//{
//	auto div = b.b / a.b;
//	return (a.a * div + b.a) / b.b;
//}
//
//template<typename A1, typename T, typename... Gs, typename B1, typename B2, 
//	typename std::enable_if_t<((expr::factor_count<Gs, B2>::value > 0) && ...), int> = 0>
//auto operator+(OpBinaryDiv<A1, OpNLVariable<T, Gs...>> const& a, OpBinaryDiv<B1, B2> const& b)
//{
//	auto div = b.b / a.b;
//	return (a.a * div + b.a) / b.b;
//}
//
//template<typename A1, typename T, typename G, typename B1, typename B2, 
//typename std::enable_if_t<(expr::factor_count<G, B2>::value > 0), int> = 0>
//auto operator-(OpBinaryDiv<A1, OpTerm<T, G>> const& a, OpBinaryDiv<B1, B2> const& b)
//{
//	auto div = b.b / a.b;
//	return (a.a * div - b.a) / b.b;
//}
//
//template<typename A1, typename T, typename... Gs, typename B1, typename B2, 
//typename std::enable_if_t<((expr::factor_count<Gs, B2>::value > 0) && ...), int> = 0>
//auto operator-(OpBinaryDiv<A1, OpNLVariable<T, Gs...>> const& a, OpBinaryDiv<B1, B2> const& b)
//{
//	auto div = b.b / a.b;
//	return (a.a * div - b.a) / b.b;
//}
//
//template<typename A1, typename A2, typename B1, typename T, typename G, 
//typename std::enable_if_t<(expr::factor_count<G, A2>::value > 0), int> = 0>
//auto operator+(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, OpTerm<T, G>> const& b)
//{
//	auto div = a.b / b.b;
//	return (a.a + b.a * div) / a.b;
//}
//
//template<typename A1, typename A2, typename B1, typename T, typename... Gs, 
//typename std::enable_if_t<((expr::factor_count<Gs, A2>::value > 0) && ...), int> = 0>
//auto operator+(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, OpNLVariable<T, Gs...>> const& b)
//{
//	auto div = a.b / b.b;
//	return (a.a + b.a * div) / a.b;
//}
//
//template<typename A1, typename A2, typename B1, typename T, typename G, 
//typename std::enable_if_t<(expr::factor_count<G, A2>::value > 0), int> = 0>
//auto operator-(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, OpTerm<T, G>> const& b)
//{
//	auto div = a.b / b.b;
//	return (a.a - b.a * div) / a.b;
//}
//
//template<typename A1, typename A2, typename B1, typename T, typename... Gs, 
//	typename std::enable_if_t<((expr::factor_count<Gs, A2>::value > 0) && ...), int> = 0>
//auto operator-(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, OpNLVariable<T, Gs...>> const& b)
//{
//	auto div = a.b / b.b;
//	return (a.a - b.a * div) / a.b;
//}
//
//
//
//template<typename A1, typename T1, typename G1, typename B1, typename T2, typename G2, 
//	typename std::enable_if_t<(expr::factor_count<G1, G2>::value > 0), int> = 0>
//auto operator+(OpBinaryDiv<A1, OpTerm<T1, G1>> const& a, OpBinaryDiv<B1, OpTerm<T2, G2>> const& b)
//{
//	auto div = b.b / a.b;
//	return symphas::internal::terminate_div((a.a * div + b.a), expr::inverse(b.b));
//}
//
//template<typename A1, typename T1, typename... G1s, typename B1, typename T2, typename G2, 
//	typename std::enable_if_t<(expr::factor_count_list<G2, G1s...>::value > 0), int> = 0>
//auto operator+(OpBinaryDiv<A1, OpNLVariable<T1, G1s...>> const& a, OpBinaryDiv<B1, OpTerm<T2, G2>> const& b)
//{
//	auto div = a.b / b.b;
//	return symphas::internal::terminate_div((a.a + b.a * div), expr::inverse(a.b));
//}
//
//template<typename A1, typename T1, typename G1, typename B1, typename T2, typename... G2s, 
//	typename std::enable_if_t<(expr::factor_count_list<G1, G2s...>::value > 0), int> = 0>
//auto operator+(OpBinaryDiv<A1, OpTerm<T1, G1>> const& a, OpBinaryDiv<B1, OpNLVariable<T2, G2s...>> const& b)
//{
//	auto div = b.b / a.b;
//	return symphas::internal::terminate_div((a.a * div + b.a), expr::inverse(b.b));
//}
//
//template<typename A1, typename T1, typename... G1s, typename B1, typename T2, typename... G2s, 
//	typename std::enable_if_t<(
//		((expr::factor_count_list<G1s, G2s...>::value >= expr::factor_count_list<G1s, G1s...>::value) && ...) &&
//		!((expr::factor_count_list<G1s, G2s...>::value == expr::factor_count_list<G1s, G1s...>::value) && ...)
//		), int> = 0>
//auto operator+(OpBinaryDiv<A1, OpNLVariable<T1, G1s...>> const& a, OpBinaryDiv<B1, OpNLVariable<T2, G2s...>> const& b)
//{
//	auto div = b.b / a.b;
//	return symphas::internal::terminate_div((a.a * div + b.a), expr::inverse(b.b));
//}
//
//template<typename A1, typename T1, typename... G1s, typename B1, typename T2, typename... G2s, 
//	typename std::enable_if_t<(
//		((expr::factor_count_list<G2s, G1s...>::value >= expr::factor_count_list<G2s, G2s...>::value) && ...) &&
//		!((expr::factor_count_list<G2s, G1s...>::value == expr::factor_count_list<G2s, G2s...>::value) && ...)
//		), int> = 0>
//auto operator+(OpBinaryDiv<A1, OpNLVariable<T1, G1s...>> const& a, OpBinaryDiv<B1, OpNLVariable<T2, G2s...>> const& b)
//{
//	auto div = a.b / b.b;
//	return symphas::internal::terminate_div((a.a + b.a * div), expr::inverse(a.b));
//}
//
//template<typename A1, typename T1, typename... G1s, typename B1, typename T2, typename... G2s,
//	typename std::enable_if_t<((expr::factor_count_list<G1s, G2s...>::value == expr::factor_count_list<G1s, G1s...>::value) && ...), int> = 0>
//auto operator+(OpBinaryDiv<A1, OpNLVariable<T1, G1s...>> const& a, OpBinaryDiv<B1, OpNLVariable<T2, G2s...>> const& b)
//{
//	auto div = a.b / b.b;
//	return symphas::internal::terminate_div((a.a + b.a * div), expr::inverse(a.b));
//}
//
//
//template<typename A1, typename T1, typename G1, typename B1, typename T2, typename G2, 
//	typename std::enable_if_t<(expr::factor_count<G1, G2>::value > 0), int> = 0>
//auto operator-(OpBinaryDiv<A1, OpTerm<T1, G1>> const& a, OpBinaryDiv<B1, OpTerm<T2, G2>> const& b)
//{
//	auto div = b.b / a.b;
//	return symphas::internal::terminate_div((a.a * div - b.a), expr::inverse(b.b));
//}
//
//template<typename A1, typename T1, typename... G1s, typename B1, typename T2, typename G2, 
//	typename std::enable_if_t<(expr::factor_count_list<G2, G1s...>::value > 0), int> = 0>
//auto operator-(OpBinaryDiv<A1, OpNLVariable<T1, G1s...>> const& a, OpBinaryDiv<B1, OpTerm<T2, G2>> const& b)
//{
//	auto div = a.b / b.b;
//	return symphas::internal::terminate_div((a.a - b.a * div), expr::inverse(a.b));
//}
//
//template<typename A1, typename T1, typename G1, typename B1, typename T2, typename... G2s, 
//	typename std::enable_if_t<(expr::factor_count_list<G1, G2s...>::value > 0), int> = 0>
//auto operator-(OpBinaryDiv<A1, OpTerm<T1, G1>> const& a, OpBinaryDiv<B1, OpNLVariable<T2, G2s...>> const& b)
//{
//	auto div = b.b / a.b;
//	return symphas::internal::terminate_div((a.a * div - b.a), expr::inverse(b.b));
//}
//
//template<typename A1, typename T1, typename... G1s, typename B1, typename T2, typename... G2s,
//	typename std::enable_if_t<(
//		((expr::factor_count_list<G1s, G2s...>::value >= expr::factor_count_list<G1s, G1s...>::value) && ...) &&
//		!((expr::factor_count_list<G1s, G2s...>::value == expr::factor_count_list<G1s, G1s...>::value) && ...)
//		), int> = 0>
//auto operator-(OpBinaryDiv<A1, OpNLVariable<T1, G1s...>> const& a, OpBinaryDiv<B1, OpNLVariable<T2, G2s...>> const& b)
//{
//	auto div = b.b / a.b;
//	return symphas::internal::terminate_div((a.a * div - b.a), expr::inverse(b.b));
//}
//
//template<typename A1, typename T1, typename... G1s, typename B1, typename T2, typename... G2s,
//	typename std::enable_if_t<(
//		((expr::factor_count_list<G2s, G1s...>::value >= expr::factor_count_list<G2s, G2s...>::value) && ...) &&
//		!((expr::factor_count_list<G2s, G1s...>::value == expr::factor_count_list<G2s, G2s...>::value) && ...)
//		), int> = 0>
//auto operator-(OpBinaryDiv<A1, OpNLVariable<T1, G1s...>> const& a, OpBinaryDiv<B1, OpNLVariable<T2, G2s...>> const& b)
//{
//	auto div = a.b / b.b;
//	return symphas::internal::terminate_div((a.a - b.a * div), expr::inverse(a.b));
//}
//
//template<typename A1, typename T1, typename... G1s, typename B1, typename T2, typename... G2s,
//	typename std::enable_if_t<((expr::factor_count_list<G1s, G2s...>::value == expr::factor_count_list<G1s, G1s...>::value) && ...), int> = 0>
//auto operator-(OpBinaryDiv<A1, OpNLVariable<T1, G1s...>> const& a, OpBinaryDiv<B1, OpNLVariable<T2, G2s...>> const& b)
//{
//	auto div = a.b / b.b;
//	return symphas::internal::terminate_div((a.a - b.a * div), expr::inverse(a.b));
//}

// ******************************************************************************************
//Rules for simplification of term expressions.
// ******************************************************************************************


namespace expr
{
	template<typename...>
	struct can_combine;

	template<typename... G1s, expr::exp_key_t... X1s>
	struct can_combine<Term<G1s, X1s>...>
	{
		template<typename...>
		struct with_check;

		template<typename... G2s, expr::exp_key_t... X2s>
		struct with_check<Term<G2s, X2s>...>
		{
			static const bool value =
				((expr::is_combinable<G2s> &&
					symphas::lib::index_of_type<Term<G2s, X2s>, Term<G1s, X1s>...> >= 0) && ...) &&
				(sizeof...(G1s) == sizeof...(G2s));
		};

		template<typename... Ts>
		static constexpr bool with = with_check<Ts...>::value;
	};
}


//! Addition of two multi variables with data that can be combined.
template<typename V1, typename... G1s, expr::exp_key_t... X1s, typename V2, typename... G2s, expr::exp_key_t... X2s,
	typename std::enable_if<expr::can_combine<Term<G1s, X1s>...>::template with<Term<G2s, X2s>...>, int>::type = 0>
auto operator+(OpTerms<V1, Term<G1s, X1s>...> const& a, OpTerms<V2, Term<G2s, X2s>...> const& b)
{
	return (a.term + b.term) * OpTerms(OpIdentity{}, symphas::internal::cast_term<OpTerms<Term<G1s, X1s>...>>::cast(a));
}

//! Addition of two multi variables with data that can be combined.
template<typename V1, typename... G1s, expr::exp_key_t... X1s, typename V2, typename... G2s, expr::exp_key_t... X2s,
	typename std::enable_if<expr::can_combine<Term<G1s, X1s>...>::template with<Term<G2s, X2s>...>, int>::type = 0>
auto operator-(OpTerms<V1, Term<G1s, X1s>...> const& a, OpTerms<V2, Term<G2s, X2s>...> const& b)
{
	return (a.term - b.term) * OpTerms(OpIdentity{}, symphas::internal::cast_term<OpTerms<Term<G1s, X1s>...>>::cast(a));
}

// ******************************************************************************************



template<typename A1, typename A2, typename E,
	typename std::enable_if_t<(expr::eval_type<A1>::rank == 0 && expr::eval_type<E>::rank == 0), int> = 0>
auto operator+(OpBinaryDiv<A1, A2> const& a, OpExpression<E> const& b)
{
	return (a.a + (*static_cast<E const*>(&b)) * a.b) / a.b;
}

template<typename A1, typename A2, typename E,
	typename std::enable_if_t<(expr::eval_type<A1>::rank == 0 && expr::eval_type<E>::rank == 0), int> = 0>
auto operator-(OpBinaryDiv<A1, A2> const& a, OpExpression<E> const& b)
{
	return (a.a - (*static_cast<E const*>(&b)) * a.b) / a.b;
}

template<typename E, typename B1, typename B2,
	typename std::enable_if_t<(expr::eval_type<E>::rank == 0 && expr::eval_type<B1>::rank == 0), int> = 0>
auto operator+(OpExpression<E> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return ((*static_cast<E const*>(&a)) * b.b + b.a) / b.b;
}

//template<typename A1, typename A2, typename A3, typename... Bs>
//auto operator-(OpBinaryDiv<OpBinaryDiv<A1, A2>, A3> const& a, OpAdd<Bs...> const& b)
//{
//	return expr::terms_after_first(a);
//}

template<typename A1, typename A2, typename... Bs,
	typename std::enable_if_t<(expr::eval_type<A1>::rank == 0 && expr::eval_type<OpAdd<Bs...>>::rank == 0), int> = 0>
auto operator+(OpBinaryDiv<A1, A2> const& a, OpAdd<Bs...> const& b)
{
	return (a.a + b * a.b) / a.b;
}

template<typename A1, typename A2, typename... Bs,
	typename std::enable_if_t<(expr::eval_type<A1>::rank == 0 && expr::eval_type<OpAdd<Bs...>>::rank == 0), int> = 0>
auto operator-(OpBinaryDiv<A1, A2> const& a, OpAdd<Bs...> const& b)
{
	return (a.a - b * a.b) / a.b;
}

template<typename... As, typename B1, typename B2,
	typename std::enable_if_t<(expr::eval_type<OpAdd<As...>>::rank == 0 && expr::eval_type<B1>::rank == 0), int> = 0>
auto operator+(OpAdd<As...> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return (a * b.b + b.a) / b.b;
}


template<typename... As, typename B1, typename B2,
	typename std::enable_if_t<(expr::eval_type<OpAdd<As...>>::rank == 0 && expr::eval_type<B1>::rank == 0), int> = 0>
auto operator-(OpAdd<As...> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return (a * b.b - b.a) / b.b;
}


template<typename A1, typename A2, typename B1, typename B2,
	typename std::enable_if_t<(expr::eval_type<A1>::rank == 0 && expr::eval_type<B1>::rank == 0), int> = 0>
auto operator-(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return (a.a * b.b - b.a * a.b) / (a.b * b.b);
}

template<typename A1, typename A2, typename B1, typename B2,
	typename std::enable_if_t<(expr::eval_type<A1>::rank == 0 && expr::eval_type<B1>::rank == 0), int> = 0>
auto operator+(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return (a.a * b.b + b.a * a.b) / (a.b * b.b);
}

template<typename A1, typename A2, typename B1, typename B2,
	typename std::enable_if_t<(expr::eval_type<A1>::rank > 0 || expr::eval_type<B1>::rank > 0), int> = 0>
	auto operator-(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return expr::make_add(a, -b);
}

template<typename A1, typename A2, typename B1, typename B2,
	typename std::enable_if_t<(expr::eval_type<A1>::rank > 0 || expr::eval_type<B1>::rank > 0), int> = 0>
	auto operator+(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return expr::make_add(a, b);
}



// two more rules are in expressiontransforms.h, when there are unbalanced factors between each denominator.

// ******************************************************************************************


template<typename E0, typename... Es>
auto OpAdd<E0, Es...>::operator-() const
{
	return expr::neg_terms(*this, std::make_index_sequence<sizeof...(Es) + 1>{});
}


template<typename E1, typename E2>
auto OpBinaryMul<E1, E2>::operator-() const
{
	return expr::make_mul(-a, b);
}

template<typename E1, typename E2>
auto OpBinaryDiv<E1, E2>::operator-() const
{
	return expr::make_div(-a, b);
}

template<typename E2>
auto operator+(scalar_t a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) + *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator+(OpExpression<E1> const& a, scalar_t b)
{
	return *static_cast<const E1*>(&a) + expr::make_literal(b);
}

template<typename E2>
auto operator-(scalar_t a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) - *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator-(OpExpression<E1> const& a, scalar_t b)
{
	return *static_cast<const E1*>(&a) - expr::make_literal(b);
}

template<typename E2>
auto operator+(int a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) + *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator+(OpExpression<E1> const& a, int b)
{
	return *static_cast<const E1*>(&a) + expr::make_literal(b);
}

template<typename E2>
auto operator-(int a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) - *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator-(OpExpression<E1> const& a, int b)
{
	return *static_cast<const E1*>(&a) - expr::make_literal(b);
}

template<typename E2>
auto operator+(scalar_t a, OpOperator<E2> const& b)
{
	return expr::make_literal(a) + *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator+(OpOperator<E1> const& a, scalar_t b)
{
	return *static_cast<const E1*>(&a) + expr::make_literal(b);
}

template<typename E2>
auto operator-(scalar_t a, OpOperator<E2> const& b)
{
	return expr::make_literal(a) - *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator-(OpOperator<E1> const& a, scalar_t b)
{
	return *static_cast<const E1*>(&a) - expr::make_literal(b);
}

template<typename E2>
auto operator+(int a, OpOperator<E2> const& b)
{
	return expr::make_literal(a) + *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator+(OpOperator<E1> const& a, int b)
{
	return *static_cast<const E1*>(&a) + expr::make_literal(b);
}

template<typename E2>
auto operator-(int a, OpOperator<E2> const& b)
{
	return expr::make_literal(a) - *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator-(OpOperator<E1> const& a, int b)
{
	return *static_cast<const E1*>(&a) - expr::make_literal(b);
}



template<typename T, size_t D, typename E2>
auto operator+(any_vector_t<T, D> const& a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) + *static_cast<const E2*>(&b);
}

template<typename T, size_t D, typename E1>
auto operator+(OpExpression<E1> const& a, any_vector_t<T, D> const& b)
{
	return *static_cast<const E1*>(&a) + expr::make_literal(b);
}


template<typename T, size_t D, typename E2>
auto operator-(any_vector_t<T, D> const& a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) - *static_cast<const E2*>(&b);
}


template<typename T, size_t D, typename E1>
auto operator-(OpExpression<E1> const& a, any_vector_t<T, D> const& b)
{
	return *static_cast<const E1*>(&a) - expr::make_literal(b);
}

template<typename T, size_t D, typename E2>
auto operator+(any_vector_t<T, D> const& a, OpOperator<E2> const& b)
{
	return expr::make_literal(a) + *static_cast<const E2*>(&b);
}

template<typename T, size_t D, typename E1>
auto operator+(OpOperator<E1> const& a, any_vector_t<T, D> const& b)
{
	return *static_cast<const E1*>(&a) + expr::make_literal(b);
}


template<typename T, size_t D, typename E2>
auto operator-(any_vector_t<T, D> const& a, OpOperator<E2> const& b)
{
	return expr::make_literal(a) - *static_cast<const E2*>(&b);
}


template<typename T, size_t D, typename E1>
auto operator-(OpOperator<E1> const& a, any_vector_t<T, D> const& b)
{
	return *static_cast<const E1*>(&a) - expr::make_literal(b);
}



template<typename E2>
auto operator*(scalar_t a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) * *static_cast<const E2*>(&b);
}

template<typename E1, typename = std::enable_if_t<!(expr::is_identity<E1> || expr::is_fraction<E1>), int>>
auto operator*(OpExpression<E1> const& a, scalar_t b)
{
	return b * *static_cast<const E1*>(&a);
}

template<typename E1>
auto operator/(OpExpression<E1> const& a, scalar_t b)
{
	return *static_cast<const E1*>(&a) / expr::make_literal(b);
}

template<typename E2>
auto operator/(scalar_t a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) / *static_cast<const E2*>(&b);
}

template<typename E2>
auto operator*(int a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) * *static_cast<const E2*>(&b);
}

template<typename E1, typename = std::enable_if_t<!(expr::is_identity<E1> || expr::is_fraction<E1>), int>>
auto operator*(OpExpression<E1> const& a, int b)
{
	return b * *static_cast<const E1*>(&a);
}

template<typename E2>
auto operator/(int a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) / *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator/(OpExpression<E1> const& a, int b)
{
	return (*static_cast<const E1*>(&a)) / expr::make_literal(b);
}


template<typename T, size_t D, typename E2>
auto operator*(any_vector_t<T, D> const& a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) * *static_cast<const E2*>(&b);
}

template<typename T, size_t D, typename E1, typename = std::enable_if_t<!(expr::is_identity<E1> || expr::is_fraction<E1>), int>>
auto operator*(OpExpression<E1> const& a, any_vector_t<T, D> const& b)
{
	return *static_cast<const E1*>(&a) * expr::make_literal(b);
}

template<typename E2>
auto operator*(scalar_t a, OpOperator<E2> const& b)
{
	return expr::make_literal(a) * *static_cast<const E2*>(&b);
}

template<typename E1, typename = std::enable_if_t<!(expr::is_identity<E1> || expr::is_fraction<E1>), int>>
auto operator*(OpOperator<E1> const& a, scalar_t b)
{
	return b * *static_cast<const E1*>(&a);
}

template<typename E1>
auto operator/(OpOperator<E1> const& a, scalar_t b)
{
	return *static_cast<const E1*>(&a) / expr::make_literal(b);
}

template<typename E2>
auto operator/(scalar_t a, OpOperator<E2> const& b)
{
	return expr::make_literal(a) / *static_cast<const E2*>(&b);
}

template<typename E2>
auto operator*(int a, OpOperator<E2> const& b)
{
	return expr::make_literal(a) * *static_cast<const E2*>(&b);
}

template<typename E1, typename = std::enable_if_t<!(expr::is_identity<E1> || expr::is_fraction<E1>), int>>
auto operator*(OpOperator<E1> const& a, int b)
{
	return b * *static_cast<const E1*>(&a);
}

template<typename E2>
auto operator/(int a, OpOperator<E2> const& b)
{
	return expr::make_literal(a) / *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator/(OpOperator<E1> const& a, int b)
{
	return (*static_cast<const E1*>(&a)) / expr::make_literal(b);
}

template<typename T, size_t D, typename E2>
auto operator*(any_vector_t<T, D> const& a, OpOperator<E2> const& b)
{
	return expr::make_literal(a) * *static_cast<const E2*>(&b);
}

template<typename T, size_t D, typename E1, typename = std::enable_if_t<!(expr::is_identity<E1> || expr::is_fraction<E1>), int>>
auto operator*(OpOperator<E1> const& a, any_vector_t<T, D> const& b)
{
	return *static_cast<const E1*>(&a) * expr::make_literal(b);
}






/*
 *
 * Tensor rules.
 *
 ******************************************************************************/



 //! Ensure a value constant is always multiplied on the left by anything else.
template<typename coeff_t, typename T, size_t... Ns, 
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t>), int> = 0>
auto operator*(coeff_t const& a, OpTensor<T, Ns...> const& tensor)
{
	return expr::make_tensor<Ns...>(a * symphas::internal::tensor_cast::cast(tensor));
}

//! Ensure a value constant is always multiplied on the left by anything else.
template<typename coeff_t, typename T, size_t... Ns,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t>), int> = 0>
auto operator*(OpTensor<T, Ns...> const& tensor, coeff_t const& a)
{
	return a * tensor;
}

 //! Ensure a value constant is always multiplied on the left by anything else.
template<typename E, typename T>
auto operator*(E const& a, OpLiteral<T> const& b)
{
	return b * a;
}

//! Ensure a value constant is always multiplied on the left by anything else.
template<typename E, typename T, size_t... Ns, typename std::enable_if_t<expr::has_coeff<E>, int> = 0>
auto operator*(OpExpression<E> const& e, OpTensor<T, Ns...> const& tensor)
{
	if constexpr (expr::eval_type<decltype(expr::coeff(std::declval<E>()))>::rank > 0)
	{
		return (expr::coeff(*static_cast<E const*>(&e)) * tensor) * (symphas::internal::tensor_cancel{} * (*static_cast<E const*>(&e)));
	}
	else
	{
		return tensor * (*static_cast<E const*>(&e));
	}
}

//! Ensure a value constant is always multiplied on the left by anything else.
template<typename E, typename T, size_t... Ns, typename std::enable_if_t<expr::has_coeff<E>, int> = 0>
auto operator*(OpOperator<E> const& e, OpTensor<T, Ns...> const& tensor)
{
	if constexpr (expr::eval_type<decltype(expr::coeff(std::declval<E>()))>::rank > 0)
	{
		return (expr::coeff(*static_cast<E const*>(&e)) * tensor) * (symphas::internal::tensor_cancel{} *(*static_cast<E const*>(&e)));
	}
	else
	{
		return tensor * (*static_cast<E const*>(&e));
	}
}

//! Ensure a value constant is always multiplied on the left by anything else.
template<typename A, typename B, typename T, size_t... Ns>
auto operator*(OpBinaryMul<A, B> const& e, OpTensor<T, Ns...> const& tensor)
{
	return e.a * (e.b * tensor);
}


//! Ensure a value constant is always multiplied on the left by anything else.
template<typename A, typename B, typename T, size_t... Ns>
auto operator*(OpTensor<T, Ns...> const& tensor, OpBinaryMul<A, B> const& e)
{
	return (tensor * e.a) * e.b;
}


template<typename T, size_t N0, size_t N1, size_t... Ns, typename S, size_t M0, size_t M1, size_t... Ms>
auto operator+(OpTensor<T, N0, N1, Ns...> const& a, OpTensor<S, M0, M1, Ms...> const& b)
{
	return expr::make_add(a, b);
}

template<typename T, size_t... Ns, typename coeff_t,
	typename std::enable_if_t<expr::is_tensor<coeff_t>, int> = 0>
auto operator+(OpTensor<T, Ns...> const& a, coeff_t const& b)
{
	return expr::make_add(a, b);
}

template<typename T, size_t... Ns, typename coeff_t,
	typename std::enable_if_t<expr::is_tensor<coeff_t>, int> = 0>
auto operator+(coeff_t const& a, OpTensor<T, Ns...> const& b)
{
	return expr::make_add(a, b);
}

template<typename T, size_t... Ns, typename coeff_t,
	typename std::enable_if_t<expr::is_tensor<coeff_t>, int> = 0>
auto operator-(OpTensor<T, Ns...> const& a, coeff_t const& b)
{
	return a + (-b);
}

template<typename T, size_t... Ns, typename coeff_t,
	typename std::enable_if_t<expr::is_tensor<coeff_t>, int> = 0>
auto operator-(coeff_t const& a, OpTensor<T, Ns...> const& b)
{
	return a + (-b);
}





template<typename T, typename coeff_t, typename = std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t>), int>>
auto operator+(OpTensor<T, 0, 1> const& a, coeff_t const& b)
{
	return symphas::internal::tensor_cast::cast(a) + b;
}

template<typename T, typename coeff_t, typename = std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t>), int>>
auto operator+(OpTensor<T, 0, 0, 1, 1> const& a, coeff_t const& b)
{
	return symphas::internal::tensor_cast::cast(a) + b;
}

template<typename T, typename coeff_t, typename = std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t>), int>>
auto operator-(OpTensor<T, 0, 1> const& a, coeff_t const& b)
{
	return symphas::internal::tensor_cast::cast(a) - b;
}

template<typename T, typename coeff_t, typename = std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t>), int>>
auto operator-(OpTensor<T, 0, 0, 1, 1> const& a, coeff_t const& b)
{
	return symphas::internal::tensor_cast::cast(a) - b;
}


template<typename T, typename coeff_t, typename = std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t>), int>>
auto operator+(coeff_t const& a, OpTensor<T, 0, 1> const& b)
{
	return a + symphas::internal::tensor_cast::cast(b);
}

template<typename T, typename coeff_t, typename = std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t>), int>>
auto operator+(coeff_t const& a, OpTensor<T, 0, 0, 1, 1> const& b)
{
	return a + symphas::internal::tensor_cast::cast(b);
}

template<typename T, typename coeff_t, typename = std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t>), int>>
auto operator-(coeff_t const& a, OpTensor<T, 0, 1> const& b)
{
	return a - symphas::internal::tensor_cast::cast(b);
}

template<typename T, typename coeff_t, typename = std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<coeff_t>), int>>
auto operator-(coeff_t const& a, OpTensor<T, 0, 0, 1, 1> const& b)
{
	return a - symphas::internal::tensor_cast::cast(b);
}

template<typename T, size_t... Ns>
auto operator*(OpTensor<T, Ns...> const&, expr::symbols::Symbol)
{
	return OpTensor<expr::symbols::Symbol, Ns...>{};
}

template<typename T, size_t... Ns>
auto operator*(expr::symbols::Symbol, OpTensor<T, Ns...> const&)
{
	return OpTensor<expr::symbols::Symbol, Ns...>{};
}


template<>
struct expr::eval_type<int> : expr::eval_type<OpLiteral<int>> {};
template<>
struct expr::eval_type<scalar_t> : expr::eval_type<OpLiteral<scalar_t>> {};
template<>
struct expr::eval_type<complex_t> : expr::eval_type<OpLiteral<complex_t>> {};
template<typename T, size_t D>
struct expr::eval_type<any_vector_t<T, D>> : expr::eval_type<OpTensor<T, 0, D>> {};


