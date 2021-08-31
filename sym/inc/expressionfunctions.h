
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

#include "expressionsprint.h"

namespace symphas::internal
{
	//! Constructs a function of the given expression
	/*!
	 * Create a function of an expression.
	 * 
	 * \param f The function which is applied.
	 */
	template<auto f>
	struct make_function
	{
		template<typename A>
		static auto get(A&& a);

		template<typename V, typename E>
		static auto get(V v, OpExpression<E> const& e);

		template<typename V, typename T, typename G>
		static auto get(V v, OpLVariable<T, G> const& e);

		template<typename V, typename G>
		static auto get_g(V v, G g);

	};
}



namespace expr
{
	//! Create a function expression of another expression.
	/*!
	 * The function is applied to the given expression through an expression
	 * that represents a function.
	 * 
	 * \param a The expression to which the function is applied.
	 */
	template<auto f, typename A>
	auto make_function(A&& a)
	{
		return symphas::internal::make_function<f>::template get(std::forward<A>(a));
	}

	//! Create a function expression of another expression.
	/*!
	 * The function is applied to the given expression through an expression
	 * that represents a function.
	 *
	 * \param v The coefficient to the function expression.
	 * \param a The expression to which the function is applied.
	 */
	template<auto f, typename V, typename A>
	auto make_function(V&& v, A&& a)
	{
		return symphas::internal::make_function<f>::template get(std::forward<V>(v), std::forward<A>(a));
	}
}




// **********************************************************************************


//! Applies a particular function to an expression.
/*!
 * The result of an expression is evaluated and then a chosen function is 
 * applied.
 * 
 * \tparam f The function applied to the result of the expression.
 * \tparam V The coefficient type of this expression object.
 * \tparam E The type of the expression the function applies to.
 */
template<auto f, typename V, typename E>
struct OpFuncApply : OpExpression<OpFuncApply<f, V, E>>
{
	//! Create a new function expression, applied to an expression.
	/*!
	 * Creates a new function expression with the given coefficient,
	 * applied to the given expression.
	 * 
	 * \param value The coefficient of the function expression.
	 * \param e The expression on which the function is applied.
	 */
	OpFuncApply(V value, E const& e) : value{ value }, e{ e } {}

	//! Create a new function expression, applied to an expression.
	/*!
	 * Creates a new function expression 
	 * applied to the given expression.
	 *
	 * \param e The expression on which the function is applied.
	 */
	OpFuncApply(E const& e) : value{ value }, e { e } {}

	inline auto eval(iter_type n) const
	{
		return value * f(e.eval(n));
	}


	auto operator-() const
	{
		return symphas::internal::make_function<f>::template get(-value, e);
	}

	size_t print(FILE* out) const
	{
		char* str = new char[e.print_length() + 1];
		e.print(str);

		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_f_op<E>::template print<f>(out, str);

		delete[] str;
		return n;
	}

	size_t print(char* out) const
	{
		char* str = new char[e.print_length() + 1];
		e.print(str);

		size_t n = expr::print_with_coeff(out, value);
		n += symphas::internal::print_f_op<E>::template print<f>(out + n, str);

		delete[] str;
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + e.print_length()
			+ symphas::internal::print_f_op<E>::template print_length<f>();
	}

	E e;			//!< The expression to which the function is applied.
	V value;		//!< Coefficient of the function expression term.
};


template<typename S1, auto f2, typename V2, typename E2>
auto operator*(OpLiteral<S1> const& a, OpFuncApply<f2, V2, E2> const& b)
{
	return symphas::internal::make_function<f2>::template get(a.value * b.value, b.e);
}



namespace symphas::internal
{

	template<auto f>
	template<typename A>
	inline auto make_function<f>::get(A&& a)
	{
		return get(OpIdentity{}, std::forward<A>(a));
	}

	template<auto f>
	template<typename V, typename E>
	inline auto make_function<f>::get(V v, OpExpression<E> const& e)
	{
		return OpFuncApply<f, V, E>(v, *static_cast<E const*>(&e));
	}

	template<auto f>
	template<typename V, typename T, typename G>
	inline auto make_function<f>::get(V v, OpLVariable<T, G> const& e)
	{
		return OpFuncApply<f, V, OpLVariable<T, G>>(v, e);
	}

	template<auto f>
	template<typename V, typename G>
	inline auto make_function<f>::get_g(V v, G g)
	{
		return OpFuncApply<f, V, OpLVariable<OpIdentity, G>>(v, g);
	}

}





//! Alias to construct a conjugate function to compute the conjugate.
template<typename E>
using OpConjugate = OpFuncApply<&symphas::math::conj<typename expr::eval_type<E>::type>, OpIdentity, E>;


//! Alias to construct a modulus function to compute the modulus.
template<typename E>
using OpModulus = OpFuncApply<&symphas::math::modulus<typename expr::eval_type<E>::type>, OpIdentity, E>;

//! Alias to construct the function to return real part.
template<typename E>
using OpRealPart = OpFuncApply<&symphas::math::real<typename expr::eval_type<E>::type>, OpIdentity, E>;

//! Alias to construct the function to return imaginary part.
template<typename E>
using OpImagPart = OpFuncApply<&symphas::math::imag<typename expr::eval_type<E>::type>, OpIdentity, E>;


//! Alias to construct the function to compute the cosine function.
template<typename E>
using OpCos = OpFuncApply<&symphas::math::cos<typename expr::eval_type<E>::type>, OpIdentity, E>;

//! Alias to construct the function to compute the sine function.
template<typename E>
using OpSin = OpFuncApply<&symphas::math::sin<typename expr::eval_type<E>::type>, OpIdentity, E>;







namespace expr
{
	template<typename E>
	OpConjugate<E> conj(OpExpression<E> const& e)
	{
		return OpConjugate<E>(OpIdentity{}, *static_cast<E const*>(&e));
	}

	template<typename E>
	OpModulus<E> modulus(OpExpression<E> const& e)
	{
		return OpModulus<E>(OpIdentity{}, *static_cast<E const*>(&e));
	}

	template<typename E>
	OpRealPart<E> real(OpExpression<E> const& e)
	{
		return OpRealPart<E>(OpIdentity{}, *static_cast<E const*>(&e));
	}

	template<typename E>
	OpImagPart<E> imag(OpExpression<E> const& e)
	{
		return OpImagPart<E>(OpIdentity{}, *static_cast<E const*>(&e));
	}

	template<typename E>
	OpCos<E> cos(OpExpression<E> const& e)
	{
		return OpCos<E>(OpIdentity{}, *static_cast<E const*>(&e));
	}

	template<typename E>
	OpSin<E> sin(OpExpression<E> const& e)
	{
		return OpSin<E>(OpIdentity{}, *static_cast<E const*>(&e));
	}

}



// ***************************************************************************************


//! Expression representing the evaluation of a function.
/*!
 * The expression does not define any data. Given a function, the expression
 * will be passed to the function or functor which is provided.
 *
 * The function needs to be compatible with at least two argument types: the
 * expression and the index. Other arguments can then be provided.
 *
 * \tparam V The coefficient term of the expression function.
 * \tparam E The type of the expression to which the function is applied.
 * \tparam F The function type or functor object type.
 * \tparam Arg0 The type of the first argument given to the functor.
 * \tparam Args... The type of the remaining arguments given to the functor.
 */
template<typename V, typename E, typename F, typename Arg0, typename... Args>
struct OpFunc : OpExpression<OpFunc<V, E, F, Arg0, Args...>>
{
	//! Create a function of an expression.
	/*!
	 * Creates a function expression of the given expression, by applying
	 * the instance of the object type `F` as a functor. The arguments
	 * to the functor are also provided.
	 *
	 * \param name The identifying name of the function.
	 * \param value The coefficient value.
	 * \param e The expression to which the functor is applied.
	 * \param f The function pointer or functor object instance.
	 * \param arg0 The first argument which is passed to the function.
	 * \param args... The remaining arguments which are passed to the function.
	 */
	OpFunc(std::string name, V value, E const& e, F f, Arg0&& arg0, Args&&... args) :
		name{ name }, value{ value }, e{ e }, f{ f }, tt{ std::forward<Arg0>(arg0), std::forward<Args>(args)... } {}

	//! Create a function of an expression.
	/*!
	 * Creates a function expression of the given expression, by applying
	 * the instance of the object type `F` as a functor. The arguments
	 * to the functor are also provided.
	 *
	 * \param name The identifying name of the function.
	 * \param value The coefficient value.
	 * \param e The expression to which the functor is applied.
	 * \param f The function pointer or functor object instance.
	 * \param tt The list of arguments to the function.
	 */
	OpFunc(std::string name, V value, E const& e, F f, std::tuple<Arg0, Args...> const& tt) :
		name{ name }, value{ value }, e{ e }, f{ f }, tt{ tt } {}

	//! Create a function of an expression.
	/*!
	 * Creates a function expression of the given expression, by applying
	 * the instance of the object type `F` as a functor. The arguments
	 * to the functor are also provided.
	 *
	 * \param value The coefficient value.
	 * \param e The expression to which the functor is applied.
	 * \param f The function pointer or functor object instance.
	 * \param arg0 The first argument which is passed to the function.
	 * \param args... The remaining arguments which are passed to the function.
	 */
	OpFunc(V value, E const& e, F f, Arg0&& arg0, Args&&... args) :
		name{ "?" }, value{ value }, e{ e }, f{ f }, tt{ std::forward<Arg0>(arg0), std::forward<Args>(args)... } {}

	//! Create a function of an expression.
	/*!
	 * Creates a function expression of the given expression, by applying
	 * the instance of the object type `F` as a functor. The arguments
	 * to the functor are also provided.
	 *
	 * \param value The coefficient value.
	 * \param e The expression to which the functor is applied.
	 * \param f The function pointer or functor object instance.
	 * \param tt The list of arguments to the function.
	 */
	OpFunc(V value, E const& e, F f, std::tuple<Arg0, Args...> const& tt) :
		name{ "?" }, value{ value }, e{ e }, f{ f }, tt{ tt } {}

	inline auto eval(iter_type n) const
	{
		return value * std::apply(f, std::tuple_cat(std::make_tuple(e, n), tt));
	}

	auto operator-() const
	{
		return OpFunc<decltype(-value), E, F, Arg0, Args...>(-value, e, f, tt);
	}

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += fprintf(out, SYEX_LAMBDA_FUNC_FMT_A, name.c_str());
		n += e.print(out);
		n += fprintf(out, SYEX_LAMBDA_FUNC_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += sprintf(out + n, SYEX_LAMBDA_FUNC_FMT_A, name.c_str());
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_LAMBDA_FUNC_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + e.print_length()
			+ name.length() + SYEX_LAMBDA_FUNC_FMT_LEN;
	}


	friend struct expr::compound_get;

	V value;

protected:
	E e;							//!< Expression to which a function is applied.
	F f;							//!< Functor which is applied to an expression.
	std::tuple<Arg0, Args...> tt;	//!< Arguments to the function.
	std::string name;				//!< Identifier of the function.
};

//! Expression representing the evaluation of a function.
/*!
 * Specialization where the function does not take any additional arguments.
 *
 * The expression does not define any data. Given a function, the expression
 * will be passed to the function or functor which is provided.
 *
 * The function needs to be compatible with at least two argument types: the
 * expression and the index. No other arguments can be provided.
 *
 * \tparam V The coefficient term of the expression function.
 * \tparam E The type of the expression to which the function is applied.
 * \tparam F The function type or functor object type.
 */
template<typename V, typename E, typename F>
struct OpFunc<V, E, F, void> : OpExpression<OpFunc<V, E, F, void>>
{
	//! Create a function of an expression.
	/*!
	 * Creates a function expression of the given expression, by applying
	 * the instance of the object type `F` as a functor.
	 *
	 * \param name The identifying name of the function.
	 * \param value The coefficient value.
	 * \param e The expression to which the functor is applied.
	 * \param f The function pointer or functor object instance.
	 */
	OpFunc(std::string name, V value, E const& e, F f) :
		name{ name }, value{ value }, e{ e }, f{ f } {}

	//! Create a function of an expression.
	/*!
	 * Creates a function expression of the given expression, by applying
	 * the instance of the object type `F` as a functor.
	 *
	 * \param name The identifying name of the function.
	 * \param value The coefficient value.
	 * \param e The expression to which the functor is applied.
	 * \param f The function pointer or functor object instance.
	 */
	OpFunc(V value, E const& e, F f) :
		name{ "?" }, value{ value }, e{ e }, f{ f } {}

	inline auto eval(iter_type n) const
	{
		return value * f(e, n);
	}

	auto operator-() const
	{
		return OpFunc<decltype(-value), E, F, void>(name, -value, e, f);
	}

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += fprintf(out, SYEX_LAMBDA_FUNC_FMT_A, name.c_str());
		n += e.print(out);
		n += fprintf(out, SYEX_LAMBDA_FUNC_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += sprintf(out + n, SYEX_LAMBDA_FUNC_FMT_A, name.c_str());
		n += e.print(out + n);
		n += sprintf(out + n, SYEX_LAMBDA_FUNC_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + e.print_length()
			+ name.length() + SYEX_LAMBDA_FUNC_FMT_LEN;
	}


	friend struct expr::compound_get;

	V value;

protected:

	E e;							//!< Expression to which a function is applied.
	F f;							//!< Functor which is applied to an expression.
	std::string name;				//!< Identifier of the function.
};


template<typename V, typename E, typename F>
OpFunc(std::string, V, E, F)->OpFunc<V, E, F, void>;
template<typename V, typename E, typename F>
OpFunc(V, E, F)->OpFunc<V, E, F, void>;



template<typename S1, typename V2, typename E2, typename F2, typename Arg, typename... Args>
auto operator*(OpLiteral<S1> const& a, OpFunc<V2, E2, F2, Arg, Args...> const& b)
{
	return OpFunc(b.name, a.value * b.value, b.e, b.f, b.tt);
}

template<typename S1, typename V2, typename E2, typename F2>
auto operator*(OpLiteral<S1> const& a, OpFunc<V2, E2, F2, void> const& b)
{
	return OpFunc(a.name, a.value * b.value, b.e, b.f);
}






