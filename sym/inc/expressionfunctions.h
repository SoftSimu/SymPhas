
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
#include "expressionaggregates.h"
#include "expressiontransforms.h"

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

		template<typename V, typename T>
		static auto get(V v, OpLiteral<T> const& e);

		template<typename V>
		static auto get(V v, OpIdentity);

		template<typename V>
		static auto get(V v, OpNegIdentity);

		template<typename V>
		static auto get(V v, OpVoid);

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

#ifdef PRINTABLE_EQUATIONS

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

#endif

	E e;			//!< The expression to which the function is applied.
	V value;		//!< Coefficient of the function expression term.
};

template<auto f, typename V, typename E>
OpFuncApply(V, E)->OpFuncApply<f, V, E>;


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
	template<typename V, typename T>
	inline auto make_function<f>::get(V v, OpLiteral<T> const& e)
	{
		return expr::make_literal(OpFuncApply<f, V, OpLiteral<T>>(v, e).eval(0));
	}

	template<auto f>
	template<typename V>
	inline auto make_function<f>::get(V v, OpIdentity)
	{
		return expr::make_literal(OpFuncApply<f, V, OpIdentity>(v, OpIdentity{}).eval(0));
	}

	template<auto f>
	template<typename V>
	inline auto make_function<f>::get(V v, OpNegIdentity)
	{
		return expr::make_literal(OpFuncApply<f, V, OpNegIdentity>(v, OpNegIdentity{}).eval(0));
	}

	template<auto f>
	template<typename V>
	inline auto make_function<f>::get(V v, OpVoid)
	{
		return expr::make_literal(OpFuncApply<f, V, OpVoid>(v, OpVoid{}).eval(0));
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
auto constexpr func_Conjugate = &symphas::math::conj<typename expr::eval_type<E>::type>;


//! Alias to construct a modulus function to compute the modulus.
template<typename E>
auto constexpr func_Modulus = &symphas::math::modulus<typename expr::eval_type<E>::type>;

//! Alias to construct the function to return real part.
template<typename E>
auto constexpr func_RealPart = &symphas::math::real<typename expr::eval_type<E>::type>;

//! Alias to construct the function to return imaginary part.
template<typename E>
auto constexpr func_ImagPart = &symphas::math::imag<typename expr::eval_type<E>::type>;


//! Alias to construct the function to compute the cosine function.
template<typename E>
auto constexpr func_Cos = &symphas::math::cos<typename expr::eval_type<E>::type>;

//! Alias to construct the function to compute the sine function.
template<typename E>
auto constexpr func_Sin = &symphas::math::sin<typename expr::eval_type<E>::type>;

//! Alias to construct the function to compute the sine function.
template<typename E>
auto constexpr func_Tan = &symphas::math::tan<typename expr::eval_type<E>::type>;

//! Alias to construct the function to compute the cosine function.
template<typename E>
auto constexpr func_Cosh = &symphas::math::cosh<typename expr::eval_type<E>::type>;

//! Alias to construct the function to compute the sine function.
template<typename E>
auto constexpr func_Sinh = &symphas::math::sinh<typename expr::eval_type<E>::type>;

//! Alias to construct the function to compute the sine function.
template<typename E>
auto constexpr func_Tanh = &symphas::math::tanh<typename expr::eval_type<E>::type>;

//! Alias to construct the function to compute the cosine function.
template<typename E>
auto constexpr func_Acos = &symphas::math::acos<typename expr::eval_type<E>::type>;

//! Alias to construct the function to compute the sine function.
template<typename E>
auto constexpr func_Asin = &symphas::math::asin<typename expr::eval_type<E>::type>;

//! Alias to construct the function to compute the sine function.
template<typename E>
auto constexpr func_Atan = &symphas::math::atan<typename expr::eval_type<E>::type>;

//! Alias to construct the function to compute the cosine function.
template<typename E>
auto constexpr func_Acosh = &symphas::math::acosh<typename expr::eval_type<E>::type>;

//! Alias to construct the function to compute the sine function.
template<typename E>
auto constexpr func_Asinh = &symphas::math::asinh<typename expr::eval_type<E>::type>;

//! Alias to construct the function to compute the sine function.
template<typename E>
auto constexpr func_Atanh = &symphas::math::atanh<typename expr::eval_type<E>::type>;


//! Alias to construct the function to compute the sine function.
template<typename E>
auto constexpr func_Sqrt = &symphas::math::sqrt<typename expr::eval_type<E>::type>;






namespace expr
{
	template<typename E>
	auto conj(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Conjugate<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto modulus(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Modulus<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto real(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_RealPart<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto imag(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_ImagPart<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto cos(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Cos<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto sin(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Sin<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto tan(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Tan<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto cosh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Cosh<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto sinh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Sinh<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto tanh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Tanh<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto acos(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Acos<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto asin(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Asin<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto atan(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Atan<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto acosh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Acosh<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto asinh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Asinh<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto atanh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Atanh<E>>::get(*static_cast<E const*>(&e));
	}


	template<typename E>
	auto sqrt(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_Sqrt<E>>::get(*static_cast<E const*>(&e));
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

#ifdef PRINTABLE_EQUATIONS

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

#endif

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


// *************************************************************************************************




namespace expr
{
	template<size_t N>
	using subvar = OpLVariable<OpIdentity, Variable<N>>;

	template<size_t N>
	subvar<N> _v = subvar<N>{ OpVoid{} };

	inline subvar<0> _x = _v<0>;
	inline subvar<1> _y = _v<1>;
	inline subvar<2> _z = _v<2>;
}


//! A function into which other expressions can be substituted when called.
/*!
 * A function into which other expressions can be substituted when called.
 * The function is created with a given number parameters, which may be zero.
 * 
 * \tparam E The type of the expression function.
 * \tparam ArgsNs... The indices of the independent variables, subvar<N>.
 */
template<typename E, size_t... ArgNs>
struct OpFuncSubstitutable
{
	OpFuncSubstitutable(E const& e) : e{ e } {}

protected:

	template<size_t... Ss>
	inline auto swap_x()
	{
		return e;
	}

	template<size_t S, size_t... Ss, typename X, typename... Xs>
	auto swap_x(X const& x, Xs const& ...xs)
	{
		return expr::transform::swap_grid<Variable<S>>(swap_x<Ss...>(xs...), x);
	}

public:

	E operator()() const
	{
		return e;
	}

	template<typename... Xs, std::enable_if_t<(sizeof...(Xs) == sizeof...(ArgNs)), int> = 0>
	auto operator()(OpExpression<Xs> const& ...xs)
	{
		return swap_x<ArgNs...>(*static_cast<const Xs*>(&xs)...);
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

	E e;				//!< The substitutable function.
};


namespace expr
{

	namespace
	{

		template<typename G>
		struct tuple_list_to_ptrs
		{
			using type = G;
		};

		template<typename... Gs>
		struct tuple_list_to_ptrs<std::tuple<Gs...>>
		{
			using type = std::tuple<Gs*...>;
		};


		template<size_t N>
		constexpr std::true_type _is_subvar(Variable<N>)
		{
			return {};
		}

		constexpr std::false_type _is_subvar(...)
		{
			return {};
		}

		template<typename T>
		constexpr auto is_subvar(T t)
		{
			return _is_subvar(t);
		}

		template<size_t N>
		constexpr std::index_sequence<N> subvar_id(Variable<N>*)
		{
			return {};
		}


		constexpr auto get_subvar_list(std::tuple<>)
		{
			return std::index_sequence<>{};
		}

		template<typename G0, typename... Gs>
		constexpr auto get_subvar_list(std::tuple<G0*, Gs*...>)
		{
			if constexpr (std::invoke_result_t<decltype(is_subvar<G0>), G0>::value)
			{
				return symphas::lib::fold_unique_ids(symphas::lib::seq_join(subvar_id((G0*)nullptr), get_subvar_list(std::tuple<Gs*...>{})));
			}
			else
			{
				return get_subvar_list(std::tuple<Gs*...>{});
			}
		}
	}

	template<typename E>
	auto get_independent_variables(OpExpression<E> const& e)
	{
		get_subvar_list(typename tuple_list_to_ptrs<typename expr::op_types<E>::type>::type{});
	}

	template<typename E, size_t... ArgNs>
	auto get_independent_variables(OpFuncSubstitutable<E, ArgNs...> const& e)
	{
		return std::index_sequence<ArgNs...>{};
	}

	namespace
	{


		template<typename E, size_t... ArgNs>
		OpFuncSubstitutable<E, ArgNs...> make_new_function(OpExpression<E> const& e, std::index_sequence<ArgNs...>)
		{
			return { *static_cast<const E*>(&e) };
		}


		//! Create a substitutable function.
		template<typename E>
		auto make_substitutable_function(OpExpression<E> const& e)
		{
			auto all_ids = get_subvar_list(typename tuple_list_to_ptrs<typename expr::op_types<E>::type>::type{});
			return make_new_function(*static_cast<const E*>(&e), all_ids);
		}

		//! Create a substitutable function.
		template<typename E, size_t... Ns>
		OpFuncSubstitutable<E, Ns...> make_substitutable_function(OpExpression<E> const& e, std::index_sequence<Ns...> ns)
		{
			return make_new_function(*static_cast<const E*>(&e), ns);
		}


		template<size_t... Ns>
		struct OpSymbolicFunction;

		template<size_t N0, size_t... Ns>
		struct OpSymbolicFunction<N0, Ns...>
		{

		protected:

			template<typename S>
			struct convert_vars;

			template<int... Ms>
			struct convert_vars<std::integer_sequence<int, Ms...>>
			{
				using type = std::index_sequence<size_t(Ms)...>;
			};

			template<typename S>
			struct has_enough_vars;

			template<size_t... Ms>
			struct has_enough_vars<std::index_sequence<Ms...>>
			{
				using check_seq = symphas::lib::seq_sub_t<
					decltype(symphas::lib::sort_ids(std::declval<std::index_sequence<Ms...>>())),
					decltype(symphas::lib::sort_ids(std::declval<std::index_sequence<N0, Ns...>>()))>;

				static const bool value = std::is_same<
					std::index_sequence<0>,
					decltype(symphas::lib::fold_unique_ids(std::declval<typename convert_vars<check_seq>::type>()))>::value;
			};

		public:

			template<typename E, 
				std::enable_if_t<has_enough_vars<decltype(
					get_subvar_list(std::declval<typename tuple_list_to_ptrs<typename expr::op_types<E>::type>::type>())
					)>::value, int> = 0>
			auto operator=(OpExpression<E> const& e)
			{
				auto all_ids = get_subvar_list(typename tuple_list_to_ptrs<typename expr::op_types<E>::type>::type{});
				return make_substitutable_function(*static_cast<const E*>(&e));
			}

			OpSymbolicFunction(OpSymbolicFunction<N0, Ns...> const&) = delete;
			OpSymbolicFunction(OpSymbolicFunction<N0, Ns...>&&) = delete;
		};

		template<>
		struct OpSymbolicFunction<>
		{
			template<typename E>
			auto operator=(OpExpression<E> const& e)
			{
				return make_substitutable_function(*static_cast<const E*>(&e));
			}

			OpSymbolicFunction(OpSymbolicFunction<> const&) = delete;
			OpSymbolicFunction(OpSymbolicFunction<>&&) = delete;
		};


		template<typename... Ts>
		struct OpSymbolicFunctionSwap;

		template<typename... Ts, size_t... Ns, typename... Gs>
		struct OpSymbolicFunctionSwap<OpLVariable<Ts, Variable<Ns, Gs>>...>
		{

		protected:

			template<size_t N, typename E>
			inline auto swap_all_grids(OpExpression<E> const& e)
			{
				return *static_cast<E const*>(&e);
			}

			template<size_t N, typename E, typename X, typename... Xs>
			auto swap_all_grids(OpExpression<E> const& e)
			{
				return expr::transform::swap_grid<X>(swap_all_grids<N, E, Xs...>(*static_cast<E const*>(&e)), _v<N - sizeof...(Xs) - 1>);
			}

		public:

			template<typename E>
			auto operator=(OpExpression<E> const& e)
			{
				return make_substitutable_function(
					swap_all_grids<sizeof...(Ts), E, OpLVariable<Ts, Variable<Ns, Gs>>...>(*static_cast<const E*>(&e)));
			}

			OpSymbolicFunctionSwap(OpSymbolicFunctionSwap<OpLVariable<Ts, Variable<Ns, Gs>>...> const&) = delete;
			OpSymbolicFunctionSwap(OpSymbolicFunctionSwap<OpLVariable<Ts, Variable<Ns, Gs>>...>&&) = delete;
		};

	}

	template<size_t... Ns>
	OpSymbolicFunction<Ns...> fn(subvar<Ns>...) { return {}; }

	OpSymbolicFunction<> fn() { return {}; }

	template<typename... Ts, size_t... Ns, typename... Gs>
	OpSymbolicFunctionSwap<OpLVariable<Ts, Variable<Ns, Gs>>...> fn(OpLVariable<Ts, Variable<Ns, Gs>>...) { return {}; }

	


}






