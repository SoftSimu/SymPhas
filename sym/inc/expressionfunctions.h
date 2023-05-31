
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

#include <tuple>

#include "expressionrules.h"

namespace expr
{
	namespace
	{
		template<size_t... Ns>
		constexpr auto multiply_fold_seq(std::index_sequence<Ns...>)
		{
			return (Ns * ... * 1);
		}
	}

	/*!
	 * Compute the factorial.
	 */
	template<size_t N, size_t M = 1>
	constexpr auto factorial()
	{
		return expr::make_integer < multiply_fold_seq(
			symphas::lib::seq_add_t<
			symphas::lib::seq_repeating_value_t<N - M, size_t, M + 1>,
			std::make_index_sequence<N - M>>{})
			> ();
	}

	template<size_t N, size_t K>
	constexpr auto choose()
	{
		if constexpr (K == N || K == 0)
		{
			return OpIdentity{};
		}
		else
		{
			using num =
				symphas::lib::seq_add_t<
				symphas::lib::seq_repeating_value_t<K, size_t, N - K + 1>,
				std::make_index_sequence<K>>;

			using den =
				symphas::lib::seq_add_t<
				symphas::lib::seq_repeating_value_t<K, size_t, 1>,
				std::make_index_sequence<K>>;

			return expr::make_fraction < multiply_fold_seq(num{}), multiply_fold_seq(den{}) > ();
		}
	}



	//! Computes the expression multiplied `N` times with itself.
	/*!
	 * Computes the expression multiplied `N` times with itself.
	 *
	 * \param e The expression to multiply with itself.
	 */
	template<size_t N, typename E, typename = std::enable_if_t<(N < MAX_EXPONENT), int>>
	auto pow(OpExpression<E> const& e);

	//! Computes the expression multiplied `N` times with itself.
	/*!
	 * Computes the expression multiplied `N` times with itself.
	 *
	 * \param e The expression to multiply with itself.
	 */
	template<size_t N, typename E, typename = std::enable_if_t<(N < MAX_EXPONENT), int>>
	auto pow(OpOperator<E> const& e);

	//! Computes the expression multiplied `N` times with itself.
	/*!
	 * Computes the expression multiplied `N` times with itself.
	 *
	 * \param e The expression to multiply with itself.
	 */
	template<size_t N, typename E0, typename... Es, typename = std::enable_if_t<(N > 1 && N < MAX_EXPONENT), int>>
	auto pow(OpAdd<E0, Es...> const& e);

	//! Computes the expression multiplied `N` times with itself.
	/*!
	 * Computes the expression multiplied `N` times with itself.
	 *
	 * \param e The expression to multiply with itself.
	 */
	template<size_t N, typename A, typename B, typename = std::enable_if_t<(N > 1 && N < MAX_EXPONENT), int>>
	auto pow(OpBinaryDiv<A, B> const& e);

	//! Computes the expression multiplied `N` times with itself.
	/*!
	 * Computes the expression multiplied `N` times with itself.
	 *
	 * \param e The expression to multiply with itself.
	 */
	template<size_t N, typename V, typename... Gs, expr::exp_key_t... Xs, typename = std::enable_if_t<(N > 1 && N < MAX_EXPONENT), int>>
	auto pow(OpTerms<V, Term<Gs, Xs>...> const& e);

	//! Computes the expression multiplied `N` times with itself.
	/*!
	 * Computes the expression multiplied `N` times with itself.
	 *
	 * \param e The expression to multiply with itself.
	 */
	template<size_t N, typename T, typename = std::enable_if_t<(N > 1 && N < MAX_EXPONENT), int>>
	auto pow(OpLiteral<T> const& e);

	//! Computes the expression multiplied `N` times with itself.
	/*!
	 * Computes the expression multiplied `N` times with itself.
	 *
	 * \param e The expression to multiply with itself.
	 */
	template<size_t N, size_t NN, size_t DD, typename = std::enable_if_t<(N > 1 && N < MAX_EXPONENT), int>>
	auto pow(OpFractionLiteral<NN, DD>);

	//! Computes the expression multiplied `N` times with itself.
	/*!
	 * Computes the expression multiplied `N` times with itself.
	 *
	 * \param e The expression to multiply with itself.
	 */
	template<size_t N, size_t NN, size_t DD, typename = std::enable_if_t<(N > 1 && N < MAX_EXPONENT), int>>
	auto pow(OpNegFractionLiteral<NN, DD>);

	//! Computes the expression multiplied `N` times with itself.
	/*!
	 * Computes the expression multiplied `N` times with itself.
	 *
	 * \param e The expression to multiply with itself.
	 */
	template<size_t N>
	auto pow(OpIdentity);

	//! Computes the expression multiplied `N` times with itself.
	/*!
	 * Computes the expression multiplied `N` times with itself.
	 *
	 * \param e The expression to multiply with itself.
	 */
	template<size_t N>
	auto pow(OpNegIdentity);


	inline auto pow(expr::symbols::Symbol)
	{
		return expr::symbols::Symbol{};
	}

	namespace
	{

		template<typename G>
		decltype(auto) pow_get(G&& g, size_t)
		{
			return std::forward<G>(g);
		}

		template<size_t N, typename... Gs, expr::exp_key_t... Xs, size_t... Is>
		auto pow_terms(OpTerms<Term<Gs, Xs>...> const& e, std::index_sequence<Is...>)
		{
			return OpTerms(OpIdentity{}, expr::get<Is>(e).template pow<N>()...);
		}


		template<typename E0, typename E1, size_t... Ns>
		auto pow_add(E0 const& a, OpVoid, std::index_sequence<Ns...>)
		{
			constexpr size_t N0 = sizeof...(Ns) - 1;
			return expr::pow<N0>(a);
		}

		template<typename E0, typename E1, typename E2, typename... Es, size_t... Ns>
		auto pow_add(E0 const& a, OpAdd<E1, E2, Es...> const& b, std::index_sequence<Ns...>)
		{
			constexpr size_t N0 = sizeof...(Ns) - 1;
			return ((expr::choose<N0, Ns>() * expr::pow<N0 - Ns>(a) * expr::pow<Ns>(b)) + ...);
		}

		template<typename E0, typename E1, size_t... Ns>
		auto pow_add(E0 const& a, E1 const& b, std::index_sequence<Ns...>)
		{
			constexpr size_t N0 = sizeof...(Ns) - 1;
			return expr::make_add((expr::choose<N0, Ns>() * expr::pow<N0 - Ns>(a) * expr::pow<Ns>(b))...);
		}

	}

	template<size_t N, typename E, typename>
	auto pow(OpExpression<E> const& e)
	{
		if constexpr (N == 0)
		{
			return OpIdentity{};
		}
		else if constexpr (N == 1)
		{
			return *static_cast<E const*>(&e);
		}
		else if constexpr (N == 2)
		{
			return (*static_cast<E const*>(&e)) * (*static_cast<E const*>(&e));
		}
		else
		{
			constexpr size_t N2 = N / 2;
			constexpr size_t N0 = N - N2 * 2;

			auto p = pow<N2>(*static_cast<E const*>(&e));
			if constexpr (N0 == 0)
			{
				return p * p;
			}
			else
			{
				return p * p * *static_cast<E const*>(&e);
			}
		}
	}

	template<size_t N, typename E, typename>
	auto pow(OpOperator<E> const& e)
	{
		if constexpr (N == 0)
		{
			return OpIdentity{};
		}
		else if constexpr (N == 1)
		{
			return *static_cast<E const*>(&e);
		}
		else if constexpr (N == 2)
		{
			return (*static_cast<E const*>(&e)) * (*static_cast<E const*>(&e));
		}
		else
		{
			constexpr size_t N2 = N / 2;
			constexpr size_t N0 = N - N2 * 2;

			auto p = pow<N2>(*static_cast<E const*>(&e));
			if constexpr (N0 == 0)
			{
				return p * p;
			}
			else
			{
				return p * p * *static_cast<E const*>(&e);
			}
		}
	}

	template<size_t N, typename E0, typename... Es, typename>
	auto pow(OpAdd<E0, Es...> const& e)
	{
		return pow_add(expr::get<0>(e), expr::terms_after_first(e), std::make_index_sequence<N + 1>{});
	}

	template<size_t N, typename A, typename B, typename>
	auto pow(OpBinaryDiv<A, B> const& e)
	{
		return expr::pow<N>(e.a) / expr::pow<N>(e.b);
	}

	template<size_t N, typename V, typename... Gs, expr::exp_key_t... Xs, typename>
	auto pow(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		return expr::pow<N>(e.term) * pow_terms<N>(expr::terms_after_first(e), std::make_index_sequence<sizeof...(Gs)>{});
	}

	template<size_t N, typename T, typename>
	auto pow(OpLiteral<T> const& e)
	{
		using symphas::math::pow;
		return expr::make_literal(pow<N>(e.value));
	}

	template<size_t N, size_t NN, size_t DD, typename>
	auto pow(OpFractionLiteral<NN, DD>)
	{
		if constexpr (N == 2)
		{
			return expr::make_fraction<NN* NN, DD* DD>();
		}
		else
		{
			constexpr size_t N2 = N / 2;
			constexpr size_t N0 = N - N2 * 2;

			auto p = pow<N2>(OpFractionLiteral<NN, DD>{});
			if constexpr (N0 == 0)
			{
				return p * p;
			}
			else
			{
				return p * p * OpFractionLiteral<NN, DD>{};
			}
		}
	}

	template<size_t N, size_t NN, size_t DD, typename>
	auto pow(OpNegFractionLiteral<NN, DD>)
	{
		if constexpr (N % 2 == 1)
		{
			return -pow<N>(OpFractionLiteral<NN, DD>{});
		}
		else
		{
			return pow<N>(OpFractionLiteral<NN, DD>{});
		}
	}

	template<size_t N>
	auto pow(OpIdentity)
	{
		return OpIdentity{};
	}

	template<size_t N>
	auto pow(OpNegIdentity)
	{
		if constexpr (N % 2 == 1)
		{
			return OpNegIdentity{};
		}
		else
		{
			return OpIdentity{};
		}
	}
	
	template<size_t N, size_t D, bool flag, typename E>
	auto pow(E const& e)
	{
		if constexpr (flag)
		{
			return expr::inverse(expr::pow<N, D, false>(e));
		}
		else
		{
			if constexpr (D == 1)
			{
				return expr::pow<N>(e);
			}
			else
			{
				return expr::make_pow<expr::Xk<N, D, false>>(e);
			}
		}
	}

	template<expr::exp_key_t X, typename E>
	auto pow(E const& e)
	{
		return pow<expr::_Xk_t<X>::N, expr::_Xk_t<X>::D, expr::_Xk_t<X>::sign>(e);
	}
}


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

		template<typename V, typename coeff_t, 
			typename = std::enable_if_t<(expr::is_coeff<coeff_t> || std::is_same<OpVoid, coeff_t>::value), int>>
		static auto get(V v, coeff_t const&);

		//template<typename V, typename G>
		//static auto get_g(V v, G g);


		// If passed an OpLiteral, uses its value rather than the object.
		template<typename V, typename E1>
		static auto get(OpLiteral<V> const& value, E1&& e1)
		{
			return get(value.value, std::forward<E1>(e1));
		}

		// If passed an OpLiteral, uses its value rather than the object.
		//template<typename V, typename E1>
		//static auto get_g(OpLiteral<V> const& value, E1&& e1)
		//{
		//	return get(value.value, std::forward<E1>(e1));
		//}

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

	//! Create a function expression of another expression.
	/*!
	 * The function is applied to the given expression through an expression
	 * that represents a function.
	 *
	 * \param v The coefficient to the function expression.
	 * \param a The expression to which the function is applied.
	 */
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto make_function(V const& v, OpExpression<E> const& e, std::string name, F f, std::tuple<Arg0, Args...> const& args)
	{
		return OpFunction<V, E, F, Arg0, Args...>(v, *static_cast<E const*>(&e), f, args);
	}

	//! Create a function expression of another expression.
	/*!
	 * The function is applied to the given expression through an expression
	 * that represents a function.
	 *
	 * \param v The coefficient to the function expression.
	 * \param a The expression to which the function is applied.
	 */
	template<typename V, typename E, typename F, typename Arg0, typename... Args>
	auto make_function(OpExpression<E> const& e, std::string name, F f, std::tuple<Arg0, Args...> const& args)
	{
		return OpFunction<OpIdentity, E, F, Arg0, Args...>(OpIdentity{}, *static_cast<E const*>(&e), f, args);
	}

	//! Create a function expression of another expression.
	/*!
	 * The function is applied to the given expression through an expression
	 * that represents a function.
	 *
	 * \param v The coefficient to the function expression.
	 * \param a The expression to which the function is applied.
	 */
	template<typename V, typename E, typename F>
	auto make_function(V const& v, OpExpression<E> const& e, std::string name, F f)
	{
		return OpFunction<V, E, F, void>(v, *static_cast<E const*>(&e), f);
	}

	//! Create a function expression of another expression.
	/*!
	 * The function is applied to the given expression through an expression
	 * that represents a function.
	 *
	 * \param v The coefficient to the function expression.
	 * \param a The expression to which the function is applied.
	 */
	template<typename V, typename E, typename F>
	auto make_function(OpExpression<E> const& e, std::string name, F f)
	{
		return OpFunction<OpIdentity, E, F, void>(OpIdentity{}, *static_cast<E const*>(&e), f);
	}


}




// **********************************************************************************

//! Applies a particular function to an expression.
/*!
 * The result of an expression is evaluated and then a chosen function is
 * applied.
 *
 * \tparam f The function applied to the result of the expression.
 * \tparam D The object which constructs the derivative of the function.
 * \tparam V The coefficient type of this expression object.
 * \tparam E The type of the expression the function applies to.
 */
template<auto f, typename V, typename E>
struct OpFunctionApply : OpExpression<OpFunctionApply<f, V, E>>
{
	OpFunctionApply() : value{ V{} }, e{} {}

	//! Create a new function expression, applied to an expression.
	/*!
	 * Creates a new function expression with the given coefficient,
	 * applied to the given expression.
	 *
	 * \param value The coefficient of the function expression.
	 * \param e The expression on which the function is applied.
	 */
	OpFunctionApply(V value, E const& e) : value{ value }, e{ e } {}

	//! Create a new function expression, applied to an expression.
	/*!
	 * Creates a new function expression
	 * applied to the given expression.
	 *
	 * \param e The expression on which the function is applied.
	 */
	OpFunctionApply(E const& e) : value{ OpIdentity{} }, e{ e } {}

	inline auto eval(iter_type n) const
	{
		if constexpr (std::is_same_v<expr::eval_type_t<E>, expr::symbols::Symbol>)
		{
			return expr::symbols::Symbol{};
		}
		else
		{
			return expr::eval(value) * f(e.eval(n));
		}
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

	V value;		//!< Coefficient of the function expression term.
	E e;			//!< The expression to which the function is applied.
};




template<typename coeff_t, auto f2, typename V2, typename E2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V2>), int> = 0>
auto operator*(coeff_t const& value, OpFunctionApply<f2, V2, E2> const& b)
{
	return symphas::internal::make_function<f2>::template get(value * b.value, b.e);
}

template<typename coeff_t, typename tensor_t, auto f2, typename E2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpFunctionApply<f2, tensor_t, E2> const& b)
{
	return (value * b.value) * symphas::internal::make_function<f2>::template get(OpIdentity{}, b.e);
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
		return OpFunctionApply<f, V, E>(v, *static_cast<E const*>(&e));
	}

	template<auto f>
	template<typename V, typename coeff_t, typename>
	inline auto make_function<f>::get(V v, coeff_t const& coeff)
	{
		return expr::make_literal(OpFunctionApply<f, V, coeff_t>(v, coeff).eval(0));
	}
}







namespace expr
{


	template<auto f, typename V, typename E>
	auto make_applied_function(V const& v, OpExpression<E> const& e)
	{
		return symphas::internal::make_function<f>::get(v, *static_cast<E const*>(&e));
	}


	namespace
	{
		//! Alias to construct a conjugate function to compute the conjugate.
		template<typename E>
		auto constexpr func_conjugate = &symphas::math::conj<typename expr::eval_type<E>::type>;

		//! Alias to construct a modulus function to compute the modulus.
		template<typename E>
		auto constexpr func_modulus = &symphas::math::modulus<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to return real part.
		template<typename E>
		auto constexpr func_realpart = &symphas::math::real<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to return imaginary part.
		template<typename E>
		auto constexpr func_imagpart = &symphas::math::imag<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the natural logarithm function.
		template<typename E>
		auto constexpr func_log = &symphas::math::log<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the cosine function.
		template<typename E>
		auto constexpr func_cos = &symphas::math::cos<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_sin = &symphas::math::sin<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_tan = &symphas::math::tan<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_csc = &symphas::math::csc<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_sec = &symphas::math::sec<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_cot = &symphas::math::cot<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the cosine function.
		template<typename E>
		auto constexpr func_cosh = &symphas::math::cosh<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_sinh = &symphas::math::sinh<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_tanh = &symphas::math::tanh<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the cosine function.
		template<typename E>
		auto constexpr func_acos = &symphas::math::acos<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_asin = &symphas::math::asin<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_atan = &symphas::math::atan<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the cosine function.
		template<typename E>
		auto constexpr func_acosh = &symphas::math::acosh<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_asinh = &symphas::math::asinh<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_atanh = &symphas::math::atanh<typename expr::eval_type<E>::type>;

		//! Alias to construct the function to compute the sine function.
		template<typename E>
		auto constexpr func_sqrt = &symphas::math::sqrt<typename expr::eval_type<E>::type>;

	}

	




	template<typename E>
	auto conj(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_conjugate<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto modulus(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_modulus<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto real(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_realpart<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto imag(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_imagpart<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto log(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_log<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto cos(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_cos<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto sin(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_sin<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto tan(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_tan<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto cosh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_cosh<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto sinh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_sinh<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto tanh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_tanh<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto acos(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_acos<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto asin(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_asin<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto atan(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_atan<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto acosh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_acosh<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto asinh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_asinh<E>>::get(*static_cast<E const*>(&e));
	}

	template<typename E>
	auto atanh(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_atanh<E>>::get(*static_cast<E const*>(&e));
	}


	template<typename E>
	auto sqrt(OpExpression<E> const& e)
	{
		return symphas::internal::make_function<func_sqrt<E>>::get(*static_cast<E const*>(&e));
	}
}


template<size_t D, typename T0, typename T1>
auto expr::make_unit_vector(T0 const& direction0, T1 const& direction1)
{
	if constexpr (D == 1)
	{
		return OpIdentity{};
	}
	else if constexpr (D == 2)
	{
		auto s = expr::sin(direction0);
		auto c = expr::cos(direction0);
		return (expr::make_column_vector<0, D>() * c + expr::make_column_vector<1, D>() * s);
	}
	else
	{
		auto s0 = expr::sin(direction0);
		auto c0 = expr::cos(direction0);
		auto s1 = expr::sin(direction1);
		auto c1 = expr::cos(direction1);
		return (expr::make_column_vector<0, D>() * c0 * s1 + expr::make_column_vector<1, D>() * s0 * s1 + expr::make_column_vector<2, D>() * c1);
	}
}

template<size_t D, typename T0>
auto expr::make_unit_vector(T0 const& direction0)
{
	return expr::make_unit_vector<D>(direction0, OpVoid{});
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
struct OpFunction : OpExpression<OpFunction<V, E, F, Arg0, Args...>>
{
	OpFunction() : name{ "" }, value{ V{} }, e{}, f{}, tt{} {}

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
	OpFunction(std::string name, V value, E const& e, F f, Arg0&& arg0, Args&&... args) :
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
	OpFunction(std::string name, V value, E const& e, F f, std::tuple<Arg0, Args...> const& tt) :
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
	OpFunction(V value, E const& e, F f, Arg0&& arg0, Args&&... args) :
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
	OpFunction(V value, E const& e, F f, std::tuple<Arg0, Args...> const& tt) :
		name{ "?" }, value{ value }, e{ e }, f{ f }, tt{ tt } {}

	inline auto eval(iter_type n) const
	{
		return expr::eval(value) * function_value(n, std::make_index_sequence<sizeof...(Args) + 1>{});
	}

	auto operator-() const
	{
		return OpFunction<decltype(-value), E, F, Arg0, Args...>(-value, e, f, tt);
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

protected:

	template<size_t... Is>
	auto function_value(iter_type n, std::index_sequence<Is...>) const
	{
		f(e.eval(n), std::get<Is>(tt)...);
	}

public:

	V value;
	F f;							//!< Functor which is applied to an expression.
	std::tuple<Arg0, Args...> tt;	//!< Arguments to the function.
	std::string name;				//!< Identifier of the function.
	E e;							//!< Expression to which a function is applied.
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
struct OpFunction<V, E, F, void> : OpExpression<OpFunction<V, E, F, void>>
{

	OpFunction() : name{ "" }, value{ V{} }, e{}, f{} {}

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
	OpFunction(std::string name, V value, E const& e, F f) :
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
	OpFunction(V value, E const& e, F f) :
		name{ "?" }, value{ value }, e{ e }, f{ f } {}

	inline auto eval(iter_type n) const
	{
		return expr::eval(value) * f(e, n);
	}

	auto operator-() const
	{
		return OpFunction<decltype(-value), E, F, void>(name, -value, e, f);
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

	V value;
	E e;							//!< Expression to which a function is applied.
	F f;							//!< Functor which is applied to an expression.
	std::string name;				//!< Identifier of the function.
};


template<typename V, typename E, typename F>
OpFunction(std::string, V, E, F)->OpFunction<V, E, F, void>;
template<typename V, typename E, typename F>
OpFunction(V, E, F)->OpFunction<V, E, F, void>;
template<typename V, typename E, typename F>
OpFunction(std::string, OpLiteral<V>, E, F)->OpFunction<V, E, F, void>;
template<typename V, typename E, typename F>
OpFunction(OpLiteral<V>, E, F)->OpFunction<V, E, F, void>;
template<typename V, typename E, typename F, typename... Args>
OpFunction(std::string, OpLiteral<V>, E, F, std::tuple<Args...>)->OpFunction<V, E, F, Args...>;
template<typename V, typename E, typename F, typename... Args>
OpFunction(OpLiteral<V>, E, F, std::tuple<Args...>)->OpFunction<V, E, F, Args...>;



template<typename coeff_t, typename V2, typename E2, typename F2, typename Arg, typename... Args,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V2>), int> = 0>
auto operator*(coeff_t const& value, OpFunction<V2, E2, F2, Arg, Args...> const& b)
{
	return OpFunction(b.name, value * b.value, b.e, b.f, b.tt);
}

template<typename coeff_t, typename tensor_t, typename E2, typename F2, typename Arg, typename... Args,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpFunction<tensor_t, E2, F2, Arg, Args...> const& b)
{
	return (value * b.value) * OpFunction(b.name, OpIdentity{}, b.e, b.f, b.tt);
}

template<typename coeff_t, typename V2, typename E2, typename F2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V2>), int> = 0>
auto operator*(coeff_t const& value, OpFunction<V2, E2, F2, void> const& b)
{
	return OpFunction(b.name, value * b.value, b.e, b.f);
}

template<typename coeff_t, typename tensor_t, typename E2, typename F2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpFunction<tensor_t, E2, F2, void> const& b)
{
	return (value * b.value) * OpFunction(b.name, OpIdentity{}, b.e, b.f);
}



template<typename V1, typename E1, typename F1, typename V2>
auto operator*(OpFunction<V2, E1, F1, void> const& a, OpFunction<V2, E1, F1, void> const& b)
{
	return (a.value + b.value) * OpFunction(b.name, OpIdentity{}, a.e, a.f);
}

template<typename V1, typename E1, typename F1, typename Arg, typename... Args, typename V2>
auto operator*(OpFunction<V1, E1, F1, Arg, Args...> const& a, OpFunction<V2, E1, F1, Arg, Args...> const& b)
{
	return (a.value + b.value) * OpFunction(b.name, OpIdentity{}, b.e, b.f, b.tt);
}






