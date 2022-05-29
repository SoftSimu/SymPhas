
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


#include "expressionaggregates.h"

//! \cond

#ifdef LATEX_PLOT
#define SYEX_EXP_FMT_A "\\exp\\left\\{"
#define SYEX_EXP_FMT_B "\\right\\}"
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
		static auto get(V, OpExpression<E> const&);

		template<typename V, typename S, typename G>
		static auto get(V, OpLVariable<S, G> const&);

		template<typename V, typename S>
		static auto get(V v, OpLiteral<S> const& e);
		template<typename V>
		static auto get(V v, OpVoid const e);
		template<typename V>
		static auto get(V v, OpIdentity const e);
		template<typename V>
		static auto get(V v, OpNegIdentity const e);

		//! Constructs the exponential using a grid instead of an expression.
		/*!
		 * Used for the exponential specialization for the OpLVariable.
		 */
		template<typename V, typename G>
		static auto get_g(V v, G g);
	};
}



namespace expr
{
	//! Create an exponential expression term.
	/*!
	 * Create an exponential expression term with the given expression,
	 * which will be put in the exponent.
	 * 
	 * \param a The term that is exponentiated.
	 */
	template<typename A>
	auto make_exponential(A&& a)
	{
		return symphas::internal::make_exponential::template get(std::forward<A>(a));
	}

	//! Create an exponential expression term.
	/*!
	 * Create an exponential expression term with the given expression,
	 * which will be put in the exponent.
	 *
	 * \param v The coefficient of the exponential.
	 * \param a The term that is exponentiated.
	 */
	template<typename V, typename A>
	auto make_exponential(V&& v, A&& a)
	{
		return symphas::internal::make_exponential::template get(std::forward<V>(v), std::forward<A>(a));
	}
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
	using G_T = typename expr::eval_type<E>::type;
	static const int G_U = expr::grid_dim<E>::dimension;
	using G = typename expr::grid_type<E>::type;

	//! Generate an exponential term with the given expression.
	/*!
	 * Generate an exponential term with the given expression.
	 * 
	 * \param value The coefficient of the exponential.
	 * \param e The expression that is exponentiated.
	 */
	OpExponential(V value, E const& e) : value{ value }, e{ e }, data{ expr::property::data_len(e) }
	{
		update();
	}

	inline auto eval(iter_type n) const
	{
		using namespace std;
		using namespace symphas::math;
		return value * exp(data[n]);
	}

	void update()
	{
		expr::result(e, data);
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

	friend struct expr::compound_get;


protected:

	Block<G_T> data;				//!< Grid storing the result of the first expression.


};



//! Exponential function of an expression.
/*!
 * The expression is natural number \f$e\f$ to the power of the given
 * expression. Specialization where a variable is exponentiated.
 *
 * \tparam V The type of the coefficient.
 * \tparam E The type of the expression which is exponentiated.
 */
template<typename V, typename S, typename G>
struct OpExponential<V, OpLVariable<S, G>> : OpExpression<OpExponential<V, OpLVariable<S, G>>>
{

	using E = OpLVariable<S, G>;

	//! Create the exponential expression.
	/*!
	 * Create the expression for the exponential of the given expression.
	 * 
	 * \param value The coefficient of the exponential.
	 * \param e The expression which is exponentiated.
	 */
	OpExponential(V value, E const& e) : value{ value }, pow{ e.value }, data{ e.data } {}

	//! Create the exponential expression.
	/*!
	 * Create the expression for the exponential of the given expression.
	 *
	 * \param value The coefficient of the exponential.
	 * \param g The data that is exponentiated.
	 */
	OpExponential(V value, G data) : value{ value }, pow{ OpIdentity{} }, data{ data } {}

	inline auto eval(iter_type n) const
	{
		using std::exp;
		using namespace symphas::math;
		return value * exp(pow * expr::BaseData<G>::get(data, n));
	}

	auto operator-() const
	{
		return symphas::internal::make_exponential::get_g(-value, data);
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += fprintf(out, SYEX_EXP_FMT_A);
		n += expr::print_with_coeff(out, pow);
		n += fprintf(out, "%s", expr::get_op_name(data));
		n += fprintf(out, SYEX_EXP_FMT_B);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = expr::print_with_coeff(out, value);
		n += sprintf(out + n, SYEX_EXP_FMT_A);
		n += expr::print_with_coeff(out + n, pow);
		n += sprintf(out + n, "%s", expr::get_op_name(data));
		n += sprintf(out + n, SYEX_EXP_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + expr::coeff_print_length(pow) 
			+ SYEX_EXP_FMT_LEN + std::strlen(expr::get_op_name(data));
	}

#endif

	V value;			//!< Value multiplying the result of this convolution.
	S pow;				//!< The coefficient of the original variable.
	G data;				//!< Grid storing the result of the first expression.
};



namespace symphas::internal
{

	template<typename A>
	inline auto make_exponential::get(A&& a)
	{
		return get(OpIdentity{}, std::forward<A>(a));
	}

	template<typename V, typename E>
	inline auto make_exponential::get(V v, OpExpression<E> const& a)
	{
		return OpExponential<V, E>(v, *static_cast<const E*>(&a));
	}

	template<typename V, typename S, typename G>
	inline auto make_exponential::get(V v, OpLVariable<S, G> const& a)
	{
		return OpExponential<V, OpLVariable<S, G>>(v, a);
	}

	template<typename V, typename G>
	inline auto make_exponential::get_g(V v, G g)
	{
		return OpExponential<V, OpLVariable<OpIdentity, G>>(v, g);
	}


	template<typename V, typename S>
	inline auto make_exponential::get(V v, OpLiteral<S> const& e)
	{
		using namespace std;
		using namespace symphas::math;
		return expr::make_literal(v * exp(e.value));
	}
	template<typename V>
	inline auto make_exponential::get(V v, OpVoid const)
	{
		return v;
	}
	template<typename V>
	inline auto make_exponential::get(V v, OpIdentity const e)
	{
		return get(v, expr::make_literal(e.eval()));
	}
	template<typename V>
	inline auto make_exponential::get(V v, OpNegIdentity const e)
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

template<typename T, typename V, typename E>
auto operator*(OpLiteral<T> const a, OpExponential<V, E> const& b)
{
	return symphas::internal::make_exponential::template get(a.value * b.value, expr::compound_get::template expr(b));
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




