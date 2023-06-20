
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
#include "expressionproperties.h"
#include "expressionaggregates.h"


namespace symphas::internal
{

	//! Implementation of convolution expression construction.
	/*!
	 * Implementation of functions which generalize and simplify the way to
	 * construct convolution expressions. Wraps the template deduction necessary
	 * to initialize a convolution expression.
	 */
	struct make_integral
	{
		template<typename V, typename E, typename T>
		static auto get(V const& v, OpExpression<E> const& e, T const& domain);

		template<typename V, typename E>
		static auto get(V const& v, OpExpression<E> const& e, symphas::grid_info const& domain);

		template<typename V, typename E, typename T>
		static auto get(V const& v, OpOperator<E> const& e, T const& domain);

		template<typename V, typename E>
		static auto get(V const& v, OpOperator<E> const& e, symphas::grid_info const& domain);
	};

}



namespace expr
{
	template<typename T, typename V, typename A>
	auto make_integral(V&& v, A&& a, T&& domain)
	{
		return symphas::internal::make_integral::get(std::forward<V>(v), std::forward<A>(a), std::forward<T>(domain));
	}

	template<typename A, typename T>
	auto make_integral(A&& a, T&& domain)
	{
		return make_integral(OpIdentity{}, std::forward<A>(a), std::forward<T>(domain));
	}

	template<size_t Z, typename V, typename A>
	auto make_integral(V&& v, A&& a)
	{
		return make_integral(std::forward<V>(v), std::forward<A>(a), Variable<Z>{});
	}

	template<size_t Z, typename A>
	auto make_integral(A&& a)
	{
		return make_integral(OpIdentity{}, std::forward<A>(a), Variable<Z>{});
	}

	template<typename V, typename A>
	auto make_domain_integral(V&& v, A&& a, symphas::grid_info const& info)
	{
		return symphas::internal::make_integral::get(std::forward<V>(v), std::forward<A>(a), info);
	}

	template<typename A>
	auto make_domain_integral(A&& a, symphas::grid_info const& info)
	{
		return make_domain_integral(OpIdentity{}, std::forward<A>(a), info);
	}
}



//! Represents the integration of an expression.
/*!
 * The variable of integration is chosen by the variable index specified
 * by the first template parameter.
 * 
 * \tparam T The index of the variable of integration.
 * \tparam V The type of the coefficient.
 * \tparam E The type of the expression which is being integrated.
 */
template<typename V, typename E, typename T>
struct OpIntegral : OpExpression<OpIntegral<V, E, T>>
{
	using result_t = expr::eval_type_t<E>;

	OpIntegral() : value{ V{} }, domain{ T{} }, data{}, e{} {}
	OpIntegral(V value, E const& e, T const& domain) : value{ value }, domain{ domain }, data{ result_t{} }, e{ e } {}

	inline auto eval(iter_type n = 0) const
	{
		return expr::eval(value) * expr::eval(data, n);
	}

	auto operator-() const
	{
		return symphas::internal::make_integral::get(-value, e, domain);
	}

	auto update()
	{
		auto len = expr::data_length(e);
		if (len > 0)
		{
			expr::prune::update(e);
			data = e.eval(0);

//#			pragma omp parallel for reduction (+:data)
			for (iter_type i = 1; i < len; ++i)
			{
				data += e.eval(i);
			}
		}
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::integral_print<T>{}(out, domain, value * e);
	}

	size_t print(char* out) const
	{
		return expr::integral_print<T>{}(out, domain, value * e);
	}

	size_t print_length() const
	{
		return expr::integral_print<T>{}(domain, value * e);
	}

#endif

	template<typename V0, typename E0, typename T0>
	friend auto const& expr::get_enclosed_expression(OpIntegral<V0, E0, T0> const&);
	template<typename V0, typename E0, typename T0>
	friend auto& expr::get_enclosed_expression(OpIntegral<V0, E0, T0>&);
	template<typename V0, typename E0, typename T0>
	friend auto const& expr::get_result_data(OpIntegral<V0, E0, T0> const&);
	template<typename V0, typename E0, typename T0>
	friend auto& expr::get_result_data(OpIntegral<V0, E0, T0>&);


	V value;						// value multiplying the result of this derivative
	symphas::grid_info domain;		//!< Information about the domain of integration. 

protected:

	result_t data;						//!< Grid storing the resulting values.
	E e;								//!< expression object specifying grid values
};


//! Represents the integration of an expression.
/*!
 * The variable of integration is chosen by the variable index specified
 * by the first template parameter.
 *
 * \tparam T The index of the variable of integration.
 * \tparam V The type of the coefficient.
 * \tparam E The type of the expression which is being integrated.
 */
template<typename V, typename E, typename T>
struct OpIntegral<V, E, expr::variational_t<T>> : OpExpression<OpIntegral<V, E, expr::variational_t<T>>>
{
	using result_t = expr::eval_type_t<E>;

	OpIntegral() : value{ V{} }, domain{ symphas::grid_info{ nullptr, 0 } }, data{}, e{} {}
	OpIntegral(V value, E const& e, symphas::grid_info const& domain) 
		: value{ value }, domain{ domain }, data{ result_t{} }, e{ e } {}

	inline auto eval(iter_type n = 0) const
	{
		return expr::eval(value) * expr::eval(data, n);
	}

	auto operator-() const
	{
		return expr::make_domain_integral(-value, e, domain);
	}

	auto update()
	{
		auto len = domain.num_points();
		if (len > 0)
		{
			auto iter_data = expr::eval_iters(e);
			iter_type* iters = iter_data.first;
			len_type n = iter_data.second;

			if (n > 0)
			{
				data = std::reduce(
#ifdef EXECUTION_HEADER_AVAILABLE
					std::execution::par_unseq,
#endif
					static_cast<const E*>(&e)->begin(iters),
					static_cast<const E*>(&e)->end(iters, n), result_t{}
				);
			}
			else
			{
				data = std::reduce(
#ifdef EXECUTION_HEADER_AVAILABLE
					std::execution::par_unseq,
#endif
					static_cast<const E*>(&e)->begin(),
					static_cast<const E*>(&e)->end(len), result_t{}
				);

			}

			double r = 1;
			for (auto width : domain.get_widths())
			{
				r *= width;
			}
			expr::result(expr::make_term(1. / r, data), data);
		}
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::integral_print<expr::variational_t<T>>{}(out, domain, value * e);
	}

	size_t print(char* out) const
	{
		return expr::integral_print<expr::variational_t<T>>{}(out, domain, value * e);
	}

	size_t print_length() const
	{
		return expr::integral_print<expr::variational_t<T>>{}(domain, value * e);
	}

#endif

	template<typename V0, typename E0, typename T0>
	friend auto const& expr::get_enclosed_expression(OpIntegral<V0, E0, T0> const&);
	template<typename V0, typename E0, typename T0>
	friend auto& expr::get_enclosed_expression(OpIntegral<V0, E0, T0>&);
	template<typename V0, typename E0, typename T0>
	friend auto const& expr::get_result_data(OpIntegral<V0, E0, T0> const&);
	template<typename V0, typename E0, typename T0>
	friend auto& expr::get_result_data(OpIntegral<V0, E0, T0>&);


	V value;						// value multiplying the result of this derivative
	symphas::grid_info domain;		//!< Information about the domain of integration. 

protected:

	result_t data;						//!< Grid storing the resulting values.
	E e;								//!< expression object specifying grid values
};


//! Represents the integration of an expression.
/*!
 * The variable of integration is chosen by the variable index specified
 * by the first template parameter.
 *
 * \tparam T The index of the variable of integration.
 * \tparam V The type of the coefficient.
 * \tparam E The type of the expression which is being integrated.
 */
template<typename V, typename E, typename T>
struct OpIntegral<V, E, SymbolicDerivative<T>> : OpExpression<OpIntegral<V, E, SymbolicDerivative<T>>>
{
	OpIntegral() : e{}, value{ V{} }, domain{ SymbolicDerivative<T>{} } {}
	OpIntegral(V value, E const& e, SymbolicDerivative<T> domain = SymbolicDerivative<T>{}) : e{ e }, value{ value } {}

	inline auto eval(iter_type n = 0) const
	{
		return expr::symbols::Symbol{};
	}

	auto operator-() const
	{
		return symphas::internal::make_integral::get(-value, e, SymbolicDerivative<T>{});
	}

	auto update() {}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::integral_print<T>{}(out, domain, value * e);
	}

	size_t print(char* out) const
	{
		return expr::integral_print<T>{}(out, domain, value * e);
	}

	size_t print_length() const
	{
		return expr::integral_print<T>{}(domain, value * e);
	}

#endif

	template<typename V0, typename E0, typename T0>
	friend auto const& expr::get_enclosed_expression(OpIntegral<V0, E0, T0> const&);
	template<typename V0, typename E0, typename T0>
	friend auto& expr::get_enclosed_expression(OpIntegral<V0, E0, T0>&);

	V value;							// value multiplying the result of this derivative
	SymbolicDerivative<T> domain;

protected:

	//result_t data;					//!< Grid storing the resulting values.
	E e;								//!< expression object specifying grid values
};


template<typename V, typename E, typename T>
using OpDomainIntegral = OpIntegral<V, E, expr::variational_t<T>>;


template<typename coeff_t, typename V2, typename E2, typename T2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V2>), int> = 0>
auto operator*(coeff_t const& value, OpIntegral<V2, E2, T2> const& b)
{
	return expr::make_integral(value * b.value, expr::get_enclosed_expression(b), b.domain);
}

template<typename coeff_t, typename tensor_t, typename E2, typename T2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t>&& expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpIntegral<tensor_t, E2, T2> const& b)
{
	return (value * b.value) * expr::make_integral(OpIdentity{}, expr::get_enclosed_expression(b), b.domain);
}


template<typename coeff_t, typename V2, typename E2, typename T2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V2>), int> = 0>
auto operator*(coeff_t const& value, OpIntegral<V2, E2, expr::variational_t<T2>> const& b)
{
	return expr::make_domain_integral(value * b.value, expr::get_enclosed_expression(b), b.domain);
}

template<typename coeff_t, typename tensor_t, typename E2, typename T2,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpIntegral<tensor_t, E2, expr::variational_t<T2>> const& b)
{
	return (value * b.value) * expr::make_domain_integral(OpIdentity{}, expr::get_enclosed_expression(b), b.domain);
}




namespace symphas::internal
{
	template<typename V, typename E, typename T>
	inline auto make_integral::get(V const& v, OpExpression<E> const& e, T const& domain)
	{
		return OpIntegral<V, E, T>(v, *static_cast<E const*>(&e), domain);
	}

	template<typename V, typename E>
	inline auto make_integral::get(V const& v, OpExpression<E> const& e, symphas::grid_info const& domain)
	{
		return OpIntegral<V, E, expr::variational_t<symphas::grid_info>>(v, *static_cast<E const*>(&e), domain);
	}

	template<typename V, typename E, typename T>
	inline auto make_integral::get(V const& v, OpOperator<E> const& e, T const& domain)
	{
		return OpIntegral<V, E, T>(v, *static_cast<E const*>(&e), domain);
	}

	template<typename V, typename E>
	inline auto make_integral::get(V const& v, OpOperator<E> const& e, symphas::grid_info const& domain)
	{
		return OpIntegral<V, E, expr::variational_t<symphas::grid_info>>(v, *static_cast<E const*>(&e), domain);
	}
}



