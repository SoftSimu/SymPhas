
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
 * PURPOSE: Forward declaration of all expression objects.
 *
 * ***************************************************************************
 */

#pragma once

#include "spslib.h"

//! Contains all elements constituting the symbolic algebra functionality.
/*!
 * Defines elements which support the symbolic algebra functionality, such
 * as creating new expression terms, manipulating and transforming expressions,
 * and printing expressions.
 */
namespace expr
{
	using exp_key_t = unsigned int;
}

namespace symphas::internal
{
	//! Convert values of the exponent into an exponent key.
	/*! 
	 * Convert values of the exponent into an exponent key. If \p b is
	 * given as zero, then it will consider a to be a key, and return the
	 * real value of the exponent from the key.
	 * 
	 * \tparam a The numerator of the exponent.
	 * \tparam b The denominator of the exponent.
	 * \tparam sign If true, then it means the exponent is negative.
	 */
	template<expr::exp_key_t a, expr::exp_key_t b = 1, bool sign = false>
	struct exp_compute_key;

	struct tensor_fold;
	struct tensor_cast;

	//! A mapping to map half complex plane (FFTW format) to sequential complex values.
	struct HCTS {};
	//! A sequential complex values to half complex plane (FFTW format).
	struct STHC {};
}

namespace expr
{
	template<exp_key_t a, exp_key_t b = 1, bool sign = false>
	using Xk_t = symphas::internal::exp_compute_key<a, b, sign>;

	template<exp_key_t a>
	using _Xk_t = symphas::internal::exp_compute_key<a, 0, false>;

	//! Exponent key helper to convert fraction into an exponent key.
	template<exp_key_t a, exp_key_t b = 1, bool sign = false>
	constexpr auto Xk = Xk_t<a, b, sign>::value;

	//! Exponent key helper to convert an exponent key into an exponent value.
	template<exp_key_t value>
	constexpr auto _Xk = Xk<value, 0, false>;


	template<exp_key_t X1, exp_key_t X2>
	constexpr int compute_XXk_N =
		((expr::_Xk_t<X1>::sign) ? -1 : 1) * int(expr::_Xk_t<X1>::N * expr::_Xk_t<X2>::D) +
		((expr::_Xk_t<X2>::sign) ? -1 : 1) * int(expr::_Xk_t<X2>::N * expr::_Xk_t<X1>::D);

	template<exp_key_t X1, exp_key_t X2>
	constexpr int compute_XXk_D = expr::_Xk_t<X1>::D * expr::_Xk_t<X2>::D;

	template<exp_key_t X1, exp_key_t X2>
	using XXk_t = Xk_t<expr::exp_key_t((compute_XXk_N<X1, X2> < 0) ? -compute_XXk_N<X1, X2> : compute_XXk_N<X1, X2>), compute_XXk_D<X1, X2>, (compute_XXk_N<X1, X2> < 0)>;

	namespace symbols
	{
		struct Symbol 
		{
			auto eval(iter_type = 0) const
			{
				return *this;
			}

			inline auto operator+(expr::symbols::Symbol)
			{
				return expr::symbols::Symbol{};
			}

			inline auto operator-(expr::symbols::Symbol)
			{
				return expr::symbols::Symbol{};
			}

			inline auto operator*(expr::symbols::Symbol)
			{
				return expr::symbols::Symbol{};
			}

			inline auto operator/(expr::symbols::Symbol)
			{
				return expr::symbols::Symbol{};
			}

			template<typename E>
			auto operator+(E)
			{
				return expr::symbols::Symbol{};
			}

			template<typename E>
			auto operator-(E)
			{
				return expr::symbols::Symbol{};
			}

			template<typename E>
			auto operator*(E)
			{
				return expr::symbols::Symbol{};
			}

			template<typename E>
			auto operator/(E)
			{
				return expr::symbols::Symbol{};
			}

			auto& operator+=(Symbol)
			{
				return *this;
			}

			auto& operator-=(Symbol)
			{
				return *this;
			}
		};


		template<typename E>
		auto operator+(E, expr::symbols::Symbol)
		{
			return expr::symbols::Symbol{};
		}

		template<typename E>
		auto operator-(E, expr::symbols::Symbol)
		{
			return expr::symbols::Symbol{};
		}

		template<typename E>
		auto operator*(E, expr::symbols::Symbol)
		{
			return expr::symbols::Symbol{};
		}

		template<typename E>
		auto operator/(E, expr::symbols::Symbol)
		{
			return expr::symbols::Symbol{};
		}

	}

	template<size_t N>
	auto pow(symbols::Symbol)
	{
		return symbols::Symbol{};
	}

}



/*!
 * \defgroup Op Symbolic Algebra Objects
 * @{
 */

template<typename T, size_t... Ns>
struct OpTensor;
struct OpIdentity;
struct OpNegIdentity;
struct OpVoid;
template<typename T>
struct OpLiteral;
template<size_t N, size_t D>
struct OpFractionLiteral;
template<size_t N, size_t D>
struct OpNegFractionLiteral;

template<typename... Es>
struct OpAdd;
template<typename E1, typename E2>
struct OpBinaryMul;
template<typename E1, typename E2>
struct OpBinaryDiv;



template<typename E>
struct OpExpression;
template<typename E>
struct OpOperator;


template<size_t Z, typename G = OpVoid>
struct Variable;
template<typename G>
struct NamedData;


template<typename G, expr::exp_key_t X>
struct Term;
template<typename... Ts>
struct OpTerms;
template<typename V, typename G>
using OpTerm = OpTerms<V, Term<G, 1>>;



template<typename G>
struct SymbolicDerivative;
template<typename Dd, typename V, typename E, typename T>
struct OpFuncDerivative;
template<size_t O, typename V, typename T>
struct OpOperatorDerivative;
template<Axis ax, size_t O, typename V, typename Sp>
struct OpOperatorDirectionalDerivative;
template<typename V, typename Sp, size_t... Os>
struct OpOperatorMixedDerivative;

template<typename V, typename E1, typename E2>
struct OpFuncConvolution;
template<size_t D>
struct GaussianSmoothing;

template<typename A1, typename A2, typename E>
struct OpCombination;
template<typename A1, typename A2>
struct OpOperatorCombination;
template<typename A1, typename A2, typename E>
struct OpChain;
template<typename A1, typename A2>
struct OpOperatorChain;
template<typename E>
struct OpOperator;
template<typename V, typename E>
struct OpExponential;
template<typename G, typename V, typename E>
struct OpMap;

template<typename V, typename E, typename F, typename Arg, typename... Args>
struct OpFunc;
template<auto f, typename V, typename E>
struct OpFuncApply;


template<size_t D, Axis ax>
struct GridAxis;


//! @}

namespace expr::prune
{

	template<typename E>
	void update(OpExpression<E>& e);
}


