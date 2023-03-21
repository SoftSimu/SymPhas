
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

#include "expressionprototypes.h"

//! Contains all elements constituting the special symbolic functionality.
/*!
 * Defines elements for symbolic manipulation that work in conjunction with regular
 * symbolic expression types.
 */


/*!
 * \addtogroup substitutables
 * @{
 */

template<typename E, size_t... ArgNs>
struct SymbolicTemplate;
template<typename E, typename... Ts>
struct SymbolicFunction;


template<typename... Ts>
struct Substitution;


template<typename Op, typename E, typename Inds>
struct SymbolicSeries;


template<typename T>
struct SymbolicData;
template<typename T>
struct SymbolicDataArray;

template<typename... Ts>
using SymbolicTuple = SymbolicDataArray<std::tuple<Ts...>>;
template<typename... Ts>
using SymbolicTermArray = SymbolicTuple<Term<Ts>...>;
template<typename T, size_t... Zs>
using SymbolicVariableArray = SymbolicTermArray<Variable<Zs, T>...>;


namespace symphas::internal
{
	enum class SeriesOp
	{
		ADD,
		SUB,
		MUL,
		DIV
	};

	template<SeriesOp Op>
	struct ReduceOp;


	template<typename V, typename sub_t, typename... Ts>
	auto make_symbolic_eval(V const& value, sub_t const& data, Ts const&... ts);

}

namespace expr
{
	template<typename T0, typename T1>
	struct series_limits;

	using sum_op = symphas::internal::ReduceOp<symphas::internal::SeriesOp::ADD>;

	namespace symbols
	{
		template<int N, int P = 0>
		struct i_;

		template<typename, typename...>
		struct v_id_type;


		template<typename I0, typename T0>
		using index_eq = symphas::lib::types_list<I0, T0>;
		template<typename I0, typename T0>
		using index_neq = symphas::lib::types_list<std::false_type, I0, T0>;
		template<typename I0, int N0>
		using index_eq_N = symphas::lib::types_list<I0, std::integer_sequence<int, N0>>;
		template<typename I0, int N0>
		using index_neq_N = symphas::lib::types_list<std::false_type, I0, std::integer_sequence<int, N0>>;

	}
}


template<typename... Ts>
struct SymbolicCase;

template<typename I, typename G, typename A, typename B>
using SymbolicTernaryCase = SymbolicCase<
	symphas::lib::types_list<expr::symbols::index_eq<I, G>, A>,
	symphas::lib::types_list<expr::symbols::index_neq<I, G>, B>>;

template<typename I, int N, typename A, typename B>
using SymbolicTernaryNCase = SymbolicCase<
	symphas::lib::types_list<expr::symbols::index_eq_N<I, N>, A>,
	symphas::lib::types_list<expr::symbols::index_neq_N<I, N>, B>>;

template<typename C, typename T = void>
using SymbolicCaseSwap = SymbolicCase<symphas::lib::types_list<C, T>>;

template<typename E, typename T, typename Seq, typename B, typename C>
using SymbolicSum = SymbolicSeries<expr::sum_op, T, symphas::lib::types_list<E, Seq, B, C>>;

template<typename V, typename E, typename T, typename Seq, typename A, typename B, typename C>
using OpSum = OpSymbolicEval<V, SymbolicSum<E, T, Seq, B, C>, A>;




namespace expr
{
	template<typename... Ts>
	auto recreate_series(Ts&&... args);
}


/*!
 * @}
 */