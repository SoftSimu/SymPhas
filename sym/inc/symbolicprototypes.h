
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


template<typename E, typename K = void>
struct SymbolicListIndex;


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
	auto make_symbolic_eval(V const& value, sub_t const& data, Ts&&... ts);

	//template<typename V, typename Op, typename E, typename Inds, typename... Ts>
	//auto make_symbolic_eval(V const& value, SymbolicSeries<Op, E, Inds> const& data, Ts const&... ts);

	//template<typename V, typename E0, typename... T0s, typename... Ts>
	//auto make_symbolic_eval(V const& value, SymbolicFunction<E0, T0s...> const& data, Ts const&... ts);

	//template<typename V, typename E, typename... Ts>
	//auto make_symbolic_eval(V const& value, OpExpression<E> const& e, Ts const&... ts);

	//template<typename V, expr::NoiseType nt, typename T, size_t D, typename E>
	//auto make_symbolic_eval(V const& value, NoiseData<nt, T, D> const& noise, OpExpression<E> const& e);

	//template<typename V, expr::NoiseType nt, typename T, size_t D, typename E, typename... Ts>
	//auto make_symbolic_eval(V const& value, NoiseData<nt, T, D> const& noise, SymbolicFunction<E, Ts...> const& f);

	template<typename E>
	constexpr bool has_index_impl = false;
}

namespace expr
{
	template<typename T0, typename T1>
	struct series_limits;

	template<typename T, typename V, typename G>
	auto limit_0(expr::series_limits<OpTerm<V, G>, T> const& limit);
	template<typename T1, typename T2>
	auto limit_0(expr::series_limits<T1, T2> const& limit);
	template<typename T, typename V, typename G>
	auto limit_1(expr::series_limits<T, OpTerm<V, G>> const& limit);
	template<typename T1, typename T2>
	auto limit_1(expr::series_limits<T1, T2> const& limit);

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

		template<typename I0, typename T0, typename T1>
		using index_iti = symphas::lib::types_list<index_eq<I0, T0>, index_eq<I0, T1>>;
		template<typename I0, typename T0, int N1>
		using index_itN = symphas::lib::types_list<index_eq<I0, T0>, index_eq_N<I0, N1>>;
		template<typename I0, int N0, typename T1>
		using index_Nti = symphas::lib::types_list<index_eq_N<I0, N0>, index_eq<I0, T1>>;
		template<typename I0, int N0, int N1>
		using index_NtN = symphas::lib::types_list<index_eq_N<I0, N0>, index_eq_N<I0, N1>>;
	}

	
	template<typename E>
	constexpr bool has_index = symphas::internal::has_index_impl<symphas::lib::types_list<E>>;

	template<typename Is, typename E>
	constexpr bool has_selected_index = false;

    template<typename C, typename T>
    struct case_entry {};

}


template<typename... Ts>
struct SymbolicCase;

template<typename I, typename G, typename A, typename B>
using SymbolicTernaryCase = SymbolicCase<
	expr::case_entry<expr::symbols::index_eq<I, G>, A>,
	expr::case_entry<expr::symbols::index_neq<I, G>, B>>;

template<typename I, int N, typename A, typename B>
using SymbolicTernaryNCase = SymbolicCase<
	expr::case_entry<expr::symbols::index_eq_N<I, N>, A>,
	expr::case_entry<expr::symbols::index_neq_N<I, N>, B>>;

template<typename C = void, typename T = void>
using SymbolicCaseSwap = SymbolicCase<expr::case_entry<C, T>>;

template<typename E, typename T, typename Seq, typename B, typename C>
using SymbolicSum = SymbolicSeries<expr::sum_op, T, symphas::lib::types_list<E, Seq, B, C>>;

template<typename V, typename E, typename T, typename Seq, typename A, typename B, typename C>
using OpSum = OpSymbolicEval<V, SymbolicSum<E, T, Seq, B, C>, A>;

template<>
struct SymbolicCase<expr::case_entry<void, void>> : expr::symbols::Symbol {};



namespace expr
{
	template<typename... Ts>
	auto recreate_series(Ts&&... args);


	template<int N, int P, typename E>
	auto make_list(expr::symbols::i_<N, P>, OpExpression<E> const& e);

	template<typename K, typename E0, typename E>
	auto make_list(OpExpression<E0> const& index, OpExpression<E> const& e);


}


namespace expr
{

	template<typename Op, typename... Is, typename E, typename = std::enable_if_t<(
		(std::is_convertible<Is, expr::symbols::Symbol>::value && ... && true)), int>>
	auto series(OpExpression<E> const& e);

	template<typename Op, typename... Is, int N, int P>
	auto series(expr::symbols::i_<N, P> const& e);

	template<typename... Is, typename E, typename = std::enable_if_t<(
		(std::is_convertible<Is, expr::symbols::Symbol>::value && ... && true)), int>>
	auto sum(OpExpression<E> const& e);

	template<typename... Is, int N, int P>
	auto sum(expr::symbols::i_<N, P> const& e);

	template<typename... Is, typename E, typename = std::enable_if_t<(
		(std::is_convertible<Is, expr::symbols::Symbol>::value && ... && true)), int>>
	auto prod(OpExpression<E> const& e);

	template<typename... Is, int N, int P>
	auto prod(expr::symbols::i_<N, P> const& e);

}

/*!
 * @}
 */
