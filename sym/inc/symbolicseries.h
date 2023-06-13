
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
 * PURPOSE: Implements substitutable functions and series, including sum
 * and product.
 *
 * ***************************************************************************
 */

#pragma once

#include "symbolicsubstitutables.h"

#ifdef LATEX_PLOT

#define SYEX_SUM_FMT_AA "\\sum_{"
#define SYEX_SUM_FMT_AB "}"
#define SYEX_SUM_FMT_BA "{"
#define SYEX_SUM_FMT_BB "}"

#else


#define SYEX_SUM_FMT_AA "SUM["
#define SYEX_SUM_FMT_AB "]"
#define SYEX_SUM_FMT_BA "("
#define SYEX_SUM_FMT_BB ")"

#endif


#define SYEX_SUM_FMT SYEX_SUM_FMT_AA "%s" SYEX_SUM_FMT_AB SYEX_SUM_FMT_BA "%s" SYEX_SUM_FMT_BB
//#define SYEX_MUL_FMT_LEN (STR_ARR_LEN(SYEX_SUM_FMT_AA SYEX_SUM_FMT_AB SYEX_SUM_FMT_BA SYEX_SUM_FMT_BB) - 1)



/*!
 * \addtogroup substitutables Symbolic Constructs and Expression Building
 * Symbolic constructs are used to simplify expression definitions and easily substitute data.
 *
 * The symbolic algebra library defines two types of symbolic algebra groups: Classes with a name
 * beginning in **Op**, which indicates that they are a symbolic algebra expression which can
 * be _evaluated_; and, Classes which begin with **Symbolic**, which are constructs that extend
 * the functionality/usability of symbolic expressions and cannot be evaluated.
 *
 * The **Symbolic** constructs include an evaluable function (SymbolicFunction), a pattern
 * or template type that allows substitution of expressions (SymbolicTemplate), and a series
 * object (SymbolicSeries).
 * @{
 */


namespace expr
{

	template<typename Op, typename... Is, typename E, typename>
	auto series(OpExpression<E> const& e);

	template<typename Op, typename... Is, int N, int P>
	auto series(expr::symbols::i_<N, P> const& e)
	{
		return series<Op, Is...>(expr::symbols::i_op_type<N, P>{});
	}

	template<typename... Is, typename E, typename>
	auto sum(OpExpression<E> const& e);

	template<typename... Is, int N, int P>
	auto sum(expr::symbols::i_<N, P> const& e)
	{
		return sum<Is...>(expr::symbols::i_op_type<N, P>{});
	}

	template<typename... Is, typename E, typename>
	auto prod(OpExpression<E> const& e);

	template<typename... Is, int N, int P>
	auto prod(expr::symbols::i_<N, P> const& e)
	{
		return prod<Is...>(expr::symbols::i_op_type<N, P>{});
	}


	inline auto series_counter(int start, int end)
	{
		return DynamicIndex(std::move(start), start, end);
	}
}


namespace symphas::internal
{
	template<>
	struct ReduceOp<SeriesOp::ADD>
	{
		template<typename... Es>
		auto operator()(Es&&... es)
		{
			return (std::forward<Es>(es) + ... + OpVoid{});
		}
	};

	template<>
	struct ReduceOp<SeriesOp::SUB>
	{
		template<typename E0, typename... Es>
		auto operator()(E0&& e0, Es&&... es)
		{
			return std::forward<E0>(e0) - (std::forward<Es>(es) + ... + OpVoid{});
		}
	};

	template<>
	struct ReduceOp<SeriesOp::MUL>
	{
		template<typename... Es>
		auto operator()(Es&&... es)
		{
			return (OpIdentity{} * ... * std::forward<Es>(es));
		}
	};

	template<>
	struct ReduceOp<SeriesOp::DIV>
	{
		template<typename... Es>
		auto operator()(Es&&... es)
		{
			return (std::forward<Es>(es) / ... / OpIdentity{});
		}
	};


}

/*!
 * \addtogroup substitutables
 * \subsection series SymbolicSeries
 *
 * A SymbolicSeries is a series associated with a specific operation (specified when the series 
 * is defined) and applied on a particular expression. It is defined using an 
 * expression that contains expr::symbols::i_ and expr::symbols::v_
 * terms. The term `i_` represents the index of the series, where indices are distinguished using
 * their first template value. The term `v_` represents a data that is substituted into the series
 * expression, and is templated using the `i_` type. Depending on the template substitutions,
 * different values will be substituted when the series is expanded.
 *
 * The series can be directly expanded into a fully constructed symbolic expression where every
 * term of the series will appear (and symbolic algebra rules have been applied), or it can be 
 * converted into an OpSymbolicEval type where the series is iteratively evaluated using a list 
 * of data. The former method is used via the member function SymbolicSeries::expand(), and the 
 * latter is used by calling the SymbolicSeries instance directly. If the series would be
 * expanded to 200 or more terms, then OpSymbolicEval should be used. 
 * 
 * \subsubsection series_construction Constructing Series
 *
 * The specific meanings of the indices are very important, and are described using definitions:
 *
 * _Definition of the **Index Offset**_
 * 
 * > When used in an expression to define a series, `i_<N, P>` \f$\rightarrow\;{i_N + P}\f$ is 
 * > the \f$i_N\f$ `N`th unique series index offset by a value `P`. 
 *
 * The uniqueness of indices
 * is relevant in iterated series defined using multiple indices. The offset allows conveniently
 * defining adding or subtracting a value to the index without directly representing it in symbolic
 * algebra.
 * 
 * _Definition of **Placeholder** and **Substitution List**_
 * 
 * A series may be expanded between two integers, defining the starting value of the index and end
 * value. The series may be defined in such a way that data or other symbolic expressions that are
 * passed to it are substituted into the defining expression of the series when it is expanded. 
 * The the object expr::symbols::v_ acts as the **placeholder** for these substitutions, and is 
 * templated on the unique index type, expr::symbols::i_. 
 * 
 * > When used in an expression to define a series, `v_<i_<N, P>>` 
 * > \f$\rightarrow\;v_{i_N + P}\f$ is the placeholder that is indexed by \f$i_N\f$ and defines
 * > that a substitution element put there is offset by
 * > \f$P\f$ from the _zero-th_ element of the substitution list. 
 * 
 * > The _zero-th_ element of a substitution list is defined by the index offset used to define the
 * > series. 
 * 
 * Example: `sum<i_<0, 2>>(e)` defines a sum series using the `0`th unique index and defines that
 * the substitution list has its zero-th element at index 2. That is, given the list `a, b, c, d`
 * substituted into this series, `c` would be substituted into an occurence of `v_<i_<0, 0>>`
 * for the first term of the series expansion.
 * 
 * > Remark: the position of the element in the
 * > substitution list that is substituted to a placeholder does not depend the on starting value
 * > of index of the expansion.
 * 
 * Example: `sum<i_<0, 0>>(e)` defines a sum series using the `0`th unique index and defines
 * first element of the substitution list as the zero-th substituted element. If this series
 * is expanded from 3 using the list `a, b, c, d`:
 * 
 *     sum<i_<0, 0>>(e).expand<3>(a, b, c, d);
 * 
 * then the first substituted element would be `a`, regardless that the starting index is 3.
 * 
 * _Usage_
 * 
 * - `sum<i_<N, P>>(f).expand<R>(...)` \f$\rightarrow\;\sum_{i_N + P = R}f\f$ is the sum
 * of the expression \f$f\f$ over the `N`th unique index, where the elipses, `...`, is the
 * substitution list.
 * 
 * \subsubsection series_examples Constructing Series
 *
 * Below are examples of using sum (a specialization of series which adds substituted expressions),
 * selected to show different cases of selecting indices and offsets. The examples assume the 
 * aliases: `using ii = expr::symbols::i_<0, 0>;` and `using v_ii = expr::symbols::v_<ii>;`
 *
 * **Summing the numbers 1 to 10**: 
 * \f[
 * \sum_{i = 0}^{10} i
 * \f]
 * 
 *     auto v = sum<ii>(ii{}).expand<1, 10>();
 * 
 * **Summing numbers starting i = -2**: 
 * \f[
 * \sum_{i = -2}^{2} i
 * \f]
 * 
 *     auto v = sum<ii>(ii{}).expand<-2, 2>();
 *
 * **Summing a list**
 * 
 * Assume that `data` is an existing variable that can be used in an expression. We will substitute
 * it into the sum \f$\sum_{i}i\,v_i\f$ where \f$i\f$ is the index and \f$v_i\f$ is the \f$i\f$th 
 * substituted function.
 * 
 *     auto op = expr::make_term(data);
 *     auto v = sum<ii>(ii{} * v_ii{}).expand<1>(op, op, op);
 *
 * This represents the sum starting at index \f$i=1\f$ up to \f$i=3\f$, and will result in `op`
 * with a coefficient of 6. _Note that the upper limit of the sum is inferred from the number of 
 * arguments and not passed explicitly._
 * 
 * **Summing a list with an offset**
 *
 * Extending the above, we can define the sum over: \f$i\,v_i\,v_{i+1}\f$. In this case it would
 * be useful to consider the example where 3 different terms are substituted:
 * 
 *     using v_iip1 = expr::symbols::v_<expr::symbols::i_<0, -1>>;
 *     auto op1 = expr::make_term<1>(data);
 *     auto op2 = expr::make_term<2>(data);
 *     auto op3 = expr::make_term<3>(data);
 *     auto v = sum<ii>(ii{} * v_ii{} * v_iip1{}).expand<1>(op1, op2, op3);
 * 
 * In this case, the sum limits will also be inferred using the offset of the substituted data. 
 * Since there is a substitution with an offset of 1, the sum will contain only 2 terms. The result
 * is equivalent to the expression: `expr::make_integer<1>()*op1*op2 + expr::make_integer<2>()*op2*op3`. 
 * 
 * \subsubsection series_iterated Iterated Series
 * 
 * Series can be defined using multiple indices, which represents an iterated series. In the
 * example of sums, defining the iterated sum: `sum<ii, jj, kk>(f)`, where `jj` and `kk` are defined
 * similarly to `ii` and **must** have a different unique number, represents:
 * \f[
 * \sum_{i} \sum_{j} \sum_{k} f
 * \f]
 * 
 * The substitution list can be different for placeholders of different index.
 * 
 * \subsubsection series_expansion Expanding Series
 * 
 * Series are expanded into an expression containing every term of the series expansion between 
 * the given limits. To expand a series, the function SymbolicSeries::expand is called with
 * specific template arguments defining the limits of the index.
 * 
 * The following are two examples using iterated series. The first example is a dynamic series
 * expansion, and the second is a direct expansion. As a reminder, a dynamic expansion generates
 * a function that evalutes the series in runtime bsaed loop. A direct expansion will write
 * out every term as a symbolic expression.
 * 
 *     auto rtsum = sum<ii, jj, kk>(ii{} + kk{})(
 *             expr::series_limits(1, 4), 
 *             expr::series_limits(ii{}, 3), 
 *             expr::series_limits(ii{}, jj{}))();
 * 
 *     auto rtsumv = sum<ii, jj, kk>(ii{} + kk{})
 *             .expand<expr::series_limits_ntn<1, 4>, 
 *                     expr::series_limits_itn<ii, 3>, 
 *                     expr::series_limits_iti<ii, jj>>();
 * 
 * The above examples also demonstrate how to write the limits of the series. Limits can be
 * defined between indices. In the second example, compile-time constant types are used to 
 * indicate the limits, but may be cumbersome to write.
 * An alternative to this method for expanding the sum directly is to use the `select` function, 
 * and then simply call the result with the substitution list or tuple:
 * 
 *     auto rtsum = sum(ii{} + jj{} * kk{} * v_jj{})
 *            .select<10, 10, 10>(
 *                    (ii{} = one) && (ii{} = expr::val<3>),  
 *                    (jj{} = ii{} + one) && (jj{} = expr::val<6>),
 *                    (kk{} = ii{} + jj{}) && (kk{} = expr::val<9>));
 *     auto expanded = rtsum(ops);
 * 
 * The number of elements iterated for each index are specified as template arguments to `select`.
 * The limits are specified through usual algebaic operations (equals, plus, etc), and the lower
 * and top limits for each index are specified for each variable by combining them with the `&&` 
 * operator. 
 * 
 * > **Important**: The order of the indices in the limit specification must correspond to the order
 * > of the indices that is desired for the summation. Additionally, even if the indices are manually 
 * > specified to `sum`, the order given in select will still be used instead.
 */



namespace expr
{
	/*
	template<>
	struct series_limits<int, int>
	{
		series_limits(int _0, int _1) : _0{ _0 }, _1{ _1 } {}
		series_limits(OpVoid, int _1) : _0{ 0 }, _1{ _1 } {}
		series_limits(int _0, OpVoid) : _0{ _0 }, _1{ 0 } {}
		series_limits(OpVoid, OpVoid) : _0{ 0 }, _1{ 0 } {}
		int _0;
		int _1;
	};

	template<typename S>
	struct series_limits<int, OpTerm<OpIdentity, S>>
	{
		series_limits(int _0, S const& _1) :
			_0{ _0 }, _1{ _1 }, _1c{ 1 } {}

		template<typename T>
		series_limits(int _0, OpTerm<T, S> const& term) :
			_0{ _0 }, _1{ expr::get<1>(term).data() }, _1c{ expr::coeff(term) } {}

		int _0;
		S _1;
		int _1c;
	};

	template<typename S>
	struct series_limits<OpTerm<OpIdentity, S>, int>
	{
		series_limits(S const& _0, int _1) :
			_0{ _0 }, _0c{ 1 }, _1{ _1 } {}

		template<typename T>
		series_limits(OpTerm<T, S> const& term, int _1) :
			_0{ expr::get<1>(term).data() }, _0c{ expr::coeff(term) }, _1{ _1 } {}

		S _0;
		int _0c;
		int _1;
	};*/

	template<typename E1, typename E2>
	struct series_limits
	{
		series_limits() = default;

		series_limits(E1 const& _0, E2 const& _1) :
			_0{ _0 }, _1{ _1 } {}

		E1 _0;
		E2 _1;
	};

	template<typename S0, typename S1>
	struct series_limits<OpTerm<OpIdentity, S0>, OpTerm<OpIdentity, S1>>
	{
		series_limits() = default;

		series_limits(S0 const& _0, S1 const& _1) :
			_0{ _0 }, _0c{ 1 },
			_1{ _1 }, _1c{ 1 } {}

		template<typename T0, typename T1>
		series_limits(
			OpTerm<T0, S0> const& term0,
			OpTerm<T1, S1> const& term1) :
			_0{ expr::get<1>(term0).data() }, _0c{ double(expr::coeff(term0)) },
			_1{ expr::get<1>(term1).data() }, _1c{ double(expr::coeff(term1)) } {}

		template<typename T1>
		series_limits(
			S0 const& _0,
			OpTerm<T1, S1> const& term1) :
			_0{ _0 }, _0c{ 1 },
			_1{ expr::get<1>(term1).data() }, _1c{ double(expr::coeff(term1)) } {}

		template<typename T0>
		series_limits(
			OpTerm<T0, S0> const& term0,
			S1 const& _1) :
			_0{ expr::get<1>(term0).data() }, _0c{ double(expr::coeff(term0)) },
			_1{ _1 }, _1c{ 1 } {}

		S0 _0;
		double _0c;
		S1 _1;
		double _1c;
	};


	template<typename E0, typename S>
	struct series_limits<E0, OpTerm<OpIdentity, S>>
	{
		series_limits() = default;

		series_limits(E0 const& _0, S const& _1) :
			_0{ _0 }, _1{ _1 }, _1c{ 1 } {}

		template<typename T>
		series_limits(E0 const& _0, OpTerm<T, S> const& term) :
			_0{ _0 }, _1{ expr::get<1>(term).data() }, _1c{ double(expr::coeff(term)) } {}

		E0 _0;
		S _1;
		double _1c;
	};

	template<typename S, typename E1>
	struct series_limits<OpTerm<OpIdentity, S>, E1>
	{
		series_limits() = default;

		series_limits(S const& _0, E1 const& _1) :
			_0{ _0 }, _0c{ 1 }, _1{ _1 } {}

		template<typename T>
		series_limits(OpTerm<T, S> const& term, E1 const& _1) :
			_0{ expr::get<1>(term).data() }, _0c{ double(expr::coeff(term)) }, _1{ _1 } {}

		S _0;
		double _0c;
		E1 _1;
	};



	template<typename T, typename U, std::enable_if_t<(is_simple_data<U> && !is_simple_data<T>), int> = 0>
	series_limits(T, U) -> series_limits<T, int>;
	template<typename T, typename U, std::enable_if_t<(is_simple_data<U> && !is_simple_data<T>), int> = 0>
	series_limits(U, T) -> series_limits<int, T>;
	template<typename U1, typename U2, std::enable_if_t<(is_simple_data<U1> && is_simple_data<U2>), int> = 0>
	series_limits(U1, U2)->series_limits<int, int>;
	template<typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(U, U) -> series_limits<int, int>;

	template<typename V, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(U, OpLiteral<V>) -> series_limits<int, int>;
	template<typename V, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(OpLiteral<V>, U) -> series_limits<int, int>;

	template<typename T, typename V, std::enable_if_t<!is_simple_data<T>, int> = 0>
	series_limits(T, OpLiteral<V>) -> series_limits<T, int>;
	template<typename T, typename V, std::enable_if_t<!is_simple_data<T>, int> = 0>
	series_limits(OpLiteral<V>, T) -> series_limits<int, T>;
	template<typename V0, typename V1>
	series_limits(OpLiteral<V0>, OpLiteral<V1>) -> series_limits<int, int>;

	template<typename T, typename S, typename V, std::enable_if_t<!is_simple_data<T>, int> = 0>
	series_limits(OpTerm<V, S>, T) -> series_limits<OpTerm<OpIdentity, S>, T>;
	template<typename T, typename S, typename V, std::enable_if_t<!is_simple_data<T>, int> = 0>
	series_limits(T, OpTerm<V, S>) -> series_limits<T, OpTerm<OpIdentity, S>>;
	template<typename T, typename S, typename V>
	series_limits(OpTerm<T, S>, OpLiteral<V>) -> series_limits<OpTerm<OpIdentity, S>, int>;
	template<typename T, typename S, typename V>
	series_limits(OpLiteral<V>, OpTerm<T, S>) -> series_limits<int, OpTerm<OpIdentity, S>>;
	template<typename T, typename S, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(OpTerm<T, S>, U) -> series_limits<OpTerm<OpIdentity, S>, int>;
	template<typename T, typename S, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(U, OpTerm<T, S>) -> series_limits<int, OpTerm<OpIdentity, S>>;
	template<typename T0, typename S0, typename T1, typename S1>
	series_limits(OpTerm<T0, S0>, OpTerm<T1, S1>) -> series_limits<OpTerm<OpIdentity, S0>, OpTerm<OpIdentity, S1>>;
	template<int N, int P, typename S, typename V>
	series_limits(OpTerm<V, S>, expr::symbols::i_<N, P>) -> series_limits<OpTerm<OpIdentity, S>, OpTerm<OpIdentity, expr::symbols::i_<N, P>>>;
	template<int N, int P, typename S, typename V>
	series_limits(expr::symbols::i_<N, P>, OpTerm<V, S>) -> series_limits<OpTerm<OpIdentity, expr::symbols::i_<N, P>>, OpTerm<OpIdentity, S>>;




	template<int N, int P, typename T, std::enable_if_t<!is_simple_data<T>, int> = 0>
	series_limits(expr::symbols::i_<N, P>, T) -> series_limits<OpTerm<OpIdentity, expr::symbols::i_<N, P>>, T>;
	template<int N, int P, typename T, std::enable_if_t<!is_simple_data<T>, int> = 0>
	series_limits(T, expr::symbols::i_<N, P>) -> series_limits<T, OpTerm<OpIdentity, expr::symbols::i_<N, P>>>;
	template<int N, int P, typename V>
	series_limits(expr::symbols::i_<N, P>, OpLiteral<V>) -> series_limits<OpTerm<OpIdentity, expr::symbols::i_<N, P>>, int>;
	template<int N, int P, typename V>
	series_limits(OpLiteral<V>, expr::symbols::i_<N, P>) -> series_limits<int, OpTerm<OpIdentity, expr::symbols::i_<N, P>>>;
	template<int N, int P, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(expr::symbols::i_<N, P>, U) -> series_limits<OpTerm<OpIdentity, expr::symbols::i_<N, P>>, int>;
	template<int N, int P, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(U, expr::symbols::i_<N, P>) -> series_limits<int, OpTerm<OpIdentity, expr::symbols::i_<N, P>>>;
	template<int N0, int P0, int N1, int P1>
	series_limits(expr::symbols::i_<N0, P0>, expr::symbols::i_<N1, P1>)
		-> series_limits<OpTerm<OpIdentity, expr::symbols::i_<N0, P0>>, OpTerm<OpIdentity, expr::symbols::i_<N1, P1>>>;

	template<size_t N, typename T, std::enable_if_t<!is_simple_data<T>, int> = 0>
	series_limits(expr::symbols::placeholder_N_symbol_<N>, T)
		-> series_limits<OpTerm<OpIdentity, expr::symbols::placeholder_N_symbol_<N>>, T>;
	template<size_t N, typename T, std::enable_if_t<!is_simple_data<T>, int> = 0>
	series_limits(T, expr::symbols::placeholder_N_symbol_<N>)
		-> series_limits<T, OpTerm<OpIdentity, expr::symbols::placeholder_N_symbol_<N>>>;
	template<size_t N, typename V>
	series_limits(expr::symbols::placeholder_N_symbol_<N>, OpLiteral<V>)
		-> series_limits<OpTerm<OpIdentity, expr::symbols::placeholder_N_symbol_<N>>, int>;
	template<size_t N, typename V>
	series_limits(OpLiteral<V>, expr::symbols::placeholder_N_symbol_<N>)
		-> series_limits<int, OpTerm<OpIdentity, expr::symbols::placeholder_N_symbol_<N>>>;
	template<size_t N, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(expr::symbols::placeholder_N_symbol_<N>, U)
		-> series_limits<OpTerm<OpIdentity, expr::symbols::placeholder_N_symbol_<N>>, int>;
	template<size_t N, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(U, expr::symbols::placeholder_N_symbol_<N>)
		-> series_limits<int, OpTerm<OpIdentity, expr::symbols::placeholder_N_symbol_<N>>>;
	template<size_t N0, size_t N1>
	series_limits(expr::symbols::placeholder_N_symbol_<N0>, expr::symbols::placeholder_N_symbol_<N1>)
		-> series_limits<OpTerm<OpIdentity, expr::symbols::placeholder_N_symbol_<N0>>, OpTerm<OpIdentity, expr::symbols::placeholder_N_symbol_<N1>>>;


	template<int N0, int P0, size_t N1>
	series_limits(expr::symbols::i_<N0, P0>, expr::symbols::placeholder_N_symbol_<N1>)
		-> series_limits<OpTerm<OpIdentity, expr::symbols::i_<N0, P0>>, OpTerm<OpIdentity, expr::symbols::placeholder_N_symbol_<N1>>>;
	template<size_t N0, int N1, int P1>
	series_limits(expr::symbols::placeholder_N_symbol_<N0>, expr::symbols::i_<N1, P1>)
		-> series_limits<OpTerm<OpIdentity, expr::symbols::placeholder_N_symbol_<N0>>, OpTerm<OpIdentity, expr::symbols::i_<N1, P1>>>;


	template<size_t N, typename T2, std::enable_if_t<!is_simple_data<T2>, int> = 0>
	series_limits(OpFractionLiteral<N, 1>, T2) -> series_limits<int, T2>;
	template<size_t N, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(OpFractionLiteral<N, 1>, U) -> series_limits<int, int>;
	template<size_t N, typename V>
	series_limits(OpFractionLiteral<N, 1>, OpLiteral<V>) -> series_limits<int, int>;
	template<typename T1, size_t N, std::enable_if_t<!is_simple_data<T1>, int> = 0>
	series_limits(T1, OpFractionLiteral<N, 1>) -> series_limits<T1, int>;
	template<size_t N, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(U, OpFractionLiteral<N, 1>) -> series_limits<int, int>;
	template<size_t N, typename V>
	series_limits(OpLiteral<V>, OpFractionLiteral<N, 1>) -> series_limits<int, int>;
	template<size_t N, typename T2, std::enable_if_t<!is_simple_data<T2>, int> = 0>
	series_limits(OpNegFractionLiteral<N, 1>, T2) -> series_limits<int, T2>;
	template<size_t N, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(OpNegFractionLiteral<N, 1>, U) -> series_limits<int, int>;
	template<size_t N, typename V>
	series_limits(OpNegFractionLiteral<N, 1>, OpLiteral<V>) -> series_limits<int, int>;
	template<typename T1, size_t N, std::enable_if_t<!is_simple_data<T1>, int> = 0>
	series_limits(T1, OpNegFractionLiteral<N, 1>) -> series_limits<T1, int>;
	template<size_t N, typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(U, OpNegFractionLiteral<N, 1>) -> series_limits<int, int>;
	template<size_t N, typename V>
	series_limits(OpLiteral<V>, OpNegFractionLiteral<N, 1>) -> series_limits<int, int>;
	template<typename T2, std::enable_if_t<!is_simple_data<T2>, int> = 0>
	series_limits(OpIdentity, T2) -> series_limits<int, T2>;
	template<typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(OpIdentity, U)->series_limits<int, int>;
	template<typename V>
	series_limits(OpIdentity, OpLiteral<V>) -> series_limits<int, int>;
	template<typename T1, std::enable_if_t<!is_simple_data<T1>, int> = 0>
	series_limits(T1, OpIdentity) -> series_limits<T1, int>;
	template<typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(U, OpIdentity)->series_limits<int, int>;
	template<typename V>
	series_limits(OpLiteral<V>, OpIdentity) -> series_limits<int, int>;
	template<typename T2, std::enable_if_t<!is_simple_data<T2>, int> = 0>
	series_limits(OpNegIdentity, T2) -> series_limits<int, T2>;
	template<typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(OpNegIdentity, U)->series_limits<int, int>;
	template<typename V>
	series_limits(OpNegIdentity, OpLiteral<V>) -> series_limits<int, int>;
	template<typename T1, std::enable_if_t<!is_simple_data<T1>, int> = 0>
	series_limits(T1, OpNegIdentity) -> series_limits<T1, int>;
	template<typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(U, OpNegIdentity)->series_limits<int, int>;
	template<typename V>
	series_limits(OpLiteral<V>, OpNegIdentity) -> series_limits<int, int>;
	template<typename T2, std::enable_if_t<!is_simple_data<T2>, int> = 0>
	series_limits(OpVoid, T2) -> series_limits<int, T2>;
	template<typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(OpVoid, U)->series_limits<int, int>;
	template<typename V>
	series_limits(OpVoid, OpLiteral<V>) -> series_limits<int, int>;
	template<typename T1, std::enable_if_t<!is_simple_data<T1>, int> = 0>
	series_limits(T1, OpVoid) -> series_limits<T1, int>;
	template<typename U, std::enable_if_t<is_simple_data<U>, int> = 0>
	series_limits(U, OpVoid)->series_limits<int, int>;
	template<typename V>
	series_limits(OpLiteral<V>, OpVoid) -> series_limits<int, int>;



	template<typename T, typename V, typename G>
	auto limit_0(expr::series_limits<OpTerm<V, G>, T> const& limit)
	{
		return limit._0c * limit._0;
	}

	template<typename T1, typename T2>
	auto limit_0(expr::series_limits<T1, T2> const& limit)
	{
		return limit._0;
	}

	template<typename T, typename V, typename G>
	auto limit_1(expr::series_limits<T, OpTerm<V, G>> const& limit)
	{
		return limit._1c * limit._1;
	}

	template<typename T1, typename T2>
	auto limit_1(expr::series_limits<T1, T2> const& limit)
	{
		return limit._1;
	}


	//! Defines the limits of the series expansion.
	/*!
	 * Represents a limit from a start and end index (end index is included) for expanding
	 * the series. This is particularly relevant when expanding a series in terms of functions,
	 * and only a specific range of the list of functions should be substituted.
	 *
	 * Instead of passing a value to the template parameters of the function SymbolicSeries::expand,
	 * this object is used. Likewise when using SymbolicSeries::select, except it is passed
	 * as a regular function argument.
	 *
	 * \tparam N0 The start index.
	 * \tparam N1 The end index, which is included in the expansion of the series.
	 */
	template<int N0, int N1>
	using series_limits_ntn = std::integer_sequence<int, N0, N1>;

	//! Index between number and index.
	template<int N0, typename I>
	using series_limits_nti = symphas::lib::types_list<std::integer_sequence<int, N0>, I>;

	//! Limit between index and number.
	template<typename I, int N1>
	using series_limits_itn = symphas::lib::types_list<I, std::integer_sequence<int, N1>>;

	//! Limit between indices.
	template<typename I0, typename I1>
	using series_limits_iti = symphas::lib::types_list<I0, I1>;


	template<typename S, int N1>
	using series_limits_stn = symphas::lib::types_list<S, std::integer_sequence<int, N1>>;

	template<typename S, typename I1>
	using series_limits_sti = symphas::lib::types_list<S, I1>;

	template<int N0, typename S>
	using series_limits_nts = symphas::lib::types_list<std::integer_sequence<int, N0>, S>;

	template<typename I0, typename S>
	using series_limits_its = symphas::lib::types_list<I0, S>;


	template<typename S0, typename S1>
	using series_limits_sts = symphas::lib::types_list<S0, S1>;

	//! Indicates to go until the end of the list; computed when series is expanded.
	template<typename I0>
	using series_limits_ite = symphas::lib::types_list<I0, void>;

	//! Indicates to go until the end of the list; computed when series is expanded.
	template<int N0>
	using series_limits_nte = symphas::lib::types_list<std::integer_sequence<int, N0>, void>;

	//! Indicates to go until the end of the list; computed when series is expanded.
	template<typename S>
	using series_limits_ste = symphas::lib::types_list<S, void>;


	template<int N, typename I, typename limits_t>
	struct update_limits_impl
	{
		using type = limits_t;
	};

	template<int N, int I0, int P0, int P1, int N0>
	struct update_limits_impl<N, expr::symbols::i_<I0, P0>, series_limits_nti<N0, expr::symbols::i_<I0, P1>>>
	{
		using type = series_limits_ntn<N0, N + P1>;
	};

	template<int N, int I0, int P0, int P1, int N1>
	struct update_limits_impl<N, expr::symbols::i_<I0, P0>, series_limits_itn<expr::symbols::i_<I0, P1>, N1>>
	{
		using type = series_limits_ntn<N + P1, N1>;
	};

	template<int N, int I0, int P0, int P1, int I1, int P2>
	struct update_limits_impl<N, expr::symbols::i_<I0, P0>, series_limits_iti<expr::symbols::i_<I0, P1>, expr::symbols::i_<I1, P2>>>
	{
		using type = series_limits_nti<N + P1, expr::symbols::i_<I1, P2>>;
	};

	template<int N, int I0, int P0, int P1, int I1, int P2>
	struct update_limits_impl<N, expr::symbols::i_<I0, P0>, series_limits_iti<expr::symbols::i_<I1, P1>, expr::symbols::i_<I0, P2>>>
	{
		using type = series_limits_itn<expr::symbols::i_<I1, P1>, N + P2>;
	};

	template<int N, int I0, int P0, int P1, int P2>
	struct update_limits_impl<N, expr::symbols::i_<I0, P0>, series_limits_iti<expr::symbols::i_<I0, P1>, expr::symbols::i_<I0, P2>>>
	{
		using type = series_limits_ntn<N + P1, N + P2>;
	};
	
	template<int N, int I0, int... P0s, typename... Es>
	auto update_limit_sum(OpAdd<Es...> const& add, symphas::lib::types_list<expr::symbols::i_<I0, P0s>...>)
	{
		return expr::transform::swap_grid<expr::symbols::i_<I0, P0s>...>(
			add, expr::make_integer<N + P0s>()...);
	}

	template<int N, int I0, typename... Es>
	auto update_limit_sum(OpAdd<Es...> const& add)
	{
		using I = expr::symbols::i_<I0, 0>;
		using indices = typename symphas::internal::select_all_i_<I, expr::op_types_t<OpAdd<Es...>>>::type;
		return update_limit_sum<N>(add, indices{});
	}


	template<typename T>
	inline constexpr int literal_value = T{};
	template<>
	inline constexpr int literal_value<OpIdentity> = 1;
	template<>
	inline constexpr int literal_value<OpNegIdentity> = -1;
	template<size_t N>
	constexpr int literal_value<OpFractionLiteral<N, 1>> = int(N);
	template<size_t N>
	constexpr int literal_value<OpNegFractionLiteral<N, 1>> = -int(N);

	template<int N, typename I, int I0, int P0, typename S, int N1>
	auto update_limit_sum_result(OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, series_limits_stn<S, N1>)
	{
		return typename update_limits_impl<N, I,
			series_limits_itn<expr::symbols::i_<I0, P0>, N1>>::type{};
	}
	template<int N, typename I, int I0, int P0, typename S, int I1, int P1>
	auto update_limit_sum_result(OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, series_limits_sti<S, expr::symbols::i_<I1, P1>>)
	{
		return typename update_limits_impl<N, I,
			series_limits_iti<expr::symbols::i_<I0, P0>, expr::symbols::i_<I1, P1>>>::type{};
	}
	template<int N, typename I, int I0, int P0, typename S, int N0>
	auto update_limit_sum_result(OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, series_limits_nts<N0, S>)
	{
		return typename update_limits_impl<N, I,
			series_limits_nti<N0, expr::symbols::i_<I0, P0>>>::type{};
	}
	template<int N, typename I, int I1, int P1, typename S, int I0, int P0>
	auto update_limit_sum_result(OpTerm<OpIdentity, expr::symbols::i_<I1, P1>>, series_limits_its<expr::symbols::i_<I0, P0>, S>)
	{
		return typename update_limits_impl<N, I,
			series_limits_iti<expr::symbols::i_<I0, P0>, expr::symbols::i_<I1, P1>>>::type{};
	}
	template<int N, typename I, typename... Es, typename S, int N1>
	auto update_limit_sum_result(OpAdd<Es...>, series_limits_stn<S, N1>)
	{
		return series_limits_stn<OpAdd<Es...>, N1>{};
	}
	template<int N, typename I, typename... Es, typename S, int I1, int P1>
	auto update_limit_sum_result(OpAdd<Es...>, series_limits_sti<S, expr::symbols::i_<I1, P1>>)
	{
		return series_limits_sti<OpAdd<Es...>, expr::symbols::i_<I1, P1>>{};
	}
	template<int N, typename I, typename... Es, typename S, int N0>
	auto update_limit_sum_result(OpAdd<Es...>, series_limits_nts<N0, S>)
	{
		return series_limits_nts<N0, OpAdd<Es...>>{};
	}
	template<int N, typename I, typename... Es, typename S, int I0, int P0>
	auto update_limit_sum_result(OpAdd<Es...>, series_limits_its<expr::symbols::i_<I0, P0>, S>)
	{
		return series_limits_its<expr::symbols::i_<I0, P0>, OpAdd<Es...>>{};
	}
	template<int N, typename I, typename T, typename S, int N1>
	auto update_limit_sum_result(T, series_limits_stn<S, N1>)
	{
		return series_limits_ntn<literal_value<T>, N1>{};
	}
	template<int N, typename I, typename T, typename S, int I1, int P1>
	auto update_limit_sum_result(T, series_limits_sti<S, expr::symbols::i_<I1, P1>>)
	{
		return typename update_limits_impl<N, I,
			series_limits_nti<literal_value<T>, expr::symbols::i_<I1, P1>>>::type{};
	}
	template<int N, typename I, typename T, typename S, int N0>
	auto update_limit_sum_result(T, series_limits_nts<N0, S>)
	{
		return series_limits_ntn<N0, literal_value<T>>{};
	}
	template<int N, typename I, typename T, typename S, int I0, int P0>
	auto update_limit_sum_result(T, series_limits_its<expr::symbols::i_<I0, P0>, S>)
	{
		return typename update_limits_impl<N, I,
			series_limits_itn<expr::symbols::i_<I0, P0>, literal_value<T>>>::type{};
	}



	template<int N, int I0, int P0, typename... Es, int N1>
	struct update_limits_impl<N, expr::symbols::i_<I0, P0>, series_limits_stn<OpAdd<Es...>, N1>>
	{
	protected:
		using result_t = decltype(update_limit_sum<N, I0>(OpAdd<Es...>{}));

	public:
		using type = decltype(update_limit_sum_result<N, expr::symbols::i_<I0, P0>>
			(result_t{}, series_limits_stn<OpAdd<Es...>, N1>{}));
	};

	template<int N, int I0, int P0, int N0, typename... Es>
	struct update_limits_impl<N, expr::symbols::i_<I0, P0>, series_limits_nts<N0, OpAdd<Es...>>>
	{
	protected:
		using result_t = decltype(update_limit_sum<N, expr::symbols::i_<I0, P0>>(OpAdd<Es...>{}));

	public:
		using type = decltype(update_limit_sum_result<N, expr::symbols::i_<I0, P0>>
			(result_t{}, series_limits_nts<N0, OpAdd<Es...>>{}));
	};

	template<int N, int I0, int P0, typename... Es, int I1, int P1>
	struct update_limits_impl<N, expr::symbols::i_<I0, P0>, series_limits_sti<OpAdd<Es...>, expr::symbols::i_<I1, P1>>>
	{
	protected:
		using result_t = decltype(update_limit_sum<N, I0>(OpAdd<Es...>{}));

	public:
		using type = decltype(update_limit_sum_result<N, expr::symbols::i_<I0, P0>>
			(result_t{}, series_limits_sti<OpAdd<Es...>, expr::symbols::i_<I1, P1>>{}));
	};

	template<int N, int I0, int P0, typename... Es, int P1>
	struct update_limits_impl<N, expr::symbols::i_<I0, P0>, series_limits_sti<OpAdd<Es...>, expr::symbols::i_<I0, P1>>>
	{
	protected:
		using result_t = decltype(update_limit_sum<N, I0>(OpAdd<Es...>{}));

	public:
		using type = decltype(update_limit_sum_result<N, expr::symbols::i_<I0, P0>>
			(result_t{}, series_limits_stn<OpAdd<Es...>, N + P1>{}));
	};

	template<int N, int I0, int P0, int I1, int P1, typename... Es>
	struct update_limits_impl<N, expr::symbols::i_<I0, P0>, series_limits_its<expr::symbols::i_<I1, P1>, OpAdd<Es...>>>
	{
	protected:
		using result_t = decltype(update_limit_sum<N, I0>(OpAdd<Es...>{}));

	public:
		using type = decltype(update_limit_sum_result<N, expr::symbols::i_<I0, P0>>
			(result_t{}, series_limits_its<expr::symbols::i_<I1, P1>, OpAdd<Es...>>{}));
	};

	template<int N, int I0, int P0, int P1, typename... Es>
	struct update_limits_impl<N, expr::symbols::i_<I0, P0>, series_limits_its<expr::symbols::i_<I0, P1>, OpAdd<Es...>>>
	{
	protected:
		using result_t = decltype(update_limit_sum<N, I0>(OpAdd<Es...>{}));

	public:
		using type = decltype(update_limit_sum_result<N, expr::symbols::i_<I0, P0>>
			(result_t{}, series_limits_nts<N + P1, OpAdd<Es...>>{}));
	};


	template<int N, typename I, typename limits_t>
	using update_limits = typename update_limits_impl<N, I, limits_t>::type;


	namespace
	{
		template<typename I0, typename I1>
		series_limits_iti<I0, I1> construct_limits()
		{
			return {};
		}

		template<typename I0, int N>
		series_limits_itn<I0, N> construct_limits()
		{
			return {};
		}

		template<int N, typename I1>
		series_limits_nti<N, I1> construct_limits()
		{
			return {};
		}

		template<int N0, int N1>
		series_limits_ntn<N0, N1> construct_limits()
		{
			return {};
		}
	}
}


// Defines cases when data are substituted into arguments.


template<typename... T1s, typename... T2s>
struct SymbolicCase<expr::case_entry<T1s, T2s>...> : expr::symbols::Symbol
{
	using parent_type = expr::symbols::Symbol;

	SymbolicCase(std::pair<T1s, T2s> const&...) : parent_type() {}
	SymbolicCase() : parent_type() {}

	std::tuple<T2s...> cases;

};

template<int N0, int P0, typename G, typename A, typename B>
struct SymbolicCase<
	expr::case_entry<expr::symbols::index_eq<expr::symbols::i_<N0, P0>, G>, A>,
	expr::case_entry<expr::symbols::index_neq<expr::symbols::i_<N0, P0>, G>, B>> : expr::symbols::Symbol
{
	using parent_type = expr::symbols::Symbol;

	SymbolicCase(expr::symbols::index_eq<expr::symbols::i_<N0, P0>, G>, A left, B right) : 
		parent_type(), cases{ left, right } {}
	SymbolicCase() : parent_type(), cases{} {}

	std::tuple<A, B> cases;
};

template<int N0, int P0, int N, typename A, typename B>
struct SymbolicCase<
	expr::case_entry<expr::symbols::index_eq_N<expr::symbols::i_<N0, P0>, N>, A>,
	expr::case_entry<expr::symbols::index_neq_N<expr::symbols::i_<N0, P0>, N>, B>> : expr::symbols::Symbol
{
	using parent_type = expr::symbols::Symbol;

	SymbolicCase(expr::symbols::index_eq_N<expr::symbols::i_<N0, P0>, N>, A left, B right) :
		parent_type(), cases{ left, right } {}
	SymbolicCase() : parent_type(), cases{} {}

	std::tuple<A, B> cases;
};

template<typename C>
struct SymbolicCase<expr::case_entry<C, void>> : expr::symbols::Symbol
{
	using parent_type = expr::symbols::Symbol;

	SymbolicCase(C = C{}) : parent_type() {}
};

template<typename... T1s, typename... T2s>
SymbolicCase(std::pair<T1s, T2s>...) -> SymbolicCase<expr::case_entry<T1s, T2s>...>;

// For SymbolicTernaryCase
template<int N0, int P0, typename G, typename A, typename B>
SymbolicCase(expr::symbols::index_eq<expr::symbols::i_<N0, P0>, G>, A, B) -> SymbolicCase<
	expr::case_entry<expr::symbols::index_eq<expr::symbols::i_<N0, P0>, G>, A>,
	expr::case_entry<expr::symbols::index_neq<expr::symbols::i_<N0, P0>, G>, B>>;
template<int N0, int P0, int N, typename A, typename B>
SymbolicCase(expr::symbols::index_eq_N<expr::symbols::i_<N0, P0>, N>, A, B) -> SymbolicCase<
	expr::case_entry<expr::symbols::index_eq_N<expr::symbols::i_<N0, P0>, N>, A>,
	expr::case_entry<expr::symbols::index_neq_N<expr::symbols::i_<N0, P0>, N>, B>>;

// For SymbolicCaseSwap
template<typename C>
SymbolicCase(C) -> SymbolicCase<expr::case_entry<C, void>>;

ALLOW_COMBINATION((typename... Ts), (SymbolicCase<Ts...>))

namespace symphas::internal
{

	template<typename limits_t>
	constexpr int series_limits_0 = 1;
	template<int N0, int N1>
	constexpr int series_limits_0<expr::series_limits_ntn<N0, N1>> = N0;
	template<typename limits_t>
	constexpr int series_limits_1 = 0;
	template<int N0, int N1>
	constexpr int series_limits_1<expr::series_limits_ntn<N0, N1>> = N1;

	template<typename limits_t>
	struct series_limits_of
	{
		static const int _0 = series_limits_0<limits_t>;
		static const int _1 = series_limits_1<limits_t>;
	};


	template<typename X>
	struct filter_i_selection
	{
		using type = X;
	};

	template<int N, int P>
	struct filter_i_selection<expr::symbols::i_<N, P>>
	{
		using type = expr::symbols::index_eq_N<expr::symbols::i_<N, P>, 1>;
	};


	using namespace expr::symbols;
	using namespace expr;

	template<typename T>
	len_type limit_length(T)
	{
		return T{};
		//static_assert(false, "limits can only be inferred from constant length " 
		//	"arrays or Substitution-type arguments");
	}

	template<typename G>
	len_type limit_length(DynamicVariable<G> const& sub)
	{
		return sub.end() - sub.start() + 1;
	}

	template<typename V, typename G>
	len_type limit_length(OpTerm<V, DynamicVariable<G>> const& sub)
	{
		auto const& [v, t] = sub;
		return limit_length(t);
	}

	template<typename T, size_t N>
	len_type limit_length(T (&)[N])
	{
		return N;
	}

	template<typename... Ts>
	len_type limit_length(std::tuple<Ts...> const& sub)
	{
		return sizeof...(Ts);
	}

	template<typename T>
	len_type limit_length(SymbolicDataArray<T> const& sub)
	{
		return sub.length();
	}

	template<int N0, int N1, typename sub_t>
	auto parse_limit(series_limits_ntn<N0, N1>, sub_t&& sub)
	{
		return expr::series_limits(N0, N1);
	}

	template<int N0, typename I, typename sub_t>
	auto parse_limit(series_limits_nti<N0, I>, sub_t&& sub)
	{
		return expr::series_limits(N0, I{});
	}

	template<typename I, int N1, typename sub_t>
	auto parse_limit(series_limits_itn<I, N1>, sub_t&& sub)
	{
		return expr::series_limits(I{}, N1);
	}

	template<typename I0, typename I1, typename sub_t>
	auto parse_limit(series_limits_iti<I0, I1>, sub_t&& sub)
	{
		return expr::series_limits(I0{}, I1{});
	}

	template<int N0, typename sub_t>
	auto parse_limit(series_limits_nte<N0>, sub_t&& sub)
	{
		return expr::series_limits(N0, limit_length(std::forward<sub_t>(sub)) + N0 - 1);
	}

	template<int N0, int P0, typename sub_t>
	auto parse_limit(series_limits_ite<expr::symbols::i_<N0, P0>>, sub_t&& sub)
	{
		return expr::series_limits(expr::symbols::i_<N0, P0>{}, limit_length(std::forward<sub_t>(sub)));
	}

	template<typename S, typename sub_t>
	auto parse_limit(series_limits_ste<S>, sub_t&& sub)
	{
		return expr::series_limits(S{}, limit_length(std::forward<sub_t>(sub)));
	}


	template<typename Seq, typename series_t, typename... Rest>
	struct expand_series;

	template<size_t N0, size_t... Ns, typename series_t, typename... limit_ts>
	struct expand_series<
		std::index_sequence<N0, Ns...>, series_t, types_list<limit_ts...>> :
		expand_series<std::index_sequence<N0, Ns...>, series_t, types_list<>, types_list<limit_ts...>> {};
	
	//! Expand the sum where all indices are set to specific values.
	template<typename... I0s, typename... limit_ts, typename E, typename Op>
	struct expand_series<
		std::index_sequence<>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>, types_list<>>
	{
		using all_types = types_list<
			expand_series<std::index_sequence<>, SymbolicSeries<Op, E, types_list<I0s...>>,
			types_list<limit_ts...>, types_list<>>>;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return SymbolicSeries<Op, E, types_list<I0s...>>(*static_cast<E const*>(&e))
				(parse_limit(limit_ts{}, std::forward<Ts>(subs))...)(std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return SymbolicSeries<Op, E, types_list<I0s...>>(*static_cast<E const*>(&e))
				.template expand<limit_ts...>(std::forward<Ts>(subs)...);
		}
	};


	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int I0, int P0, int I1, int P1, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<index_neq<i_<I1, P1>, i_<I0, P0>>, Xs...>>
	{
		using second_limit_t = std::conditional_t<
			(N0 > 0),
			expr::series_limits_itn<i_<I0, P0 + 1>, N0>,
			expr::series_limits_ite<i_<I0, P0 + 1>>>;

		using first_t = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>, 
			types_list<limit_ts..., expr::series_limits_nti<1, i_<I0, P0 - 1>>>, types_list<Xs...>>;
		using second_t = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>,
			types_list<limit_ts..., second_limit_t>, types_list<Xs...>>;

		using all_types = expand_types_list<typename first_t::all_types, typename second_t::all_types>;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return Op{}(
				first_t{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...), 
				second_t{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...));
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return Op{}.expand(
				first_t{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...), 
				second_t{}.expand(*static_cast<E const*>(&e), std::forward<Ts>(subs)...));
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int I0, int P0, int I1, int P1, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<index_eq<i_<I1, P1>, i_<I0, P0>>, Xs...>>
	{
		using second_limit_t = std::conditional_t<
			(N0 > 0),
			expr::series_limits_itn<i_<I0, P0>, N0>,
			expr::series_limits_ite<i_<I0, P0>>>;

		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>, 
			types_list<limit_ts..., second_limit_t>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};


	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int I0, int P0, int M0, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<index_neq_N<i_<I0, P0>, M0>, Xs...>>
	{
		using second_limit_t = std::conditional_t<
			(N0 > 0),
			expr::series_limits_ntn<M0 + 1, N0>,
			expr::series_limits_nte<M0 + 1>>;

		using first_t = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I0, P0>>>,
			types_list<limit_ts..., expr::series_limits_ntn<1, M0 - 1>>, types_list<Xs...>>;
		using second_t = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I0, P0>>>,
			types_list<limit_ts..., second_limit_t>, types_list<Xs...>>;

		using all_types = expand_types_list<typename first_t::all_types, typename second_t::all_types>;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return Op{}(
				first_t{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...),
				second_t{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...));
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return Op{}.expand(
				first_t{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...),
				second_t{}.expand(*static_cast<E const*>(&e), std::forward<Ts>(subs)...));
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int I0, int P0, int M0, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<index_eq_N<i_<I0, P0>, M0>, Xs...>>
	{
		using second_limit_t = std::conditional_t<
			(N0 > 0),
			expr::series_limits_ntn<M0, N0 + M0 - 1>,
			expr::series_limits_nte<M0>>;

		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I0, P0>>>,
			types_list<limit_ts..., second_limit_t>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int I0, int P0, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<i_<I0, P0>, Xs...>>
	{
		using second_limit_t = std::conditional_t<
			(N0 > 0),
			expr::series_limits_ntn<1, N0>,
			expr::series_limits_nte<1>>;

		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I0, P0>>>,
			types_list<limit_ts..., second_limit_t>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int I0, int P0, typename... Es, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<index_eq<i_<I0, P0>, OpAdd<Es...>>, Xs...>>
	{
		using second_limit_t = std::conditional_t<
			(N0 > 0),
			expr::series_limits_stn<OpAdd<Es...>, N0>,
			expr::series_limits_ste<OpAdd<Es...>>>;

		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I0, P0>>>,
			types_list<limit_ts..., second_limit_t>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int I0, int P0, typename... Es, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<index_neq<i_<I0, P0>, OpAdd<Es...>>, Xs...>>
	{
		using second_limit_t = std::conditional_t<
			(N0 > 0),
			expr::series_limits_stn<OpAdd<OpIdentity, Es...>, N0>,
			expr::series_limits_ste<OpAdd<OpIdentity, Es...>>>;

		using first_t = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I0, P0>>>,
			types_list<limit_ts..., expr::series_limits_nts<1, OpAdd<Es..., OpNegIdentity>>>, types_list<Xs...>>;
		using second_t = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I0, P0 + 1>>>,
			types_list<limit_ts..., second_limit_t>, types_list<Xs...>>;

		using all_types = expand_types_list<typename first_t::all_types, typename second_t::all_types>;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return Op{}(
				first_t{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...),
				second_t{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...));
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return Op{}.expand(
				first_t{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...),
				second_t{}.expand(*static_cast<E const*>(&e), std::forward<Ts>(subs)...));
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int I00, int P00, int I01, int P01, int I1, int P1, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<symphas::lib::types_list<index_eq<i_<I1, P1>, i_<I00, P00>>, index_eq<i_<I1, P1>, i_<I01, P01>>>, Xs...>>
	{
		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>,
			types_list<limit_ts..., expr::series_limits_iti<i_<I00, P00>, i_<I01, P01>>>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int I00, int P00, int M1, int I1, int P1, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<symphas::lib::types_list<index_eq<i_<I1, P1>, i_<I00, P00>>, index_eq_N<i_<I1, P1>, M1>>, Xs...>>
	{
		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>,
			types_list<limit_ts..., expr::series_limits_itn<i_<I00, P00>, M1>>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int M0, int I01, int P01, int I1, int P1, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<symphas::lib::types_list<index_eq_N<i_<I1, P1>, M0>, index_eq<i_<I1, P1>, i_<I01, P01>>>, Xs...>>
	{
		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>,
			types_list<limit_ts..., expr::series_limits_nti<M0, i_<I01, P01>>>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int M0, int M1, int I1, int P1, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<symphas::lib::types_list<index_eq_N<i_<I1, P1>, M0>, index_eq_N<i_<I1, P1>, M1>>, Xs...>>
	{
		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>,
			types_list<limit_ts..., expr::series_limits_ntn<M0, M1>>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		typename... Es, int I01, int P01, int I1, int P1, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<symphas::lib::types_list<index_eq<i_<I1, P1>, OpAdd<Es...>>, index_eq<i_<I1, P1>, i_<I01, P01>>>, Xs...>>
	{
		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>,
			types_list<limit_ts..., expr::series_limits_sti<OpAdd<Es...>, i_<I01, P01>>>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int I00, int P00, typename... Es, int I1, int P1, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<symphas::lib::types_list<index_eq<i_<I1, P1>, i_<I00, P00>>, index_eq<i_<I1, P1>, OpAdd<Es...>>>, Xs...>>
	{
		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>,
			types_list<limit_ts..., expr::series_limits_its<i_<I00, P00>, OpAdd<Es...>>>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		typename... Es, int M1, int I1, int P1, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<symphas::lib::types_list<index_eq<i_<I1, P1>, OpAdd<Es...>>, index_eq_N<i_<I1, P1>, M1>>, Xs...>>
	{
		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>,
			types_list<limit_ts..., expr::series_limits_stn<OpAdd<Es...>, M1>>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		int M0, typename... Es, int I1, int P1, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<symphas::lib::types_list<index_eq_N<i_<I1, P1>, M0>, index_eq<i_<I1, P1>, OpAdd<Es...>>>, Xs...>>
	{
		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>,
			types_list<limit_ts..., expr::series_limits_nts<M0, OpAdd<Es...>>>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};

	//! Expand the sum where all indices are set to specific values.
	template<size_t N0, size_t... Ns, typename... I0s, typename... limit_ts, typename E, typename Op,
		typename... E0s, typename... E1s, int I1, int P1, typename... Xs>
	struct expand_series<
		std::index_sequence<N0, Ns...>, SymbolicSeries<Op, E, types_list<I0s...>>,
		types_list<limit_ts...>,
		types_list<symphas::lib::types_list<index_eq<i_<I1, P1>, OpAdd<E0s...>>, index_eq<i_<I1, P1>, OpAdd<E1s...>>>, Xs...>>
	{
		using type = expand_series<std::index_sequence<Ns...>, SymbolicSeries<Op, E, types_list<I0s..., i_<I1, P1>>>,
			types_list<limit_ts..., expr::series_limits_sts<OpAdd<E0s...>, OpAdd<E1s...>>>, types_list<Xs...>>;

		using all_types = typename type::all_types;

		template<typename... Ts>
		auto operator()(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}

		template<typename... Ts>
		auto expand(OpExpression<E> const& e, Ts&&... subs) const
		{
			return type{}(*static_cast<E const*>(&e), std::forward<Ts>(subs)...);
		}
	};




	template<typename E, typename T1, typename T2>
	auto normalize_indices(OpExpression<E> const& e, DynamicIndex const& index, 
		symphas::lib::types_list<>, expr::series_limits<T1, T2> const& limit)
	{
		return *static_cast<E const*>(&e);
	}

	template<typename E, int I0, int P00, int... P0s, typename T1, typename T2>
	auto normalize_indices(OpExpression<E> const& e, DynamicIndex const& index,
		symphas::lib::types_list<expr::symbols::i_<I0, P00>, expr::symbols::i_<I0, P0s>...>,
		expr::series_limits<T1, T2> const& limit)
	{
		auto e0 = expr::transform::swap_grid<expr::symbols::i_<I0, P00>, expr::symbols::i_<I0, P0s>...>(
			*static_cast<E const*>(&e), 
			(expr::symbols::i_<I0, 0>{} + expr::make_integer<P00>()), (expr::symbols::i_<I0, 0>{} + expr::make_integer<P0s>())...);
		return expr::transform::swap_grid<OpCoeffSwap<expr::symbols::i_<I0, 0>>>(e0, index);
	}

	//template<int P0, typename E>
	//auto normalize_placeholders(OpExpression<E> const& e, symphas::lib::types_list<>)
	//{
	//	return *static_cast<E const*>(&e);
	//}

	//template<int P0, typename E, int I0, int P00, int... P0s>
	//auto normalize_placeholders(OpExpression<E> const& e,
	//	symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I0, P00>>, expr::symbols::v_id_type<expr::symbols::i_<I0, P0s>>...>)
	//{
	//	return expr::transform::swap_grid<expr::symbols::v_id_type<expr::symbols::i_<I0, P00>>, expr::symbols::v_id_type<expr::symbols::i_<I0, P0s>>...>
	//		(*static_cast<E const*>(&e), expr::symbols::v_id_type<expr::symbols::i_<I0, P00 + P0>>{}, expr::symbols::v_id_type<expr::symbols::i_<I0, P0s + P0>>{}...);
	//}

	template<typename v_types, typename i_types, int I0, int P0, typename E, typename T1, typename T2>
	auto normalize(OpExpression<E> const& e, DynamicIndex *index, 
		expr::symbols::i_<I0, P0>, expr::series_limits<T1, T2> const& limit)
	{
		//return normalize_indices(normalize_placeholders<P0>(*static_cast<E const*>(&e), v_types{}), index, i_types{}, limit);
		return normalize_indices(*static_cast<E const*>(&e), *index, i_types{}, limit);
	}

	template<typename v_types, typename E>
	auto normalize(OpExpression<E> const& e, DynamicIndex* index, symphas::lib::types_list<>, std::tuple<> const& limits)
	{
		return *static_cast<E const*>(&e);
	}

	template<typename v_types, typename i_types0, typename... i_types, int I0, int P0, int... I0s, int... P0s, typename E, typename... T1s, typename... T2s>
	auto normalize(OpExpression<E> const& e, DynamicIndex* index,
		symphas::lib::types_list<expr::symbols::i_<I0, P0>, expr::symbols::i_<I0s, P0s>...>,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	{
		auto vs = select_v_i_<expr::symbols::i_<I0, 0>, v_types>{};
		//auto e0 = normalize_indices(normalize_placeholders<P0>(*static_cast<E const*>(&e), vs), i_types0{}, std::get<0>(limits));
		auto e0 = normalize_indices(*static_cast<E const*>(&e), *index, i_types0{}, std::get<0>(limits));
		return normalize<v_types, i_types...>(e0, index + 1, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{}, symphas::lib::get_tuple_ge<1>(limits));
	}

	//template<typename v_types, typename... i_types, int... I0s, int... P0s, typename E, typename... T1s, typename... T2s>
	//auto normalize(OpExpression<E> const& e, DynamicIndex* index, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, std::tuple<expr::series_limits<T1s, T2s>...> const& limits)
	//{
	//	return normalize<v_types, i_types...>(*static_cast<E const*>(&e), index, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>{}, limits);
	//}

	template<typename Op, typename E, typename S, size_t... Ns, 
		typename... T1s, typename... T2s, int... I0s, int... P0s, typename... Is>
	auto construct_series(
		OpExpression<E> const& e,
		const DynamicIndex(&index)[sizeof...(I0s)],
		SymbolicTemplate<S, Ns...> const& tmpl,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>,
		expr::series_limits<T1s, T2s> const&... limits);

	template<typename Op, typename E, typename... Ts, 
		typename... T1s, typename... T2s, int... I0s, int... P0s, typename... Is>
	auto construct_series(
		OpExpression<E> const& e,
		const DynamicIndex(&index)[sizeof...(I0s)],
		Substitution<SymbolicDataArray<Ts>...> const& substitution,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>,
		expr::series_limits<T1s, T2s> const&... limits);


	template<int... I0s, typename... Ts, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<int, T2> const& limit);

	template<int... I0s, typename... Ts, typename V, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpBinaryMul<V, DynamicIndex>, T2> const& limit);

	template<int... I0s, typename... Ts, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<DynamicIndex, T2> const& limit);

	template<int... I0s, typename... Ts, typename T2, typename G>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, G>, T2> const& limit);

	template<int... I0s, typename... Ts, int I0, int P0, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, T2> const& limit);

	template<int... I0s, typename... Ts, typename... E0s, typename T2, size_t... Ns>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit, std::index_sequence<Ns...>);

	template<int... I0s, typename... Ts, typename... E0s, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit);

	template<int... I0s, typename... Ts, typename T2, typename A, typename B>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<std::pair<A, B>, T2> const& limit);

	template<int... I0s, typename... Ts, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, int> const& limit);

	template<int... I0s, typename... Ts, typename V, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpBinaryMul<V, DynamicIndex>> const& limit);

	template<int... I0s, typename... Ts, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, DynamicIndex> const& limit);

	template<int... I0s, typename... Ts, typename T1, typename G>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, G>> const& limit);

	template<int... I0s, typename... Ts, typename T1, int I1, int P1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, expr::symbols::i_<I1, P1>>> const& limit);

	template<int... I0s, typename... Ts, typename T1, typename... E0s, size_t... Ns>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit, std::index_sequence<Ns...>);

	template<int... I0s, typename... Ts, typename T1, typename... E0s>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit);

	template<int... I0s, typename... Ts, typename T1, typename A, typename B>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, std::pair<A, B>> const& limit);

	template<int... I0s, typename... Ts, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<int, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename V, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpBinaryMul<V, DynamicIndex>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<DynamicIndex, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename T2, typename G>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, G>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, int I0, int P0, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename... E0s, typename T2, size_t... Ns>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)],
		std::index_sequence<Ns...>);

	template<int... I0s, typename... Ts, typename... E0s, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename T2, typename A, typename B>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<std::pair<A, B>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, int> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename V, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpBinaryMul<V, DynamicIndex>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, DynamicIndex> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename T1, typename G>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, G>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename T1, int I1, int P1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, expr::symbols::i_<I1, P1>>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename T1, typename... E0s, size_t... Ns>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)],
		std::index_sequence<Ns...>);

	template<int... I0s, typename... Ts, typename T1, typename... E0s>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename T1, typename A, typename B>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, std::pair<A, B>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);


	template<int... I0s, typename... Ts, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<int, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, size_t N0, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<expr::symbols::placeholder_N_<N0>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename V, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpBinaryMul<V, DynamicIndex>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<DynamicIndex, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename G, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, G>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, int I0, int P0, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<int... I0s, typename... Ts, typename... E0s, size_t... Ns, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)],
		std::index_sequence<Ns...>);

	template<int... I0s, typename... Ts, typename... E0s, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)]);

	template<bool left_flag>
	inline auto limit_side(expr::symbols::Symbol const&, expr::symbols::Symbol const&)
	{
		return "df";
		//static_assert(false, "limits of a series must evaluate to an integer if a series is to be computed");
	}

	template<bool left_flag>
	inline auto limit_side(expr::symbols::Symbol const&, int value)
	{
		return "df";
		//static_assert(false, "limits of a series must evaluate to an integer if a series is to be computed");
	}

	template<bool left_flag>
	inline auto limit_side(int value, expr::symbols::Symbol const&)
	{
		return "df";
		//static_assert(false, "limits of a series must evaluate to an integer if a series is to be computed");
	}

	template<bool left_flag>
	inline auto limit_side(int value1, int value2)
	{
		return (left_flag) ? std::max(value1, value2) : std::min(value1, value2);
	}

	template<bool left_flag, typename E1, typename E2>
	inline auto limit_side(OpExpression<E1> const& left, OpExpression<E2> const& right)
	{
		return limit_side<left_flag>(expr::eval(*static_cast<E1 const*>(&left)), expr::eval(*static_cast<E2 const*>(&right)));
	}


	inline auto check_limit_range(expr::symbols::Symbol)
	{
		return 0;
	}

	inline auto check_limit_range(int range)
	{
		return range;
	}

	template<typename E>
	auto check_limit_range(OpExpression<E> const& e)
	{
		return check_limit_range(static_cast<E const*>(&e)->eval());
	}

	template<int... I0s, typename... Ts, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<int, T2> const& limit)
	{
		return limit._0;
	}

	template<int... I0s, typename... Ts, typename V, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpBinaryMul<V, DynamicIndex>, T2> const& limit)
	{
		return limit._0.a.eval() * limit._0.b.start();
	}

	template<int... I0s, typename... Ts, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<DynamicIndex, T2> const& limit)
	{
		return limit._0.start();
	}

	template<int... I0s, typename... Ts, typename T2, typename G>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, G>, T2> const& limit)
	{
		return expr::eval(limit._0c) * expr::eval(limit._0);
	}

	template<int... I0s, typename... Ts, int I0, int P0, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, T2> const& limit)
	{
		static const int N = symphas::lib::index_of_value<int, I0, I0s...>;
		if constexpr (N < 0)
		{
			return 0;
		}
		else
		{
			auto start = limit_start(std::integer_sequence<int, I0s...>{}, limits, std::get<size_t(N)>(limits));
			return limit._0c * start;
		}
	}

	template<int... I0s, typename... Ts, typename... E0s, typename T2, size_t... Ns>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit, std::index_sequence<Ns...>)
	{
		return (limit_start(std::integer_sequence<int, I0s...>{}, limits,
			expr::series_limits(expr::get<Ns>(limit._0), limit._1))
			+ ...);
	}

	template<int... I0s, typename... Ts, typename... E0s, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit)
	{
		return limit_start(std::integer_sequence<int, I0s...>{}, limits,
			limit, std::make_index_sequence<sizeof...(E0s)>{});
	}

	template<int... I0s, typename... Ts, typename T2, typename A, typename B>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<std::pair<A, B>, T2> const& limit)
	{
		auto [a, b] = limit._0;
		auto la = limit_start(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(a, limit._1));
		auto lb = limit_start(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(b, limit._1));
		return limit_side<false>(la, lb);
	}

	template<int... I0s, typename... Ts, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, int> const& limit)
	{
		return limit._1;
	}

	template<int... I0s, typename... Ts, typename V, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpBinaryMul<V, DynamicIndex>> const& limit)
	{
		return limit._1.a.eval() * limit._1.b.start();
	}

	template<int... I0s, typename... Ts, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, DynamicIndex> const& limit)
	{
		return limit._1.start();
	}

	template<int... I0s, typename... Ts, typename T1, typename G>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, G>> const& limit)
	{
		return expr::eval(limit._1c) * expr::eval(limit._1);
	}

	template<int... I0s, typename... Ts, typename T1, int I1, int P1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, expr::symbols::i_<I1, P1>>> const& limit)
	{
		static const int N = symphas::lib::index_of_value<int, I1, I0s...>;
		if constexpr (N < 0)
		{
			return 0;
		}
		else
		{
			auto end = limit_end(std::integer_sequence<int, I0s...>{}, limits, std::get<size_t(N)>(limits));
			return limit._1c * end;
		}
	}

	template<int... I0s, typename... Ts, typename T1, typename... E0s, size_t... Ns>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit, std::index_sequence<Ns...>)
	{
		return (limit_end(std::integer_sequence<int, I0s...>{}, limits,
			expr::series_limits(limit._0, expr::get<Ns>(limit._1)))
			+ ...);
	}

	template<int... I0s, typename... Ts, typename T1, typename... E0s>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit)
	{
		return limit_end(std::integer_sequence<int, I0s...>{}, limits,
			limit, std::make_index_sequence<sizeof...(E0s)>{});
	}

	template<int... I0s, typename... Ts, typename T1, typename A, typename B>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, std::pair<A, B>> const& limit)
	{
		auto [a, b] = limit._1;
		auto la = limit_end(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(limit._0, a));
		auto lb = limit_end(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(limit._0, b));
		return limit_side<true>(la, lb);
	}

	//template<int... I0s, typename... Ts, typename T1, typename T2>
	//auto limit_range(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<T1, T2> const& limit)
	//{
	//	return check_limit_range((limit_end(std::integer_sequence<int, I0s...>{}, limits, limit)
	//		- limit_start(std::integer_sequence<int, I0s...>{}, limits, limit) + 1));
	//}


	template<int... I0s, typename... Ts, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<int, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return limit._0;
	}

	template<int... I0s, typename... Ts, typename V, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpBinaryMul<V, DynamicIndex>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return limit._0.a.eval() * limit._0.b.index();
	}

	template<int... I0s, typename... Ts, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<DynamicIndex, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return limit._0.index();
	}

	template<int... I0s, typename... Ts, typename T2, typename G>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, G>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return expr::eval(limit._0c) * expr::eval(limit._0);
	}

	template<int... I0s, typename... Ts, int I0, int P0, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		static const int N = symphas::lib::index_of_value<int, I0, I0s...>;
		if constexpr (N < 0)
		{
			return limit._0c * P0;
		}
		else
		{
			auto start = limit_start(std::integer_sequence<int, I0s...>{}, limits, std::get<size_t(N)>(limits));
			return limit._0c * (start + offsets[N] + P0);
		}
	}

	template<int... I0s, typename... Ts, typename... E0s, typename T2, size_t... Ns>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)],
		std::index_sequence<Ns...>)
	{
		return (limit_start(std::integer_sequence<int, I0s...>{}, limits,
			expr::series_limits(expr::get<Ns>(limit._0), limit._1), offsets)
			+ ...);
	}

	template<int... I0s, typename... Ts, typename... E0s, typename T2>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return limit_start(std::integer_sequence<int, I0s...>{}, limits,
			limit, offsets, std::make_index_sequence<sizeof...(E0s)>{});
	}

	template<int... I0s, typename... Ts, typename T2, typename A, typename B>
	auto limit_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<std::pair<A, B>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		auto [a, b] = limit._0;
		auto la = limit_start(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(a, limit._1), offsets);
		auto lb = limit_start(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(b, limit._1), offsets);
		return limit_side<true>(la, lb);
	}

	template<int... I0s, typename... Ts, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, int> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return limit._1;
	}

	template<int... I0s, typename... Ts, typename V, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpBinaryMul<V, DynamicIndex>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return limit._1.a.eval() * limit._1.b.index();
	}

	template<int... I0s, typename... Ts, typename T1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, DynamicIndex> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return limit._1.index();
	}

	template<int... I0s, typename... Ts, typename T1, typename G>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, G>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return expr::eval(limit._1c) * expr::eval(limit._1);
	}

	template<int... I0s, typename... Ts, typename T1, int I1, int P1>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, expr::symbols::i_<I1, P1>>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		static const int N = symphas::lib::index_of_value<int, I1, I0s...>;
		if constexpr (N < 0)
		{
			return limit._1c * P1;
		}
		else
		{
			auto start = limit_start(std::integer_sequence<int, I0s...>{}, limits, std::get<size_t(N)>(limits));
			return limit._1c * (start + offsets[N] + P1);
		}
	}

	template<int... I0s, typename... Ts, typename T1, typename... E0s, size_t... Ns>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)],
		std::index_sequence<Ns...>)
	{
		return (limit_end(std::integer_sequence<int, I0s...>{}, limits,
			expr::series_limits(limit._0, expr::get<Ns>(limit._1)), offsets) + ...);
	}

	template<int... I0s, typename... Ts, typename T1, typename... E0s>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return limit_end(std::integer_sequence<int, I0s...>{}, limits,
			limit, offsets, std::make_index_sequence<sizeof...(E0s)>{});
	}

	template<int... I0s, typename... Ts, typename T1, typename A, typename B>
	auto limit_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, std::pair<A, B>> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		auto [a, b] = limit._1;
		auto la = limit_end(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(limit._0, a), offsets);
		auto lb = limit_end(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(limit._0, b), offsets);
		return limit_side<false>(la, lb);
	}

	template<int... I0s, typename... Ts, typename T1, typename T2>
	auto limit_range(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return check_limit_range((limit_end(std::integer_sequence<int, I0s...>{}, limits, limit, offsets)
			- limit_start(std::integer_sequence<int, I0s...>{}, limits, limit, offsets) + 1));
	}




	template<bool left_flag>
	inline auto limit_side_dimension(expr::symbols::Symbol const&, expr::symbols::Symbol const&)
	{
		return 0;
	}

	template<bool left_flag>
	inline auto limit_side_dimension(expr::symbols::Symbol const&, int value)
	{
		return value;
	}

	template<bool left_flag>
	inline auto limit_side_dimension(int value, expr::symbols::Symbol const&)
	{
		return value;
	}

	template<bool left_flag>
	inline auto limit_side_dimension(int value1, int value2)
	{
		return (left_flag) ? std::max(value1, value2) : std::min(value1, value2);
	}

	template<bool left_flag, typename E1, typename E2>
	inline auto limit_side_dimension(OpExpression<E1> const& left, OpExpression<E2> const& right)
	{
		return limit_side_dimension<left_flag>(expr::eval(*static_cast<E1 const*>(&left)), expr::eval(*static_cast<E2 const*>(&right)));
	}


	template<int... I0s, typename... Ts, typename T1, typename T2>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, T2> const& limit);
	template<int... I0s, typename... Ts, typename V, typename T2>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpBinaryMul<V, DynamicIndex>, T2> const& limit);
	template<int... I0s, typename... Ts, typename T2>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<DynamicIndex, T2> const& limit);
	template<int... I0s, typename... Ts, typename T2, typename G>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, G>, T2> const& limit);
	template<int... I0s, typename... Ts, int I0, int P0, typename T2>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, T2> const& limit);
	template<int... I0s, typename... Ts, typename... E0s, typename T2, size_t... Ns>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit, std::index_sequence<Ns...>);
	template<int... I0s, typename... Ts, typename... E0s, typename T2>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit);
	template<int... I0s, typename... Ts, typename T2, typename A, typename B>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<std::pair<A, B>, T2> const& limit);
	template<int... I0s, typename... Ts, typename T1, typename T2>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, T2> const& limit);
	template<int... I0s, typename... Ts, typename T1, typename V>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpBinaryMul<V, DynamicIndex>> const& limit);
	template<int... I0s, typename... Ts, typename T1>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, DynamicIndex> const& limit);
	template<int... I0s, typename... Ts, typename T1, typename G>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, G>> const& limit);
	template<int... I0s, typename... Ts, typename T1, int I1, int P1>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, expr::symbols::i_<I1, P1>>> const& limit);
	template<int... I0s, typename... Ts, typename T1, typename... E0s, size_t... Ns>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit, std::index_sequence<Ns...>);
	template<int... I0s, typename... Ts, typename T1, typename... E0s>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit);
	template<int... I0s, typename... Ts, typename T1, typename A, typename B>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, std::pair<A, B>> const& limit);

	template<int... I0s, typename... Ts, typename T1, typename T2>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, T2> const& limit)
	{
		return limit._0;
	}

	template<int... I0s, typename... Ts, typename V, typename T2>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpBinaryMul<V, DynamicIndex>, T2> const& limit)
	{
		return limit._0.a * expr::make_literal(limit._0.b.start());
	}

	template<int... I0s, typename... Ts, typename T2>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<DynamicIndex, T2> const& limit)
	{
		return expr::make_literal(limit._0.start());
	}

	template<int... I0s, typename... Ts, typename T2, typename G>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, G>, T2> const& limit)
	{
		return expr::eval(limit._0c) * expr::eval(limit._0);
	}

	template<int... I0s, typename... Ts, int I0, int P0, typename T2>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, T2> const& limit)
	{
		static const int N = symphas::lib::index_of_value<int, I0, I0s...>;
		if constexpr (N < 0)
		{
			return expr::make_literal(limit._0c) * val<P0>;
		}
		else
		{
			auto start = limit_dimension_start(std::integer_sequence<int, I0s...>{}, limits, std::get<size_t(N)>(limits));
			return expr::make_literal(limit._0c) * (start + val<P0>);
		}
	}

	template<int... I0s, typename... Ts, typename... E0s, typename T2, size_t... Ns>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit, std::index_sequence<Ns...>)
	{
		return (limit_dimension_start(std::integer_sequence<int, I0s...>{}, limits,
			expr::series_limits(expr::get<Ns>(limit._0), expr::limit_1(limit)))
			+ ...);
	}

	template<int... I0s, typename... Ts, typename... E0s, typename T2>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit)
	{
		return limit_dimension_start(std::integer_sequence<int, I0s...>{}, limits,
			limit, std::make_index_sequence<sizeof...(E0s)>{});
	}

	template<int... I0s, typename... Ts, typename T2, typename A, typename B>
	auto limit_dimension_start(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<std::pair<A, B>, T2> const& limit)
	{
		auto [a, b] = limit._0;
		auto la = limit_dimension_start(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(a, expr::limit_1(limit)));
		auto lb = limit_dimension_start(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(b, expr::limit_1(limit)));
		return limit_side_dimension<false>(la, lb);
	}

	template<int... I0s, typename... Ts, typename T1, typename T2>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, T2> const& limit)
	{
		return limit._1;
	}

	template<int... I0s, typename... Ts, typename T1, typename V>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpBinaryMul<V, DynamicIndex>> const& limit)
	{
		return limit._1.a * expr::make_literal(limit._1.b.end());
	}

	template<int... I0s, typename... Ts, typename T1>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, DynamicIndex> const& limit)
	{
		return expr::make_literal(limit._1.end());
	}

	template<int... I0s, typename... Ts, typename T1, typename G>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, G>> const& limit)
	{
		return expr::eval(limit._1c) * expr::eval(limit._1);
	}

	template<int... I0s, typename... Ts, typename T1, int I1, int P1>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpTerm<OpIdentity, expr::symbols::i_<I1, P1>>> const& limit)
	{
		static const int N = symphas::lib::index_of_value<int, I1, I0s...>;
		if constexpr (N < 0)
		{
			return expr::make_literal(limit._1c) * val<P1>;
		}
		else
		{
			auto end = limit_dimension_end(std::integer_sequence<int, I0s...>{}, limits, std::get<size_t(N)>(limits));
			return expr::make_literal(limit._1c) * (end + val<P1>);
		}
	}

	template<int... I0s, typename... Ts, typename T1, typename... E0s, size_t... Ns>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit, std::index_sequence<Ns...>)
	{
		return (limit_dimension_end(std::integer_sequence<int, I0s...>{}, limits,
			expr::series_limits(expr::limit_0(limit), expr::get<Ns>(expr::limit_1(limit))))
			+ ...);
	}

	template<int... I0s, typename... Ts, typename T1, typename... E0s>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, OpAdd<E0s...>> const& limit)
	{
		return limit_dimension_end(std::integer_sequence<int, I0s...>{}, limits,
			limit, std::make_index_sequence<sizeof...(E0s)>{});
	}

	template<int... I0s, typename... Ts, typename T1, typename A, typename B>
	auto limit_dimension_end(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, std::pair<A, B>> const& limit)
	{
		auto [a, b] = limit._1;
		auto la = limit_dimension_end(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(expr::limit_0(limit), a));
		auto lb = limit_dimension_end(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(expr::limit_0(limit), b));
		return limit_side_dimension<true>(expr::eval(la), expr::eval(lb));
	}

	template<int... I0s, typename... Ts, typename T1, typename T2>
	auto limit_dimension(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, T2> const& limit)
	{
		return check_limit_range((limit_dimension_end(std::integer_sequence<int, I0s...>{}, limits, limit)
			- limit_dimension_start(std::integer_sequence<int, I0s...>{}, limits, limit) + expr::symbols::one));
	}



	//template<size_t N, int... I0s, typename... Ts, typename T2>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<int, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)]);
	//template<size_t N, int... I0s, typename... Ts, typename T2>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<DynamicIndex, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)]);
	//template<size_t N, int... I0s, typename... Ts, typename T2, typename G>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<OpTerm<OpIdentity, G>, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)]);
	//template<size_t N, int... I0s, typename... Ts, int I0, int P0, typename T2>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)]);
	//template<size_t N, int... I0s, typename... Ts, typename... E0s, typename T2, size_t... Ns>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<OpAdd<E0s...>, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)],
	//	std::index_sequence<Ns...>);
	//template<size_t N, int... I0s, typename... Ts, typename... E0s, typename T2>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<OpAdd<E0s...>, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)]);
	//template<size_t N, int... I0s, typename... Ts, typename T2, typename A, typename B>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<std::pair<A, B>, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)]);

	//template<size_t N, int... I0s, typename... Ts, typename T2>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<int, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)])
	//{
	//	return limit._0 + offsets[N];
	//}

	//template<size_t N, int... I0s, typename... Ts, typename T2>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<DynamicIndex, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)])
	//{
	//	return limit._0.start() + offsets[N];
	//}

	//template<size_t N, int... I0s, typename... Ts, typename T2, typename G>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<OpTerm<OpIdentity, G>, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)])
	//{
	//	return expr::make_term(limit._0c, limit._0).eval();
	//}

	//template<size_t N, int... I0s, typename... Ts, int I0, int P0, typename T2>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)])
	//{
	//	static const int M = symphas::lib::index_of_value<int, I0, I0s...>;
	//	if constexpr (M < 0)
	//	{
	//		return offsets[N];
	//	}
	//	else
	//	{
	//		iter_type _offsets[sizeof...(I0s)]{};
	//		std::copy(offsets, offsets + sizeof...(I0s), _offsets);
	//		_offsets[M] = 0;

	//		auto start = limit_current<M>(std::integer_sequence<int, I0s...>{}, limits, std::get<size_t(M)>(limits), offsets);
	//		return limit._0c * start + offsets[N];
	//	}
	//}

	//template<size_t N, int... I0s, typename... Ts, typename... E0s, typename T2, size_t... Ns>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<OpAdd<E0s...>, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)],
	//	std::index_sequence<Ns...>)
	//{
	//	iter_type _offsets[sizeof...(I0s)]{};
	//	std::copy(offsets, offsets + sizeof...(I0s), _offsets);
	//	_offsets[N] = 0;

	//	return (limit_current<N>(std::integer_sequence<int, I0s...>{}, limits,
	//		expr::series_limits(expr::get<Ns>(limit._0), limit._1), _offsets)
	//		+ ...) + offsets[N];
	//}

	//template<size_t N, int... I0s, typename... Ts, typename... E0s, typename T2>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<OpAdd<E0s...>, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)])
	//{
	//	return limit_current<N>(std::integer_sequence<int, I0s...>{}, limits,
	//		limit, offsets, std::make_index_sequence<sizeof...(E0s)>{});
	//}

	//template<size_t N, int... I0s, typename... Ts, typename T2, typename A, typename B>
	//auto limit_current(
	//	std::integer_sequence<int, I0s...>,
	//	std::tuple<Ts...> limits,
	//	expr::series_limits<std::pair<A, B>, T2> const& limit,
	//	iter_type(&offsets)[sizeof...(I0s)])
	//{
	//	auto [a, b] = limit._0;
	//	auto la = limit_current<N>(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(a, limit._1));
	//	auto lb = limit_current<N>(std::integer_sequence<int, I0s...>{}, limits, expr::series_limits(b, limit._1));
	//	return limit_side<false>(la, lb);
	//}



	template<int... I0s, typename... Ts, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<int, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return 0;
	}

	template<int... I0s, typename... Ts, size_t N0, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<expr::symbols::placeholder_N_<N0>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return 0;
	}

	template<int... I0s, typename... Ts, typename V, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpBinaryMul<V, DynamicIndex>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return expr::limit_0(limit);// limit._0.a.eval()* limit._0.b.index();
	}

	template<int... I0s, typename... Ts, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<DynamicIndex, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return expr::limit_0(limit);//limit._0.index();
	}

	template<int... I0s, typename... Ts, typename G, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, G>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return 0;
	}

	template<int... I0s, typename... Ts, int I0, int P0, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpTerm<OpIdentity, expr::symbols::i_<I0, P0>>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		static const int N = symphas::lib::index_of_value<int, I0, I0s...>;
		if constexpr (N < 0)
		{
			return limit._0c * P0;
		}
		else
		{
			return limit._0c * (offsets[N] + P0);
		}
	}

	template<int... I0s, typename... Ts, typename... E0s, size_t... Ns, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)],
		std::index_sequence<Ns...>)
	{
		return (compute_offset(
			std::integer_sequence<int, I0s...>{}, limits,
			expr::series_limits(expr::get<Ns>(limit._0), limit._1), offsets)
			+ ...);
	}

	template<int... I0s, typename... Ts, typename... E0s, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<OpAdd<E0s...>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return compute_offset(
			std::integer_sequence<int, I0s...>{}, limits,
			limit, offsets, std::make_index_sequence<sizeof...(E0s)>{});
	}

	template<int... I0s, typename... Ts, typename A, typename B, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<std::pair<A, B>, T2> const& limit,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		auto [a, b] = expr::limit_0(limit);
		auto offset0 = compute_offset(
			std::integer_sequence<int, I0s...>{}, limits, 
			expr::series_limits(a, expr::limit_1(limit)), offsets);
		auto offset1 = compute_offset(
			std::integer_sequence<int, I0s...>{}, limits, 
			expr::series_limits(b, expr::limit_1(limit)), offsets);
		return limit_side<true>(expr::eval(offset0), expr::eval(offset1));
	}

	template<int... I0s, typename... Ts, typename T1, typename T2>
	auto compute_offset(
		std::integer_sequence<int, I0s...>,
		std::tuple<Ts...> limits,
		expr::series_limits<T1, T2> const& limit)
	{
		len_type offsets[sizeof...(I0s)]{};
		return compute_offset(
			std::integer_sequence<int, I0s...>{}, limits,
			limit, offsets);
	}



	template<int N, typename... Gs, int I0, int P0, size_t L>
	auto select_data_of_index(
		SymbolicDataArray<std::tuple<Gs...>> const& substitution0,
		expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>,
		iter_type(&offsets)[L])
	{
		return std::make_tuple((substitution0.data[P0 + offsets[N]]));
	}

	template<int N, typename G, int I0, int P0, size_t L>
	auto select_data_of_index(
		SymbolicDataArray<G> const& substitution0,
		expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>,
		iter_type(&offsets)[L])
	{
		return std::make_tuple(NamedData(SymbolicData(&(substitution0.data[P0 + offsets[N]]), false), ""));
	}

	template<int N, typename G, int I0, int P0, size_t L>
	auto select_data_of_index(
		SymbolicDataArray<NamedData<G>> const& substitution0,
		expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>,
		iter_type(&offsets)[L])
	{
		return std::make_tuple(SymbolicData(&(substitution0.data[P0 + offsets[N]]), false));
	}

	template<int N, int I0, int P0, size_t L>
	auto select_data_of_index(
		SymbolicDataArray<expr::symbols::Symbol> const& substitution0,
		expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>,
		iter_type(&offsets)[L])
	{
		return std::make_tuple();
	}

	template<int N, int I0, int P0, size_t L>
	auto select_data_of_index(
		Substitution<> const& substitution,
		expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>,
		iter_type(&offsets)[L])
	{
		return std::make_tuple();
	}

	template<int N, typename... Gs, int I0, int P0, size_t L>
	auto select_data_of_index(
		Substitution<SymbolicDataArray<Gs>...> const& substitution,
		expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>,
		iter_type(&offsets)[L])
	{
		return select_data_of_index<N>(
			std::get<size_t(N)>(substitution),
			expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>{},
			offsets);
	}


	template<int N, typename... Gs, int I0, int P0>
	auto select_data_of_index(
		SymbolicDataArray<std::tuple<Gs...>> const& substitution0,
		expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>)
	{
		return std::make_tuple(expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>{});
	}

	template<int N, typename G, int I0, int P0>
	auto select_data_of_index(
		SymbolicDataArray<G> const& substitution0,
		expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>)
	{
		return std::make_tuple(expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>{});
	}

	template<int N, int I0, int P0>
	auto select_data_of_index(
		SymbolicDataArray<expr::symbols::Symbol> const& substitution0,
		expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>)
	{
		return std::make_tuple();
	}

	template<int N, int I0, int P0>
	auto select_data_of_index(
		Substitution<> const& substitution,
		expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>)
	{
		return std::make_tuple();
	}

	template<int N, typename... Gs, int I0, int P0>
	auto select_data_of_index(
		Substitution<SymbolicDataArray<Gs>...> const& substitution,
		expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>)
	{
		return select_data_of_index<N>(
			std::get<size_t(N)>(substitution),
			expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>{});
	}




	template<int... I0s, int... P0s>
	auto args_for_indices(
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>)
	{
		return std::make_tuple(expr::symbols::i_<I0s, 0>{}...);
	}

	template<size_t... Ms, size_t N>
	auto args_for_indices(
		iter_type(&offsets)[N],
		iter_type(&values)[N],
		std::index_sequence<Ms...>)
	{
		return std::make_tuple((values[Ms])...);
	}

	template<typename... T0s, int... I0s, int... P0s, int... I00s, int... P00s>
	auto args_for_placeholders(
		Substitution<SymbolicDataArray<T0s>...> const& substitution,
		symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I00s, P00s>>...>,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		iter_type(&offsets)[sizeof...(I0s)])
	{
		return std::tuple_cat(
			select_data_of_index<symphas::lib::index_of_value<int, I00s, I0s...>>
			(substitution, expr::symbols::v_id_type<expr::symbols::i_<I00s, P00s>>{}, offsets)...);
	}

	template<size_t D, int I0, int P0>
	auto as_grid_data(std::tuple<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>> const&)
	{
		return GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I0, P0>>, D>{};
	}

	template<size_t D, typename T0>
	auto as_grid_data(std::tuple<T0> const& data)
	{
		return std::get<0>(data);
	}

	template<typename... T0s, int... I00s, int... P00s, size_t... Ds, int... I0s, int... P0s>
	auto placeholder_list(
		Substitution<SymbolicDataArray<T0s>...> const& substitution,
		symphas::lib::types_list<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I00s, P00s>>, Ds>...>,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>)
	{
		using id_types = symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>;

		return std::tuple_cat(
			std::make_tuple(as_grid_data<Ds>(select_data_of_index<symphas::lib::index_of_value<int, I00s, I0s...>>
				(substitution, expr::symbols::v_id_type<expr::symbols::i_<I00s, P00s>>{}))...),
			args_for_indices(id_types{}));
	}



	//! Generates both the substitution list and expanded arguments with offsets.
	/*!
	 * The given data elements that should be substituted are passed to this function. The
	 * data elements will be used to generate a Substitution list, and that substitution list
	 * will be used to generate a list of SymbolicDataArray instances. Each of those
	 * SymbolicDataArray instances is connected with the Substitution, including with a necessary
	 * offset. The SymbolicDataArray instances will be generated and put in a tuple in the same
	 * order that they need to be passed into the function defining the series.
	 */
	template<typename... T0s, typename... Is, int... I0s, int... P0s, size_t N = sizeof...(I0s)>
	auto expanded_arg_list(
		Substitution<SymbolicDataArray<T0s>...> const& substitution,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>)
	{
		using v_id_types = symphas::lib::types_list<expr::symbols::v_id_type<Is>...>;
		using id_types = symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>;

		iter_type offsets[N]{};
		return std::tuple_cat(
			args_for_placeholders(substitution, v_id_types{}, id_types{}, offsets),
			args_for_indices(offsets, offsets, std::make_index_sequence<N>{}));
	}

	template<typename... T0s, typename... Is, int... I0s, int... P0s>
	auto expanded_arg_list(
		Substitution<SymbolicDataArray<T0s>...> const& substitution,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		iter_type(&offsets)[sizeof...(I0s)],
		iter_type(&values)[sizeof...(I0s)])
	{
		using v_id_types = symphas::lib::types_list<expr::symbols::v_id_type<Is>...>;
		using id_types = symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>;

		return std::tuple_cat(
			args_for_placeholders(substitution, v_id_types{}, id_types{}, offsets),
			args_for_indices(offsets, values, std::make_index_sequence<sizeof...(I0s)>{}));
	}

	template<typename E0, typename... T0s>
	static auto _get_function_list(SymbolicFunction<E0, T0s...>)
	{
		return SymbolicFunctionArray<E0, T0s...>(0);
	}

	template<typename S, size_t... Ns, typename... T0s>
	static auto _get_function_list(SymbolicTemplate<S, Ns...> const& tmpl, std::tuple<T0s...>)
	{
		return _get_function_list(expr::function_of(expr::symbols::arg<Ns, T0s>...) = tmpl);
	}

	template<typename Op, typename E, typename... Ts, int... I0s, int... P0s,
		typename... T1s, typename... T2s>
	auto _get_function_list(
		SymbolicSeries<Op, E, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>> const& series,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		Substitution<SymbolicDataArray<Ts>...> const& substitution)
	{
		if constexpr (sizeof...(Ts) == 0)
		{
			return 0;
		}
		else
		{
			using id_types = symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>;
			using all_var_types = typename expr::op_types<E>::type;
			using v_id_types = symphas::lib::expand_types_list<
				symphas::internal::select_v_i_<expr::symbols::i_<I0s, P0s>, all_var_types>...>;

			DynamicIndex index[sizeof...(I0s)];
			auto tmpl = series.get_template(&index[0], limits, v_id_types{});
			auto args = expanded_arg_list(substitution, v_id_types{}, id_types{});
			return _get_function_list(tmpl, args);
		}
	}

	template<typename A, typename B, typename C>
	auto get_function_list(
		A const& series,
		B const& limits,
		C const& substitution)
	{
		return _get_function_list(series, limits, substitution);
	}
}


template<typename Op, int I0, int P0>
struct SymbolicSeries<Op, OpVoid, symphas::lib::types_list<expr::symbols::i_<I0, P0>>>
{
	SymbolicSeries() = default;
	SymbolicSeries(OpVoid) {}

	template<typename... Ts>
	auto expand(Ts&&... subs) const
	{
		return OpVoid{};
	}

	template<typename... Ts>
	auto operator()(Ts&&... subs) const
	{
		return OpVoid{};
	}
};

//! Defines a series of an expression with indices and substitutable expressions.
/*!
 * A series is defined using an expression that contains expr::symbols::i_ and expr::symbols::v_
 * terms. The term i_ represents the index of the series, where indices are distinguished using
 * their first template value. The term v_ represents a data that is substituted into the series
 * expression, and is templated using the i_ type. Depending on the template substitutions,
 * different values will be substituted when the series is expanded.
 *
 * The series can be fully expanded into a symbolic expression, or it can be converted into an
 * OpSymbolicEval type where the series is iteratively evaluated using a list of data. The latter
 * should be the method used when the substituted data list is long. I.e. if the series would be
 * expanded to 200 or more terms, then OpSymbolicEval should be used.
 */
template<typename Op, typename E, int I0, int P0>
struct SymbolicSeries<Op, E, symphas::lib::types_list<expr::symbols::i_<I0, P0>>>
{
protected:

	using v_types = symphas::internal::select_v_i_<expr::symbols::i_<I0, P0>, expr::op_types_t<E>>;
	using all_indices_of_id = symphas::internal::select_all_i_<expr::symbols::i_<I0, P0>, expr::op_types_t<E>>;

	
	template<typename limits_t>
	using limits_of = symphas::internal::series_limits_of<limits_t>;

	static const int v_type_count = symphas::lib::types_list_size<v_types>::value;

public:

	SymbolicSeries() = default;
	SymbolicSeries(E const& e) : e{ e } {}

	template<size_t N, typename X0>
	auto select(X0) const
	{
		return SymbolicSeries<Op, E, 
			symphas::lib::types_list<std::index_sequence<N>,
				typename symphas::internal::filter_i_selection<X0>::type>>(e);
	}

	template<typename X0>
	auto select(X0) const
	{
		return select<0>(X0{});
	}

	//! Expand the sum between the given indices.
	/*!
	 * Expand the sum between the given indices, where `I0` is the first index, and `I1` is the
	 * final index. The final index is included in the sum (**not** C-style, where the last index of
	 * a loop is omitted).
	 *
	 * This expand function can only be used if there are no substitutable terms
	 * corresponding to this index.
	 */
	template<int N0, int N1, size_t N = v_type_count,
		typename std::enable_if_t<(N == 0), int> = 0>
	auto expand() const
	{
		if constexpr (N1 < N0)
		{
			return Op{}();
		}
		else
		{
			return expand_integer<N0>(
				all_indices_of_id{},
				std::make_integer_sequence<int, N1 - N0 + 1>{});
		}
	}

	template<typename limits_t, size_t N = v_type_count,
		int N0 = limits_of<limits_t>::_0, int N1 = limits_of<limits_t>::_1,
		typename std::enable_if_t<(N == 0 && N1 >= N0), int> = 0>
	auto expand() const
	{
		return expand<N0, N1>();
	}


	//! Expand the sum for all substitutions.
	/*!
	 * Expand the sum, starting at the given index. The index is then incremented for each term
	 * that is provided. The provided terms are substituted into the corresponding
	 * expr::symbols::v_ variables.
	 */
	template<int N0, int N1, typename... Es>
	auto expand(std::tuple<Es...> const& subs) const
	{
		if constexpr (N1 >= N0)
		{
			return expand_function<N0>(std::make_integer_sequence<int, N1 - N0 + 1>{}, all_indices_of_id{}, subs);
		}
		else
		{
			return Op{}();
		}
	}

	//! Expand the sum for all substitutions.
	/*!
	 * Expand the sum, starting at the given index. The index is then incremented for each term
	 * that is provided. The provided terms are substituted into the corresponding
	 * expr::symbols::v_ variables.
	 */
	template<int N0, int N1, typename... Es>
	auto expand(Es const& ...subs) const
	{
		return expand<N0, N1>(std::make_tuple(subs...));
	}

	//! Expand the sum for all substitutions.
	/*!
	 * Expand the sum, starting at the given index. The index is then incremented for each term
	 * that is provided. The provided terms are substituted into the corresponding
	 * expr::symbols::v_ variables.
	 */
	template<int N0, typename... Es>
	auto expand(std::tuple<Es...> const& subs) const
	{
		return expand<N0, int(N0 + sizeof...(Es) - 1)>(subs);
	}

	//! Expand the sum for all substitutions.
	/*!
	 * Expand the sum, starting at the given index. The index is then incremented for each term
	 * that is provided. The provided terms are substituted into the corresponding
	 * expr::symbols::v_ variables.
	 */
	template<int N0, typename... Es>
	auto expand(Es const& ...subs) const
	{
		return expand<N0, int(N0 + sizeof...(Es) - 1)>(subs...);
	}
	
	//! Expand the sum for all substitutions.
	/*!
	 * Expand the sum, starting at the given index. The index is then incremented for each term
	 * that is provided. The provided terms are substituted into the corresponding
	 * expr::symbols::v_ variables.
	 */
	template<typename limits_t, typename... Es,
		int N0 = limits_of<limits_t>::_0, int N1 = limits_of<limits_t>::_1,
		typename std::enable_if_t<(N1 >= N0), int> = 0>
	auto expand(std::tuple<Es...> const& subs) const
	{
		return expand<N0, N1>(subs);
	}

	//! Expand the sum for all substitutions.
	/*!
	 * Expand the sum, starting at the given index. The index is then incremented for each term
	 * that is provided. The provided terms are substituted into the corresponding
	 * expr::symbols::v_ variables.
	 */
	template<typename limits_t, typename... Es,
		int N0 = limits_of<limits_t>::_0, int N1 = limits_of<limits_t>::_1,
		typename std::enable_if_t<(N1 >= N0), int> = 0>
	auto expand(Es const& ...subs) const
	{
		return expand<N0, N1>(subs...);
	}

	template<typename T1, typename T2>
	auto operator()(expr::series_limits<T1, T2> const& limit) const
	{
		DynamicIndex index[1];
		return symphas::internal::construct_series<Op>(
			e, index, get_template(&index[0], limit, v_types{}),
			symphas::lib::types_list<expr::symbols::i_<I0, 0>>{},
			v_types{},
			limit);
	}

	template<typename T1, typename T2, size_t... Ns>
	auto operator()(std::tuple<expr::series_limits<T1, T2>> const& limit) const
	{
		return this->operator()(std::get<0>(limit));
	}

	E e;				//!< The expression which is summed.

	template<typename T1, typename T2, int... P0s>
	auto get_template(
		DynamicIndex* index,
		expr::series_limits<T1, T2> const& limit,
		symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I0, P0s>>...>) const
	{
		auto e0 = symphas::internal::normalize<
			symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I0, P0s>>...>, 
			all_indices_of_id>
			(e, index, expr::symbols::i_<I0, P0>{}, limit);
		return (expr::template_of(expr::symbols::v_id_type<expr::symbols::i_<I0, P0s>>{}..., expr::symbols::i_<I0, 0>{}) = e);
	}


	template<typename T1, typename T2, int... P0s>
	auto get_template(
		DynamicIndex *index,
		std::tuple<expr::series_limits<T1, T2>> const& limit,
		symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I0, P0s>>...>) const
	{
		return get_template(index, std::get<0>(limit), symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I0, P0s>>...>{});
	}

protected:

	template<int N0, int I00, int... IJ0s>
	auto swap_ii(symphas::lib::types_list<expr::symbols::i_<I00, IJ0s>...>) const
	{
		return expr::transform::swap_grid<expr::symbols::i_<I00, IJ0s>...>(e, expr::make_integer<N0 + IJ0s>()...);
	}

	template<int N0, int... IJ0s, int... Ns>
	auto expand_integer(
		symphas::lib::types_list<expr::symbols::i_<I0, IJ0s>...>,	// the indices to swap with integers
		std::integer_sequence<int, Ns...>							// integers to substitute
	) const
	{
		return Op{}(swap_ii<N0 + Ns - P0>(symphas::lib::types_list<expr::symbols::i_<I0, IJ0s>...>{})...);
	}

	template<size_t N, typename Ex, typename... Es>
	auto substitute_functions(
		symphas::lib::types_list<>,
		Ex const& ex,
		std::tuple<Es...> const& subs) const
	{
		return ex;
	}

	template<size_t N, typename Ex, typename... Es, int P00, int... P00s>
	auto substitute_functions(
		symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I0, P00>>, expr::symbols::v_id_type<expr::symbols::i_<I0, P00s>>...>,
		Ex const& ex, 
		std::tuple<Es...> const& subs) const
	{
		if constexpr (((int(N) + P00 < 0 || int(N) + P00 >= sizeof...(Es)) || ... || (int(N) + P00s < 0 || int(N) + P00s >= sizeof...(Es))))
		{
			return Op{}();
		}
		else
		{
			return expr::transform::swap_grid<
				expr::symbols::v_id_type<expr::symbols::i_<I0, P00>>, expr::symbols::v_id_type<expr::symbols::i_<I0, P00s>>...>
				(ex, std::get<N + P00>(subs), std::get<N + P00s>(subs)...);
		}
	}

	template<int N0, int... IJ0s, typename... Es, int... Ns>
	auto expand_function(
		std::integer_sequence<int, Ns...>,
		symphas::lib::types_list<expr::symbols::i_<I0, IJ0s>...>, 
		std::tuple<Es...> const& subs) const
	{
		if constexpr (sizeof...(IJ0s) > 0)
		{
			return Op{}(
				substitute_functions<size_t(Ns + P0)>(
					v_types{},
					swap_ii<N0 + Ns>(symphas::lib::types_list<expr::symbols::i_<I0, IJ0s>...>{}),
					subs)...);
		}
		else
		{
			return Op{}(
				substitute_functions<size_t(Ns + P0)>(v_types{}, e, subs)...);
		}
	}

	template<int N0>
	auto expand() const
	{
		return Op{}();
	}

};


template<typename Op, int I0, int P0, int I1, int P1, int... Is, int... Ps>
struct SymbolicSeries<Op, OpVoid, symphas::lib::types_list<
	expr::symbols::i_<I0, P0>, expr::symbols::i_<I1, P1>, expr::symbols::i_<Is, Ps>...
	>>
{
	SymbolicSeries() = default;
	SymbolicSeries(OpVoid) {}

	template<typename... Ts>
	auto expand(Ts&&... subs) const
	{
		return OpVoid{};
	}

	template<typename... Ts>
	auto operator()(Ts&&... subs) const
	{
		return OpVoid{};
	}

};

template<typename Op, typename E, int I0, int P0, int I1, int P1, int... Is, int... Ps>
struct SymbolicSeries<Op, E, symphas::lib::types_list<
	expr::symbols::i_<I0, P0>, expr::symbols::i_<I1, P1>, expr::symbols::i_<Is, Ps>...
	>>
{

protected:

	using all_v_types = symphas::internal::select_v_<expr::op_types_t<E>>;
	using all_v_nested_types = symphas::internal::select_v_nested_<
		symphas::lib::types_list<expr::symbols::i_<I0, P0>, expr::symbols::i_<I1, P1>, expr::symbols::i_<Is, Ps>...>, all_v_types>;
	using v_types = symphas::lib::expand_types_list<
		symphas::internal::select_v_i_<expr::symbols::i_<I0, P0>, all_v_types>,
		symphas::internal::select_v_i_<expr::symbols::i_<I1, P1>, all_v_types>,
		symphas::internal::select_v_i_<expr::symbols::i_<Is, Ps>, all_v_types>...
	>;

	//static constexpr int v_type_count = symphas::lib::types_list_size<
	//	symphas::lib::expand_types_list<
	//		//typename symphas::internal::select_v_<expr::symbols::i_<I0, P0>, all_v_types>::type,
	//		typename symphas::internal::select_v_<expr::symbols::i_<I1, P1>, all_v_types>::type,
	//		typename symphas::internal::select_v_<expr::symbols::i_<Is, Ps>, all_v_types>::type>...>::value;

	//static constexpr int all_v_type_count = symphas::lib::types_list_size<all_v_types>::value;
	static constexpr int nested_v_type_count = symphas::lib::types_list_size<all_v_nested_types>::value;

	//template<typename V>
	//static constexpr int index_of_type = symphas::lib::index_of_type<V, all_v_types>;

	template<typename II>
	using all_indices_of_id = symphas::internal::select_all_i_<II, expr::op_types_t<E>>;

	template<typename limits_t>
	using limits_of = symphas::internal::series_limits_of<limits_t>;

public:

	SymbolicSeries() = default;
	SymbolicSeries(E const& e) : e{ e } {}

	template<size_t... Ns, typename... Xs,
		typename = std::enable_if_t<(sizeof...(Xs) == sizeof...(Is) + 2 && sizeof...(Ns) == sizeof...(Xs)), int>>
	auto select(Xs...)
	{
		return SymbolicSeries<Op, E, 
			symphas::lib::types_list<std::index_sequence<Ns...>,
			typename symphas::internal::filter_i_selection<Xs>::type...>>(e);
	}

	template<typename... Xs>
	auto select(Xs...)
	{
		return select<symphas::internal::typei<0, Xs>...>(Xs{}...);
	}

	//! Expand the sum between the given indices.
	/*!
	 * Expand the sum between the given indices, where `I0` is the first index, and `I1` is the
	 * final index. The final index is included in the sum (**not** C-style, where the last index of
	 * a loop is omitted).
	 *
	 * This expand function can only be used if there are no substitutable terms
	 * corresponding to this index.
	 */
	template<int N0, int N1>
	auto expand_1() const
	{
		auto expr = SymbolicSeries<Op, E, symphas::lib::types_list<expr::symbols::i_<I0, P0>>>(e).template expand<N0, N1>();

		using I_rest = symphas::lib::types_list<expr::symbols::i_<I1, P1>, expr::symbols::i_<Is, Ps>...>;
		return SymbolicSeries<Op, decltype(expr), I_rest>(expr);
	}

	template<typename limits_t,
		int N0 = limits_of<limits_t>::_0, int N1 = limits_of<limits_t>::_1>
	auto expand_1() const
	{
		return expand_1<N0, N1>();
	}


	template<typename... T1s, typename... T2s, typename = std::enable_if_t<(sizeof...(Is) + 2 == sizeof...(T1s)), int>>
	auto operator()(expr::series_limits<T1s, T2s> const&... limits) const
	{
		using i_types = symphas::lib::types_list<expr::symbols::i_<I0, P0>, expr::symbols::i_<I1, P1>, expr::symbols::i_<Is, Ps>...>;

		DynamicIndex index[sizeof...(Is) + 2];
		return symphas::internal::construct_series<Op>(
			e, index, get_template(&index[0], std::make_tuple(limits...), v_types{}),
			i_types{},
			v_types{},
			limits...);
	}

	template<typename... T1s, typename... T2s, size_t... Ns>
	auto operator()(std::tuple<expr::series_limits<T1s, T2s>...> const& limits, std::index_sequence<Ns...>) const
	{
		return this->operator()(std::get<Ns>(limits)...);
	}

	template<typename... T1s, typename... T2s>
	auto operator()(std::tuple<expr::series_limits<T1s, T2s>...> const& limits) const
	{
		return this->operator()(limits, std::make_index_sequence<sizeof...(T1s)>{});
	}

	//! Start the iterated sum expansion using a single list of functions.
	/*!
	 * Expand the sum, starting at the given index. The rest of the indices have offsets that
	 * indicate the starting position of the substitutions in the list. The offsets are not based
	 * on the primary index, can only be positive, and do not loop around the list.
	 * 
	 * The list of functions is non-nested; it is assumed that the iterated sum has expressions 
	 * that are indexed by only a single value, e.g. $\psi_i$.
	 * 
	 * \param subs The list of expressions that is substituted, only a single list
	 * since it is assumed that the iterated sum has expressions that are indexed by only
	 * a single value, e.g. $\psi_i$.
	 * 
	 * \tparam I The starting value of the first iteration index.
	 * \tparam ISs... The starting value of the other indices is set to the current running index
	 * of the specified index.
	 */
	template<int N0, int N1, typename... limits_ts, typename... Ts, size_t N = nested_v_type_count,
		typename std::enable_if_t<(N == 0 && sizeof...(limits_ts) == sizeof...(Is) + 1), int> = 0>
	auto expand(std::tuple<Ts...> const& subs) const
	{
		return expand_iterated<N0>(subs, std::make_integer_sequence<int, N1 - N0 + 1>{}, limits_ts{}...);
	}

	template<typename limits_t, typename... limits_ts, typename... Ts,
		size_t N = nested_v_type_count, int N0 = limits_of<limits_t>::_0, int N1 = limits_of<limits_t>::_1,
		typename std::enable_if_t<(N == 0 && N1 >= N0 && sizeof...(limits_ts) == sizeof...(Is) + 1), int> = 0>
	auto expand(std::tuple<Ts...> const& subs) const
	{
		return expand<N0, N1, limits_ts...>(subs);
	}

	template<int N0, int... Ns, typename... Ts, size_t N = nested_v_type_count,
		typename std::enable_if_t<(N == 0 && sizeof...(Ns) == sizeof...(Is) + 1), int> = 0>
	auto expand(std::tuple<Ts...> const& subs) const
	{
		constexpr int M = sizeof...(Ts);
		return expand<expr::series_limits_ntn<N0, N0 + M>, expr::series_limits_ntn<Ns, Ns + M>...>(subs);
	}

	//! Start the iterated sum expansion using a single list of functions.
	template<int N0, int N1, typename... limits_ts, typename... Ts, size_t N = nested_v_type_count,
		typename std::enable_if_t<(N == 0 && sizeof...(limits_ts) == sizeof...(Is) + 1), int> = 0>
	auto expand(Ts const&... subs) const
	{
		return expand<N0, N1, limits_ts...>(std::make_tuple(subs...));
	}

	//! Start the iterated sum expansion using a single list of functions.
	template<typename limits_t, typename... limits_ts, typename... Ts,
		size_t N = nested_v_type_count,	int N0 = limits_of<limits_t>::_0, int N1 = limits_of<limits_t>::_1,
		typename std::enable_if_t<(N == 0 && sizeof...(limits_ts) == sizeof...(Is) + 1), int> = 0>
	auto expand(Ts const&... subs) const
	{
		return expand<N0, N1, limits_ts...>(std::make_tuple(subs...));
	}

	template<int N0, int... Ns, typename... Ts, size_t N = nested_v_type_count,
		typename std::enable_if_t<(N == 0 && sizeof...(Ns) == sizeof...(Is) + 1), int> = 0>
	auto expand(Ts const&... subs) const
	{
		constexpr int M = sizeof...(Ts);
		return expand<expr::series_limits_ntn<N0, N0 + M>, expr::series_limits_ntn<Ns, Ns + M>...>(std::make_tuple(subs...));
	}


	//! Start the iterated sum expansion using a single list of functions.
	/*!
	 * Expand the sum, starting at the given index. The rest of the indices have offsets that
	 * indicate the starting position of the substitutions in the list. The offsets are not based
	 * on the primary index, can only be positive, and do not loop around the list.
	 *
	 * The number of tuples provided to this expansion function is the same as the number of
	 * indices; a different expansion list will be used for placeholders of each index. The lists
	 * will be substituted into the placeholder corresponding with the index that matches the
	 * position 
	 *
	 * \tparam I The starting value of the first iteration index.
	 * \tparam ISs... The starting value of the other indices is set to the current running index
	 * of the specified index.
	 */
	template<int N0, int N1, typename... limits_ts, typename T00, typename... T0s, 
		typename T01, typename... T1s, typename... TTs, size_t N = nested_v_type_count,
		typename std::enable_if_t<(N == 0 && sizeof...(limits_ts) == sizeof...(Is) + 1), int> = 0>
	auto expand(std::tuple<T00, T0s...> const& subs0, std::tuple<T01, T1s...> const& subs1, TTs const&... subss) const
	{
		return expand_iterated<N0>(std::make_tuple(subs0, subs1, subss...), std::make_integer_sequence<int, N1 - N0 + 1>{}, limits_ts{}...);
	}

	template<typename limits_t, typename... limits_ts, typename T00, typename... T0s,
		typename T01, typename... T1s, typename... TTs, size_t N = nested_v_type_count,
		int N0 = limits_of<limits_t>::_0, int N1 = limits_of<limits_t>::_1,
		typename std::enable_if_t<(N == 0 && N1 >= N0 && sizeof...(limits_ts) == sizeof...(Is) + 1), int> = 0>
	auto expand(std::tuple<T00, T0s...> const& subs0, std::tuple<T01, T1s...> const& subs1, TTs const&... subss) const
	{
		return expand<N0, N1, limits_ts...>(subs0, subs1, subss...);
	}

	template<int N0, int... Ns, typename T00, typename... T0s,
		typename T01, typename... T1s, typename... TTs, size_t N = nested_v_type_count,
		typename std::enable_if_t<(N == 0 && sizeof...(Ns) == sizeof...(Is) + 1), int> = 0>
	auto expand(std::tuple<T00, T0s...> const& subs0, std::tuple<T01, T1s...> const& subs1, TTs const&... subss) const
	{
		constexpr int M = sizeof...(T0s) + 1;
		return expand<expr::series_limits_ntn<N0, N0 + M>, expr::series_limits_ntn<Ns, Ns + M>...>(subs0, subs1, subss...);
	}

	template<typename... T1s, typename... T2s, int... I0s, int... P0s>
	auto get_template(
		DynamicIndex* index,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>...>) const
	{
		using namespace expr::symbols;
		using I_list = symphas::lib::types_list<expr::symbols::i_<I0, P0>, expr::symbols::i_<I1, P1>, expr::symbols::i_<Is, Ps>...>;

		auto e0 = symphas::internal::normalize<
			v_types, 
			all_indices_of_id<i_<I0, 0>>,
			all_indices_of_id<i_<I1, 0>>,
			all_indices_of_id<i_<Is, 0>>...>(
			e, index, I_list{}, limits);

		return (expr::template_of(
			expr::symbols::v_id_type<expr::symbols::i_<I0s, P0s>>{}..., 
			expr::symbols::i_<I0, 0>{},
			expr::symbols::i_<I1, 0>{},
			expr::symbols::i_<Is, 0>{}...) = e0);

	}

protected:

	template<int N0, int N, typename... Ts, int M0, int M1, typename... limits_ts>
	auto construct_iterated(std::tuple<Ts...> const& subs,
		expr::series_limits_ntn<M0, M1>, limits_ts...) const
	{
		constexpr int NN = N0 + N; // Current value of the major index.

		auto sum_major = expr::series<Op, expr::symbols::i_<I0, P0 + N>>(e).template expand<NN, NN>(subs);
		if constexpr (M0 <= M1)
		{
			return expr::series<Op, expr::symbols::i_<I1, P1 + N>, expr::symbols::i_<Is, Ps>...>(sum_major)
				.template expand<expr::series_limits_ntn<M0, M1>, limits_ts...>(subs);
		}
		else
		{
			return Op{}();
		}
	}

	template<int N0, int N, typename... Ts, typename... limits_ts>
	auto expand_iterated(std::tuple<Ts...> const& subs, limits_ts...) const
	{
		return construct_iterated<N0, N>(subs, expr::update_limits<N0 + N, expr::symbols::i_<I0, P0>, limits_ts>{}...); // All indices to be substituted, Is
	}

	//! Start the iterated sum expansion.
	/*!
	 * Start the iterated sum expansion.
	 * 
	 * \param subs The list of expressions that is substituted, only a single list
	 * since it is assumed that the iterated sum has expressions that are indexed by only
	 * a single value, e.g. $\psi_i$. 
	 * 
	 * \tparam Ns... The starting index ID of an index from Is...
	 * \tparam Ps... The starting index value of an index from Is...
	 */
	template<int N0, typename... Ts, int... Ns, typename... limits_ts>
	auto expand_iterated(std::tuple<Ts...> const& subs, 
		std::integer_sequence<int, Ns...>,
		limits_ts...) const
	{
		return Op{}(expand_iterated<N0, Ns>(subs, limits_ts{}...)...); // All indices to be substituted, Is
	}


	template<typename... Ts, typename... Rest>
	auto construct_iterated(std::tuple<std::tuple<Ts...>, Rest...> const& subs) const
	{
		if constexpr (sizeof...(Is) == 0)
		{
			return std::get<1>(subs);
		}
		else
		{
			return symphas::lib::get_tuple_ge<1>(subs);
		}
	}

	template<typename... Ts>
	auto next_sub(std::tuple<std::tuple<Ts...>> const& subs)
	{
		return std::get<0>(subs);
	}

	template<typename... T0s, typename... T1s, typename... Rest>
	auto next_sub(std::tuple<std::tuple<T0s...>, std::tuple<T1s...>, Rest...> const& subs)
	{
		return symphas::lib::get_tuple_ge<1>(subs);
	}

	template<int N0, int N, typename... Ts, typename... Rest, int M0, int M1, typename... limits_ts>
	auto construct_iterated(std::tuple<std::tuple<Ts...>, Rest...> const& subs,
		expr::series_limits_ntn<M0, M1>, limits_ts...) const
	{
		constexpr int NN = N0 + N; // Current value of the major index.
		auto sum_major = expr::series<Op, expr::symbols::i_<I0, P0 + N>>(e).template expand<NN, NN>(std::get<0>(subs));
		if constexpr (M0 <= M1)
		{
			return expr::series<Op, expr::symbols::i_<I1, P1 + N>, expr::symbols::i_<Is, Ps>...>(sum_major)
				.template expand<expr::series_limits_ntn<M0, M1>, limits_ts...>(next_sub(subs));
		}
		else
		{
			return Op{}();
		}
	}

	template<int N0, int N, typename... Ts, typename... Rest, typename... limits_ts>
	auto expand_iterated(
		std::tuple<std::tuple<Ts...>, Rest...> const& subs, 
		limits_ts...) const
	{
		return construct_iterated<N0, N>(subs, expr::update_limits<N0 + N, expr::symbols::i_<I0, P0>, limits_ts>{}...); // All indices to be substituted, Is
	}

	//! Start the iterated sum expansion.
	/*!
	 * Start the iterated sum expansion.
	 *
	 * \param subs The list of expressions that is substituted, only a single list
	 * since it is assumed that the iterated sum has expressions that are indexed by only
	 * a single value, e.g. $\psi_i$.
	 *
	 * \tparam Ns... The starting index ID of an index from Is...
	 * \tparam Ps... The starting index value of an index from Is...
	 */
	template<int N0, typename... Ts, typename... Rest, int... Ns, typename... range_ts>
	auto expand_iterated(std::tuple<std::tuple<Ts...>, Rest...> const& subs,
		std::integer_sequence<int, Ns...>,
		range_ts...) const
	{
		return Op{}(expand_iterated<N0, Ns>(subs, range_ts{}...)...); // All indices to be substituted, Is
	}

public:

	E e;				//!< The expression which is summed.

};


template<typename Op, typename E>
struct SymbolicSeries<Op, E, symphas::lib::types_list<>>
{

	using i_types = symphas::internal::select_unique_i_<expr::op_types_t<E>>;

	SymbolicSeries() = default;
	SymbolicSeries(E const& e) : e{ e } {}
	
	template<size_t... Ns, typename... Xs, typename = std::enable_if_t<(sizeof...(Ns) == sizeof...(Xs)), int>>
	auto select(Xs...) const
	{
		return SymbolicSeries<Op, E, 
			symphas::lib::types_list<std::index_sequence<Ns...>,
				typename symphas::internal::filter_i_selection<Xs>::type...>>(e);
	}

	template<typename... Xs>
	auto select(Xs...) const
	{
		return select<symphas::internal::typei<0, Xs>...>(Xs{}...);
	}

	template<typename... T1s, typename... T2s>
	auto operator()(expr::series_limits<T1s, T2s> const& ...limits) const
	{
		return series_with_sorted(i_types{}, std::make_tuple(limits...),
			std::make_index_sequence<symphas::lib::types_list_size<i_types>::value>{});
	}

	template<typename... T1s, typename... T2s>
	auto operator()(std::tuple<expr::series_limits<T1s, T2s>...> const& limits) const
	{
		return construct_with_limits(limits, std::make_index_sequence<sizeof...(T1s)>{});
	}

	template<typename... T1s, typename... T2s, typename... Ts>
	auto operator()(Substitution<Ts...> const& substitution, expr::series_limits<T1s, T2s> const& ...limits) const
	{
		return series_with_sorted(i_types{}, substitution, std::make_tuple(limits...), 
			std::make_index_sequence<symphas::lib::types_list_size<i_types>::value>{});
	}

	template<typename... T1s, typename... T2s, typename... Ts>
	auto operator()(std::tuple<expr::series_limits<T1s, T2s>...> const& limits, Substitution<Ts...> const& substitution) const
	{
		return construct_with_limits(limits, substitution, std::make_index_sequence<sizeof...(T1s)>{});
	}

protected:

	template<int N0, int... Ns, typename... T1s, typename... T2s, size_t... Is>
	auto series_with_sorted(
		symphas::lib::types_list<expr::symbols::i_<N0, 0>, expr::symbols::i_<Ns, 0>...>,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		std::index_sequence<Is...>
	) const
	{
		return expr::series<Op, expr::symbols::i_<N0, 0>, expr::symbols::i_<Ns, 0>...>(e)(limits);
	}

	template<typename... T1s, typename... T2s, size_t... Is>
	auto series_with_sorted(
		symphas::lib::types_list<>,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		std::index_sequence<Is...>
	) const
	{
		return e;
	}

	template<typename... T1s, typename... T2s, size_t... Is>
	auto construct_with_limits(
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits, 
		std::index_sequence<Is...>) const
	{
		return series_with_sorted(i_types{}, limits,
			std::make_index_sequence<symphas::lib::types_list_size<i_types>::value>{});
	}

	template<int N0, int... Ns, typename... T1s, typename... T2s, typename... Ts, size_t... Is>
	auto series_with_sorted(
		symphas::lib::types_list<expr::symbols::i_<N0, 0>, expr::symbols::i_<Ns, 0>...>,
		Substitution<Ts...> const& substitution,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		std::index_sequence<Is...>
	) const
	{
		return expr::series<Op, expr::symbols::i_<N0, 0>, expr::symbols::i_<Ns, 0>...>(e)(std::get<Is>(limits)...)(substitution);
	}

	template<typename... T1s, typename... T2s, typename... Ts, size_t... Is>
	auto series_with_sorted(
		symphas::lib::types_list<>, 
		Substitution<Ts...> const& substitution,
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		std::index_sequence<Is...>
	) const
	{
		return e;
	}

	template<typename... T1s, typename... T2s, size_t... Is, typename... Ts>
	auto construct_with_limits(
		std::tuple<expr::series_limits<T1s, T2s>...> const& limits,
		Substitution<Ts...> const& substitution,
		std::index_sequence<Is...>) const
	{
		return series_with_sorted(i_types{}, substitution, limits,
			std::make_index_sequence<symphas::lib::types_list_size<i_types>::value>{});
	}


public:

	E e;
};


template<typename Op, size_t N0, size_t... Ns, typename X0, typename... Xs>
struct SymbolicSeries<Op, OpVoid, symphas::lib::types_list<std::index_sequence<N0, Ns...>, X0, Xs...>>
{
	SymbolicSeries() = default;
	SymbolicSeries(OpVoid) {}

	template<typename... Ts>
	auto expand(Ts&&... subs) const
	{
		return OpVoid{};
	}

	template<typename... Ts>
	auto operator()(Ts&&... subs) const
	{
		return OpVoid{};
	}

};

template<typename Op, typename E, size_t N0, size_t... Ns, typename X0, typename... Xs>
struct SymbolicSeries<Op, E, symphas::lib::types_list<std::index_sequence<N0, Ns...>, X0, Xs...>>
{
	SymbolicSeries() = default;
	SymbolicSeries(E const& e) : e{ e } {}

	template<typename... call_ts, typename... Ts>
	auto expand(symphas::lib::types_list<call_ts...> calls, Ts&&... subs) const
	{
		return Op{}(call_ts{}.expand(e, std::forward<Ts>(subs)...)...);
	}

	template<typename... Ts>
	auto expand(Ts&&... subs) const
	{
		using call_list = typename symphas::internal::expand_series<
			std::index_sequence<N0, Ns...>,
			SymbolicSeries<Op, E, symphas::lib::types_list<>>,
			symphas::lib::types_list<X0, Xs...>>::all_types;
		return this->operator()(call_list{}, std::forward<Ts>(subs)...);
	}

	template<typename... call_ts, typename... Ts>
	auto operator()(symphas::lib::types_list<call_ts...> calls, Ts&&... subs) const
	{
		return Op{}(call_ts{}(e, std::forward<Ts>(subs)...)...);
	}

	template<typename... Ts>
	auto operator()(Ts&&... subs) const
	{
		using call_list = typename symphas::internal::expand_series<
			std::index_sequence<N0, Ns...>,
			SymbolicSeries<Op, E, symphas::lib::types_list<>>,
			symphas::lib::types_list<X0, Xs...>>::all_types;
		return this->operator()(call_list{}, std::forward<Ts>(subs)...);
	}

	E e;
};



template<typename Op, typename S, size_t... Ns, int... I0s, int... P0s, 
	typename... T1s, typename... T2s, typename... Is>
struct SymbolicSeries<Op, SymbolicTemplate<S, Ns...>, 
	symphas::lib::types_list<OpVoid,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>>
{
	SymbolicSeries() = default;
	SymbolicSeries(OpVoid, SymbolicTemplate<S, Ns...> const& tmpl, 
		const DynamicIndex(&index)[sizeof...(I0s)], expr::series_limits<T1s, T2s> const&... limits) {}

	template<typename... Ts>
	auto operator()(Ts&&... subs) const
	{
		return OpVoid{};
	}
};

template<typename Op, typename E, int... I0s, int... P0s, typename... T1s, typename... T2s, typename... Is>
struct SymbolicSeries<Op, void,
	symphas::lib::types_list<E,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
		symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>
	>>
{
	SymbolicSeries() = default;
	SymbolicSeries(
		E const& e,
		expr::series_limits<T1s, T2s> const&... limits) :
			e{ e }, limits{ limits... },
			dims{ symphas::internal::limit_dimension(std::integer_sequence<int, I0s...>{}, this->limits, limits)... }
	{}

	E e;
	std::tuple<expr::series_limits<T1s, T2s>...> limits;
	len_type dims[sizeof...(I0s)];
};


template<typename Op, typename E, typename S, size_t... Ns, int... I0s, int... P0s, 
	typename... T1s, typename... T2s, typename... Is>
struct SymbolicSeries<Op, SymbolicTemplate<S, Ns...>, 
	symphas::lib::types_list<E,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>
	>> :
	SymbolicSeries<Op, void,
		symphas::lib::types_list<E,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
			symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
			symphas::lib::types_list<expr::symbols::v_id_type<Is>...>
	>>
{

	using parent_type = SymbolicSeries<Op, void,
		symphas::lib::types_list<E,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>, 
			symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
			symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>>;
	using parent_type::e;
	using parent_type::limits;
	using parent_type::dims;


	SymbolicSeries() = default;
	SymbolicSeries(E const& e, SymbolicTemplate<S, Ns...> const& tmpl, 
		const DynamicIndex(&index)[sizeof...(I0s)], expr::series_limits<T1s, T2s> const&... limits)
		: parent_type(e, limits...), index{}, tmpl{ tmpl } 
	{
		std::copy(index, index + sizeof...(I0s), this->index);
	}


	template<typename... Ts>
	auto operator()(Ts&& ...datas) const
	{
		return substitute_args(std::forward<Ts>(datas)...);
	}


protected:

	template<typename... Ts>
	auto substitute_args(Substitution<SymbolicDataArray<Ts>...> const& substitution) const
	{
		using v_id_types = symphas::lib::types_list<expr::symbols::v_id_type<Is>...>;
		using id_types = symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>;

		auto args = symphas::internal::expanded_arg_list(substitution, v_id_types{}, id_types{});
		return make_eval(substitution, args);
	}

	template<typename... Ts>
	auto substitute_args(Ts const& ...datas) const
	{
		using placeholder_mask_t = std::integer_sequence<bool,
			(symphas::lib::types_list_size<
				symphas::internal::select_v_i_<
					expr::symbols::i_<I0s, 0>,
					symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>
				>::value > 0)...
			>;

		auto substitution = get_substitution(
			std::make_index_sequence<sizeof...(I0s)>{},
			placeholder_mask_t{},
			datas...);

		return substitute_args(substitution);
	}

	template<typename... Ts, size_t... Ms, bool... Bs>
	auto get_substitution(std::index_sequence<Ms...>, std::integer_sequence<bool, Bs...>, Ts const& ...datas) const
	{
		return get_substitution(
			std::index_sequence<Ms...>{},
			build_arg_list<0>(
				std::index_sequence<Ms...>{}, 
				std::integer_sequence<bool, Bs...>{},
				datas...));
	}

	template<typename... Ts, size_t... Ms>
	auto get_substitution(std::index_sequence<Ms...>, std::tuple<SymbolicDataArray<Ts>...> const& datas) const
	{
		return expr::arg_list(std::get<Ms>(datas)...);
	}

	template<size_t N>
	auto build_arg_list(std::index_sequence<>, std::integer_sequence<bool>) const
	{
		return std::make_tuple();
	}

	template<size_t N, size_t M0, size_t... Ms, bool B0, bool... Bs>
	auto build_arg_list(std::index_sequence<M0, Ms...>, std::integer_sequence<bool, B0, Bs...>) const
	{
		using variable_list = symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I0s, 0>>...>;

		return std::make_tuple(SymbolicDataArray<symphas::lib::type_at_index<M0, symphas::lib::unroll_types_list<variable_list>>>(),
			SymbolicDataArray<symphas::lib::type_at_index<Ms, symphas::lib::unroll_types_list<variable_list>>>()...);
	}

	template<size_t N, typename T0, typename... Ts>
	auto build_arg_list(std::index_sequence<>, std::integer_sequence<bool>, T0&& data0, Ts&& ...datas) const 
	{
		return std::make_tuple();
	}

	template<size_t N, typename T0, typename... Ts, size_t M0, size_t... Ms, bool B0, bool... Bs>
	auto build_arg_list(std::index_sequence<M0, Ms...>, std::integer_sequence<bool, B0, Bs...>, T0&& data0, Ts&& ...datas) const
	{
		if constexpr (B0)
		{
			auto&& next_list = build_arg_list<N + 1>(
				std::index_sequence<Ms...>{},
				std::integer_sequence<bool, Bs...>{},
				std::forward<Ts>(datas)...);

			return std::tuple_cat(
				std::make_tuple(symphas::internal::as_array_arg(std::forward<T0>(data0), dims[M0])),
				next_list);
		}
		else
		{
			auto&& next_list = build_arg_list<N + 1>(
				std::index_sequence<Ms...>{},
				std::integer_sequence<bool, Bs...>{},
				std::forward<T0>(data0), std::forward<Ts>(datas)...);

			using curr_i_type = expr::symbols::i_<symphas::lib::seq_index_value<N, std::integer_sequence<int, I0s...>>::value, 0>;
			return std::tuple_cat(
				std::make_tuple(SymbolicDataArray<expr::symbols::v_id_type<curr_i_type>>()),
				next_list);
		}
	}

	template<typename sub_t, typename... Ts, size_t... Ms, size_t... Ls>
	auto make_eval(sub_t const& substitution, std::tuple<Ts...> const& args, 
		std::index_sequence<Ms...>, std::index_sequence<Ls...>) const
	{
		using v_id_types = symphas::lib::types_list<expr::symbols::v_id_type<Is>...>;
		using id_types = symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>;

		auto f = to_function(std::get<Ms>(args)...);
		
		auto series = symphas::internal::construct_series<Op>(
			e, index, substitution, 
			id_types{}, 
			v_id_types{}, 
			std::get<Ls>(limits)...);
		
		return symphas::internal::make_symbolic_eval(
			OpIdentity{}, series, f);
	}

	template<typename sub_t, typename... Ts>
	auto make_eval(sub_t const& substitution, std::tuple<Ts...> const& args) const
	{
		return make_eval(substitution, args, 
			std::make_index_sequence<sizeof...(Ts)>{}, std::make_index_sequence<sizeof...(T1s)>{});
	}

	template<size_t N>
	auto update_name(int const& arg) const
	{
		return arg;
	}

#ifdef PRINTABLE_EQUATIONS
	template<size_t N, typename T>
	auto update_name(T const& arg) const
	{
		using v_named_t = symphas::lib::type_at_index<N, expr::symbols::v_id_type<Is>...>;
		NamedData arg0(T(arg), expr::get_op_name(v_named_t{}));
		return arg0;
	}

	template<size_t N, typename T>
	auto update_name(NamedData<T> const& arg) const
	{
		//constexpr int I = symphas::lib::seq_index_value<N, std::integer_sequence<int, I0s...>>::value;
		using v_named_t = symphas::lib::type_at_index<N, expr::symbols::v_id_type<Is>...>;
		NamedData arg0(T(arg), expr::get_op_name(v_named_t{}));
		return arg0;
	}
#endif

	template<typename E0, typename... Ts>
	SymbolicFunction<E0, symphas::lib::types_list<Ts...>> _to_function(
		SymbolicFunction<E0, Ts...> const& f) const
	{
		return f;
	}

	template<typename... Ts>
	auto to_function(Ts const&... args) const
	{
		auto f = (expr::function_of(expr::symbols::arg<Ns, Ts>...) = tmpl);
		f.set_data(update_name<symphas::lib::index_of_value<size_t, Ns, Ns...>>(args)...);
		return _to_function(f);
	}
	
	DynamicIndex index[sizeof...(I0s)];
	SymbolicTemplate<S, Ns...> tmpl;
};



template<typename Op, typename... Ts, typename E, int... I0s, int... P0s,
	typename... T1s, typename... T2s, typename... Is>
struct SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
	symphas::lib::types_list<E,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>
	>> :
	SymbolicSeries<Op, void,
		symphas::lib::types_list<E,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
			symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
			symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>>
{

	using parent_type = SymbolicSeries<Op, void,
		symphas::lib::types_list<E,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>>;
	using parent_type::e;
	using parent_type::limits;
	using parent_type::dims;

	using this_type = SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>,
		symphas::lib::types_list<E,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>
		>>;

	SymbolicSeries() = default;
	SymbolicSeries(
		E const& e,
		Substitution<SymbolicDataArray<Ts>...> const& substitution,
		const DynamicIndex(&index)[sizeof...(I0s)],
		expr::series_limits<T1s, T2s> const&... limits) :
		parent_type(e, limits...), substitution{ substitution },
		persistent{ sizeof...(I0s), get_length(std::make_index_sequence<sizeof...(I0s)>{}) }
	{
		std::copy(index, index + sizeof...(I0s), this->index);
	}


protected:

	template<int I0, int... P00s>
	auto offsets_for_index(
		symphas::lib::types_list<expr::symbols::v_id_type<expr::symbols::i_<I0, P00s>>...>) const
	{
		return std::integer_sequence<int, P00s...>{};
	}


	auto offsets_for_index(
		symphas::lib::types_list<>) const
	{
		return std::integer_sequence<int>{};
	}

	template<size_t... Ms>
	auto index_values(
		iter_type(&values)[sizeof...(I0s)],
		iter_type(&offsets)[sizeof...(I0s)],
		std::index_sequence<Ms...>) const
	{
		((values[Ms] = offsets[Ms] + symphas::internal::check_limit_range(
			symphas::internal::limit_start(std::integer_sequence<int, I0s...>{}, limits, std::get<Ms>(limits))
		)), ...);
	}

	template<typename E0, typename... T0s>
	auto eval_series(SymbolicFunctionArray<E0, T0s...> const& persistent, iter_type n) const
	{
		using namespace symphas::internal;

		if (persistent.len > 0)
		{
			auto result = persistent[0][n];

			for (iter_type i = 1; i < persistent.len; ++i)
			{
				result = Op{}(result, persistent[i][n]);
			}
			return result;
		}
		else
		{
			return Op{}(typename SymbolicFunctionArray<E0, T0s...>::eval_type{});
		}
	}

	template<typename E0, typename... T0s>
	auto eval_series(SymbolicFunctionArray<E0, T0s...> const& persistent) const
	{
		using namespace symphas::internal;
		if (persistent.len > 0)
		{
			auto result = persistent[0]();
			expr::result(Op{}(), result);
			auto op = expr::make_term(result);

			for (iter_type i = 0; i < persistent.len; ++i)
			{
				expr::result(Op{}(op, expr::make_term(persistent[i]())), result);
			}
			return result;
		}
		else
		{
			return Op{}(typename SymbolicFunctionArray<E0, T0s...>::eval_type{});
		}
	}

	//! Evaluates for only one element.
	template<size_t N, typename E0, typename... T0s>
	auto build_series_evals(SymbolicFunction<E0, T0s...>& f, iter_type(&offsets)[sizeof...(I0s)],
		SymbolicFunction<E0, T0s...>** fs, iter_type** set_offsets) const
	{
		using namespace symphas::internal;
		using std::swap;

		using id_types = symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>;
		using v_id_types = symphas::lib::types_list<expr::symbols::v_id_type<Is>...>;

		len_type pos = 0;

		if constexpr (N >= sizeof...(I0s))
		{
			return pos;
		}
		else
		{
			auto& index_mutable = const_cast<DynamicIndex(&)[sizeof...(I0s)]>(index);


			iter_type values[sizeof...(I0s)]{};
			index_values(values, offsets, std::make_index_sequence<sizeof...(I0s)>{});

			auto lim = std::get<N>(limits);
			auto offset = compute_offset(std::integer_sequence<int, I0s...>{}, limits, lim, offsets)
				+ symphas::lib::seq_index_value<N, std::integer_sequence<int, P0s...>>::value;

			offsets[N] = expr::eval(offset);
			auto range = limit_range(std::integer_sequence<int, I0s...>{}, limits, lim, offsets)/* - offset*/;

			for (iter_type i = 0; i < range; ++i)
			{
				if constexpr (N == sizeof...(I0s) - 1)
				{
					for (iter_type n = 0; n < sizeof...(I0s); ++n)
					{
						index_mutable[n] = (int)offsets[n];
					}

					auto args0 = symphas::internal::expanded_arg_list(substitution, v_id_types{}, id_types{}, offsets, values);
					std::copy(offsets, offsets + sizeof...(I0s), set_offsets[pos]);
					
					fs[pos]->set_data_tuple(args0);
					expr::transform::fix_coeffs(fs[pos]->e);

					++pos;
				}
				else
				{
					pos += build_series_evals<N + 1>(f, offsets, fs + pos, set_offsets);
				}
				offsets[N] += 1;
				values[N] += 1;
			}
			offsets[N] = expr::eval(offset);
		}
		return pos;
	}

	template<size_t... Ns>
	len_type get_length(std::index_sequence<Ns...>) const
	{
		return (dims[Ns] * ...);
	}

public:


	template<typename E0, typename... T0s>
	void update_function_list(SymbolicFunction<E0, T0s...>& f, 
		SymbolicFunctionArray<E0, T0s...>& persistent) const
	{
		iter_type offsets[sizeof...(I0s)]{};
		auto len = build_series_evals<0>(f, offsets, persistent.data, persistent.offsets);
		persistent.len = len;
	}

	template<typename E0, typename... T0s>
	void update(SymbolicFunction<E0, T0s...>& f) 
	{
		persistent.init(f);
		update_function_list(f, persistent);
	}

	template<typename... T0s>
	auto grid_symbols(symphas::lib::types_list<T0s...>) const
	{
		return symphas::lib::types_list<GridSymbol<expr::symbols::v_id_type<Is>, expr::grid_dim<T0s>::value>...>{};
	}


	template<int... I00s, int... P00s, size_t... Ds, int... I01s, int... P01s>
	auto substitute_placeholders(
		symphas::lib::types_list<GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I00s, P00s>>, Ds>...>,
		symphas::lib::types_list<expr::symbols::i_<I01s, P01s>...>) const
	{
		if constexpr (sizeof...(I00s) > 0)
		{
			return expr::apply_operators(
				expr::transform::swap_grid<expr::symbols::v_id_type<expr::symbols::i_<I00s, P00s>>...>
				(e, GridSymbol<expr::symbols::v_id_type<expr::symbols::i_<I00s, P00s>>, Ds>{}...));
		}
		else
		{
			return e;
		}
	}

	template<typename E0, typename... T0s>
	auto substitute_placeholders(SymbolicFunction<E0, T0s...> const& f) const
	{
		using id_types = symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>;
		using t_list = symphas::lib::types_before_index<sizeof...(T0s) - sizeof...(I0s), T0s...>;

		auto vs = grid_symbols(t_list{});
		auto is = id_types{};

		return substitute_placeholders(vs, is);

		//return f[symphas::internal::placeholder_list(substitution, grid_symbols(t_list{}), id_types{})];
	}

	auto eval() const
	{
		return eval_series(persistent);
	}

	auto eval(iter_type n) const
	{
		return eval_series(persistent, n);
	}

	Substitution<SymbolicDataArray<Ts>...> substitution;

	using func_list_t = std::invoke_result_t<
		decltype(&symphas::internal::get_function_list<
			SymbolicSeries<Op, E, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>,
			std::tuple<expr::series_limits<T1s, T2s>...>,
			Substitution<SymbolicDataArray<Ts>...>>),
		SymbolicSeries<Op, E, symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>>,
		std::tuple<expr::series_limits<T1s, T2s>...>,
		Substitution<SymbolicDataArray<Ts>...>>;

	func_list_t persistent;
	DynamicIndex index[sizeof...(I0s)];
};




namespace expr
{

	//! Defines a series using the given operation iterated over the specified indices.
	/*!
	 * A series is repeated over a given index, which is of type
	 * expr::symbols::i_, with an `int` template parameter indicating the "value" of the index, 
	 * where multiple indices are supported in order to allow iterated series.
	 *
	 * A series is expanded by calling its member function SymbolicSeries::expand.
	 * Terms or other expressions can also be substituted into series when they are expanded.
	 *
	 * A series is written like a regular expression, with the additional use of expr::symbols::i_ 
	 * to substitute the current series index, as well of the expr::symbols::v_, 
	 * which takes the series's corresponding index, expr::symbols::i_, in order to allow terms or 
	 * expressions to be substituted into the expansion of the series.
	 *
	 * Multiple indices can be defined, which indicates multiple sums.
	 * 
	 * See \ref symbolicseries for more information.
	 */
	template<typename Op, typename... Is, typename E, typename>
	auto series(OpExpression<E> const& e)
	{
		return SymbolicSeries<Op, E, symphas::lib::types_list<Is...>>(*static_cast<E const*>(&e));
	}


	//! Defines a summation of an expression with indices and substitutable expressions.
	/*!
	 * See series().
	 */
	template<typename... Is, typename E, typename>
	auto sum(OpExpression<E> const& e)
	{
		using Op = symphas::internal::ReduceOp<symphas::internal::SeriesOp::ADD>;
		return series<Op, Is...>(*static_cast<E const*>(&e));
	}


	//! Defines a product of an expression with indices and substitutable expressions.
	/*!
	 * See series().
	 */
	template<typename... Is, typename E, typename>
	auto prod(OpExpression<E> const& e)
	{
		using Op = symphas::internal::ReduceOp<symphas::internal::SeriesOp::MUL>;
		return series<Op, Is...>(*static_cast<E const*>(&e));
	}



}


namespace symphas::internal
{

	template<typename Op, typename E, typename S, size_t... Ns, 
		typename... T1s, typename... T2s, int... I0s, int... P0s, typename... Is>
	auto construct_series(
		OpExpression<E> const& e,
		const DynamicIndex(&index)[sizeof...(I0s)],
		SymbolicTemplate<S, Ns...> const& tmpl,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>,
		expr::series_limits<T1s, T2s> const&... limits)
	{
		using construct_series_list_t = symphas::lib::types_list<E,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
			symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
			symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>;
		return SymbolicSeries<Op, SymbolicTemplate<S, Ns...>, construct_series_list_t>(*static_cast<E const*>(&e), tmpl, index, limits...);
	}

	template<typename Op, typename E, typename... Ts, 
		typename... T1s, typename... T2s, int... I0s, int... P0s, typename... Is>
	auto construct_series(
		OpExpression<E> const& e,
		const DynamicIndex(&index)[sizeof...(I0s)],
		Substitution<SymbolicDataArray<Ts>...> const& substitution,
		symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
		symphas::lib::types_list<expr::symbols::v_id_type<Is>...>,
		expr::series_limits<T1s, T2s> const&... limits)
	{
		using construct_series_list_t = symphas::lib::types_list<E,
			symphas::lib::types_list<expr::symbols::i_<I0s, P0s>...>,
			symphas::lib::types_list<expr::series_limits<T1s, T2s>...>,
			symphas::lib::types_list<expr::symbols::v_id_type<Is>...>>;
		return SymbolicSeries<Op, Substitution<SymbolicDataArray<Ts>...>, construct_series_list_t>(*static_cast<E const*>(&e), substitution, index, limits...);
	}
}




 //! @}




