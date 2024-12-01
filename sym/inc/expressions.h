
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
 * PURPOSE: Defines the foundational symbolic algebra elements, such as the
 * base expression object and the binary operations.
 *
 * ***************************************************************************
 */

#pragma once

#include <iostream>

#include "expressionproperties.h"
#include "expressionsprint.h"


// ******************************************************************************************

//! \cond



#ifdef LATEX_PLOT

#define SYEX_MUL_SEP_OP " "
#define SYEX_DIV_SEP_OP "}{"
#define SYEX_NLV_SEP_OP " "

//#define SYEX_MUL_FMT_AA "\\left["
#define SYEX_MUL_FMT_AA ""
#define SYEX_MUL_FMT_AB "\\left("
#define SYEX_MUL_FMT_BA "\\right)"
//#define SYEX_MUL_FMT_BB "\\right]"
#define SYEX_MUL_FMT_BB ""

#define SYEX_DIV_SEP SYEX_DIV_SEP_OP
#define SYEX_DIV_FMT_A "\\frac{"
#define SYEX_DIV_FMT_B "}" 

#define SYEX_SCIENTIFIC_A SYEX_MUL_FMT_AB
#define SYEX_SCIENTIFIC_B SYEX_MUL_FMT_BA
#define SYEX_SCIENTIFIC_VALUE_A "\\times 10^{"
#define SYEX_SCIENTIFIC_VALUE_B "}"


#else

#define SYEX_MUL_SEP_OP "*"
#define SYEX_DIV_SEP_OP "/"
#define SYEX_NLV_SEP_OP "*"

//#define SYEX_MUL_FMT_AA "["
#define SYEX_MUL_FMT_AA ""
#define SYEX_MUL_FMT_AB "("
#define SYEX_MUL_FMT_BA ")"
//#define SYEX_MUL_FMT_BB "]"
#define SYEX_MUL_FMT_BB ""

#define SYEX_DIV_SEP  ") " SYEX_DIV_SEP_OP " ("
#define SYEX_DIV_FMT_A SYEX_MUL_FMT_AB
#define SYEX_DIV_FMT_B SYEX_MUL_FMT_BA

#define SYEX_SCIENTIFIC_A SYEX_MUL_FMT_AB
#define SYEX_SCIENTIFIC_B SYEX_MUL_FMT_BA
#define SYEX_SCIENTIFIC_VALUE_A "E"
#define SYEX_SCIENTIFIC_VALUE_B ""

#endif

#define SYEX_MUL_SEP  SYEX_MUL_FMT_BA SYEX_MUL_SEP_OP SYEX_MUL_FMT_AB
#define SYEX_MUL_SEP_A  SYEX_MUL_FMT_BA SYEX_MUL_SEP_OP
#define SYEX_MUL_SEP_B  SYEX_MUL_SEP_OP SYEX_MUL_FMT_AB
#define SYEX_MUL_SEP_AB  SYEX_MUL_SEP_OP

#define SYEX_MUL_FMT_A SYEX_MUL_FMT_AA SYEX_MUL_FMT_AB
#define SYEX_MUL_FMT_B SYEX_MUL_FMT_BA SYEX_MUL_FMT_BB

#define SYEX_MUL_FMT SYEX_MUL_FMT_A "%s" SYEX_MUL_SEP "%s" SYEX_MUL_FMT_B
#define SYEX_MUL_FMT_LEN (STR_ARR_LEN(SYEX_MUL_FMT_A SYEX_MUL_SEP SYEX_MUL_FMT_B) - 1)

#define SYEX_DIV_FMT SYEX_DIV_FMT_A "%s" SYEX_DIV_SEP "%s" SYEX_DIV_FMT_B
#define SYEX_DIV_FMT_LEN (STR_ARR_LEN(SYEX_DIV_FMT_A SYEX_DIV_SEP SYEX_DIV_FMT_B) - 1)


#define SYEX_FRA_SEP SYEX_DIV_SEP_OP
#define SYEX_FRA_FMT_A SYEX_DIV_FMT_A
#define SYEX_FRA_FMT_B SYEX_MUL_FMT_B 

#define SYEX_FRA_FMT SYEX_DIV_FMT_A "%zd" SYEX_DIV_SEP_OP "%zd" SYEX_DIV_FMT_B
#define SYEX_FRA_FMT_LEN (STR_ARR_LEN(SYEX_DIV_FMT_A SYEX_DIV_SEP_OP SYEX_DIV_FMT_B) - 1)

#define SYEX_LAMBDA_FUNC_FMT_A "%s("
#define SYEX_LAMBDA_FUNC_FMT_B ")"
#define SYEX_LAMBDA_FUNC_FMT SYEX_DIV_FMT_A SYEX_LAMBDA_FUNC_FMT_A "%s" SYEX_LAMBDA_FUNC_FMT_B
#define SYEX_LAMBDA_FUNC_FMT_LEN (STR_ARR_LEN(SYEX_LAMBDA_FUNC_FMT_A SYEX_LAMBDA_FUNC_FMT_B) - 3)

#ifdef LATEX_PLOT
#define SYEX_POW_SEP_A "^{"
#define SYEX_POW_SEP_B "}"
#define SYEX_POW_DIV_SEP "/"
#else
#define SYEX_POW_SEP_A "^"
#define SYEX_POW_SEP_B ""
#define SYEX_POW_DIV_SEP "/"
#endif



//! \endcond


//! The maximum supported power that an expression can be evaluated at.
#define MAX_EXPONENT 20


//! The display precision of values which appear in expressions.
/*!
 * Defines the number of digits which appear in floating point numbers
 * for values from expressions that are printed to the screen.
 */
#define EXPR_VALUE_DISPLAY_PRECISION 3


// ******************************************************************************************

/*
 * evaluate the expression in every index in the array
 */

namespace expr
{

#define PARALLELIZATION_CUTOFF_COUNT 1000

	struct forward_value
	{
		template<typename T>
		decltype(auto) operator()(T&& value)
		{
			return std::forward<T>(value);
		}
	};

	//! Evaluate the expression into the underlying data member.
	/*!
	 * The expression must be iterable over the entire given length.
	 * 
	 * \param e Expression that is evaluated.
	 * \param data The array containing the result of the expression.
	 * \param len The number of elements in the array.
	 */
	template<typename E, typename assign_type>
	void result(OpEvaluable<E> const& e, assign_type&& data, len_type len)
	{
		symphas::data_iterator it(std::forward<assign_type>(data));
		
#ifdef EXECUTION_HEADER_AVAILABLE
		if (params::parallelization)
			std::transform(
				std::execution::par_unseq,
				static_cast<const E*>(&e)->begin(),
				static_cast<const E*>(&e)->end(len), it,
				forward_value{});
		else
#endif
			std::transform(
				static_cast<const E*>(&e)->begin(),
				static_cast<const E*>(&e)->end(len), it,
				forward_value{});
	}

	template<typename E, typename assign_type, size_t D>
	void result(OpEvaluable<E> const& e, assign_type&& data, grid::region_interval<D> const& interval)
	{
		symphas::data_iterator_group it(std::forward<assign_type>(data), interval);

#ifdef EXECUTION_HEADER_AVAILABLE
		if (params::parallelization)
			std::transform(
				std::execution::par_unseq,
				static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
				static_cast<const E*>(&e)->end(symphas::it_grp, interval), it,
				forward_value{});
		else
#endif
			std::transform(
				static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
				static_cast<const E*>(&e)->end(symphas::it_grp, interval), it,
				forward_value{});

	}

	template<typename E, typename assign_type, size_t D>
	void result(OpEvaluable<E> const& e, assign_type&& data, grid::region_interval_multiple<D> const& regions)
	{
		for (grid::region_interval<D> region : regions)
		{
			result(*static_cast<E const*>(&e), std::forward<assign_type>(data), region);
		}
	}


	template<typename E, typename assign_type>
	void result(OpEvaluable<E> const& e, assign_type&& data, grid::region_interval<0> const& interval)
	{
		auto data_region = expr::iterable_domain(std::forward<assign_type>(data));
		if (grid::length(data_region) > 1)
		{
			result(*static_cast<E const*>(&e), std::forward<assign_type>(data), data_region);
		}
		else
		{
			result(*static_cast<E const*>(&e), std::forward<assign_type>(data), 1);
		}
	}

	template<typename E, typename assign_type>
	void result(OpEvaluable<E> const& e, assign_type&& data, grid::region_empty) {}

	template<typename E, typename assign_type>
	void result(OpEvaluable<E> const& e, assign_type&& data)
	{
		result(*static_cast<E const*>(&e), std::forward<assign_type>(data), expr::iterable_domain(*static_cast<E const*>(&e)));
	}


	template<typename L, typename R>
	void result(std::pair<L, R> const& evaluate)
	{
		auto&& [lhs, rhs] = evaluate;
		result(rhs, lhs);
	}

	template<typename G, typename R>
	void result(std::pair<OpTerm<OpIdentity, G>, R> const& evaluate)
	{
		auto&& [lhs, rhs] = evaluate;
		result(rhs, BaseData<G>::get(expr::get<1>(lhs).data()));
	}

	//! Add the result of the expression into the underlying data member.
	/*!
	 * The expression is evaluated and the result is added to the existing
	 * values in the data array.
	 * 
	 * \param e Expression that is evaluated.
	 * \param data The array of data.
	 * \param len The length of the array.
	 */
	template<typename E, typename assign_type>
	void result_accumulate(OpEvaluable<E> const& e, assign_type&& data, len_type len)
	{
		symphas::data_iterator it(std::forward<assign_type>(data));

#ifdef EXECUTION_HEADER_AVAILABLE
		if (params::parallelization)
			std::transform(
				std::execution::par_unseq,
				static_cast<const E*>(&e)->begin(), 
				static_cast<const E*>(&e)->end(len), it, it,
				[](auto expr_value, auto data_value) { return data_value + expr_value; });
		else
#endif
			std::transform(
				static_cast<const E*>(&e)->begin(),
				static_cast<const E*>(&e)->end(len), it, it,
				[] (auto expr_value, auto data_value) { return data_value + expr_value; });
	}

	template<typename E, typename assign_type, size_t D>
	void result_accumulate(OpEvaluable<E> const& e, assign_type&& data, grid::region_interval<D> const& interval)
	{
		symphas::data_iterator_group it(std::forward<assign_type>(data), interval);

		if (grid::length<D>(interval) <= PARALLELIZATION_CUTOFF_COUNT)
		{
			auto start = static_cast<const E*>(&e)->begin(symphas::it_grp, interval);
			auto end = static_cast<const E*>(&e)->end(symphas::it_grp, interval);

			for (auto eit = start; eit < end; ++eit, ++it)
			{
				*it = *it + *eit;
			}
		}
		else
		{
#ifdef EXECUTION_HEADER_AVAILABLE
			if (params::parallelization)
				std::transform(
					std::execution::par_unseq,
					static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
					static_cast<const E*>(&e)->end(symphas::it_grp, interval), it, it,
					[] (auto expr_value, auto data_value) { return data_value + expr_value; });
			else
#endif
				std::transform(
					static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
					static_cast<const E*>(&e)->end(symphas::it_grp, interval), it, it,
					[] (auto expr_value, auto data_value) { return data_value + expr_value; });
		}

	}

	template<typename V, typename... Gs, expr::exp_key_t... Xs, typename assign_type, size_t D>
	void result_accumulate(OpTerms<V, Term<Gs, Xs>...> const& e, assign_type&& data, grid::region_interval<D> const& interval)
	{
		symphas::data_iterator_group it(std::forward<assign_type>(data), interval);

		auto start = e.begin(symphas::it_grp, interval);
		auto end = e.end(symphas::it_grp, interval);

		for (auto eit = start; eit < end; ++eit, ++it)
		{
			*it = *it + *eit;
		}
	}

	template<typename assign_type>
	void result_accumulate(OpVoid, assign_type&& data) {}
	template<typename assign_type>
	void result_accumulate(OpVoid, assign_type&& data, len_type) {}
	template<typename assign_type>
	void result_accumulate(OpVoid, assign_type&& data, grid::region_interval<0>) {}
	template<typename assign_type, size_t D>
	void result_accumulate(OpVoid, assign_type&& data, grid::region_interval<D>) {}
	template<typename assign_type, size_t D>
	void result_accumulate(OpVoid, assign_type&& data, grid::region_interval_multiple<D>) {}

	template<typename V, typename... Gs, expr::exp_key_t... Xs, typename assign_type, size_t D>
	void result_accumulate(OpTerms<V, Term<Gs, Xs>...> const& e, assign_type&& data, grid::region_interval_multiple<D> const& regions)
	{
		for (grid::region_interval<D> region : regions)
		{
			result_accumulate(e, std::forward<assign_type>(data), region);
		}
	}

	template<typename E, typename assign_type, size_t D>
	void result_accumulate(OpEvaluable<E> const& e, assign_type&& data, grid::region_interval_multiple<D> const& regions)
	{
		for (grid::region_interval<D> region : regions)
		{
			result_accumulate(*static_cast<E const*>(&e), std::forward<assign_type>(data), region);
		}
	}

	template<typename E, typename assign_type>
	void result_accumulate(OpEvaluable<E> const& e, assign_type&& data, grid::region_interval<0> const& interval)
	{
		auto data_region = expr::iterable_domain(std::forward<assign_type>(data));
		if (grid::length(data_region) > 1)
		{
			result_accumulate(*static_cast<E const*>(&e), std::forward<assign_type>(data), data_region);
		}
		else
		{
			result_accumulate(*static_cast<E const*>(&e), std::forward<assign_type>(data), 1);
		}
	}

	template<typename V, typename... Gs, expr::exp_key_t... Xs, typename assign_type>
	void result_accumulate(OpTerms<V, Term<Gs, Xs>...> const& e, assign_type&& data, grid::region_interval<0> const& interval)
	{
		auto data_region = expr::iterable_domain(std::forward<assign_type>(data));
		if (grid::length(data_region) > 1)
		{
			result_accumulate(e, std::forward<assign_type>(data), data_region);
		}
		else
		{
			result_accumulate(e, std::forward<assign_type>(data), 1);
		}
	}

	template<typename E, typename assign_type>
	void result_accumulate(OpEvaluable<E> const& e, assign_type&& data, grid::region_empty) {}

	template<typename E, typename assign_type>
	void result_accumulate(OpEvaluable<E> const& e, assign_type&& data)
	{
		result_accumulate(*static_cast<E const*>(&e), std::forward<assign_type>(data), expr::iterable_domain(*static_cast<E const*>(&e)));
	}


	//! Evaluate the expression into the underlying data member.
	/*!
	 * The expression must be iterable over the entire given length.
	 *
	 * \param e Expression that is evaluated.
	 * \param data The array containing the result of the expression.
	 * \param len The number of elements in the array.
	 */
	template<typename E>
	auto result_sum(OpEvaluable<E> const& e, len_type len)
	{
#ifdef EXECUTION_HEADER_AVAILABLE
		if (params::parallelization)
			return std::reduce(
				std::execution::par_unseq,
				static_cast<const E*>(&e)->begin(),
				static_cast<const E*>(&e)->end(len));
		else
#endif
			return std::reduce(
				static_cast<const E*>(&e)->begin(),
				static_cast<const E*>(&e)->end(len));
	}

	inline auto result_sum(OpVoid) { return OpVoid{}; }
	inline auto result_sum(OpVoid, len_type) { return OpVoid{}; }
	inline auto result_sum(OpVoid, grid::region_interval<0>) { return OpVoid{}; }
	template<size_t D> inline auto result_sum(OpVoid, grid::region_interval<D>) { return OpVoid{}; }
	template<size_t D> inline auto result_sum(OpVoid, grid::region_interval_multiple<D>) { return OpVoid{}; }

	template<typename E>
	auto result_sum(OpEvaluable<E> const& e, grid::region_interval<0> const& interval)
	{
		return result_sum(*static_cast<E const*>(&e), 1);
	}

	template<typename V, typename... Gs, expr::exp_key_t... Xs, size_t D>
	auto result_sum(OpTerms<V, Term<Gs, Xs>...> const& e, grid::region_interval<D> const& interval)
	{
		auto start = e.begin(symphas::it_grp, interval);
		auto end = e.end(symphas::it_grp, interval);

		auto reduce = expr::eval_type_t<OpTerms<V, Term<Gs, Xs>...>>{};
		for (auto it = start; it < end; ++it)
		{
			reduce = reduce + *it;
		}
		return reduce;
	}

	template<typename V, typename... Gs, expr::exp_key_t... Xs, size_t D>
	auto result_sum(OpTerms<V, Term<Gs, Xs>...> const& e, grid::region_interval_multiple<D> const& regions)
	{
		expr::eval_type_t<OpTerms<V, Term<Gs, Xs>...>> sum{};
		for (grid::region_interval<D> region : regions)
		{
			sum += result_sum(e, region);
		}
		return sum;
	}

	template<typename E, size_t D>
	auto result_sum(OpEvaluable<E> const& e, grid::region_interval<D> const& interval)
	{
		if (grid::length<D>(interval) <= PARALLELIZATION_CUTOFF_COUNT)
		{
			auto start = static_cast<const E*>(&e)->begin(symphas::it_grp, interval);
			auto end = static_cast<const E*>(&e)->end(symphas::it_grp, interval);

			auto reduce = expr::eval_type_t<E>{};
			for (auto it = start; it < end; ++it)
			{
				reduce = reduce + *it;
			}
			return reduce;
		}
		else
		{
#ifdef EXECUTION_HEADER_AVAILABLE
			if (params::parallelization)
				return std::reduce(
					std::execution::par_unseq,
					static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
					static_cast<const E*>(&e)->end(symphas::it_grp, interval));
			else
#endif
				return std::reduce(
					static_cast<const E*>(&e)->begin(symphas::it_grp, interval),
					static_cast<const E*>(&e)->end(symphas::it_grp, interval));
		}
	}

	template<typename E, size_t D>
	auto result_sum(OpEvaluable<E> const& e, grid::region_interval_multiple<D> const& regions)
	{
		expr::eval_type_t<E> sum{};
		for (grid::region_interval<D> region : regions)
		{
			sum += result_sum(*static_cast<E const*>(&e), region);
		}
		return sum;
	}

	template<typename E>
	auto result_sum(OpEvaluable<E> const& e)
	{
		TIME_THIS_EXPRESSION_LIFETIME(iterable_domain, auto r = expr::iterable_domain(*static_cast<E const*>(&e));)
		return result_sum(*static_cast<E const*>(&e), expr::iterable_domain(*static_cast<E const*>(&e)));
	}




	struct matches_series;
	struct matches_mul;
	struct matches_div;
	struct matches_term;
	struct matches_integral;
	struct matches_derivative;
	struct matches_operator;

	template<typename matches_t>
	struct matching_in_mul;


	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename E, typename assign_type>
	void result_of_matching(OpEvaluable<E> const& e, assign_type&& data);

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename... Es, typename assign_type>
	void result_of_matching(OpAdd<Es...> const& e, assign_type&& data);

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename E, typename assign_type, typename region_type>
	void result_of_matching(OpEvaluable<E> const& e, assign_type&& data, region_type&& region);

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename... Es, typename assign_type, typename region_type>
	void result_of_matching(OpAdd<Es...> const& e, assign_type&& data, region_type&& region);

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename E, typename region_type>
	auto result_sum_of_matching(OpEvaluable<E> const& e, region_type&& region);

	//! Accumulates the result if the given expression matches the condition.
	template<typename condition_t, typename... condition_ts, typename... Es, typename region_type>
	auto result_sum_of_matching(OpAdd<Es...> const& e, region_type&& region);


	template<typename condition_t, typename... condition_ts, typename E, typename assign_type, typename E0 = OpVoid>
	void result_by_group(OpEvaluable<E> const& e, assign_type&& data, E0 const& init = OpVoid{});

	template<typename condition_t, typename... condition_ts, typename E, typename assign_type, typename region_type, typename E0 = OpVoid>
	void result_by_group(OpEvaluable<E> const& e, assign_type&& data, region_type&& region, E0 const& init = OpVoid{});


	template<typename condition_t>
	struct result_by_term_apply;

	template<typename condition_t, typename... condition_ts, typename... Es, typename assign_type>
	void result_by_term(OpAdd<Es...> const& e, assign_type&& data);
	template<typename condition_t, typename... condition_ts, typename E, typename assign_type>
	void result_by_term(OpEvaluable<E> const& e, assign_type&& data);

	template<typename condition_t>
	struct result_sum_by_term_apply;


	template<typename condition_t, typename... condition_ts, typename... Es>
	auto result_sum_by_term(OpAdd<Es...> const& e);
	template<typename condition_t, typename... condition_ts, typename E>
	auto result_sum_by_term(OpEvaluable<E> const& e);


}




namespace symphas::internal
{
	struct tensor_cast
	{
		template<typename T, size_t... Ns>
		static auto const& cast(OpTensor<T, Ns...> const& tensor)
		{
			return tensor.value;
		}

		template<typename T, size_t... Ns>
		static auto& cast(OpTensor<T, Ns...>& tensor)
		{
			return tensor.value;
		}
	};

	template<typename T, size_t N0, size_t N1, size_t NA, size_t NB, size_t... Ns, 
		typename S, size_t M0, size_t M1, size_t MA, size_t MB, size_t... Ms>
	auto tensor_multiply(OpTensor<T, N0, N1, NA, NB, Ns...>, OpTensor<S, M0, M1, MA, MB, Ms...>) = delete;

	// Matrix multiplication
	template<typename T, typename S, size_t N0, size_t N1, size_t N2, size_t N3, size_t NA, size_t NB, size_t NC>
	auto tensor_multiply(OpTensor<T, N0, N1, NA, NB> const& a, OpTensor<S, N2, N3, NB, NC> const& b)
	{
		if constexpr (N1 == N2)
		{
			T at = tensor_cast::cast(a);
			S bt = tensor_cast::cast(b);
			return expr::make_tensor<N0, N3, NA, NC>(at * bt);
		}
		else
		{
			return OpVoid{};
		}
	}

	// Matrix multiplication
	template<typename T, typename S, size_t N0, size_t N1, size_t NA, size_t NB, 
		typename std::enable_if_t<(NA > 0), int> = 0>
	auto tensor_multiply(OpTensor<T, N0, NA> const& a, OpTensor<S, 0, N1, 1, NB> const& b)
	{
		T at = tensor_cast::cast(a);
		S bt = tensor_cast::cast(b);
		return expr::make_tensor<N0, N1, NA, NB>(at * bt);
	}

	// Matrix multiplication
	template<typename T, typename S, size_t N0, size_t N1, size_t N2, size_t NA, size_t NB>
	auto tensor_multiply(OpTensor<T, N0, N1, NA, NB> const& a, OpTensor<S, N2, NB> const& b)
	{
		if constexpr (N1 == N2)
		{
			T at = tensor_cast::cast(a);
			S bt = tensor_cast::cast(b);
			
			if constexpr (NA == 1)
			{
				return expr::make_literal(at * bt);
			}
			else
			{
				return expr::make_tensor<N0, 0, NA, 1>(at * bt);
			}
		}
		else
		{
			return OpVoid{};
		}
	}

	template<typename T, size_t... Ms>
	auto tensor_multiply(OpTensor<OpIdentity, 0, 0> const& a, OpTensor<T, Ms...> const& b)
	{
		return OpIdentity{};
	}

	template<typename T, size_t... Ms>
	auto tensor_multiply(OpTensor<T, Ms...> const& a, OpTensor<OpIdentity, 0, 0> const& b)
	{
		return OpIdentity{};
	}

	template<typename T, typename S, size_t... Ms>
	auto tensor_multiply(OpTensor<T, Ms...> const& a, OpTensor<S, 0, 1> const& b)
	{
		auto v = symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
		return expr::make_tensor<Ms...>(v);
	}

	template<typename T, typename S, size_t... Ms>
	auto tensor_multiply(OpTensor<T, 0, 0, 1, 1> const& a, OpTensor<S, Ms...> const& b)
	{
		auto v = symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
		return expr::make_tensor<Ms...>(v);
	}

	template<typename T, typename S, size_t N0, size_t N>
	auto tensor_multiply(OpTensor<T, 0, 0, 1, 1> const& a, OpTensor<S, 0, N0, 1, N> const& b)
	{
		auto v = symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
		return expr::make_tensor<0, N0, 1, N>(v);
	}


	template<typename S>
	auto tensor_multiply(OpTensor<OpIdentity, 0, 0> const& a, OpTensor<S, 0, 1> const& b)
	{
		return OpIdentity{};
	}

	template<typename S>
	auto tensor_multiply(OpTensor<S, 0, 0, 1, 1> const& a, OpTensor<OpIdentity, 0, 0> const& b)
	{
		return OpIdentity{};
	}




	template<typename T, typename S>
	auto tensor_multiply(OpTensor<T, 0, 1> const& a, OpTensor<S, 0, 1> const& b)
	{
		return symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
	}

	template<typename T, typename S>
	auto tensor_multiply(OpTensor<T, 0, 0, 1, 1> const& a, OpTensor<S, 0, 1> const& b)
	{
		return symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
	}

	template<typename T, typename S>
	auto tensor_multiply(OpTensor<S, 0, 1> const& a, OpTensor<T, 0, 0, 1, 1> const& b)
	{
		auto v = symphas::internal::tensor_cast::cast(a) * symphas::internal::tensor_cast::cast(b);
		return expr::make_tensor<0, 0, 1, 1>(v);
	}

	template<typename T, typename S, size_t N, size_t M, size_t D>
	auto tensor_multiply(OpTensor<S, N, D> const& a, OpTensor<T, M, D> const& b)
	{
		return expr::make_row_vector<N, D>(S(a)) * b;
	}

	template<size_t N0, size_t N1, size_t NA, size_t NB, typename T>
	auto tensor_as_coeff(T const& value)
	{
		//if constexpr (N0 == 0 && N1 == 0 && NA == 1 && NB == 1)
		//{
		//	return value.eval();
		//}
		//else
		{
			using elem_type = std::invoke_result_t<decltype(&T::eval), T, iter_type>;
			any_matrix_t<elem_type, NA, NB> matrix;
			matrix[N0][N1] = value.eval();
			return matrix;
		}
	}

	template<size_t N0, size_t NA, typename T>
	auto tensor_as_coeff(T const& value)
	{
		if constexpr (N0 == 0 && NA == 0)
		{
			return value;
		}
		else
		{
			using elem_type = std::invoke_result_t<decltype(&T::eval), T, iter_type>;
			any_vector_t<elem_type, NA> vector;
			vector[N0] = value.eval();
			return vector;
		}
	}

	template<typename>
	struct tensor_to_vector;

	template<typename T, size_t N1, size_t NA>
	struct tensor_to_vector<OpTensor<T, N1, NA>>
	{
		using type = any_vector_t<std::invoke_result_t<decltype(&T::eval), T, iter_type>, NA>;
	};

	template<typename T, size_t N1, size_t N2, size_t NA, size_t NB>
	struct tensor_to_vector<OpTensor<T, N1, N2, NA, NB>>
	{
		using type = any_vector_t<any_vector_t<std::invoke_result_t<decltype(&T::eval), T, iter_type>, NB>, NA>;
	};


	template<typename T>
	using tensor_to_vector_t = typename tensor_to_vector<T>::type;

	template<size_t N1, size_t NA, typename T, typename S>
	void set_vector(any_vector_t<T, NA>& vector, S const& value)
	{
		vector[N1] = value;
	}

	template<size_t N1, size_t N2, size_t NA, size_t NB, typename T, typename S>
	void set_vector(any_vector_t<any_vector_t<T, NB>, NA>& vector, S const& value)
	{
		vector[N1][N2] = value;
	}
}



//! A value of a tensor, acts as a coefficient.
/*! 
 * A tensor turns the expression for which it is a coefficient into an entry of a
 * tensor. 
 * 
 * A tensor of rank 1 is a column vector, and of rank 2 is a matrix. A row vector is a
 * tensor of rank 2 and has only 1 row. Multiplication only for rank 1 and 2 vectors is
 * supported, and there are no tensor products.
 */
template<typename T, size_t... Ns>
struct OpTensor : OpExpression<OpTensor<T, Ns...>>
{

	explicit OpTensor(T const& entry) : value{ entry } {}
	explicit OpTensor(T&& entry) noexcept : value{ std::move(entry) } {}
	constexpr OpTensor() : value{ T{} } {}

	auto eval(iter_type n = 0) const
	{
		return symphas::internal::tensor_as_coeff<Ns...>(cast());
	}

	template<typename S, size_t... Ms>
	auto operator*(OpTensor<S, Ms...> const& a) const
	{
		return symphas::internal::tensor_multiply(*this, a);
	}

	template<typename S>
	auto operator-(OpTensor<S, Ns...> const& a) const
	{
		return expr::make_tensor<Ns...>(cast() - symphas::internal::tensor_cast::cast(a));
	}

	template<typename S>
	auto operator+(OpTensor<S, Ns...> const& a) const
	{
		return expr::make_tensor<Ns...>(cast() + symphas::internal::tensor_cast::cast(a));
	}

	auto operator-() const
	{
		return expr::make_tensor<Ns...>(-cast());
	}

	operator symphas::internal::tensor_to_vector_t<OpTensor<T, Ns...>>() const
	{
		symphas::internal::tensor_to_vector_t<OpTensor<T, Ns...>> vector;
		symphas::internal::set_vector<Ns...>(vector, symphas::internal::tensor_cast::cast(*this));
		return vector;
	}

	explicit operator T const& () const
    {
        return cast();
    }

	explicit operator T& ()
    {
        return cast();
    }

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::print_with_coeff(out, "", *this);
	}

	size_t print(char* out) const
	{
		return expr::print_with_coeff(out, "", *this);
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(*this);
	}

#endif

	friend struct symphas::internal::tensor_cast;

protected:

	const T& cast() const
	{
		return symphas::internal::tensor_cast::cast(*this);
	}

	T value;
};

namespace expr
{
	template<typename T, size_t... Ns>
	decltype(auto) coeff(OpTensor<T, Ns...> const& tensor)
	{
		return tensor;
		//return symphas::internal::tensor_cast::cast(tensor);
	}
}

namespace symphas::internal
{
	using tensor_cancel = OpTensor<OpIdentity, 0, 0>;
}


template<typename coeff_t, typename T, size_t... Ns,
	typename std::enable_if_t<(expr::is_identity<coeff_t> && expr::is_fraction<coeff_t>), int> = 0>
auto operator*(coeff_t const& value, OpTensor<T, Ns...> const& tensor)
{
	return expr::make_tensor<Ns...>(value * symphas::internal::tensor_cast::cast(tensor));
}

template<typename T1, typename T2, size_t... Ns>
auto operator*(OpLiteral<T1> const& value, OpTensor<T2, Ns...> const& tensor)
{
	return expr::make_tensor<Ns...>(value * symphas::internal::tensor_cast::cast(tensor));
}

template<typename T1, typename T2, size_t... Ns>
auto operator*(OpTensor<T1, Ns...> const& tensor, OpLiteral<T2> const& value)
{
	return expr::make_tensor<Ns...>(symphas::internal::tensor_cast::cast(tensor) * value);
}


template<typename T1, typename T2, size_t N0, size_t N1, size_t NB>
auto operator*(any_vector_t<T1, 1> const& value, OpTensor<T2, 0, N1, 1, NB> const& tensor)
{
	return expr::make_tensor<0, N1, 1, NB>(value[0] * symphas::internal::tensor_cast::cast(tensor));
}

template<typename T1, typename T2, size_t N0, size_t N1, size_t NB>
auto operator*(any_vector_t<T1, 2> const& value, OpTensor<T2, 0, N1, 1, NB> const& tensor)
{
	auto v = (*static_cast<T2 const*>(&tensor));
	return expr::make_tensor<0, N1, 2, NB>(value[0] * v)
		+ expr::make_tensor<1, N1, 2, NB>(value[1] * v);
}

template<typename T1, typename T2, size_t N0, size_t N1, size_t NB>
auto operator*(any_vector_t<T1, 3> const& value, OpTensor<T2, 0, N1, 1, NB> const& tensor)
{
	auto v = (*static_cast<T2 const*>(&tensor));
	return expr::make_tensor<0, N1, 3, NB>(value[0] * v)
		+ expr::make_tensor<1, N1, 3, NB>(value[1] * v)
		+ expr::make_tensor<2, N1, 3, NB>(value[2] * v);
}


template<typename T1, typename T2, size_t N0, size_t NA>
auto operator*(OpTensor<T1, N0, NA> const& tensor, any_vector_t<T2, 1> const& value)
{
	return expr::make_tensor<N0, NA>(value[0] * (*static_cast<T1 const*>(&tensor)));
}

template<typename T1, typename T2, size_t N0, size_t N1, size_t NA>
auto operator*(OpTensor<T1, N0, N1, NA, 2> const& tensor, any_vector_t<T2, 2> const& value)
{
	return expr::make_tensor<N0, NA>(value[N1] * (*static_cast<T1 const*>(&tensor)));
}

template<typename T1, typename T2, size_t N0, size_t N1, size_t NA>
auto operator*(OpTensor<T1, N0, N1, NA, 3> const& tensor, any_vector_t<T2, 3> const& value)
{
	return expr::make_tensor<N0, NA>(value[N1] * (*static_cast<T1 const*>(&tensor)));
}

template<typename T, size_t D>
auto operator*(OpTensor<OpIdentity, 0, 0> const&, OpLiteral<any_vector_t<T, D>> const&)
{
	return OpIdentity{};
}

template<typename T, size_t D>
auto operator*(OpTensor<OpIdentity, 0, 0> const&, any_vector_t<T, D> const&)
{
	return OpIdentity{};
}



//! Representation of a constant.
/*!
 * Stores a single value of any type. 
 * 
 * \tparam T The type of constant being stored.
 */
template<typename T>
struct OpLiteral : OpExpression<OpLiteral<T>>
{
	T value;

	OpLiteral(T value) : value{ value } {}
	constexpr OpLiteral(OpIdentity);
	constexpr OpLiteral(OpNegIdentity);
	constexpr OpLiteral(OpVoid);
	constexpr OpLiteral() : value{ T{} } {}

	inline T eval(iter_type = 0) const
	{
		return value;
	}

	operator const T&() const
	{
		return value;
	}

	operator T&()
	{
		return value;
	}


	//auto operator^(size_t exp) const
	//{
	//	return expr::make_literal(std::pow(value, exp));
	//}

	template<typename S>
	auto operator*(OpLiteral<S> const& other) const
	{
		return expr::make_literal(value * other.value);
	}

	auto operator*(OpLiteral<T> const& other) const
	{
		return expr::make_literal(value * other.value);
	}

	template<typename S>
	auto operator-(OpLiteral<S> const& a) const
	{
		auto v = expr::make_literal(value - a.value);
		return v;
	}

	template<typename S>
	auto operator*(OpLiteral<S>&& other) const
	{
		return expr::make_literal(value * other.value);
	}

	auto operator*(OpLiteral<T>&& other) const
	{
		return expr::make_literal(value * other.value);
	}

	template<typename S>
	auto operator-(OpLiteral<S>&& a) const
	{
		auto v = expr::make_literal(value - a.value);
		return v;
	}


	auto operator-() const
	{
		return expr::make_literal(-value);
	}

	template<typename S>
	auto operator+(OpLiteral<S> const& a) const
	{
		auto v = expr::make_literal(value + a.value);
		return v;
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const;
	size_t print(char* out) const;
	size_t print_length() const;

#endif

};


//! Representation of a constant.
/*!
 * Stores a single value of any type.
 *
 * \tparam T The type of constant being stored.
 */

template<typename T, size_t D>
struct OpLiteral<any_vector_t<T, D>> : OpExpression<OpLiteral<any_vector_t<T, D>>>
{
	any_vector_t<T, D> value;

	OpLiteral(any_vector_t<T, D> value) : value{ value } {}
	constexpr OpLiteral() : value{ any_vector_t<T, D>{} } {}

	inline T eval(iter_type = 0) const
	{
		return value;
	}

	operator const any_vector_t<T, D>& () const
	{
		return value;
	}

	operator any_vector_t<T, D>& ()
	{
		return value;
	}


	//auto operator^(size_t exp) const
	//{
	//	return expr::make_literal(std::pow(value, exp));
	//}

	template<typename S>
	auto operator*(OpLiteral<S> const& other) const
	{
		return expr::make_literal(value * other.value);
	}

	auto operator*(OpLiteral<T> const& other) const
	{
		return expr::make_literal(value * other.value);
	}

	template<typename S>
	auto operator-(OpLiteral<S> const& a) const
	{
		auto v = expr::make_literal(value - a.value);
		return v;
	}

	template<typename S>
	auto operator*(OpLiteral<S>&& other) const
	{
		return expr::make_literal(value * other.value);
	}

	auto operator*(OpLiteral<T>&& other) const
	{
		return expr::make_literal(value * other.value);
	}

	template<typename S>
	auto operator-(OpLiteral<S>&& a) const
	{
		auto v = expr::make_literal(value - a.value);
		return v;
	}


	auto operator-() const
	{
		return expr::make_literal(-value);
	}

	template<typename S>
	auto operator+(OpLiteral<S> const& a) const
	{
		auto v = expr::make_literal(value + a.value);
		return v;
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::print_with_coeff(out, value);
	}

	size_t print(char* out) const
	{
		return expr::print_with_coeff(out, value);
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value);
	}

#endif

};


#ifdef PRINTABLE_EQUATIONS

template<>
inline size_t OpLiteral<int>::print(FILE* out) const
{
	return fprintf(out, "%d", value);
}

template<>
inline size_t OpLiteral<double>::print(FILE* out) const
{
	if (int(value) == value)
	{
		return OpLiteral<int>(int(value)).print(out);
	}
	else
	{
		int p = symphas::lib::exponent10(value);
		if (std::abs(p) >= EXPR_VALUE_DISPLAY_PRECISION)
		{
			double b = symphas::lib::base10(value);
			return fprintf(out, SYEX_SCIENTIFIC_A "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf"
				SYEX_SCIENTIFIC_VALUE_A "%d" SYEX_SCIENTIFIC_VALUE_B SYEX_SCIENTIFIC_B, b, p);
		}
		else
		{
			return fprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf", value);
		}
	}
}

template<>
inline size_t OpLiteral<complex_t>::print(FILE* out) const
{
	return fprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf + %." STR(EXPR_VALUE_DISPLAY_PRECISION) "lfi", real(value), imag(value));
}

template<>
inline size_t OpLiteral<int>::print(char* out) const
{
	return sprintf(out, "%d", value);
}

template<>
inline size_t OpLiteral<double>::print(char* out) const
{
	if (int(value) == value)
	{
		return OpLiteral<int>(int(value)).print(out);
	}
	else
	{
		int p = symphas::lib::exponent10(value);
		if (std::abs(p) >= EXPR_VALUE_DISPLAY_PRECISION)
		{
			double b = symphas::lib::base10(value);
			return sprintf(out, SYEX_SCIENTIFIC_A "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf"
				SYEX_SCIENTIFIC_VALUE_A "%d" SYEX_SCIENTIFIC_VALUE_B SYEX_SCIENTIFIC_B, b, p);
		}
		else
		{
			return sprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf", value);
		}
	}
}

template<>
inline size_t OpLiteral<complex_t>::print(char* out) const
{
	if (real(value) == 0)
	{
		if (imag(value) < 0)
		{
			return sprintf(out, "-%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lfi", std::abs(imag(value)));
		}
		else if (imag(value) == 0)
		{
			return sprintf(out, "0");
		}
		else
		{
			return sprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lfi", imag(value));
		}
	}
	else
	{
		if (imag(value) < 0)
		{
			return sprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf - %." STR(EXPR_VALUE_DISPLAY_PRECISION) "lfi", real(value), std::abs(imag(value)));
		}
		else if (imag(value) == 0)
		{
			return sprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf", real(value)); //1
		}
		else
		{
			return sprintf(out, "%." STR(EXPR_VALUE_DISPLAY_PRECISION) "lf + %." STR(EXPR_VALUE_DISPLAY_PRECISION) "lfi", real(value), imag(value));
		}
	}
}

template<>
inline size_t OpLiteral<int>::print_length() const
{
	size_t len = symphas::lib::num_digits(std::abs(value));
	if (value < 0)
	{
		return len + 1;
	}
	else
	{
		return len;
	}
}

template<>
inline size_t OpLiteral<double>::print_length() const
{
	if (int(value) == value)
	{
		return OpLiteral<int>(int(value)).print_length();
	}
	else
	{
		int p = symphas::lib::exponent10(value);
		if (std::abs(p) >= EXPR_VALUE_DISPLAY_PRECISION)
		{
			size_t len = symphas::lib::num_digits(p);
			return 2 + EXPR_VALUE_DISPLAY_PRECISION + len 
				+ STR_ARR_LEN(SYEX_SCIENTIFIC_VALUE_A SYEX_SCIENTIFIC_VALUE_B)
				+ STR_ARR_LEN(SYEX_SCIENTIFIC_A SYEX_SCIENTIFIC_B)
				+ ((value < 0) ? 1 : 0) + ((p < 0) ? 1 : 0);
		}
		else
		{
			size_t len = symphas::lib::num_digits(static_cast<int>(std::abs(value)));
			return len + 1 + EXPR_VALUE_DISPLAY_PRECISION + ((value < 0) ? 1 : 0);
		}
	}
}

template<>
inline size_t OpLiteral<complex_t>::print_length() const
{
	size_t len = 0;
	len += symphas::lib::num_digits(static_cast<int>(std::abs(real(value))));
	len += symphas::lib::num_digits(static_cast<int>(std::abs(imag(value))));

	if (real(value) == 0)
	{
		if (imag(value) < 0)
		{
			len += 3 + EXPR_VALUE_DISPLAY_PRECISION;
		}
		else if (imag(value) == 0)
		{
			len += 1;
		}
		else
		{
			len += 2 + EXPR_VALUE_DISPLAY_PRECISION;
		}
	}
	else
	{
		if (imag(value) == 0)
		{
			len += 1 + EXPR_VALUE_DISPLAY_PRECISION;
		}
		else
		{
			len += 6 + EXPR_VALUE_DISPLAY_PRECISION * 2;
		}
	}


	if (real(value) < 0)
	{
		return len + 1;
	}
	else
	{
		return len;
	}
}

#endif


template<size_t N, size_t D>
struct OpNegFractionLiteral;
template<size_t N, size_t D>
struct OpFractionLiteral;

template<size_t N, size_t D>
struct OpFractionLiteral : OpExpression<OpFractionLiteral<N, D>>
{
	constexpr inline scalar_t eval(iter_type = 0) const
	{
		return static_cast<scalar_t>(N) / D;
	}

	constexpr operator const scalar_t() const
	{
		return eval();
	}

	constexpr auto operator-() const
	{
		return OpNegFractionLiteral<N, D>{};
	}
	
	constexpr auto operator--(int) const
	{
		return expr::numeric_range(OpFractionLiteral<N, D>{});
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return fprintf(out, SYEX_FRA_FMT, N, D);
	}

	size_t print(char* out) const
	{
		return sprintf(out, SYEX_FRA_FMT, N, D);
	}

	size_t print_length() const
	{
		size_t len1 = symphas::lib::num_digits<N>();
		size_t len2 = symphas::lib::num_digits<N>();
		return len1 + len2 + SYEX_FRA_FMT_LEN;
	}

#endif

};

template<size_t N>
struct OpFractionLiteral<N, 1> : OpExpression<OpFractionLiteral<N, 1>>
{
	constexpr inline scalar_t eval(iter_type = 0) const
	{
		return static_cast<scalar_t>(N);
	}

	constexpr operator const scalar_t() const
	{
		return eval();
	}

	constexpr auto operator-() const
	{
		return OpNegFractionLiteral<N, 1>{};
	}

	constexpr auto operator--(int) const
	{
		return expr::numeric_range(OpFractionLiteral<N, 1>{});
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return fprintf(out, "%zd", N);
	}

	size_t print(char* out) const
	{
		return sprintf(out, "%zd", N);
	}

	size_t print_length() const
	{
		return symphas::lib::num_digits<N>();
	}

#endif

};


template<size_t N, size_t D>
struct OpNegFractionLiteral : OpExpression<OpNegFractionLiteral<N, D>>
{
	constexpr inline scalar_t eval(iter_type = 0) const
	{
		return -static_cast<scalar_t>(N) / D;
	}

	constexpr operator const scalar_t() const
	{
		return eval();
	}

	constexpr auto operator-() const
	{
		return OpFractionLiteral<N, D>{};
	}

	constexpr auto operator--(int) const
	{
		return expr::numeric_range(OpNegFractionLiteral<N, D>{});
	}


#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = 0;
		n += fprintf(out, "-");
		n += OpFractionLiteral<N, D>{}.print(out);
		return n;
	}

	size_t print(char* out) const
	{
		size_t n = 0;
		n += sprintf(out, "-");
		n += OpFractionLiteral<N, D>{}.print(out + n);
		return n;
	}

	size_t print_length() const
	{
		return OpFractionLiteral<N, D>{}.print_length() + 1;
	}

#endif

};


template<typename T, size_t N, size_t D>
auto operator*(OpLiteral<T> const& a, OpFractionLiteral<N, D>)
{
	return expr::make_literal(a.value * OpFractionLiteral<N, D>{}.eval());
}

template<typename T, size_t N, size_t D>
auto operator*(OpLiteral<T> const& a, OpNegFractionLiteral<N, D>)
{
	return expr::make_literal(a.value * OpNegFractionLiteral<N, D>{}.eval());
}

// ******************************************************************************************


//! Specialized version which returns multiplicative identity.
template<>
inline decltype(auto) expr::make_literal<OpIdentity>(OpIdentity const&)
{
	return OpIdentity{};
}

//! Specialized version which returns negative multiplicative identity.
template<>
inline decltype(auto) expr::make_literal<OpNegIdentity>(OpNegIdentity const&)
{
	return OpNegIdentity{};
}

//! Specialized version which returns additive identity.
template<>
inline decltype(auto) expr::make_literal<OpVoid>(OpVoid const&)
{
	return OpVoid{};
}


template<typename T>
constexpr OpLiteral<T>::OpLiteral(OpIdentity) : value{ symphas::lib::get_identity<T>() } {}
template<typename T>
constexpr OpLiteral<T>::OpLiteral(OpNegIdentity) : value{ -symphas::lib::get_identity<T>() } {}
template<typename T>
constexpr OpLiteral<T>::OpLiteral(OpVoid) : value{ 0 } {}


namespace expr
{
	namespace symbols
	{
		//! Common fractions.
		constexpr OpFractionLiteral<1, 2> _2{};
		constexpr OpFractionLiteral<1, 4> _4{};

		//! Pi to 7th order
		constexpr OpFractionLiteral<355, 113> Pi{};
	}
}

// ******************************************************************************************

namespace symphas::internal
{
	using namespace symphas::lib;


	template<int... Is>
	constexpr static bool all_positive(std::integer_sequence<int, Is...>)
	{
		return ((Is > 0) && ...) || (sizeof...(Is) == 1 && ((Is == 0) && ...));
	}

	template<size_t... Ns>
	struct construct_tensor
	{

		static constexpr size_t N = sizeof...(Ns);
		using pos_t = seq_join_t<types_before_index<N / 2, std::index_sequence<Ns>...>>;
		using dim_t = seq_join_t<types_after_at_index<N / 2, std::index_sequence<Ns>...>>;

		template<typename T>
		auto operator()(T const& v)
		{
			static_assert(symphas::internal::all_positive(seq_sub_t<dim_t, pos_t>{}),
				"incorrect tensor arguments");
			return OpTensor<T, Ns...>{ v };
		}
	};

	template<size_t N0, size_t NA>
	struct construct_tensor<N0, 0, NA, 1> : construct_tensor<N0, NA> {};


	template<>
	struct construct_tensor<0, 0, 1, 1>
	{
		template<typename T>
		auto operator()(T const& v)
		{
			return OpTensor<T, 0, 0, 1, 1>{ v };
		}
	};
}


template<size_t N0, size_t N1, size_t... Ns, typename T>
constexpr auto expr::make_tensor(T const& v)
{
	return symphas::internal::construct_tensor<N0, N1, Ns...>{}(expr::make_literal(v));
}

template<size_t N0, size_t N1, size_t... Ns>
constexpr auto expr::make_tensor()
{
	return expr::make_tensor<N0, N1, Ns...>(OpIdentity{});
}

template<size_t I, size_t N, typename T>
constexpr auto expr::make_column_vector(T&& v)
{
	return make_tensor<I, N>(std::forward<T>(v));
}

template<size_t I, size_t N>
constexpr auto expr::make_column_vector()
{
	return make_column_vector<I, N>(OpIdentity{});
}

template<size_t I, size_t N, typename T>
constexpr auto expr::make_row_vector(T&& v)
{
	return make_tensor<0, I, 1, N>(std::forward<T>(v));
}

template<size_t I, size_t N>
constexpr auto expr::make_row_vector()
{
	return make_row_vector<I, N>(OpIdentity{});
}


template<size_t R, size_t R0>
auto expr::make_filled_column_vector()
{
	if constexpr (R0 >= R)
	{
		return OpIdentity{};
	}
	else
	{
		auto c = expr::make_column_vector<R0, R>();
		if constexpr (R0 == R - 1)
		{
			return c;
		}
		else
		{
			return c + make_filled_column_vector<R, R0 + 1>();
		}
	}
}

template<size_t R, size_t R0>
auto expr::make_filled_row_vector()
{
	if constexpr (R0 >= R)
	{
		return OpIdentity{};
	}
	else
	{
		auto c = expr::make_row_vector<R0, R>();
		if constexpr (R0 == R - 1)
		{
			return c;
		}
		else
		{
			return c + make_filled_row_vector<R, R0 + 1>();
		}
	}
}

template<typename T>
constexpr decltype(auto) expr::make_literal(T const& v)
{
	if constexpr (expr::is_fraction<T>)
	{
		return T{};
	}
	else if constexpr (expr::is_tensor<T>)
	{
		return v;
	}
	else if constexpr (expr::is_arr_coeff<T>)
	{
		return v;
	}
	else
	{
		return OpLiteral{ v };
	}
}

template<size_t N, size_t D>
constexpr auto expr::make_fraction()
{
	static_assert(D != 0, "dividing by zero");
	if constexpr (N == 0)
	{
		return OpVoid{};
	}
	else if constexpr (N == D)
	{
		return OpIdentity{};
	}
	else
	{
		constexpr size_t GCD = GCD_of<N, D>;
		return OpFractionLiteral<N / GCD, D / GCD>{};
	}
}

template<int I>
constexpr auto expr::make_integer()
{
	if constexpr (I < 0)
	{
		return -expr::make_fraction<static_cast<size_t>(-I), 1>();
	}
	else
	{
		return expr::make_fraction<static_cast<size_t>(I), 1>();
	}
}



// ******************************************************************************************


struct DynamicIndex : OpExpression<DynamicIndex>
{
	DynamicIndex(iter_type&& data, len_type start_index, len_type end_index) :
		data{ &data0 }, data0{ data }, start_index{ start_index }, end_index{ end_index } {}
	DynamicIndex(iter_type&& data, len_type start_index) : 
		data{ &data0 }, data0{ data }, start_index{ start_index }, end_index{ data } {}
	DynamicIndex(iter_type&& data) : DynamicIndex(std::move(data), len_type(data)) {}

	DynamicIndex(iter_type& data, len_type start_index, len_type end_index) :
		data{ &data }, data0{ data }, start_index{ start_index }, end_index{ end_index } {}
	DynamicIndex(iter_type& data) : DynamicIndex(data, data, data) {}

	DynamicIndex(iter_type* data, len_type start_index, len_type end_index) : DynamicIndex(*data, start_index, end_index) {}
	DynamicIndex(iter_type* data) : DynamicIndex(*data) {}

	DynamicIndex() : DynamicIndex(0) {}

	DynamicIndex(DynamicIndex const& other) :
		data{ (other.is_local()) ? &data0 : other.data },
		data0{ other.data0 },
		start_index{ other.start_index }, end_index{ other.end_index } {}

	DynamicIndex(DynamicIndex&& other) : DynamicIndex()
	{
		swap(*this, other);
	}

	DynamicIndex& operator=(DynamicIndex other)
	{
		swap(*this, other);
		return *this;
	}

	DynamicIndex operator=(iter_type index)
	{
		DynamicIndex other(*this);
		*other.data = index;
		return other;
	}

	friend void swap(DynamicIndex& first, DynamicIndex& second)
	{
		using std::swap;
		bool first_local = first.is_local();
		bool second_local = second.is_local();
		
		iter_type* second_data_ptr = second.data;
		second.data = (first_local) ? &second.data0 : first.data;
		first.data = (second_local) ? &first.data0 : second_data_ptr;

		swap(first.data0, second.data0);
		swap(first.start_index, second.start_index);
		swap(first.end_index, second.end_index);
	}

	auto eval(iter_type = 0) const
	{
		return index();
	}

	auto operator-() const;

#ifdef PRINTABLE_EQUATIONS
	size_t print(FILE* out) const
	{
		if (is_constant())
		{
			return expr::make_literal(index()).print(out);
		}
		else
		{
			return fprintf(out, "%s", expr::get_op_name(*this));
		}
	}

	size_t print(char* out) const
	{
		if (is_constant())
		{
			return expr::make_literal(index()).print(out);
		}
		else
		{
			return sprintf(out, "%s", expr::get_op_name(*this));
		}
	}

	size_t print_length() const
	{
		if (is_constant())
		{
			return expr::make_literal(index()).print_length();
		}
		else
		{
			return std::strlen(expr::get_op_name(*this));
		}
	}
#endif

	iter_type index() const
	{
		return *data;
	}

	len_type start() const
	{
		return start_index;
	}

	len_type end() const
	{
		return end_index;
	}

	len_type range() const
	{
		return end() - start() + 1;
	}

	explicit operator iter_type() const
	{
		return index();
	}

	bool is_constant() const
	{
		return (*data == end_index && (start_index == end_index || start_index == 0));
	}

	bool is_local() const
	{
		return (data == &data0);
	}

	void set_data(DynamicIndex const& other)
	{
		data = other.data;
		data0 = other.data0;
	}

	void fix(iter_type index)
	{
		DynamicIndex fixed(std::move(index));
		swap(*this, fixed);
	}

	void fix(DynamicIndex fixed)
	{
		swap(*this, fixed);
	}

	void fix(DynamicIndexSet const& fixed);

	bool same(DynamicIndex const& other)
	{
		return data == other.data;
	}

	bool same(iter_type *other)
	{
		return data == other;
	}

	friend struct DynamicIndexSet;
	friend auto operator==(DynamicIndex const& index, DynamicIndex& other);

protected:

	iter_type* data;	            //!< The reference to the index.
	iter_type data0;	            //!< The value of the index when it is local.
	len_type start_index;			//!< The smallest value the index can attain.
	len_type end_index;				//!< The maximum value the index can attain.
};

struct DynamicIndexSet
{

	DynamicIndexSet() : index{ 0 }, value{ 0 } {}
	DynamicIndexSet(DynamicIndex const& index, iter_type value) : index{ index.data }, value{ std::move(value) } {}
	DynamicIndexSet(DynamicIndex const& index, iter_type *value) : index{ index.data }, value{ *value } {}

	operator bool()
	{
		return *index == value.index();
	}

	iter_type *index;
	DynamicIndex value;
};

inline auto operator==(DynamicIndex const& index, iter_type value)
{
	return DynamicIndexSet(index, value);
}

inline auto operator==(DynamicIndex const& index, iter_type *value)
{
	return DynamicIndexSet(index, value);
}

inline auto operator==(DynamicIndex const& index, DynamicIndex& other)
{
	return (other.is_local()) ? (index == other.data0) : (index == other.data);
}

inline void DynamicIndex::fix(DynamicIndexSet const& fixed)
{
	if (same(fixed.index))
	{
		fix(fixed.value);
	}
}

namespace expr
{
	namespace
	{
		auto get_limit_data(DynamicIndex const& index)
		{
			if (index.is_constant())
			{
				return limit_data(index.index(), index.index());
			}
			else
			{
				return limit_data(0, 0, true);
			}
		}

		template<typename V>
		auto get_limit_data(OpBinaryMul<V, DynamicIndex> const& mul)
		{
			auto index = mul.b;
			if (index.is_constant())
			{
				auto c = mul.a.eval();
				return limit_data<OpVoid>(c * index.index(), c * index.index());
			}
			else
			{
				return limit_data<OpVoid>(0, 0, true);
			}
		}
	}
}

namespace symphas::internal
{
	template<typename T>
	struct coeff_data
	{
		coeff_data(const T* data, len_type len, len_type stride = 1) :
			data{ (len > 0) ? new T[len]{} : nullptr }, len{ len }
		{
			for (iter_type i = 0; i < len; ++i)
			{
				this->data[i] = data[i * stride];
			}
		}

		coeff_data(len_type len) : data{ (len > 0) ? new T[len]{} : nullptr }, len{ len } {}

		constexpr coeff_data() : coeff_data(nullptr, 0) {}

		coeff_data(coeff_data<T> const& other) : coeff_data(other.data, other.len) {}
		coeff_data(coeff_data<T>&& other) : coeff_data()
		{
			swap(*this, other);
		}

		auto operator=(coeff_data<T> other)
		{
			swap(*this, other);
		}

		friend void swap(coeff_data<T>& first, coeff_data<T>& second)
		{
			using std::swap;
			swap(first.data, second.data);
			swap(first.len, second.len);
		}

		const auto& operator[](iter_type i) const
		{
			return data[i];
		}

		auto& operator[](iter_type i)
		{
			return data[i];
		}

#ifdef PRINTABLE_EQUATIONS
		size_t print(FILE* out) const
		{
			return expr::make_literal(data[0]).print(out);
		}

		size_t print(char* out) const
		{
			return expr::make_literal(data[0]).print(out);
		}

		size_t print_length() const
		{
			return expr::make_literal(data[0]).print_length();
		}
#endif

		T* data;
		len_type len;

		~coeff_data()
		{
			delete[] data;
		}
	};
}


template<typename T, typename I>
struct OpCoeff;

template<typename T>
struct OpCoeff<T, void> : OpExpression<OpCoeff<T, void>>
{
	OpCoeff(const T* data, len_type len, len_type stride = 1) : data{ data, len, stride } {}
	OpCoeff(len_type len) : data{ len } {}
	OpCoeff(symphas::internal::coeff_data<T> const& data) : data{ data } {}
	OpCoeff(symphas::internal::coeff_data<T>&& data) : data{ std::move(data) } {}

	constexpr OpCoeff() : OpCoeff(nullptr, 0) {}

	auto eval(iter_type n = 0) const
	{
		if (data.len > 0)
		{
			return data[0];
		}
		else
		{
			return T{};
		}
	}

	auto operator-() const
	{
		OpCoeff neg(*this);
		for (iter_type i = 0; i < data.len; ++i)
		{
			neg.data[i] = -neg.data[i];
		}
		return neg;
	}

	template<int N, int P>
	auto operator()(expr::symbols::i_<N, P>) const
	{
		return OpCoeff<T, expr::symbols::i_<N, P>>(data);
	}

	auto operator()(DynamicIndex const& index) const
	{
		return OpCoeff<T, DynamicIndex>(index, data);
	}

	template<size_t N0>
	auto operator()(expr::symbols::placeholder_N_symbol_<N0>) const
	{
		return OpCoeff<T, expr::symbols::placeholder_N_symbol_<N0>>(data);
	}

	auto operator[](iter_type i) const
	{
		return expr::make_literal(data[i]);
	}

	template<int N, int P>
	auto operator[](expr::symbols::i_<N, P>) const
	{
		return this->operator()(expr::symbols::i_<N, P>{});
	}

	auto operator[](DynamicIndex const& index) const
	{
		return this->operator()(index);
	}

	template<size_t N0>
	auto operator[](expr::symbols::placeholder_N_symbol_<N0>) const
	{
		return this->operator()(expr::symbols::placeholder_N_symbol_<N0>{});
	}


#ifdef PRINTABLE_EQUATIONS
	size_t print(FILE* out) const
	{
		return data.print(out);
	}

	size_t print(char* out) const
	{
		return data.print(out);
	}

	size_t print_length() const
	{
		return data.print_length();
	}
#endif

	symphas::internal::coeff_data<T> data;
};

template<typename T, typename I>
struct OpCoeff : OpExpression<OpCoeff<T, I>>
{
	OpCoeff(const T* data, len_type len, len_type stride = 1) : data{ data, len, stride } {}
	OpCoeff(len_type len) : data{ len } {}
	OpCoeff(symphas::internal::coeff_data<T> const& data) : data{ data } {}
	OpCoeff(symphas::internal::coeff_data<T>&& data) : data{ std::move(data) } {}

	constexpr OpCoeff() : OpCoeff(nullptr, 0) {}

	auto eval(iter_type n = 0) const
	{
		if (data.len > 0)
		{
			return data[0];
		}
		else
		{
			return T{};
		}
	}

	auto operator-() const
	{
		OpCoeff neg(*this);
		for (iter_type i = 0; i < data.len; ++i)
		{
			neg.data[i] = -neg.data[i];
		}
		return neg;
	}

	auto operator()(DynamicIndex const& index) const
	{
		return OpCoeff<T, DynamicIndex>(index, data);
	}

	template<size_t N0>
	auto operator()(expr::symbols::placeholder_N_symbol_<N0>) const
	{
		return OpCoeff<T, expr::symbols::placeholder_N_symbol_<N0>>(data);
	}

	template<typename E>
	auto operator()(OpExpression<E> const& e) const
	{
		return *this * (*static_cast<E const*>(&e));
	}

	template<typename E>
	auto operator()(OpOperator<E> const& e) const
	{
		return *this * (*static_cast<E const*>(&e));
	}

	auto operator[](iter_type i) const
	{
		return expr::make_literal(data[i]);
	}

#ifdef PRINTABLE_EQUATIONS
	size_t print(FILE* out) const
	{
		return data.print(out);
	}

	size_t print(char* out) const
	{
		return data.print(out);
	}

	size_t print_length() const
	{
		return data.print_length();
	}
#endif

	symphas::internal::coeff_data<T> data;
};


template<typename T>
struct OpCoeff<T, DynamicIndex> : 
	OpExpression<OpCoeff<T, DynamicIndex>>
{
	OpCoeff(DynamicIndex const& index, len_type len) : data{ len }, index{ index } {}
	OpCoeff(DynamicIndex const& index, symphas::internal::coeff_data<T> const& data) : data{ data }, index{index} {}
	OpCoeff(DynamicIndex const& index, symphas::internal::coeff_data<T>&& data) : data{ std::move(data) }, index{ index } {}

	constexpr OpCoeff() : OpCoeff({}, 0) {}

	auto eval(iter_type n = 0) const
	{
		if ((int)index < data.len)
		{
			return data[(int)index];
		}
		else
		{
			return T{};
		}
	}

	auto operator-() const
	{
		OpCoeff neg(index, data);
		for (iter_type i = 0; i < data.len; ++i)
		{
			neg.data[i] = -neg.data[i];
		}
		return neg;
	}

	auto operator[](iter_type i) const
	{
		return expr::make_literal(data[i]);
	}

	auto& fix()
	{
		index = DynamicIndex(index.index());
		return *this;
	}

	auto& fix(DynamicIndexSet const& set)
	{
		index.fix(set);
		return *this;
	}

#ifdef PRINTABLE_EQUATIONS
	size_t print(FILE* out) const
	{
		return data.print(out);
	}

	size_t print(char* out) const
	{
		return data.print(out);
	}

	size_t print_length() const
	{
		return data.print_length();
	}
#endif

	symphas::internal::coeff_data<T> data;
	DynamicIndex index;
};


template<typename I>
struct OpCoeff<void, I> {};

template<>
struct OpCoeff<void, DynamicIndex> {};


//template<>
//struct OpCoeff<void, void>
//{
//	OpCoeff(iter_type i) : i{ i } {}
//	int i;
//};
//
//template<typename T>
//OpCoeff(T*, len_type) -> OpCoeff<T, void>;
//OpCoeff(iter_type)->OpCoeff<void, void>;

namespace expr
{
	template<typename T>
	auto make_coeff(T* data, len_type len, len_type stride)
	{
		return OpCoeff<T, void>(data, len, stride);
	}

	template<typename T, typename I, typename T0, typename I0>
	auto init_coeff_from(OpCoeff<T0, I0> const& coeff)
	{
		if constexpr (std::is_same_v<I, DynamicIndex>)
		{
			if constexpr (std::is_same_v<I0, DynamicIndex>)
			{
				return OpCoeff<T, I>(coeff.index, coeff.data.len);
			}
			else
			{
				// nothing here since cannot initialize a
				// coefficient with a dynamic index without an
				// existing dynamic index
			}
		}
		else
		{
			return OpCoeff<T, I>(coeff.data.len);
		}
	}


	//! Specialization based on SymbolID.
	template<>
	struct SymbolID<DynamicIndex>
	{
		static const char* get(DynamicIndex const& data)
		{
			return SymbolID<expr::symbols::placeholder_N_symbol>::get(expr::symbols::placeholder_N_symbol{});
		}
	};
}


//
//template<typename E1, typename E2,
//	typename std::enable_if_t<(expr::is_identity<E1> && expr::is_identity<E2>), int> = 0>
//constexpr bool operator==(E1 const& a, E2 const& b)
//{
//	return a.eval() == b.eval();
//}
//
//template<typename A, typename B,
//	typename std::enable_if_t<(expr::is_identity<A> && expr::is_identity<B>), int> = 0>
//constexpr bool operator<(A const& a, B const& b)
//{
//	return a.eval() < b.eval();
//}
//
//template<typename E1, typename E2, 
//	typename std::enable_if_t<expr::identity_comparator_enabled<E1, E2>, int> = 0>
//inline bool operator==(E1 const& a, E2 const& b)
//{
//	return OpLiteral<identity_eval_t>(a).eval() == OpLiteral<identity_eval_t>(b).eval();
//}
//
//template<typename E1, typename E2,
//	typename std::enable_if_t<expr::identity_comparator_enabled<E1, E2>, int> = 0>
//inline bool operator<(E1 const& a, E2 const& b)
//{
//	return OpLiteral<identity_eval_t>(a).eval() < OpLiteral<identity_eval_t>(b).eval();
//}
//
//
//template<typename E1, typename E2,
//	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
//inline bool operator!=(E1 const& a, E2 const& b)
//{
//	return !(a == b);
//}
//
//template<typename E1, typename E2,
//	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
//inline bool operator<=(E1 const& a, E2 const& b)
//{
//	return (a < b || a == b);
//}
//
//template<typename E1, typename E2,
//	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
//inline bool operator>(E1 const& a, E2 const& b)
//{
//	return !(a < b || a == b);
//}
//
//template<typename E1, typename E2,
//	typename std::enable_if_t<(expr::identity_comparator_enabled<E1, E2> || (expr::is_identity<E1> && expr::is_identity<E2>)), int> = 0>
//inline bool operator>=(E1 const& a, E2 const& b)
//{
//	return !(a < b);
//}

namespace expr
{
	//! Constructs the addition expression.
	/*!
	 * Directly constructs the addition expression between two
	 * expressions without applying any rules.
	 *
	 * \param a The first term in the addition.
	 * \param b The second term in the addition.
	 * \param es The rest of the terms in the addition.
	 */
	OpVoid make_add();
	OpVoid make_add(OpVoid, OpVoid);
	template<typename E>
	auto make_add(E&& a);
	template<typename E>
	auto make_add(E&& e, OpVoid);
	template<typename E>
	auto make_add(OpVoid, E&& e);
	template<typename... E0s, typename E0>
	auto make_add(OpAdd<E0s...> const& add, E0 const& e0);
	template<typename... E0s>
	auto make_add(OpAdd<E0s...> const& add, OpVoid);
	template<typename... E0s, typename E0>
	auto make_add(E0 const& e0, OpAdd<E0s...> const& add);
	template<typename... E0s>
	auto make_add(OpVoid, OpAdd<E0s...> const& add);
	template<typename A, typename B>
	auto make_add(OpExpression<A> const& a, OpExpression<B> const& b);
	template<typename A, typename B>
	auto make_add(OpOperator<A> const& a, OpExpression<B> const& b);
	template<typename A, typename B>
	auto make_add(OpExpression<A> const& a, OpOperator<B> const& b);
	template<typename E1, typename E2, typename E3, typename... Es>
	auto make_add(E1&& a, E2&& b, E3&& c, Es&& ...es);


	template<typename... As, typename... Bs, size_t... Is>
	auto make_add(OpAdd<As...> const& a, OpAdd<Bs...> const& b, std::index_sequence<Is...>)
	{
		return make_add(a, expr::get<Is>(b)...);
	}

	template<typename... As, typename... Bs>
	auto make_add(OpAdd<As...> const& a, OpAdd<Bs...> const& b)
	{
		if constexpr (sizeof...(As) > sizeof...(Bs))
		{
			return make_add(b, a, std::make_index_sequence<sizeof...(As)>{});
		}
		else
		{
			return make_add(a, b, std::make_index_sequence<sizeof...(Bs)>{});
		}
	}
}


template<typename T, size_t D, typename S>
bool operator>(any_vector_t<T, D> const& value, S&&)
{
	return true;
}

template<typename T, size_t D, typename S>
bool operator<(any_vector_t<T, D> const& value, S&&)
{
	return false;
}

template<typename... Es>
struct OpAddList;


template<>
struct OpAddList<>
{
	inline auto _eval(iter_type) const
	{
		return OpVoid{};
	}

	bool coeff_sign() const
	{
		return false;
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print_length() const
	{
		return 0;
	}

#endif

};

namespace symphas::internal
{
	using expr::has_coeff;

#ifdef PRINTABLE_EQUATIONS

	template<typename E>
	void print_one(E const& e, FILE* out, size_t& n)
	{
		if (expr::eval(expr::coeff(e)) > 0)
		{
			n += print_sep(out, SYEX_ADD_SEP);
			n += e.print(out);
		}
		else
		{
			n += print_sep(out, SYEX_SUB_SEP);
			n += (-e).print(out);
		}
	}

	template<typename E>
	void print_one(E const& e, char* out, size_t& n)
	{
		if (expr::eval(expr::coeff(e)) > 0)
		{
			n += print_sep(out + n, SYEX_ADD_SEP);
			n += e.print(out + n);
		}
		else
		{
			n += print_sep(out + n, SYEX_SUB_SEP);
			n += (-e).print(out + n);
		}
	}

	//! Specialization based on symphas::internal::binary_print().
	template<typename... Es, typename os_type, size_t... Is>
	size_t binary_print(OpAddList<Es...> const& e, os_type* out, std::index_sequence<Is...>)
	{
		size_t n = 0;
		(print_one(expr::get<Is + 1>(e), out, n), ...);
		return n;
	}


#endif

	template<size_t N, typename E>
	struct Nth_type_of_add;

	template<size_t N, typename E0, typename... Es>
	struct Nth_type_of_add<N, OpAdd<E0, Es...>>
	{
		using type = typename Nth_type_of_add<N - 1, OpAdd<Es...>>::type;
	};

	template<size_t N>
	struct Nth_type_of_add<N, OpAdd<>>
	{
		using type = OpAddList<>;
	};

	template<typename E0, typename... Es>
	struct Nth_type_of_add<0, OpAdd<E0, Es...>>
	{
		using type = OpAddList<E0, Es...>;
	};

	template<size_t N, typename E>
	using Nth_type_of_add_t = typename Nth_type_of_add<N, E>::type;

	template<typename cast_type>
	struct cast_add;

	template<typename... E0s>
	struct cast_add<OpAdd<E0s...>>
	{
        using cast_type = OpAdd<E0s...>;

		template<typename... Es>
		static cast_type cast(OpAdd<Es...> const& adds)
		{
			return *static_cast<OpAddList<E0s...> const*>(&adds);
		}

		template<typename... Es>
		static cast_type cast(OpAdd<Es...>& adds)
		{
			return *static_cast<OpAddList<E0s...>*>(&adds);
		}
	};

	template<typename... E0s>
	struct cast_add<OpAddList<E0s...>>
	{
        using cast_type = OpAddList<E0s...>;

		template<typename... Es>
		static const cast_type& cast(OpAdd<Es...> const& adds)
		{
			return *static_cast<cast_type const*>(&adds);
		}

		template<typename... Es>
		static cast_type& cast(OpAdd<Es...>& adds)
		{
			return *static_cast<cast_type*>(&adds);
		}

		template<typename... Es>
		static const cast_type& cast(OpAddList<Es...> const& adds)
		{
			return *static_cast<cast_type const*>(&adds);
		}

		template<typename... Es>
		static cast_type& cast(OpAddList<Es...>& adds)
		{
			return *static_cast<cast_type*>(&adds);
		}
	};

	template<typename... E0s>
	struct cast_add<symphas::lib::types_list<E0s...>>
	{
		template<typename... Es>
		static OpAdd<E0s...> cast(OpAdd<Es...> const& adds)
		{
			return *static_cast<OpAddList<E0s...> const*>(&adds);
		}

		template<typename... Es>
		static OpAdd<E0s...> cast(OpAdd<Es...>& adds)
		{
			return *static_cast<OpAddList<E0s...>*>(&adds);
		}

		template<typename... Es>
		static OpAddList<E0s...> const& cast(OpAddList<Es...> const& adds)
		{
			return *static_cast<OpAddList<E0s...> const*>(&adds);
		}

		template<typename... Es>
		static OpAddList<E0s...>& cast(OpAddList<Es...>& adds)
		{
			return *static_cast<OpAddList<E0s...>*>(&adds);
		}
	};

	template<>
	struct cast_add<void>
	{
		template<typename... Es, size_t... Is>
		static decltype(auto) cast(OpAdd<Es...> const& adds, std::index_sequence<Is...>)
		{
			return expr::make_add(expr::get<Is>(adds)...);
		}

		template<typename... Es, size_t... Is>
		static decltype(auto) cast(OpAdd<Es...>& adds, std::index_sequence<Is...>)
		{
			return expr::make_add(expr::get<Is>(adds)...);
		}
	};

	inline auto to_add(OpAddList<>)
	{
		return OpVoid{};
	}

	template<typename T>
	auto to_add(T const& e)
	{
		return OpAdd(e);
	}

	template<typename T>
	auto to_add(T& e)
	{
		return OpAdd(e);
	}

	inline auto to_add(OpVoid)
	{
		return OpVoid{};
	}

}


namespace expr
{

	template<size_t N, typename... Es>
	const auto& get(OpAdd<Es...> const& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_add<N, OpAdd<Es...>>::type;
		return symphas::internal::cast_add<Nth_type>::cast(e).term;
	}

	template<size_t N, typename... Es>
	auto& get(OpAdd<Es...>& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_add<N, OpAdd<Es...>>::type;
		return symphas::internal::cast_add<Nth_type>::cast(e).term;
	}

	template<size_t N, typename... Es>
	const auto& get(OpAddList<Es...> const& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_add<N, OpAdd<Es...>>::type;
		return symphas::internal::cast_add<Nth_type>::cast(e).term;
	}

	template<size_t N, typename... Es>
	auto& get(OpAddList<Es...>& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_add<N, OpAdd<Es...>>::type;
		return symphas::internal::cast_add<Nth_type>::cast(e).term;
	}

	template<size_t N, typename... Es>
	decltype(auto) terms_before_n(OpAdd<Es...> const& e)
	{
		if constexpr (N == 0)
		{
			return OpVoid{};
		}
		else if constexpr (N == 1)
		{
			return expr::get<0>(e);
		}
		else
		{
			return symphas::internal::cast_add<void>::cast(e, std::make_index_sequence<N>{});
		}
	}

	template<size_t N, typename... Es>
	decltype(auto) terms_before_n(OpAdd<Es...>& e)
	{
		if constexpr (N == 0)
		{
			return OpVoid{};
		}
		else if constexpr (N == 1)
		{
			return expr::get<0>(e);
		}
		else
		{
			return symphas::internal::cast_add<void>::cast(e, std::make_index_sequence<N>{});
		}
	}


	template<size_t N, typename... Es>
	decltype(auto) terms_after_n(OpAddList<Es...> const& e)
	{
		using ts = symphas::lib::types_after_at_index<N + 1, Es...>;
		if constexpr (N + 1 == sizeof...(Es) - 1)
		{
			return symphas::internal::cast_add<ts>::cast(e).term;
		}
		else if constexpr (N + 1 < sizeof...(Es))
		{
			return symphas::internal::cast_add<ts>::cast(e);
		}
		else
		{
			return OpVoid{};
		}
	}

	template<size_t N, typename... Es>
	decltype(auto) terms_after_n(OpAddList<Es...>& e)
	{
		using ts = symphas::lib::types_after_at_index<N + 1, Es...>;
		if constexpr (N + 1 == sizeof...(Es) - 1)
		{
			return symphas::internal::cast_add<ts>::cast(e).term;
		}
		else if constexpr (N + 1 < sizeof...(Es))
		{
			return symphas::internal::cast_add<ts>::cast(e);
		}
		else
		{
			return OpAddList<>{};
		}
	}

	template<size_t N, typename... Es>
	decltype(auto) terms_after_n(OpAdd<Es...> const& e)
	{
		return symphas::internal::to_add(terms_after_n<N>(*static_cast<OpAddList<Es...> const*>(&e)));
	}

	template<size_t N, typename... Es>
	decltype(auto) terms_after_n(OpAdd<Es...>& e)
	{
		return symphas::internal::to_add(terms_after_n<N>(*static_cast<OpAddList<Es...>*>(&e)));
	}

	template<typename E0, typename E1, typename E2, typename... Es>
	const auto& terms_after_first(OpAdd<E0, E1, E2, Es...> const& e)
	{
		return symphas::internal::cast_add<OpAddList<E1, E2, Es...>>::cast(e);
	}

	template<typename E0, typename E1>
	const auto& terms_after_first(OpAdd<E0, E1> const& e)
	{
		return symphas::internal::cast_add<OpAddList<E1>>::cast(e).term;
	}

	template<typename E0, typename E1, typename E2, typename... Es>
	auto& terms_after_first(OpAdd<E0, E1, E2, Es...>& e)
	{
		return symphas::internal::cast_add<OpAddList<E1, E2, Es...>>::cast(e);
	}

	template<typename E0, typename E1>
	auto& terms_after_first(OpAdd<E0, E1>& e)
	{
		return symphas::internal::cast_add<OpAddList<E1>>::cast(e).term;
	}

	template<typename... Es>
	decltype(auto) terms_after_first(OpAddList<Es...> const& e)
	{
		return terms_after_n<0>(e);
	}

	template<typename... Es>
	decltype(auto) terms_after_first(OpAddList<Es...>& e)
	{
		return terms_after_n<0>(e);
	}
}


template<typename E0, typename... Es>
struct OpAddList<E0, Es...> : OpAddList<Es...>
{
	using parent_type = OpAddList<Es...>;

	OpAddList() : parent_type(), term{ } {}
	OpAddList(E0 const& e, Es const&... es) : parent_type(es...), term{ e } {}
	OpAddList(E0 const& e, OpAdd<Es...> const& rest) : parent_type(rest), term{ e } {}
    //OpAddList(OpAdd<E0, Es...> const& list) : 
    //    parent_type(*static_cast<OpAddList<Es...> const*>(&list)), term{ list.term } {}
    //OpAddList(OpAddList<E0, Es...> const& list) : 
    //    parent_type(*static_cast<OpAddList<Es...> const*>(&list)), term{ list.term } {}

	template<size_t... Is>
	auto _eval(std::index_sequence<Is...>, iter_type n = 0) const
	{
		return (static_cast<symphas::internal::Nth_type_of_add_t<Is, OpAdd<E0, Es...>> const&>(*this).term.eval(n) + ...);
	}

    auto _eval(iter_type n = 0) const
    {
		return _eval(std::make_index_sequence<sizeof...(Es) + 1>{}, n);
        //return term.eval(n) + parent_type::_eval(n);
    }

	bool coeff_sign() const
	{
		return (expr::coeff(term).eval() < 0);
	}

    E0 term;

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = term.print(out);
		return n + symphas::internal::binary_print(*this, out, std::make_index_sequence<sizeof...(Es)>{});
	}

	size_t print(char* out) const
	{
		size_t n = term.print(out);
		return n + symphas::internal::binary_print(*this, out + n, std::make_index_sequence<sizeof...(Es)>{});
	}

	size_t print_length() const
	{
		size_t n = term.print_length() + SYEX_BINARY_FMT_LEN;
		if (OpAddList<Es...>::coeff_sign())
		{
			n -= 1;
		}
		return n + parent_type::print_length();
	}

#endif

	template<typename cast_type>
	friend struct symphas::internal::cast_add;
};

//! Binary expression, the addition of two terms.
/*!
 * Binary addition between two expressions.
 * 
 * \tparam E1 The type of the left hand side expression.
 * \tparam E2 The type of the right hand side expression.
 */
template<typename E0, typename... Es>
struct OpAdd<E0, Es...> : OpExpression<OpAdd<E0, Es...>>, OpAddList<E0, Es...>
{
protected:
	using parent_type = OpAddList<E0, Es...>;

public:

    using parent_type::term;
    using parent_type::coeff_sign;
#ifdef PRINTABLE_EQUATIONS
	using parent_type::print;
	using parent_type::print_length;
#endif

	//! Create the binary addition expression between two expressions.
	/*!
	 * Create the binary addition expression between two expressions.
	 * 
	 * \param a The expression on the left hand side of the addition operator.
	 * \param b The expression on the right hand side of the addition operator.
	 */
	OpAdd(E0 const& e, Es const&... es) : parent_type(e, es...) {}
	OpAdd(E0 const& e, OpAdd<Es...> const& rest) : parent_type(e, rest) {}
	OpAdd(OpAdd<Es...> const& rest, E0 const& e) : parent_type(e, rest) {}
    OpAdd(OpAddList<E0, Es...> const& list) : parent_type(list) {}
    
	OpAdd() : parent_type() {}

	inline auto eval(iter_type n = 0) const
	{
		return parent_type::_eval(n);
	}

	auto operator-() const;



};

template<>
struct OpAdd<> {};

template<typename E0, typename... Es>
OpAdd(E0, OpAdd<Es...>)->OpAdd<E0, Es...>;
template<typename E0, typename... Es>
OpAdd(E0, Es...)->OpAdd<E0, Es...>;
template<typename... Es>
OpAdd(OpAddList<Es...>) -> OpAdd<Es...>;

inline OpVoid expr::make_add()
{
	return OpVoid{};
}

inline OpVoid expr::make_add(OpVoid, OpVoid)
{
	return OpVoid{};
}

template<typename E>
auto expr::make_add(E&& a)
{
	return std::forward<E>(a);
}

template<typename E>
auto expr::make_add(OpVoid, E&& e)
{
	return std::forward<E>(e);
}

template<typename E>
auto expr::make_add(E&& e, OpVoid)
{
	return std::forward<E>(e);
}

template<typename... E0s, typename E0>
auto expr::make_add(OpAdd<E0s...> const& add, E0 const& e0)
{
	return OpAdd<E0, E0s...>(add, e0);
}

template<typename... E0s>
auto expr::make_add(OpAdd<E0s...> const& add, OpVoid)
{
	return add;
}

template<typename... E0s, typename E0>
auto expr::make_add(E0 const& e0, OpAdd<E0s...> const& add)
{
	return OpAdd<E0, E0s...>(add, e0);
}

template<typename... E0s>
auto expr::make_add(OpVoid, OpAdd<E0s...> const& add)
{
	return add;
}

template<typename A, typename B>
auto expr::make_add(OpExpression<A> const& a, OpExpression<B> const& b)
{
	return OpAdd<A, B>(*static_cast<A const*>(&a), *static_cast<B const*>(&b));
}

template<typename A, typename B>
auto expr::make_add(OpOperator<A> const& a, OpExpression<B> const& b)
{
	return OpAdd<A, B>(*static_cast<A const*>(&a), *static_cast<B const*>(&b));
}

template<typename A, typename B>
auto expr::make_add(OpExpression<A> const& a, OpOperator<B> const& b)
{
	return OpAdd<A, B>(*static_cast<A const*>(&a), *static_cast<B const*>(&b));
}

template<typename E1, typename E2, typename E3, typename... Es>
auto expr::make_add(E1&& a, E2&& b, E3&& c, Es&& ...es)
{
	return make_add(make_add(std::forward<E1>(a), std::forward<E2>(b)), std::forward<E3>(c), std::forward<Es>(es)...);
}


template<typename E1, typename E2>
auto operator+(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_add(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
}

template<typename E1, typename E2>
auto operator-(OpExpression<E1> const& a, OpExpression<E2> const& b)
{
	return expr::make_add(*static_cast<E1 const*>(&a), -*static_cast<E2 const*>(&b));
}


namespace symphas::internal
{
	using expr::has_coeff;
	using expr::has_pmi_coeff;
	using expr::has_nmi_coeff;

#ifdef PRINTABLE_EQUATIONS

	inline size_t mul_print_left(FILE* out, bool neg)
	{
		if (neg)
		{
			return fprintf(out, "%s", SYEX_MUL_FMT_A);
		}
		else
		{
			return fprintf(out, "%s", SYEX_MUL_FMT_AA);
		}
	}

	inline size_t mul_print_sep(FILE* out, bool neg1, bool neg2)
	{
		if (neg1 && neg2)
		{
			return fprintf(out, SYEX_MUL_SEP);
		}
		else if (neg1)
		{
			return fprintf(out, SYEX_MUL_SEP_A);
		}
		else if (neg2)
		{
			return fprintf(out, SYEX_MUL_SEP_B);
		}
		else
		{
			return fprintf(out, SYEX_MUL_SEP_AB);
		}
	}

	inline size_t mul_print_right(FILE* out, bool neg)
	{
		if (neg)
		{
			return fprintf(out, "%s", SYEX_MUL_FMT_B);
		}
		else
		{
			return fprintf(out, "%s", SYEX_MUL_FMT_BB);
		}
	}

	inline size_t mul_print_left(char* out, bool neg)
	{
		if (neg)
		{
			return sprintf(out, "%s", SYEX_MUL_FMT_A);
		}
		else
		{
			return sprintf(out, "%s", SYEX_MUL_FMT_AA);
		}
	}

	inline size_t mul_print_sep(char* out, bool neg1, bool neg2)
	{
		if (neg1 && neg2)
		{
			return sprintf(out, SYEX_MUL_SEP);
		}
		else if (neg1)
		{
			return sprintf(out, SYEX_MUL_SEP_A);
		}
		else if (neg2)
		{
			return sprintf(out, SYEX_MUL_SEP_B);
		}
		else
		{
			return sprintf(out, SYEX_MUL_SEP_AB);
		}
	}

	inline size_t mul_print_right(char* out, bool neg)
	{
		if (neg)
		{
			return sprintf(out, "%s", SYEX_MUL_FMT_B);
		}
		else
		{
			return sprintf(out, "%s", SYEX_MUL_FMT_BB);
		}
	}

	template<typename E1, typename E2>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = expr::eval(expr::coeff(*static_cast<E1 const*>(&a))) < 0 || expr::is_add<E1>;
		bool neg2 = expr::eval(expr::coeff(*static_cast<E2 const*>(&b))) < 0 || expr::is_add<E2>;

		n += mul_print_left(out, neg1);
		n += static_cast<E1 const*>(&a)->print(out);
		n += mul_print_sep(out, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out);
		n += mul_print_right(out, neg2);
		return n;
	}

	template<typename E1, typename E2>
	size_t mul_print(FILE* out, OpOperator<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = false;
		bool neg2 = expr::eval(expr::coeff(*static_cast<E2 const*>(&b))) < 0 || expr::is_add<E2>;

		n += mul_print_left(out, neg1);
		n += static_cast<E1 const*>(&a)->print(out);
		n += mul_print_sep(out, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out);
		n += mul_print_right(out, neg2);
		return n;
	}

	template<typename E1, typename E2>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpOperator<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = expr::eval(expr::coeff(*static_cast<E1 const*>(&a))) < 0 || expr::is_add<E1>;
		bool neg2 = false;

		n += mul_print_left(out, neg1);
		n += static_cast<E1 const*>(&a)->print(out);
		n += mul_print_sep(out, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out);
		n += mul_print_right(out, neg2);
		return n;
	}

	template<typename E1, typename E2>
	size_t mul_print(FILE* out, OpOperator<E1> const& a, OpOperator<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = false;
		bool neg2 = false;

		n += mul_print_left(out, neg1);
		n += static_cast<E1 const*>(&a)->print(out);
		n += mul_print_sep(out, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out);
		n += mul_print_right(out, neg2);
		return n;
	}


	template<typename E1, typename E2>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = expr::eval(expr::coeff(*static_cast<E1 const*>(&a))) < 0 || expr::is_add<E1>;
		bool neg2 = expr::eval(expr::coeff(*static_cast<E2 const*>(&b))) < 0 || expr::is_add<E2>;

		n += mul_print_left(out + n, neg1);
		n += static_cast<E1 const*>(&a)->print(out + n);
		n += mul_print_sep(out + n, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out + n);
		n += mul_print_right(out + n, neg2);
		return n;
	}

	template<typename E1, typename E2>
	size_t mul_print(char* out, OpOperator<E1> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = false;
		bool neg2 = expr::eval(expr::coeff(*static_cast<E2 const*>(&b))) < 0 || expr::is_add<E2>;

		n += mul_print_left(out + n, neg1);
		n += static_cast<E1 const*>(&a)->print(out + n);
		n += mul_print_sep(out + n, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out + n);
		n += mul_print_right(out + n, neg2);
		return n;
	}

	template<typename E1, typename E2>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpOperator<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = expr::eval(expr::coeff(*static_cast<E1 const*>(&a))) < 0 || expr::is_add<E1>;
		bool neg2 = false;

		n += mul_print_left(out + n, neg1);
		n += static_cast<E1 const*>(&a)->print(out + n);
		n += mul_print_sep(out + n, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out + n);
		n += mul_print_right(out + n, neg2);
		return n;
	}

	template<typename E1, typename E2>
	size_t mul_print(char* out, OpOperator<E1> const& a, OpOperator<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = false;
		bool neg2 = false;

		n += mul_print_left(out + n, neg1);
		n += static_cast<E1 const*>(&a)->print(out + n);
		n += mul_print_sep(out + n, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out + n);
		n += mul_print_right(out + n, neg2);
		return n;
	}

	template<typename V, typename E1, typename T, typename E2>
	size_t mul_print(FILE* out, OpIntegral<V, E1, T> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = expr::eval(expr::coeff(a)) < 0;
		bool neg2 = true;

		n += mul_print_left(out, neg1);
		n += a.print(out);
		n += mul_print_sep(out, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out);
		n += mul_print_right(out, neg2);
		return n;
	}

	template<typename V, typename E1, typename T, typename E2>
	size_t mul_print(FILE* out, OpIntegral<V, E1, T> const& a, OpOperator<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = expr::eval(expr::coeff(a)) < 0;
		bool neg2 = false;

		n += mul_print_left(out, neg1);
		n += a.print(out);
		n += mul_print_sep(out, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out);
		n += mul_print_right(out, neg2);
		return n;
	}

	template<typename E>
	size_t mul_print(FILE* out, OpIdentity, OpExpression<E> const& b)
	{
		size_t n = 0;
		n += mul_print_left(out, true);
		n += static_cast<E const*>(&b)->print(out);
		n += mul_print_right(out, true);
		return n;
	}

	template<typename E>
	size_t mul_print(FILE* out, OpIdentity, OpOperator<E> const& b)
	{
		size_t n = 0;
		n += mul_print_left(out, true);
		n += static_cast<E const*>(&b)->print(out);
		n += mul_print_right(out, true);
		return n;
	}

	template<typename E>
	size_t mul_print(FILE* out, OpNegIdentity, OpExpression<E> const& b)
	{
		size_t n = 0;
		n += mul_print_left(out, true);
		n += (-*static_cast<E const*>(&b)).print(out);
		n += mul_print_right(out, true);
		return n;
	}

	template<typename E>
	size_t mul_print(FILE* out, OpNegIdentity, OpOperator<E> const& b)
	{
		size_t n = 0;
		n += mul_print_left(out, true);
		n += (-*static_cast<E const*>(&b)).print(out);
		n += mul_print_right(out, true);
		return n;
	}

	template<typename E>
	size_t mul_print(FILE* out, OpOperator<E> const& a, OpIdentity)
	{
		size_t n = 0;
		n += mul_print_left(out, true);
		n += static_cast<E const*>(&a)->print(out);
		n += mul_print_right(out, true);
		return n;
	}

	template<typename E>
	size_t mul_print(FILE* out, OpExpression<E> const& a, OpIdentity)
	{
		size_t n = 0;
		n += mul_print_left(out, true);
		n += static_cast<E const*>(&a)->print(out);
		n += mul_print_right(out, true);
		return n;
	}

	template<typename V, typename E1, typename T, typename E2>
	size_t mul_print(char* out, OpIntegral<V, E1, T> const& a, OpExpression<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = expr::eval(expr::coeff(a)) < 0;
		bool neg2 = true;

		n += mul_print_left(out + n, neg1);
		n += a.print(out + n);
		n += mul_print_sep(out + n, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out + n);
		n += mul_print_right(out + n, neg2);
		return n;
	}

	template<typename V, typename E1, typename T, typename E2>
	size_t mul_print(char* out, OpIntegral<V, E1, T> const& a, OpOperator<E2> const& b)
	{
		size_t n = 0;
		bool neg1 = expr::eval(expr::coeff(a)) < 0;
		bool neg2 = false;

		n += mul_print_left(out + n, neg1);
		n += a.print(out + n);
		n += mul_print_sep(out + n, neg1, neg2);
		n += static_cast<E2 const*>(&b)->print(out + n);
		n += mul_print_right(out + n, neg2);
		return n;
	}

	template<typename V, typename E1, typename T, typename E2>
	size_t mul_print(FILE* out, OpExpression<E1> const& a, OpIntegral<V, E2, T> const& b)
	{
		bool neg1 = expr::eval(expr::coeff(*static_cast<E1 const*>(&a))) < 0;
		bool neg2 = false;

		size_t n = 0;
		n += fprintf(out, "%s", SYEX_MUL_FMT_AA);
		n += static_cast<E1 const*>(&a)->print(out);
		n += fprintf(out, "%s", SYEX_MUL_SEP_OP);
		n += b.print(out);
		n += fprintf(out, "%s", SYEX_MUL_FMT_BB);
		return n;
	}

	template<typename V, typename E1, typename T, typename E2>
	size_t mul_print(FILE* out, OpOperator<E1> const& a, OpIntegral<V, E2, T> const& b)
	{
		bool neg1 = false;
		bool neg2 = false;

		size_t n = 0;
		n += fprintf(out, "%s", SYEX_MUL_FMT_AA);
		n += static_cast<E1 const*>(&a)->print(out);
		n += fprintf(out, "%s", SYEX_MUL_SEP_OP);
		n += b.print(out);
		n += fprintf(out, "%s", SYEX_MUL_FMT_BB);
		return n;
	}

	template<typename V, typename E1, typename T, typename E2>
	size_t mul_print(char* out, OpExpression<E1> const& a, OpIntegral<V, E2, T> const& b)
	{
		size_t n = 0;
		bool neg1 = expr::eval(expr::coeff(*static_cast<E1 const*>(&a))) < 0;
		bool neg2 = false;

		n += mul_print_left(out + n, neg1);
		n += static_cast<E1 const*>(&a)->print(out);
		n += mul_print_sep(out + n, neg1, neg2);
		n += b.print(out + n);
		n += mul_print_right(out + n, neg2);
		return n;
	}

	template<typename V, typename E1, typename T, typename E2>
	size_t mul_print(char* out, OpOperator<E1> const& a, OpIntegral<V, E2, T> const& b)
	{
		size_t n = 0;
		bool neg1 = false;
		bool neg2 = false;

		n += mul_print_left(out + n, neg1);
		n += static_cast<E1 const*>(&a)->print(out);
		n += mul_print_sep(out + n, neg1, neg2);
		n += b.print(out + n);
		n += mul_print_right(out + n, neg2);
		return n;
	}

	template<typename V1, typename E1, typename T1, typename V2, typename E2, typename T2>
	size_t mul_print(FILE* out, OpIntegral<V1, E1, T1> const& a, OpIntegral<V2, E2, T2> const& b)
	{
		bool neg1 = false;
		bool neg2 = false;

		size_t n = 0;
		n += fprintf(out, "%s", SYEX_MUL_FMT_AA);
		n += a.print(out);
		n += fprintf(out, "%s", SYEX_MUL_SEP_OP);
		n += b.print(out);
		n += fprintf(out, "%s", SYEX_MUL_FMT_BB);
		return n;
	}

	template<typename V1, typename E1, typename T1, typename V2, typename E2, typename T2>
	size_t mul_print(char* out, OpIntegral<V1, E1, T1> const& a, OpIntegral<V2, E2, T2> const& b)
	{
		size_t n = 0;
		bool neg1 = false;
		bool neg2 = false;

		n += mul_print_left(out + n, neg1);
		n += a.print(out + n);
		n += mul_print_sep(out + n, neg1, neg2);
		n += b.print(out + n);
		n += mul_print_right(out + n, neg2);
		return n;
	}

	template<typename E>
	size_t mul_print(char* out, OpIdentity, OpExpression<E> const& b)
	{
		size_t n = 0;
		n += mul_print_left(out + n, true);
		n += static_cast<E const*>(&b)->print(out + n);
		n += mul_print_right(out + n, true);
		return n;
	}

	template<typename E>
	size_t mul_print(char* out, OpIdentity, OpOperator<E> const& b)
	{
		size_t n = 0;
		n += mul_print_left(out + n, true);
		n += static_cast<E const*>(&b)->print(out + n);
		n += mul_print_right(out + n, true);
		return n;
	}

	template<typename E>
	size_t mul_print(char* out, OpNegIdentity, OpExpression<E> const& b)
	{
		size_t n = 0;
		n += mul_print_left(out + n, true);
		n += (-*static_cast<E const*>(&b)).print(out + n);
		n += mul_print_right(out + n, true);
		return n;
	}

	template<typename E>
	size_t mul_print(char* out, OpNegIdentity, OpOperator<E> const& b)
	{
		size_t n = 0;
		n += mul_print_left(out + n, true);
		n += (-*static_cast<E const*>(&b)).print(out + n);
		n += mul_print_right(out + n, true);
		return n;
	}

	template<typename E>
	size_t mul_print(char* out, OpExpression<E> const& a, OpIdentity)
	{
		size_t n = 0;
		n += mul_print_left(out + n, true);
		n += static_cast<E const*>(&a)->print(out + n);
		n += mul_print_right(out + n, true);
		return n;
	}

	template<typename E>
	size_t mul_print(char* out, OpOperator<E> const& a, OpIdentity)
	{
		size_t n = 0;
		n += mul_print_left(out + n, true);
		n += static_cast<E const*>(&a)->print(out + n);
		n += mul_print_right(out + n, true);
		return n;
	}

	template<typename coeff_t, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
	size_t mul_print(char* out, coeff_t const& a, DynamicIndex const& b)
	{
		size_t n = 0;
		n += expr::print_with_coeff(out, a);
		n += b.print(out + n);
		return n;
	}

	template<typename coeff_t, std::enable_if_t<expr::is_coeff<coeff_t>, int> = 0>
	size_t mul_print(FILE* out, coeff_t const& a, DynamicIndex const& b)
	{
		size_t n = 0;
		n += expr::print_with_coeff(out, a);
		n += b.print(out);
		return n;
	}

#endif

}

//! Binary expression, the multiplication of two terms.
template<typename E1, typename E2>
struct OpBinaryMul : OpExpression<OpBinaryMul<E1, E2>>
{
	OpBinaryMul(E1 const& a, E2 const& b) : a{ a }, b{ b } {}

	OpBinaryMul() : a{ E1{} }, b{ E2{} } {}

	inline auto eval(iter_type n = 0) const
	{
		return a.eval(n) * b.eval(n);
	}

	auto operator-() const;

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return symphas::internal::mul_print(out, a, b);
	}

	size_t print(char* out) const
	{
		return symphas::internal::mul_print(out, a, b);
	}

	size_t print_length() const
	{
		return a.print_length() + b.print_length() + SYEX_MUL_FMT_LEN;
	}

#endif

	E1 a;		//!< Left hand side of the binary operator.
	E2 b;		//!< Right hand side of the binary operator.

};

namespace expr
{
	//! Constructs the binary multiplication expression.
	/*!
	 * Directly constructs the binary multiplication expression between two
	 * expressions without applying any rules.
	 *
	 * \param a The left hand side expression.
	 * \param b The right hand side expression.
	 */
	template<typename E1, typename E2>
	auto make_mul(E1 const& a, E2 const& b)
	{
		return OpBinaryMul<E1, E2>(a, b);
	}

	template<typename E1, typename E2>
	auto dot(OpExpression<E1> const& a, OpExpression<E2> const& b);

	template<typename E1, typename E2>
	auto dot(OpExpression<E1> const& a, OpOperator<E2> const& b);

	template<typename E1, typename E2>
	auto dot(OpOperator<E1> const& a, OpExpression<E2> const& b);

	template<typename E1, typename E2>
	auto dot(OpOperator<E1> const& a, OpOperator<E2> const& b);
}


inline auto DynamicIndex::operator-() const
{
	return expr::make_mul(OpNegIdentity{}, (*this));
}



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

template<typename A, typename B, typename E2,
	typename std::enable_if_t<(expr::eval_type<OpBinaryDiv<A, B>>::rank == 0 || expr::eval_type<E2>::rank == 0), int> = 0>
auto operator*(OpBinaryDiv<A, B> const& a, OpOperator<E2> const& b)
{
	return OpOperatorChain(a, *static_cast<const E2*>(&b));
}

template<typename A, typename B, typename E2,
	typename std::enable_if_t<(expr::eval_type<OpBinaryMul<A, B>>::rank == 0 || expr::eval_type<E2>::rank == 0), int> = 0>
auto operator*(OpBinaryMul<A, B> const& a, OpOperator<E2> const& b)
{
	return OpOperatorChain(a, *static_cast<const E2*>(&b));
}

//! Distributing an RHS expression between operands in addition expression.
template<typename... As, typename B2>
auto operator*(OpAdd<As...> const& a, OpOperator<B2> const& b)
{
	return OpOperatorChain(a, *static_cast<const B2*>(&b));
}

template<typename E1, typename E2,
	typename std::enable_if_t<(expr::eval_type<E1>::rank > 0 && expr::eval_type<E2>::rank > 0), int> = 0>
auto operator*(OpExpression<E1> const& a, OpOperator<E2> const& b)
{
	return expr::dot(*static_cast<const E1*>(&a), *static_cast<const E2*>(&b));
}

template<typename A, typename B, typename E,
	typename std::enable_if_t<(expr::eval_type<OpBinaryMul<A, B>>::rank > 0 && expr::eval_type<E>::rank > 0), int> = 0>
auto operator*(OpBinaryMul<A, B> const& a, OpOperator<E> const& b)
{
	return expr::dot(a, *static_cast<const E*>(&b));
}

template<typename A, typename B, typename E,
	typename std::enable_if_t<(expr::eval_type<OpBinaryMul<A, B>>::rank > 0 && expr::eval_type<E>::rank > 0), int> = 0>
auto operator*(OpBinaryDiv<A, B> const& a, OpOperator<E> const& b)
{
	return expr::dot(a, *static_cast<const E*>(&b));
}


//! Binary expression, the division of two terms.
template<typename E1, typename E2>
struct OpBinaryDiv : OpExpression<OpBinaryDiv<E1, E2>>
{
	OpBinaryDiv(E1 const& a, E2 const& b) : a{ a }, b{ b } {}

	//template<typename AA = E1, typename BB = E2,
	//	typename = std::enable_if_t<(std::is_default_constructible<AA>::value&& std::is_default_constructible<BB>::value), int>>
	OpBinaryDiv() : a{ E1{} }, b{ E2{} } {}

	inline auto eval(iter_type n) const
	{
		return a.eval(n) / b.eval(n);
	}

	auto operator-() const;

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = 0;
		n += fprintf(out, "%s", SYEX_DIV_FMT_A);
		n += a.print(out);
		n += fprintf(out, "%s", SYEX_DIV_SEP);
		n += b.print(out);
		n += fprintf(out, "%s", SYEX_DIV_FMT_B);
		return n;

	}

	size_t print(char* out) const
	{
		size_t n = 0;
		n += sprintf(out + n, "%s", SYEX_DIV_FMT_A);
		n += a.print(out + n);
		n += sprintf(out + n, "%s", SYEX_DIV_SEP);
		n += b.print(out + n);
		n += sprintf(out + n, "%s", SYEX_DIV_FMT_B);
		return n;
	}

	size_t print_length() const
	{
		return a.print_length() + b.print_length() + SYEX_DIV_FMT_LEN;
	}

#endif

	E1 a;		//!< Left hand side of the binary operator.
	E2 b;		//!< Right hand side of the binary operator.
};

namespace symphas::internal
{

	template<typename>
	struct coeff_factor_numerator
	{
		static const size_t value = 1;
	};

	template<size_t N, size_t D>
	struct coeff_factor_numerator<OpFractionLiteral<N, D>>
	{
		static const size_t value = N;
	};

	template<typename>
	struct coeff_factor_denominator
	{
		static const size_t value = 1;
	};

	template<size_t N, size_t D>
	struct coeff_factor_denominator<OpFractionLiteral<N, D>>
	{
		static const size_t value = D;
	};

	template<typename E>
	struct coeff_factor_impl
	{
		using type = typename coeff_factor_impl<expr::coeff_t<E>>::type;
	};

	template<typename T, size_t... Ns>
	struct coeff_factor_impl<OpTensor<T, Ns...>>
	{
		using type = typename coeff_factor_impl<T>::type;
	};

	template<typename T>
	struct coeff_factor_impl<OpLiteral<T>>
	{
		using type = OpIdentity;
	};

	template<>
	struct coeff_factor_impl<OpNegIdentity>
	{
		using type = OpIdentity;
	};

	template<>
	struct coeff_factor_impl<OpIdentity>
	{
		using type = OpIdentity;
	};

	template<size_t N, size_t D>
	struct coeff_factor_impl<OpNegFractionLiteral<N, D>>
	{
		using type = OpFractionLiteral<N, D>;
	};

	template<size_t N, size_t D>
	struct coeff_factor_impl<OpFractionLiteral<N, D>>
	{
		using type = OpFractionLiteral<N, D>;
	};

	template<typename... As>
	struct coeff_factor_impl<OpAdd<As...>>
	{
		using type = decltype(expr::make_fraction<
			GCD_of<coeff_factor_numerator<typename coeff_factor_impl<As>::type>::value...>,
			GCD_of<coeff_factor_denominator<typename coeff_factor_impl<As>::type>::value...>>());
	};

	template<>
	struct coeff_factor_impl<OpVoid>
	{
		using type = OpVoid;
	};


	template<typename E>
	using coeff_factor = typename coeff_factor_impl<E>::type;
}

namespace expr
{

	//! Constructs the binary division expression.
	/*!
	 * Directly constructs the binary division expression between two
	 * expressions without applying any rules.
	 *
	 * \param a The left hand side expression.
	 * \param b The right hand side expression.
	 */
	template<typename E1, typename E2, size_t R1 = expr::eval_type<E1>::rank, size_t R2 = expr::eval_type<E2>::rank,
		typename std::enable_if_t<(R1 == 0 && R2 == 0), int> = 0>
	auto make_div(OpEvaluable<E1> const& a, OpEvaluable<E2> const& b)
	{
		auto fr1 = symphas::internal::coeff_factor<E1>{};
		auto fr2 = symphas::internal::coeff_factor<E2>{};

		auto fr = (fr1 / fr2);
		constexpr int N = symphas::internal::coeff_factor_numerator<decltype(fr)>::value;
		constexpr int D = symphas::internal::coeff_factor_denominator<decltype(fr)>::value;

		auto _fr1 = expr::make_integer<N>() / fr1;
		auto _fr2 = expr::make_integer<D>() / fr2;


		return OpBinaryDiv(_fr1 * (*static_cast<const E1*>(&a)), _fr2 * (*static_cast<const E2*>(&b)));
	}

	template<typename E1, typename E2, size_t R1 = expr::eval_type<E1>::rank, size_t R2 = expr::eval_type<E2>::rank,
		typename std::enable_if_t<(R1 > 0 || R2 > 0), int> = 0>
	auto make_div(OpEvaluable<E1> const& a, OpEvaluable<E2> const& b)
	{
		return OpBinaryDiv((*static_cast<const E1*>(&a)), (*static_cast<const E2*>(&b)));
	}

	template<typename E2>
	auto make_div(OpVoid, OpExpression<E2> const& b)
	{
		return OpVoid{};
	}
}

namespace expr
{

	template<size_t N0, size_t N1, size_t... Ns>
	constexpr auto sym_T = expr::make_tensor<N0, N1, Ns...>();
	template<size_t I, size_t N>
	constexpr auto sym_R = expr::make_row_vector<I, N>(OpIdentity{});
	template<size_t I, size_t N>
	constexpr auto sym_C = expr::make_column_vector<I, N>(OpIdentity{});
}


namespace expr::symbols
{
	namespace
	{

		template<char... n>
		struct char_to_lit
		{
			constexpr auto operator()()
			{
				return expr::make_literal(symphas::lib::char_to_number<n...>{}());
			}
		};
		
		template<char... n>
		struct char_to_val
		{
			template<char... ln, size_t... Is>
			constexpr auto left(symphas::lib::char_list<ln...>, std::index_sequence<Is...>)
			{
				using namespace symphas::lib;
				return (expr::make_integer<get_value<ln>() * fixed_pow<10, sizeof...(Is) - 1 - Is>>() + ... + OpVoid{});
			}

			template<char... ln>
			constexpr auto left(symphas::lib::char_list<ln...>)
			{
				return left(symphas::lib::char_list<ln...>{}, std::make_index_sequence<sizeof...(ln)>{});
			}

			template<char... rn, size_t... Is>
			constexpr auto right(symphas::lib::char_list<'.', rn...>, std::index_sequence<Is...>)
			{
				using namespace symphas::lib;
				return (expr::make_fraction<get_value<rn>(), fixed_pow<10, Is + 1>>() + ... + OpVoid{});
			}

			template<char... rn>
			constexpr auto right(symphas::lib::char_list<'.', rn...>)
			{
				return right(symphas::lib::char_list<'.', rn...>{}, std::make_index_sequence<sizeof...(rn)>{});
			}

			constexpr auto right(symphas::lib::char_list<>)
			{
				return OpVoid{};
			}

			constexpr auto operator()()
			{
				using namespace symphas::lib;
				return left(typename filter_decimal_left<symphas::lib::char_list<>, symphas::lib::char_list<n...>>::type{})
					+ right(typename filter_decimal_right<symphas::lib::char_list<>, symphas::lib::char_list<n...>>::type{});
			}
		};
	}

	template<char... n>
	constexpr auto operator ""_n()
	{
		return char_to_val<n...>{}();
	}

	template<char... n>
	constexpr auto operator ""_c()
	{
		return char_to_lit<n...>{}();
	}
}

#undef SYEX_BINARY_FMT
#undef SYEX_BINARY_FMT_LEN
#undef SYEX_ADD_SEP
#undef SYEX_SUB_SEP
#undef SYEX_MUL_FMT
#undef SYEX_MUL_FMT_LATEX
//#undef SYEX_MUL_FMT_LEN
#undef SYEX_MUL_FMT_LATEX_LEN
#undef SYEX_DIV_FMT
#undef SYEX_DIV_FMT_LATEX
#undef SYEX_DIV_FMT_LEN
#undef SYEX_DIV_FMT_LATEX_LEN



