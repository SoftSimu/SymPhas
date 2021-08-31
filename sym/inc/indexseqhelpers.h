
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
 * PURPOSE: Used in manipulating sequences of numbers, introduced in
 * particular to help with manipulating lists of variable indices.
 *
 * ***************************************************************************
 */

#pragma once


namespace symphas::lib
{

	// **************************************************************************************


	//! Swap two numbers based on their order.
	template<size_t Q0, size_t Q1, typename std::enable_if_t<(Q0 < Q1), int> = 0>
	constexpr auto swap_ids()
	{
		return std::index_sequence<Q0, Q1>{};
	}

	//! Swap two numbers based on their order.
	template<size_t Q0, size_t Q1, typename std::enable_if_t<(Q0 >= Q1), int> = 0>
	constexpr auto swap_ids()
	{
		return std::index_sequence<Q1, Q0>{};
	}


	//! Sort numbers in an `std::index_sequence`.
	/*!
	 * Implementation of the bubble sort algorithm of `std::index_sequence`.
	 * Given a single index_sequence, a new sequence is returned containing a
	 * sorted arrangement of the input values.
	 */
	template<size_t... Qs>
	constexpr auto sort_ids(std::index_sequence<Qs...>);


	template<size_t... Qs, size_t Qmax>
	constexpr auto sort_ids_depth(std::index_sequence<Qs...>, std::index_sequence<Qmax>)
	{
		return seq_join(sort_ids(std::index_sequence<Qs...>{}), std::index_sequence<Qmax>{});
	}


	//! Specialization based on symphas::lib::sort_ids().
	template<size_t... Qs, size_t Qmax>
	constexpr auto sort_ids(std::index_sequence<Qs...>, std::index_sequence<Qmax>, std::index_sequence<>)
	{
		return sort_ids_depth(std::index_sequence<Qs...>{}, std::index_sequence<Qmax>{});
	}

	//! Specialization based on symphas::lib::sort_ids().
	template<size_t... Ys, size_t Qmax, size_t Q0, size_t... Qs, typename std::enable_if_t<(Qmax > Q0), int> = 0>
	constexpr auto sort_ids(std::index_sequence<Ys...>, std::index_sequence<Qmax>, std::index_sequence<Q0, Qs...>)
	{
		return sort_ids(std::index_sequence<Ys..., Q0>{}, std::index_sequence<Qmax>{}, std::index_sequence<Qs...>{});
	}

	//! Specialization based on symphas::lib::sort_ids().
	template<size_t... Ys, size_t Qmax, size_t Q0, size_t... Qs,
		typename std::enable_if_t<(Qmax <= Q0), int> = 0>
	constexpr auto sort_ids(std::index_sequence<Ys...>, std::index_sequence<Qmax>, std::index_sequence<Q0, Qs...>)
	{
		return sort_ids(std::index_sequence<Ys..., Qmax>{}, std::index_sequence<Q0>{}, std::index_sequence<Qs...>{});
	}

	//! Specialization based on symphas::lib::sort_ids().
	template<size_t Qmax, size_t Q0, size_t... Qs, 
		typename std::enable_if_t<(Qmax > Q0), int> = 0>
	constexpr auto sort_ids(std::index_sequence<Qmax>, std::index_sequence<Q0, Qs...>)
	{
		return sort_ids(std::index_sequence<Q0>{}, std::index_sequence<Qmax>{}, std::index_sequence<Qs...>{});
	}

	template<size_t Qmax, size_t Q0, size_t... Qs, 
		typename std::enable_if_t<(Qmax <= Q0), int> = 0>
	constexpr auto sort_ids(std::index_sequence<Qmax>, std::index_sequence<Q0, Qs...>)
	{
		return sort_ids(std::index_sequence<Qmax>{}, std::index_sequence<Q0>{}, std::index_sequence<Qs...>{});
	}


	//! Specialization based on symphas::lib::sort_ids().
	template<size_t Q0, size_t... Qs, 
		typename std::enable_if_t<(sizeof...(Qs) > 1), int> = 0>
	constexpr auto sort_ids(std::index_sequence<Q0, Qs...>)
	{
		return sort_ids(std::index_sequence<Q0>{}, std::index_sequence<Qs...>{});
	}

	//! Specialization based on symphas::lib::sort_ids().
	template<size_t Q0, size_t... Qs, 
		typename std::enable_if_t<(sizeof...(Qs) == 1), int> = 0>
	constexpr auto sort_ids(std::index_sequence<Q0, Qs...>)
	{
		return swap_ids<Q0, Qs...>();
	}

	template<size_t Q0, size_t... Qs, 
		typename std::enable_if_t<(sizeof...(Qs) == 0), int> = 0>
	constexpr auto sort_ids(std::index_sequence<Q0, Qs...>)
	{
		return std::index_sequence<Q0>{};
	}

	//! See symphas::lib::sort_ids().
	inline constexpr auto sort_ids(std::index_sequence<>)
	{
		return std::index_sequence<>{};
	}


	// **************************************************************************************


	inline constexpr auto fold_unique_ids_apply(std::index_sequence<>)
	{
		return std::index_sequence<>{};
	}

	template<size_t Q0>
	constexpr auto fold_unique_ids_apply(std::index_sequence<Q0>)
	{
		return std::index_sequence<Q0>{};
	}

	template<size_t Q0, size_t Q1, size_t... Qs, 
		typename std::enable_if_t<(Q0 == Q1), int> = 0>
	constexpr auto fold_unique_ids_apply(std::index_sequence<Q0, Q1, Qs...>)
	{
		return fold_unique_ids_apply(std::index_sequence<Q1, Qs...>{});
	}

	template<size_t Q0, size_t Q1, size_t... Qs, 
		typename std::enable_if_t<(Q0 != Q1), int> = 0>
	constexpr auto fold_unique_ids_apply(std::index_sequence<Q0, Q1, Qs...>)
	{
		return seq_join(std::index_sequence<Q0>{}, fold_unique_ids_apply(std::index_sequence<Q1, Qs...>{}));
	}



	//! Generate a list of unique values from a given sequence.
	template<size_t... Qs>
	constexpr auto fold_unique_ids(std::index_sequence<Qs...>)
	{
		return fold_unique_ids_apply(sort_ids(std::index_sequence<Qs...>{}));
	}



	// **************************************************************************************



	//! Return only values less than the prescribed limit.
	/*!
	 * This algorithm expects a sorted list
	 * Returns the list of values that are strictly less than `Zmax`.
	 * 
	 * \tparam Zmax The limit values are compared to.
	 */
	template<size_t Zmax>
	auto split_under(std::index_sequence<>)
	{
		return std::index_sequence<>{};
	}


	//! Specialization based on symphas::lib::split_under().
	template<size_t Zmax, size_t Q0, size_t... Qs, 
		typename std::enable_if_t<(Zmax > Q0), int> = 0>
	auto split_under(std::index_sequence<Q0, Qs...>)
	{
		return seq_join(std::index_sequence<Q0>{}, split_under<Zmax>(std::index_sequence<Qs...>{}));
	}

	//! Specialization based on symphas::lib::split_under().
	template<size_t Zmax, size_t Q0, size_t... Qs, 
		typename std::enable_if_t<(Zmax <= Q0), int> = 0>
	auto split_under(std::index_sequence<Q0, Qs...>)
	{
		return std::index_sequence<>{};
	}


	//! Return only values less than the prescribed limit.
	/*!
	 * This algorithm expects a sorted list
	 * Returns the list of values that are greater than or equal to `Zmax`.
	 * 
	 * \tparam Zmax The limit values are compared to.
	 */
	template<size_t Zmax>
	auto split_at_above(std::index_sequence<>)
	{
		return std::index_sequence<>{};
	}

	//! Specialization based on symphas::lib::split_at_above().
	template<size_t Zmin, size_t Q0, size_t... Qs, 
		typename std::enable_if_t<(Zmin > Q0), int> = 0>
	auto split_at_above(std::index_sequence<Q0, Qs...>)
	{
		return split_at_above<Zmin>(std::index_sequence<Qs...>{});
	}

	//! Specialization based on symphas::lib::split_at_above().
	template<size_t Zmin, size_t Q0, size_t... Qs,
		typename std::enable_if_t<(Zmin <= Q0), int> = 0>
	auto split_at_above(std::index_sequence<Q0, Qs...>)
	{
		return std::index_sequence<Q0, Qs...>{};
	}

}



