
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
 * MODULE:  sol
 * PURPOSE: Manages a group of provisional systems, used by the phase field
 * model in order to manage the provisional variables.
 *
 * ***************************************************************************
 */

#pragma once

#include <tuple>

#include "provisionalsystem.h"
#include "solversystem.h"


#ifdef CUSTOM_SYS_HEADER
#include CUSTOM_SYS_HEADER
#endif

//! Manages the provisional systems of the given type.
/*!
 * Manages the provisional systems of the given type.
 * 
 * \tparam ProvisionalSystemApplied The type of the provisional system
 * managed by this object, which is determined by the solver used for the
 * phase field model.
 * \tparam D The dimension of the provisional systems.
 * \tparam P... The types of the provisional variables.
 */
template<template<typename, size_t> typename ProvisionalSystemApplied, size_t D, typename... P>
struct ProvisionalSystemGroup
{
	ProvisionalSystemGroup(const symphas::interval_data_type vdata, const symphas::b_data_type bdata) :
		_s{ construct_systems(vdata, bdata) } {}

	template<size_t I>
	auto& grid()
	{
		return std::get<I>(_s).as_grid();
	}

	template<size_t I>
	const auto& grid() const
	{
		return std::get<I>(_s).as_grid();
	}

	// updates the systems by directly calling their update functions
	void update_systems(iter_type index, double time)
	{
		update_systems(index, time, std::make_index_sequence<sizeof...(P)>{});
	}

	//! Container of the provisional data.
	std::tuple<ProvisionalSystemApplied<P, D>...> _s;

protected:


	template<size_t... Is>
	void update_systems(iter_type index, double time, std::index_sequence<Is...>)
	{
		((..., std::get<Is>(_s).update(index, time)));
	}


	std::tuple<ProvisionalSystemApplied<P, D>...> construct_systems(const symphas::interval_data_type vdata, const symphas::b_data_type bdata)
	{
		return std::make_tuple((ProvisionalSystemApplied<P, D>(vdata, bdata))...);
	}


};


