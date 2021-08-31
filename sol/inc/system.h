
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
 * PURPOSE: Implementation of all required elements to manage and represent
 * a phase field system.
 *
 * ***************************************************************************
 */

#pragma once

#include <memory>

#include "initialconditions.h"
#include "systemlib.h"



// ****************************************************************************************************************


 //! Representation of a phase field and associated characteristics.
 /*!
  * Represents a phase field system and manages the data thereof, based on the
  * information given by ::SystemInfo.
  *
  * When the phase field is created, it will automatically populate the phase
  * field values from the given initial conditions.
  *
  * The phase field is based on a PersistenSystemData, meaning it has the ability
  * to persist its data to a file if the **io** module is enabled.
  *
  * \tparam T The type of the phase field.
  * \tparam D The dimension of the phase field.
  */
template<template<typename, size_t> typename G, typename T, size_t D>
struct PhaseFieldSystem : PersistentSystemData<G<T, D>>
{
	using parent_type = PersistentSystemData<G<T, D>>;
	using parent_type::info;
	using parent_type::data_len;

	PhaseFieldSystem(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, size_t id = 0);
	PhaseFieldSystem(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const&, size_t id = 0) :
		PhaseFieldSystem(tdata, vdata, id) {}

	//! There is no update mechanism required for the base object.
	inline void update(iter_type = 0, double = 0) const {}


protected:

	PhaseFieldSystem();
};




// **********************************************************************************************


template<template<typename, size_t> typename G, typename T, size_t D>
PhaseFieldSystem<G, T, D>::PhaseFieldSystem() : parent_type{} {}

template<template<typename, size_t> typename G, typename T, size_t D>
PhaseFieldSystem<G, T, D>::PhaseFieldSystem(
	symphas::init_data_type const& tdata, 
	symphas::interval_data_type const& vdata, 
	size_t id) :
	parent_type{ vdata, id }
{
	symphas::internal::populate_tdata(tdata, static_cast<G<T, D>&>(*this), &info, id);
}




template<typename T, size_t D>
using System = PhaseFieldSystem<Grid, T, D>;





