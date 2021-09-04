
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
 * PURPOSE: Defines a system using a grid with boundaries.
 *
 * ***************************************************************************
 */

#pragma once

#include "system.h"
#include "boundarygroup.h"
#include "expressionlogic.h"


//! Representation of a phase field and associated characteristics.
/*!
 * Specialization based on a phase field system in which boundaries must be
 * managed. See ::PhaseFieldSystem.
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
template<typename T, size_t D>
struct PhaseFieldSystem<BoundaryGrid, T, D> : PersistentSystemData<BoundaryGrid<T, D>>, BoundaryGroup<T, D>
{
	using parent_type = PersistentSystemData<BoundaryGrid<T, D>>;
	using parent_type::info;



	//! Given the system information, generate the phase field.
	/*!
	 * Given the system information, generate the phase field. This includes
	 * populating the values of the field from the initial conditions.
	 * This also includes setting the boundary conditions of the system.
	 * 
	 * \param tdata The information about the initial conditions of the
	 * phase field system.
	 * \param vdata The interval data about the system.
	 * \param bdata The information about the boundary conditions of the phase
	 * field system.
	 * \param id A special identifier number. This does not have to be unique, 
	 * but should be unique between different systems in the same model.
	 */
	PhaseFieldSystem(symphas::init_data_type const& tdata, symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata, size_t id = 0);

	//! Update the system to prepare for the next solver iteration.
	/*!
	 * Update the system to prepare for the next solver iteration. The values
	 * on the boundaries are updated based on the solution index and solution
	 * time. The way the boundaries are updated is based on the boundary
	 * types, and is implemented higher in the class heirarchy.
	 * 
	 * \param index The solution index.
	 * \param time The solution time.
	 */
	void update(iter_type index = 0, double time = 0);


protected:

	PhaseFieldSystem();

	symphas::interval_data_type get_extended_intervals(symphas::interval_data_type vdata)
	{
		if (params::extend_boundary)
		{
			for (auto& [_, interval] : vdata)
			{
				interval.set_interval_count(
					interval.left(), interval.right(),
					interval.count() + 2 * THICKNESS);
			}
		}
		return vdata;
	}

};





template<typename T, size_t D>
void PhaseFieldSystem<BoundaryGrid, T, D>::update(iter_type index, double time)
{
	BoundaryGroup<T, D>::template update_boundaries(*this, index, time);
}




template<typename T, size_t D>
PhaseFieldSystem<BoundaryGrid, T, D>::PhaseFieldSystem() : parent_type{}, BoundaryGroup<T, D>{} {}

template<typename T, size_t D>
PhaseFieldSystem<BoundaryGrid, T, D>::PhaseFieldSystem(
	symphas::init_data_type const& tdata, 
	symphas::interval_data_type const& vdata, 
	symphas::b_data_type const& bdata, size_t id) :
	parent_type{ get_extended_intervals(vdata), id }, BoundaryGroup<T, D>{ info.intervals, bdata }
{
	symphas::internal::populate_tdata(tdata, *static_cast<Grid<T, D>*>(this), &info, id);
}

template<typename T, size_t D>
using BoundarySystem = PhaseFieldSystem<BoundaryGrid, T, D>;





namespace expr
{


#ifdef MULTITHREAD

	//! Specialization based on expr::result_interior(OpExpression<E> const&, T*, len_type(&)[3]).
	template<typename T, size_t D, typename E>
	void result_interior(OpExpression<E> const& e, BoundaryGrid<T, D>& grid)
	{
		expr::result_interior(*static_cast<const E*>(&e), grid.values, grid.inners, grid.len_inner);
	}


#else

	//! Specialization based on expr::result_interior(OpExpression<E> const&, T*, len_type(&)[3]).
	template<typename T, size_t D, typename E>
	void result_interior(OpExpression<E> const& e, BoundaryGrid<T, D>& grid)
	{
		expr::result_interior(*static_cast<const E*>(&e), grid.values, grid.dims);
	}

#endif

}

DEFINE_SYMBOL_ID((typename T, size_t D), (BoundaryGrid<T, D>), return data.values)



