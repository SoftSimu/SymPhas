
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


#ifdef SOL_EXPORTS
#define DLLSOL DLLEXPORT
#else
#define DLLSOL DLLIMPORT
#endif

namespace params
{
	//! Sets the resize factor for regional grids.
	/*! 
	 * Sets the factor by which the minimum iterable area found for a field is expanded to become 
	 * the regional grid.
	 */
	DLLSOL extern double regional_resize_factor;

	//! Sets the resize time for regional grids.
	/*!
	 * Sets the time increment at which the regional grid dimensions and position is reevalauted.
	 */
	DLLSOL extern double regional_resize_time;

	//! Defines that the original shape found for the regional grid is fixed.
	/*!
	 * Defines that the original shape, meaning the precise dimensions of the grid, which is 
	 * determined from the initial state of the field, will remain fixed throughout the simulation.
	 * The position will still be reevaluated at the time interval provided by regional_resize_time.
	 */
	DLLSOL extern bool regional_resize_is_fixed;

	//! This value is added to the minimum value to get the cutoff.
	/*!
	 * The regional grid area is determined by evaluating which elements exceed the cutoff. The
	 * cutoff is determined by adding this value to the minimum value of the field.
	 */
	DLLSOL extern double regional_resize_cutoff_eps;



	//! This enumeration is used as `ParamNamesSol << bool`.
	enum ParamNamesSol
	{
		REG_RESIZE_FACTOR,	//!< See params::regional_resize_factor.
		REG_RESIZE_TIME,	//!< See params::regional_resize_time.
		REG_RESIZE_FIXED,	//!< See params::regional_resize_is_fixed.
		REG_RESIZE_CUTOFF	//!< See params::regional_resize_cutoff_eps.
	};

	struct sol_param_flag
	{
		ParamNamesSol param;
		bool value;
	};

	struct sol_param_value
	{
		ParamNamesSol param;
		double value;
	};

	inline sol_param_flag operator<<(ParamNamesSol param, bool value)
	{
		return { param, value };
	}

	inline sol_param_value operator<<(ParamNamesSol param, double value)
	{
		return { param, value };
	}

	inline void set_param(ParamNamesSol param, bool value)
	{
		switch (param)
		{
		case ParamNamesSol::REG_RESIZE_FIXED:
			regional_resize_is_fixed = value;
			break;
		default:
			break;
		}
	}

	inline void set_param(ParamNamesSol param, double value)
	{
		switch (param)
		{
		case ParamNamesSol::REG_RESIZE_FACTOR:
			regional_resize_factor = value;
			break;
		case ParamNamesSol::REG_RESIZE_TIME:
			regional_resize_time = value;
			break;
		case ParamNamesSol::REG_RESIZE_CUTOFF:
			regional_resize_cutoff_eps = value;
			break;
		default:
			break;
		}
	}

	inline void set_param(sol_param_flag const& input)
	{
		auto [param, value] = input;
		set_param(param, value);
	}

	inline void set_param(sol_param_value const& input)
	{
		auto [param, value] = input;
		set_param(param, value);
	}
}

#define REGIONAL_GRID_RESIZE_FACTOR params::regional_resize_factor
#define REGIONAL_GRID_RESIZE_TIME params::regional_resize_time
#define REGIONAL_GRID_RESIZE_IS_FIXED params::regional_resize_is_fixed
#define REGIONAL_GRID_CUTOFF_EPS params::regional_resize_cutoff_eps
#define REGIONAL_GRID_RELATIVE_DIMS_EPS 0.05

bool add_solution_params(param_map_type& param_map);

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

	//! Implemented for parallelization routines.
	template<typename... Ts>
	static void synchronize(Ts&&...) {}

protected:

	PhaseFieldSystem();

	symphas::interval_data_type get_extended_intervals(symphas::interval_data_type vdata)
	{
		if (params::extend_boundary)
		{
			for (auto& [_, interval] : vdata)
			{
				interval.set_count(interval.get_count() + 2 * BOUNDARY_DEPTH);
				interval.interval_to_domain();
			}
		}
		return vdata;
	}

};

template<typename T, size_t D>
void PhaseFieldSystem<BoundaryGrid, T, D>::update(iter_type index, double time)
{
	BoundaryGroup<T, D>::update_boundaries(*this, index, time);
}

template<typename T>
struct regional_system_info_type
{
	iter_type next_resize;
	double resize_delta;
	bool fixed_resize;
	T cutoff;

	regional_system_info_type() :
		next_resize{ 0 }, resize_delta{ REGIONAL_GRID_RESIZE_TIME }, fixed_resize{ REGIONAL_GRID_RESIZE_IS_FIXED }, cutoff{} {}
	regional_system_info_type(iter_type next_resize, double resize_delta, bool fixed_resize, T cutoff) :
		next_resize{ next_resize }, resize_delta{ resize_delta }, fixed_resize{ fixed_resize }, cutoff{ cutoff } {}
};


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
struct PhaseFieldSystem<RegionalGrid, T, D> : PersistentSystemData<RegionalGrid<T, D>>, BoundaryGroup<T, D>
{
	using parent_type = PersistentSystemData<RegionalGrid<T, D>>;
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

	//! Implemented for parallelization routines.
	template<typename... Ts>
	static void synchronize(Ts&&...) {}

protected:

	PhaseFieldSystem();

	symphas::interval_data_type get_extended_intervals(symphas::interval_data_type vdata)
	{
		if (params::extend_boundary)
		{
			for (auto& [_, interval] : vdata)
			{
				interval.set_count(interval.get_count() + 2 * BOUNDARY_DEPTH);
				interval.interval_to_domain();
			}
		}
		return vdata;
	}


	regional_system_info_type<T> regional_info;
};

#ifdef USING_MPI

template<typename T, size_t D>
struct PhaseFieldSystem<RegionalGridMPI, T, D> : PersistentSystemData<RegionalGridMPI<T, D>>, BoundaryGroup<T, D>
{
	using parent_type = PersistentSystemData<RegionalGridMPI<T, D>>;
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

	//! Implemented for parallelization routines.
	template<typename... Ts>
	static void synchronize(Ts&&...);

protected:

	PhaseFieldSystem();

	symphas::interval_data_type get_extended_intervals(symphas::interval_data_type vdata)
	{
		if (params::extend_boundary)
		{
			for (auto& [_, interval] : vdata)
			{
				interval.set_count(interval.get_count() + 2 * BOUNDARY_DEPTH);
				interval.interval_to_domain();
			}
		}
		return vdata;
	}


	regional_system_info_type<T> regional_info;
public:
	symphas::multi_thr_info_type thr_info;
};

#endif

namespace grid
{
	template<typename T, size_t D>
	auto min_value(RegionalGrid<T, D> const& grid)
	{
		symphas::data_iterator it(grid.values);
		auto e = std::min_element(it, it + grid.region.len);
		return std::min(*e, grid.empty);
	}

	template<typename T, size_t D>
	auto min_value(RegionalGrid<any_vector_t<T, D>, D> const& grid)
	{
		symphas::data_iterator_region it(grid, grid::get_iterable_domain(grid));

		auto e = std::min_element(it, it + grid.region.len, 
			[] (auto a, auto b) {
				using std::abs;
				using symphas::lib::abs;
				return abs(any_vector_t<T, D>(a)) < abs(any_vector_t<T, D>(b));
			});

		using std::abs;
		using symphas::math::abs;
		using symphas::lib::abs;
		return (abs(any_vector_t<T, D>(*e)) < abs(grid.empty)) ? *e : grid.empty;
	}
}

namespace symphas::internal
{
	template<typename T>
	void set_value_for_resize(T& cutoff, T min, T eps = {})
	{
		cutoff = min + eps;
	}

	template<typename T, size_t D>
	void set_value_for_resize(any_vector_t<T, D>& cutoff, any_vector_t<T, D> const& min, T eps = {})
	{
		for (iter_type i = 0; i < D; ++i)
		{
			cutoff[i] = min[i] + eps;
		}
	}

	template<typename T, size_t D>
	void set_value_for_resize(any_vector_t<T, D>& value, multi_value<D, T> const& min, T eps = {})
	{
		for (iter_type i = 0; i < D; ++i)
		{
			value[i] = min[i] + eps;
		}
	}


	template<typename T, size_t D>
	void set_value_for_resize(T(&value)[D], multi_value<D, T> const& min, T eps = {})
	{
		for (iter_type i = 0; i < D; ++i)
		{
			value[i] = min[i] + eps;
		}
	}

	template<typename T, size_t D>
	void update_regional_system(regional_system_info_type<T>& regional_info, RegionalGrid<T, D>& grid, symphas::grid_info& info, double time)
	{
		if (regional_info.next_resize == 0)
		{
			auto min0 = grid::min_value(grid);
			set_value_for_resize(regional_info.cutoff, min0, REGIONAL_GRID_CUTOFF_EPS);
			set_value_for_resize(grid.empty, min0);

			grid::resize_adjust_region(grid, regional_info.cutoff, REGIONAL_GRID_RESIZE_FACTOR);

			for (auto& [axis, interval] : info)
			{
				iter_type i = symphas::axis_to_index(axis);
				double offset = grid.region.boundary_size * interval.width();
				interval.set_interval(
					grid.region.origin[i] * interval.width() + offset,
					(grid.region.origin[i] + grid.region.dims[i] - 1) * interval.width() - offset);
			}
			++regional_info.next_resize;
		}
		else if (regional_info.next_resize > 0)
		{
			iter_type current_resize = iter_type(time / regional_info.resize_delta);
			if (current_resize >= regional_info.next_resize)
			{
				if (regional_info.fixed_resize)
				{
					grid::adjust_region(grid, regional_info.cutoff);
				}
				else
				{
					grid::resize_adjust_region(grid, regional_info.cutoff, REGIONAL_GRID_RESIZE_FACTOR, REGIONAL_GRID_RELATIVE_DIMS_EPS);
				}
				regional_info.next_resize = current_resize + 1;
			}
		}
	}
}

template<typename T, size_t D>
void PhaseFieldSystem<RegionalGrid, T, D>::update(iter_type index, double time)
{
	symphas::internal::update_regional_system(regional_info, *this, info, time);
	BoundaryGroup<T, D>::update_boundaries(*this, index, time);
}

#ifdef USING_MPI

template<typename T, size_t D>
void PhaseFieldSystem<RegionalGridMPI, T, D>::update(iter_type index, double time)
{
	if (thr_info.index_in_node())
	{
		symphas::internal::update_regional_system(regional_info, *this, info, time);
		BoundaryGroup<T, D>::update_boundaries(*this, index, time);
	}
}

#endif


template<typename T, size_t D>
PhaseFieldSystem<BoundaryGrid, T, D>::PhaseFieldSystem() : parent_type{}, BoundaryGroup<T, D>{} {}

template<typename T, size_t D>
PhaseFieldSystem<BoundaryGrid, T, D>::PhaseFieldSystem(
	symphas::init_data_type const& tdata, 
	symphas::interval_data_type const& vdata, 
	symphas::b_data_type const& bdata, size_t id) :
	parent_type{ get_extended_intervals(vdata), id }, BoundaryGroup<T, D>{ info.intervals, bdata }
{
	grid::region_interval<D> region(parent_type::dims, BOUNDARY_DEPTH);
	symphas::internal::populate_tdata(tdata, *static_cast<Grid<T, D>*>(this), &info, region, id);
}

template<typename T, size_t D>
using BoundarySystem = PhaseFieldSystem<BoundaryGrid, T, D>;


template<typename T, size_t D>
PhaseFieldSystem<RegionalGrid, T, D>::PhaseFieldSystem() : 
	parent_type{}, BoundaryGroup<T, D>{},
	regional_info{}
{}

template<typename T, size_t D>
PhaseFieldSystem<RegionalGrid, T, D>::PhaseFieldSystem(
	symphas::init_data_type const& tdata,
	symphas::interval_data_type const& vdata,
	symphas::b_data_type const& bdata, size_t id) :
	parent_type{ get_extended_intervals(vdata), id }, BoundaryGroup<T, D>{ info.intervals, bdata }, 
	regional_info{}
{
	grid::region_interval<D> region(RegionalGrid<T, D>::region.dims, RegionalGrid<T, D>::region.boundary_size);
	symphas::internal::populate_tdata(tdata, *this, &info, region, id);

	if (grid::has_subdomain(vdata))
	{
		grid::resize_adjust_region(*this, vdata);
		regional_info.next_resize = -1;
	}
}

template<typename T, size_t D>
using RegionalSystem = PhaseFieldSystem<RegionalGrid, T, D>;

#ifdef USING_MPI

template<typename T, size_t D>
PhaseFieldSystem<RegionalGridMPI, T, D>::PhaseFieldSystem() :
	parent_type{}, BoundaryGroup<T, D>{},
	regional_info{}
{}

template<typename T, size_t D>
PhaseFieldSystem<RegionalGridMPI, T, D>::PhaseFieldSystem(
	symphas::init_data_type const& tdata,
	symphas::interval_data_type const& vdata,
	symphas::b_data_type const& bdata, size_t id) :
	parent_type{ get_extended_intervals(vdata), tdata.thr_info }, BoundaryGroup<T, D>{ info.intervals, bdata },
	regional_info{}, thr_info{ tdata.thr_info }
{
	if (thr_info.index_in_node())
	{
		grid::region_interval<D> region(RegionalGridMPI<T, D>::region.dims, RegionalGridMPI<T, D>::region.boundary_size);
		symphas::internal::populate_tdata(tdata, *this, &info, region, id);

		if (grid::has_subdomain(vdata))
		{
			grid::resize_adjust_region(*this, vdata);
			regional_info.next_resize = -1;
		}
	}
}

template<typename T, size_t D>
using RegionalSystemMPI = PhaseFieldSystem<RegionalGridMPI, T, D>;

DEFINE_SYMBOL_ID((typename T, size_t D), (RegionalGridMPI<T, D>), return data.values)

#endif



DEFINE_SYMBOL_ID((typename T, size_t D), (BoundaryGrid<T, D>), return data.values)
DEFINE_SYMBOL_ID((typename T, size_t D), (RegionalGrid<T, D>), return data.values)
DEFINE_BASE_DATA_INHERITED((template<typename, size_t> typename grid_t, typename T, size_t D), (PhaseFieldSystem<grid_t, T, D>), (grid_t<T, D>))



namespace symphas::internal
{
	template<template<typename, size_t> typename G, typename T, size_t D>
	struct data_value_type<PhaseFieldSystem<G, T, D>> : data_value_type<G<T, D>> {};
}