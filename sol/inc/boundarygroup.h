
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
 * PURPOSE: Defines a group of boundaries, which is used for a whole
 * system in applying the boundary conditions to a grid.
 *
 * ***************************************************************************
 */

#include "boundaryupdate.h"



//! Manages the boundaries across the whole `D`-dimensional system.
/*!
 * The boundary group manages and controls the update process of the boundaries
 * in a `D`-dimensional system. The boundaries have corresponding dimension
 * `D-1`.
 * 
 * \tparam T The system value type.
 * \tparam D The dimension of the system.
 */
template<typename T, size_t D>
struct BoundaryGroup
{
	BoundaryType types[D * 2];						//!< The system boundary types.
	grid::Boundary<T, D - 1>* boundaries[D * 2];	//!< Boundaries, indexing based on symphas::index_to_side().

	//! Creates boundaries from the boundary and interval data.
	/*!
	 * Creates boundaries from the boundary and interval data.
	 * 
	 * \param vdata The interval data used to initialize the boundaries. The
	 * range of the boundaries will be taken from this information.
	 * \param bdata The information about the system boundaries which is used to
	 * initialize each boundary in the group.
	 */
	BoundaryGroup(symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata);
	BoundaryGroup(BoundaryGroup<T, D> const&);
	BoundaryGroup(BoundaryGroup<T, D>&&) noexcept;
	BoundaryGroup<T, D>& operator=(BoundaryGroup<T, D>);


	//! Update the boundaries of a grid.
	/*!
	 * Iterate over the sides and update all the boundaries
	 * in this object based on the type of the boundary.
	 * 
	 * \param grid The grid to which the boundaries apply.
	 * \param index The current solution index.
	 * \param time The current time of the solution data.
	 */
	template<size_t I = D * 2 - 1>
	void update_boundaries(Grid<T, D> &grid, iter_type index, double time);


	friend void swap(BoundaryGroup<T, D>& first, BoundaryGroup<T, D>& second)
	{
		using std::swap;
		swap(first.types, second.types);
		swap(first.boundaries, second.boundaries);
	}

	~BoundaryGroup();

protected:

	BoundaryGroup();
};





template<typename T, size_t D>
BoundaryGroup<T, D>::BoundaryGroup(symphas::interval_data_type const& vdata, symphas::b_data_type const& bdata)
{
	for (auto&& [side, data] : bdata)
	{
		types[symphas::side_to_index(side)] = data.type;
	}

	for (auto&& [side, data] : bdata)
	{
		boundaries[symphas::side_to_index(side)] = symphas::internal::new_boundary<T, D - 1>(data);
	}

	symphas::boundary_setup<D - 1>{}(boundaries, bdata, vdata);
}



template<typename T, size_t D>
BoundaryGroup<T, D>::BoundaryGroup() : types{}, boundaries{} 
{
	for (iter_type i = 0; i < D * 2; ++i)
	{
		types[i] = BoundaryType::DEFAULT;
		boundaries[i] = nullptr;
	}
}


template<typename T, size_t D>
BoundaryGroup<T, D>::BoundaryGroup(BoundaryGroup<T, D> const& other)
{
	for (iter_type i = 0; i < D * 2; ++i)
	{
		types[i] = other.types[i];
		boundaries[i] = other.boundaries[i]->new_copy();
	}
}

template<typename T, size_t D>
BoundaryGroup<T, D>::BoundaryGroup(BoundaryGroup<T, D>&& other) noexcept : BoundaryGroup()
{
	swap(*this, other);
}

template<typename T, size_t D>
BoundaryGroup<T, D>& BoundaryGroup<T, D>::operator=(BoundaryGroup<T, D> other)
{
	swap(*this, other);
	return *this;
}



template<typename T, size_t D>
BoundaryGroup<T, D>::~BoundaryGroup()
{
	for (iter_type i = 0; i < D * 2; ++i)
	{
		delete boundaries[i];
	}
}


namespace symphas::internal
{
	template<size_t I>
	struct update_boundary_call
	{
		template<typename B, typename T, size_t D>
		void operator()(B* boundaries, BoundaryType types[D * 2], Grid<T, D>& grid, iter_type, double time)
		{
			switch (types[I])
			{
			case BoundaryType::PERIODIC:
				symphas::internal::update_boundary<BoundaryType::PERIODIC, symphas::index_to_side(I), D - 1>{}(boundaries[I], grid);
				break;
			case BoundaryType::DEFAULT:
				symphas::internal::update_boundary<BoundaryType::DEFAULT, symphas::index_to_side(I), D - 1>{}(boundaries[I], grid, time);
				break;
			case BoundaryType::OPEN:
				throw;
			}
		}
	};

	template<size_t I>
	struct update_boundary_recurse;

	template<>
	struct update_boundary_recurse<0>
	{
		template<typename B, typename T, size_t D>
		void operator()(B* boundaries, BoundaryType types[D * 2], Grid<T, D>& grid, iter_type index, double time)
		{
			update_boundary_call<0>{}(boundaries, types, grid, index, time);
		}
	};

	template<size_t I>
	struct update_boundary_recurse
	{
		template<typename B, typename T, size_t D>
		void operator()(B* boundaries, BoundaryType types[D * 2], Grid<T, D>& grid, iter_type index, double time)
		{
			update_boundary_call<I>{}(boundaries, types, grid, index, time);
			update_boundary_recurse<I - 1>{}(boundaries, types, grid, index, time);
		}
	};
}


template<typename T, size_t D>
template<size_t I>
void BoundaryGroup<T, D>::update_boundaries(Grid<T, D> &grid, iter_type index, double time)
{
	symphas::internal::update_boundary_recurse<I>{}(boundaries, types, grid, index, time);
}





