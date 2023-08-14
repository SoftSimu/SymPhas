
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
 * MODULE:  datatypes
 * PURPOSE: Defines how to update the boundaries, depending on the type
 * and dimension.
 *
 * ***************************************************************************
 */

#pragma once

#include "boundaryupdatedata.h"


namespace symphas::internal
{
	//! Update a boundary for a particular grid.
	/*! 
	 * The boundary update utility with
	 * primary goal to immediately choose the update routine in order to 
	 * avoid branching and virtual calls.
	 * 
	 * The grid and boundary do not have to be necessarily related.
	 */
	template<BoundaryType type, Side side, size_t D>
	struct update_boundary
	{
		//! Updates boundary with time index.
		/*!
		 * The type is given boundary is identified and the update algorithm
		 * is used to update values in grid. The time is passed in cases where
		 * the boundary requires the time to generate values.
		 * 
		 * \param boundary The boundary which is used to update the grid.
		 * \param grid The grid where values are updated.
		 * \param time The current solution time.
		 */
		template<typename T>
		void operator()(const grid::Boundary<T, D>* boundary, Grid<T, D + 1> &grid, double time);

		//! Updates boundary.
		/*!
		 * The type is given boundary is identified and the update algorithm
		 * is used to update values in grid. 
		 *
		 * \param boundary The boundary which is used to update the grid.
		 * \param grid The grid where values are updated.
		 */
		template<typename T>
		void operator()(const grid::Boundary<T, D>* boundary, Grid<T, D + 1> &grid);
		
		//! Updates boundary with time index.
		/*!
		 * The type is given boundary is identified and the update algorithm
		 * is used to update values in grid. The time is passed in cases where
		 * the boundary requires the time to generate values.
		 *
		 * \param boundary The boundary which is used to update the grid.
		 * \param grid The grid where values are updated.
		 * \param time The current solution time.
		 */
		template<typename T>
		void operator()(const grid::Boundary<T, D>* boundary, RegionalGrid<T, D + 1>& grid, double time);

		//! Updates boundary.
		/*!
		 * The type is given boundary is identified and the update algorithm
		 * is used to update values in grid.
		 *
		 * \param boundary The boundary which is used to update the grid.
		 * \param grid The grid where values are updated.
		 */
		template<typename T>
		void operator()(const grid::Boundary<T, D>* boundary, RegionalGrid<T, D + 1>& grid);
	};
}


// *********************************************************************
/* an implementation of all BOUNDARY ITERATION ALGORITHMS
 */


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::FRONT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	// edges

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_TOP({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_BOTTOM({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BACK, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	// edges

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_TOP({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_BOTTOM({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID3_LEFT({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	// edges

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_LEFT({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_LEFT({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_LEFT_TOP({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = -grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_LEFT_BOTTOM({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	// corners

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_LEFT_TOP_FRONT({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_LEFT_BOTTOM_FRONT({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_TOP_LEFT({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_BOTTOM_LEFT({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::RIGHT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID3_RIGHT({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	// edges

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_RIGHT({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_RIGHT({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_RIGHT_TOP({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = -grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_RIGHT_BOTTOM({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);


	// corners

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_TOP_RIGHT({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_BOTTOM_RIGHT({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_RIGHT_TOP_BACK({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_RIGHT_BOTTOM_BACK({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::TOP, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_TOP({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BOTTOM, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BOTTOM({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}




template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::FRONT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::BACK, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_ALL({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::LEFT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID3_LEFT_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::RIGHT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID3_RIGHT_ALL({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::TOP, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_TOP_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::BOTTOM, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BOTTOM_ALL({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}





template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::FRONT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
    iter_type offset;
	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_3A({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_LEFT_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_RIGHT_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::BACK, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_3A({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_RIGHT_ALL({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_LEFT_ALL({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::LEFT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID3_LEFT_3A({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::RIGHT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID3_RIGHT_3A({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::TOP, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_TOP_3A({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_LEFT_TOP_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_RIGHT_TOP_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::BOTTOM, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BOTTOM_3A({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = -grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_LEFT_BOTTOM_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = -grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_RIGHT_BOTTOM_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}






template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::FRONT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
    iter_type offset;
	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_3AA({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_TOP_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_FRONT_BOTTOM_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::BACK, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_3AA({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) - grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_TOP_ALL({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

	offset = grid.dims[0] * grid.dims[1] * (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) + grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BACK_BOTTOM_ALL({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::LEFT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID3_LEFT_3AA({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::RIGHT, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID3_RIGHT_3AA({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::TOP, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_TOP_3AA({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);


}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::BOTTOM, 2>::operator()(const grid::Boundary<T, 2>*, Grid<T, 3>& grid)
{
	iter_type offset;
	offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID3_BOTTOM_3AA({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1], grid.dims[2]);

}





/* 2 dimensional boundaries
 *
 */

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 1>::operator()(const grid::Boundary<T, 1>*, Grid<T, 2>& grid)
{
	iter_type offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;

	ITER_GRID2_LEFT({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) * grid.dims[0] + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID2_LEFT_TOP({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0]);

	offset = -(grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) * grid.dims[0] + (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID2_BOTTOM_LEFT({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::RIGHT, 1>::operator()(const grid::Boundary<T, 1>*, Grid<T, 2>& grid)
{
	iter_type offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID2_RIGHT({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1]);

	offset = -(grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) * grid.dims[0] - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID2_RIGHT_BOTTOM({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);

	offset = (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) * grid.dims[0] - (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID2_TOP_RIGHT({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0]);
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::TOP, 1>::operator()(const grid::Boundary<T, 1>*, Grid<T, 2>& grid)
{
	iter_type offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID2_TOP({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0]);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BOTTOM, 1>::operator()(const grid::Boundary<T, 1>*, Grid<T, 2>& grid)
{
	iter_type offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID2_BOTTOM({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1]);
}





template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::LEFT, 1>::operator()(const grid::Boundary<T, 1>*, Grid<T, 2>& grid)
{
	iter_type offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID2_LEFT_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0], grid.dims[1]);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::RIGHT, 1>::operator()(const grid::Boundary<T, 1>*, Grid<T, 2>& grid)
{
	iter_type offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID2_RIGHT_ALL({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1]);
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::TOP, 1>::operator()(const grid::Boundary<T, 1>*, Grid<T, 2>& grid)
{
	iter_type offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID2_TOP_ALL({
		grid[INDEX] = grid[INDEX + offset];
		}, grid.dims[0]);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::BOTTOM, 1>::operator()(const grid::Boundary<T, 1>*, Grid<T, 2>& grid)
{
	iter_type offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
	ITER_GRID2_BOTTOM_ALL({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0], grid.dims[1]);
}



// one dimensional boundaries

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 0>::operator()(const grid::Boundary<T, 0>*, Grid<T, 1>& grid)
{
	iter_type offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID1_LEFT({
		grid[INDEX] = grid[INDEX + offset];
		});
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::RIGHT, 0>::operator()(const grid::Boundary<T, 0>*, Grid<T, 1>& grid)
{
	iter_type offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
	ITER_GRID1_RIGHT({
		grid[INDEX] = grid[INDEX - offset];
		}, grid.dims[0]);
}





template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::FRONT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::FRONT, Side::FRONT, Side::FRONT>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::TOP, Side::FRONT>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::TOP, Side::FRONT>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::BOTTOM, Side::FRONT>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::BOTTOM, Side::FRONT>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BACK, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::BACK, Side::BACK, Side::BACK>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::TOP, Side::BACK>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::TOP, Side::BACK>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::BOTTOM, Side::BACK>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::BOTTOM, Side::BACK>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT, Side::LEFT>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT, Side::TOP>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT, Side::BOTTOM>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT, Side::BACK>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT, Side::FRONT>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::RIGHT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT, Side::RIGHT>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT, Side::TOP>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT, Side::BOTTOM>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT, Side::BACK>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT, Side::FRONT>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::TOP, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::TOP, Side::TOP>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::TOP, Side::FRONT>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::TOP, Side::BACK>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BOTTOM, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM, Side::BOTTOM>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM, Side::FRONT>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM, Side::BACK>{}, grid);
}





template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3YZ, Side::FRONT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::FRONT, Side::FRONT>{}, grid);
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3YZ, Side::BACK, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::BACK, Side::BACK>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3YZ, Side::LEFT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3YZ, Side::RIGHT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3YZ, Side::TOP, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3YZ, Side::BOTTOM, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, grid);
}



template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::FRONT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::FRONT>{}, grid);
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::BACK, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::BACK>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::LEFT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::RIGHT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::TOP, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::BOTTOM, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{

}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::FRONT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::BACK, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::LEFT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::RIGHT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::TOP, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::TOP>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::BOTTOM, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM>{}, grid);

}



template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::FRONT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::FRONT>{}, grid);
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::BACK, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::BACK>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::LEFT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, grid);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::RIGHT, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::TOP, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::TOP>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::BOTTOM, 2>::operator()(const grid::Boundary<T, 2>*, RegionalGrid<T, 3>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM>{}, grid);

}



template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 1>::operator()(const grid::Boundary<T, 1>*, RegionalGrid<T, 2>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::TOP>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::BOTTOM>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::RIGHT, 1>::operator()(const grid::Boundary<T, 1>*, RegionalGrid<T, 2>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::TOP>{}, grid);
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::BOTTOM>{}, grid);
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::TOP, 1>::operator()(const grid::Boundary<T, 1>*, RegionalGrid<T, 2>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::TOP>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BOTTOM, 1>::operator()(const grid::Boundary<T, 1>*, RegionalGrid<T, 2>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM>{}, grid);
}





template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::LEFT, 1>::operator()(const grid::Boundary<T, 1>*, RegionalGrid<T, 2>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::RIGHT, 1>::operator()(const grid::Boundary<T, 1>*, RegionalGrid<T, 2>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, grid);
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::TOP, 1>::operator()(const grid::Boundary<T, 1>*, RegionalGrid<T, 2>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::BOTTOM, 1>::operator()(const grid::Boundary<T, 1>*, RegionalGrid<T, 2>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM>{}, grid);
}



// one dimensional boundaries

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 0>::operator()(const grid::Boundary<T, 0>*, RegionalGrid<T, 1>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, grid);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::RIGHT, 0>::operator()(const grid::Boundary<T, 0>*, RegionalGrid<T, 1>& grid)
{
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, grid);
}







// *********************************************************************
/* DEFAULT BOUNDARY ALGORITHMS
 */

 // 3 dimensions

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::FRONT, 2>::operator()(const grid::Boundary<T, 2>* b, Grid<T, 3>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 2, BoundaryType::DEFAULT> const*>(b);

	const double* x = bd->v;
	const double* y = bd->v + 2;
	double h[2];

	// go backwards or forwards in iteration depending on the interval
	int fx = (x[0] < x[1]) ? 1 : -1;
	int fy = (y[0] < y[1]) ? 1 : -1;
	h[0] = bd->h[0] * fx;
	h[1] = bd->h[1] * fy;

	iter_type L = grid.dims[0];
	iter_type M = grid.dims[1];



	// iter_j == (ENTRY % li)
	// iter_k = ((ENTRY / li) % mi)
	ITER_GRID3_FRONT({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0] + iter_j * h[1], time);
		}, L, M);

	// edges

	ITER_GRID3_FRONT_TOP({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0], time);
		}, L, M);

	ITER_GRID3_FRONT_BOTTOM({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[1], time);
		}, L, M);

	ITER_GRID3_FRONT_LEFT({
		bd->update(grid[INDEX], x[0], y[0] + iter_j * h[1], time);
		}, L, M);

	ITER_GRID3_FRONT_RIGHT({
		bd->update(grid[INDEX], x[1], y[0] + iter_j * h[1], time);
		}, L, M);

	// corners

	ITER_GRID3_LEFT_TOP_FRONT({
		bd->update(grid[INDEX], x[0], y[0], time);
		}, L, M);

	ITER_GRID3_FRONT_TOP_RIGHT({
		bd->update(grid[INDEX], x[1], y[0], time);
		}, L, M);

	ITER_GRID3_LEFT_BOTTOM_FRONT({
		bd->update(grid[INDEX], x[0], y[1], time);
		}, L, M);

	ITER_GRID3_FRONT_BOTTOM_RIGHT({
		bd->update(grid[INDEX], x[1], y[1], time);
		}, L, M);


}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::BACK, 2>::operator()(const grid::Boundary<T, 2>* b, Grid<T, 3>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 2, BoundaryType::DEFAULT> const*>(b);

	const double* x = bd->v;
	const double* y = bd->v + 2;
	double h[2];

	// go backwards or forwards in iteration depending on the interval
	int fx = (x[0] < x[1]) ? 1 : -1;
	int fy = (y[0] < y[1]) ? 1 : -1;
	h[0] = bd->h[0] * fx;
	h[1] = bd->h[1] * fy;

	iter_type L = grid.dims[0];
	iter_type M = grid.dims[1];
	iter_type N = grid.dims[2];

	// iter_j == (ENTRY % li)
	// iter_k = ((ENTRY / li) % mi)
	ITER_GRID3_BACK({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0] + iter_j * h[1], time);
		}, L, M, N);

	// edges

	ITER_GRID3_BACK_TOP({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0], time);
		}, L, M, N);

	ITER_GRID3_BACK_BOTTOM({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[1], time);
		}, L, M, N);

	ITER_GRID3_BACK_LEFT({
		bd->update(grid[INDEX], x[0], y[0] + iter_j * h[1], time);
		}, L, M, N);

	ITER_GRID3_BACK_RIGHT({
		bd->update(grid[INDEX], x[1], y[0] + iter_j * h[1], time);
		}, L, M, N);

	// corners

	ITER_GRID3_BACK_TOP_LEFT({
		bd->update(grid[INDEX], x[0], y[0], time);
		}, L, M, N);

	ITER_GRID3_RIGHT_TOP_BACK({
		bd->update(grid[INDEX], x[1], y[0], time);
		}, L, M, N);

	ITER_GRID3_BACK_BOTTOM_LEFT({
		bd->update(grid[INDEX], x[0], y[1], time);
		}, L, M, N);

	ITER_GRID3_RIGHT_BOTTOM_BACK({
		bd->update(grid[INDEX], x[1], y[1], time);
		}, L, M, N);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 2>::operator()(const grid::Boundary<T, 2>* b, Grid<T, 3>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 2, BoundaryType::DEFAULT> const*>(b);

	const double* x = bd->v;
	const double* y = bd->v + 2;
	double h[2];

	// go backwards or forwards in iteration depending on the interval
	int fx = (x[0] < x[1]) ? 1 : -1;
	int fy = (y[0] < y[1]) ? 1 : -1;
	h[0] = bd->h[0] * fx;
	h[1] = bd->h[1] * fy;

	iter_type L = grid.dims[0];
	iter_type M = grid.dims[1];
	iter_type N = grid.dims[2];

	// iter_i == (ENTRY % ni)
	// iter_j == ((ENTRY / ni) % mi)
	ITER_GRID3_LEFT({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0] + iter_j * h[1], time);
		}, L, M, N);

	// edges

	ITER_GRID3_LEFT_TOP({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0], time);
		}, L, M, N);

	ITER_GRID3_LEFT_BOTTOM({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[1], time);
		}, L, M, N);

	ITER_GRID3_FRONT_LEFT({
		bd->update(grid[INDEX], x[0], y[0] + iter_j * h[1], time);
		}, L, M);

	ITER_GRID3_BACK_LEFT({
		bd->update(grid[INDEX], x[1], y[0] + iter_j * h[1], time);
		}, L, M, N);

	// corners

	ITER_GRID3_LEFT_TOP_FRONT({
		bd->update(grid[INDEX], x[0], y[0], time);
		}, L, M);

	ITER_GRID3_BACK_TOP_LEFT({
		bd->update(grid[INDEX], x[1], y[0], time);
		}, L, M, N);

	ITER_GRID3_LEFT_BOTTOM_FRONT({
		bd->update(grid[INDEX], x[0], y[1], time);
		}, L, M);

	ITER_GRID3_BACK_BOTTOM_LEFT({
		bd->update(grid[INDEX], x[1], y[1], time);
		}, L, M, N);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 2>::operator()(const grid::Boundary<T, 2>* b, Grid<T, 3>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 2, BoundaryType::DEFAULT> const*>(b);

	const double* x = bd->v;
	const double* y = bd->v + 2;
	double h[2];

	// go backwards or forwards in iteration depending on the interval
	int fx = (x[0] < x[1]) ? 1 : -1;
	int fy = (y[0] < y[1]) ? 1 : -1;
	h[0] = bd->h[0] * fx;
	h[1] = bd->h[1] * fy;

	iter_type L = grid.dims[0];
	iter_type M = grid.dims[1];
	iter_type N = grid.dims[2];

	// iter_i == (ENTRY % ni)
	// iter_j == ((ENTRY / ni) % mi)
	ITER_GRID3_RIGHT({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0] + iter_j * h[1], time);
		}, L, M, N);


	// edges

	ITER_GRID3_RIGHT_TOP({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0], time);
		}, L, M, N);

	ITER_GRID3_RIGHT_BOTTOM({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[1], time);
		}, L, M, N);

	ITER_GRID3_FRONT_RIGHT({
		bd->update(grid[INDEX], x[0], y[0] + iter_j * h[1], time);
		}, L, M);

	ITER_GRID3_BACK_RIGHT({
		bd->update(grid[INDEX], x[1], y[0] + iter_j * h[1], time);
		}, L, M, N);


	// corners

	ITER_GRID3_FRONT_TOP_RIGHT({
		bd->update(grid[INDEX], x[0], y[0], time);
		}, L, M);

	ITER_GRID3_RIGHT_TOP_BACK({
		bd->update(grid[INDEX], x[1], y[0], time);
		}, L, M, N);

	ITER_GRID3_FRONT_BOTTOM_RIGHT({
		bd->update(grid[INDEX], x[0], y[1], time);
		}, L, M);

	ITER_GRID3_RIGHT_BOTTOM_BACK({
		bd->update(grid[INDEX], x[1], y[1], time);
		}, L, M, N);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::TOP, 2>::operator()(const grid::Boundary<T, 2>* b, Grid<T, 3>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 2, BoundaryType::DEFAULT> const*>(b);

	const double* x = bd->v;
	const double* y = bd->v + 2;
	double h[2];

	// go backwards or forwards in iteration depending on the interval
	int fx = (x[0] < x[1]) ? 1 : -1;
	int fy = (y[0] < y[1]) ? 1 : -1;
	h[0] = bd->h[0] * fx;
	h[1] = bd->h[1] * fy;

	iter_type L = grid.dims[0];
	iter_type M = grid.dims[1];
	iter_type N = grid.dims[2];

	// iter_k == (ENTRY % ni)
	// iter_i == ((ENTRY / li) % ni)
	ITER_GRID3_TOP({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0] + iter_j * h[1], time);
		}, L, M, N);


	// edges

	ITER_GRID3_FRONT_TOP({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0], time);
		}, L, M);

	ITER_GRID3_BACK_TOP({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[1], time);
		}, L, M, N);

	ITER_GRID3_LEFT_TOP({
		bd->update(grid[INDEX], x[0], y[0] + iter_j * h[1], time);
		}, L, M, N);

	ITER_GRID3_RIGHT_TOP({
		bd->update(grid[INDEX], x[1], y[0] + iter_j * h[1], time);
		}, L, M, N);

	// corners

	ITER_GRID3_LEFT_TOP_FRONT({
		bd->update(grid[INDEX], x[0], y[0], time);
		}, L, M);

	ITER_GRID3_FRONT_TOP_RIGHT({
		bd->update(grid[INDEX], x[1], y[0], time);
		}, L, M);

	ITER_GRID3_BACK_TOP_LEFT({
		bd->update(grid[INDEX], x[0], y[1], time);
		}, L, M, N);

	ITER_GRID3_RIGHT_TOP_BACK({
		bd->update(grid[INDEX], x[1], y[1], time);
		}, L, M, N);

}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::BOTTOM, 2>::operator()(const grid::Boundary<T, 2>* b, Grid<T, 3>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 2, BoundaryType::DEFAULT> const*>(b);

	const double* x = bd->v;
	const double* y = bd->v + 2;
	double h[2];

	// go backwards or forwards in iteration depending on the interval
	int fx = (x[0] < x[1]) ? 1 : -1;
	int fy = (y[0] < y[1]) ? 1 : -1;
	h[0] = bd->h[0] * fx;
	h[1] = bd->h[1] * fy;

	iter_type L = grid.dims[0];
	iter_type M = grid.dims[1];
	iter_type N = grid.dims[2];

	// iter_k == (ENTRY % ni)
	// iter_i == ((ENTRY / li) % ni)
	ITER_GRID3_BOTTOM({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0] + iter_j * h[1], time);
		}, L, M, N);

	// edges

	ITER_GRID3_FRONT_BOTTOM({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[0], time);
		}, L, M);

	ITER_GRID3_BACK_BOTTOM({
		bd->update(grid[INDEX], x[0] + iter_i * h[0], y[1], time);
		}, L, M, N);

	ITER_GRID3_LEFT_BOTTOM({
		bd->update(grid[INDEX], x[0], y[0] + iter_j * h[1], time);
		}, L, M, N);

	ITER_GRID3_RIGHT_BOTTOM({
		bd->update(grid[INDEX], x[1], y[0] + iter_j * h[1], time);
		}, L, M, N);

	// corners

	ITER_GRID3_LEFT_BOTTOM_FRONT({
		bd->update(grid[INDEX], x[0], y[0], time);
		}, L, M);

	ITER_GRID3_FRONT_BOTTOM_RIGHT({
		bd->update(grid[INDEX], x[1], y[0], time);
		}, L, M);

	ITER_GRID3_BACK_BOTTOM_LEFT({
		bd->update(grid[INDEX], x[0], y[1], time);
		}, L, M, N);

	ITER_GRID3_RIGHT_BOTTOM_BACK({
		bd->update(grid[INDEX], x[1], y[1], time);
		}, L, M, N);
}

// 2 dimension

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 1>::operator()(const grid::Boundary<T, 1>* b, Grid<T, 2>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 1, BoundaryType::DEFAULT> const*>(b);

	const double v0 = bd->v[0];
	const double v1 = bd->v[1];
	double h;

	// go backwards or forwards in iteration depending on the interval
	int fx = (v0 < v1) ? 1 : -1;
	h = bd->h * fx;

	iter_type L = grid.dims[0];
	iter_type M = grid.dims[1];

	ITER_GRID2_LEFT({
		bd->update(grid[INDEX], v0 + iter_i * h, 0, time);
		}, L, M);

	// corners

	ITER_GRID2_LEFT_TOP({
		bd->update(grid[INDEX], 0, v0, time);
		}, L);
	ITER_GRID2_BOTTOM_LEFT({
		bd->update(grid[INDEX], 0, v1, time);
		}, L, M);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 1>::operator()(const grid::Boundary<T, 1>* b, Grid<T, 2>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 1, BoundaryType::DEFAULT> const*>(b);

	const double v0 = bd->v[0];
	const double v1 = bd->v[1];
	double h;

	// go backwards or forwards in iteration depending on the interval
	int fx = (v0 < v1) ? 1 : -1;
	h = bd->h * fx;

	iter_type L = grid.dims[0];
	iter_type M = grid.dims[1];

	ITER_GRID2_RIGHT({
		bd->update(grid[INDEX], v0 + iter_i * h, 0, time);
		}, L, M);

	// corners

	ITER_GRID2_TOP_RIGHT({
		bd->update(grid[INDEX], 0, v0, time);
		}, L);
	ITER_GRID2_RIGHT_BOTTOM({
		bd->update(grid[INDEX], 0, v1, time);
		}, L, M);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::TOP, 1>::operator()(const grid::Boundary<T, 1>* b, Grid<T, 2>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 1, BoundaryType::DEFAULT> const*>(b);

	const double v0 = bd->v[0];
	const double v1 = bd->v[1];
	double h;

	// go backwards or forwards in iteration depending on the interval
	int fx = (v0 < v1) ? 1 : -1;
	h = bd->h * fx;

	iter_type L = grid.dims[0];

	ITER_GRID2_TOP({
		bd->update(grid[INDEX], v0 + iter_i * h, 0, time);
		}, L);

	// corners

	ITER_GRID2_LEFT_TOP({
		bd->update(grid[INDEX], v0, 0, time);
		}, L);
	ITER_GRID2_TOP_RIGHT({
		bd->update(grid[INDEX], v1, 0, time);
		}, L);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::BOTTOM, 1>::operator()(const grid::Boundary<T, 1>* b, Grid<T, 2>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 1, BoundaryType::DEFAULT> const*>(b);

	const double v0 = bd->v[0];
	const double v1 = bd->v[1];
	double h;

	// go backwards or forwards in iteration depending on the interval
	int fx = (v0 < v1) ? 1 : -1;
	h = bd->h * fx;

	iter_type L = grid.dims[0];
	iter_type M = grid.dims[1];

	ITER_GRID2_BOTTOM({
		bd->update(grid[INDEX], v0 + iter_i * h, 0, time);
		}, L, M);

	// corners

	ITER_GRID2_BOTTOM_LEFT({
		bd->update(grid[INDEX], v0, 0, time);
		}, L, M);
	ITER_GRID2_RIGHT_BOTTOM({
		bd->update(grid[INDEX], v1, 0, time);
		}, L, M);
}

// 1 dimension

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 0>::operator()(const grid::Boundary<T, 0>* b, Grid<T, 1>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 0, BoundaryType::DEFAULT> const*>(b);
	const double v = bd->v;

	ITER_GRID1_LEFT({
		bd->update(grid[INDEX], v, 0, time);
		});
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 0>::operator()(const grid::Boundary<T, 0>* b, Grid<T, 1>& grid, double time)
{
	auto* bd = static_cast<grid::BoundaryApplied<T, 0, BoundaryType::DEFAULT> const*>(b);
	const double v = bd->v;
	iter_type L = grid.dims[0];

	ITER_GRID1_RIGHT({
		bd->update(grid[INDEX], v, 0, time);
		}, L);
}


template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::FRONT, 2>::operator()(const grid::Boundary<T, 2>* b, RegionalGrid<T, 3>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::FRONT, Side::FRONT>{}, b, grid, time);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::BACK, 2>::operator()(const grid::Boundary<T, 2>* b, RegionalGrid<T, 3>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::BACK, Side::BACK>{}, b, grid, time);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 2>::operator()(const grid::Boundary<T, 2>* b, RegionalGrid<T, 3>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, b, grid, time);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 2>::operator()(const grid::Boundary<T, 2>* b, RegionalGrid<T, 3>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, b, grid, time);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::TOP, 2>::operator()(const grid::Boundary<T, 2>* b, RegionalGrid<T, 3>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::TOP, Side::TOP>{}, b, grid, time);
	regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::BACK>{}, b, grid, time);
	regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::FRONT>{}, b, grid, time);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::BOTTOM, 2>::operator()(const grid::Boundary<T, 2>* b, RegionalGrid<T, 3>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM, Side::BOTTOM>{}, b, grid, time);
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM, Side::BACK>{}, b, grid, time);
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM, Side::FRONT>{}, b, grid, time);
}

// 2 dimension

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 1>::operator()(const grid::Boundary<T, 1>* b, RegionalGrid<T, 2>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT>{}, b, grid, time);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 1>::operator()(const grid::Boundary<T, 1>* b, RegionalGrid<T, 2>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT>{}, b, grid, time);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::TOP, 1>::operator()(const grid::Boundary<T, 1>* b, RegionalGrid<T, 2>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::TOP>{}, b, grid, time);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::BOTTOM, 1>::operator()(const grid::Boundary<T, 1>* b, RegionalGrid<T, 2>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM>{}, b, grid, time);
}

// 1 dimension

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 0>::operator()(const grid::Boundary<T, 0>* b, RegionalGrid<T, 1>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, b, grid, time);
}

template<>
template<typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 0>::operator()(const grid::Boundary<T, 0>* b, RegionalGrid<T, 1>& grid, double time)
{
	regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, b, grid, time);
}

