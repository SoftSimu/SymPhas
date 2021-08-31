
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
 * PURPOSE: Defines functions for using grids, such as scaling grid values.
 *
 * ***************************************************************************
 */

#pragma once

#include "grid.h"
#include "dft.h"


namespace symphas::dft
{
	//! Specialization based on symphas::dft::scale() using ::Grid type.
	/*!
	 * Scales the values of the grid type based on the number of elements
	 * in the grid.
	 * 
	 * \param grid The grid containing the values to be scaled.
	 */
	template<typename T>
	void scale(Grid<T, 1>& grid)
	{
		scale(grid.values, grid.dims[0]);
	}

	//! Specialization based on symphas::dft::scale()
	/*!
	 * Scales the values of the grid type based on the number of elements
	 * in the grid.
	 *
	 * \param grid The grid containing the values to be scaled.
	 */
	template<typename T>
	void scale(Grid<T, 2>& grid)
	{
		scale(grid.values, grid.dims[0], grid.dims[1]);
	}

	//! Specialization based on symphas::dft::scale()
	/*!
	 * Scales the values of the grid type based on the number of elements
	 * in the grid.
	 *
	 * \param grid The grid containing the values to be scaled.
	 */
	template<typename T>
	void scale(Grid<T, 3>& grid)
	{
		scale(grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
	}
}


/*! 
 * \addtogroup grid
 * @{
 */


namespace grid
{
	//! Creates a sub-grid centered around the given index.
	/*!
	 * Given a grid, the get_subgrid function will copy a small part of the Grid "src" into
	 * another grid and then return that grid.
	 *
	 * The part that is copied is centered on the raw index n, and the given len_type
	 * array "extent" indicates how far in each axis to copy around the given point.
	 *
	 * The copy algorithm uses periodic boundaries, so sub-grids where the source point
	 * is near the edges will be continued onto the opposite side.
	 *
	 * One restriction is that the individual extent dimensions can't be larger than the
	 * grid dimensions.
	 * 
	 * \param src The grid from which is taken the subset (subgrid).
	 * \param extent The lengths beyond the center point from which to take
	 * the subgrid.
	 * \param n The index at which the subgrid is centered.
	 */
	template<template<typename, size_t> typename GridType, typename T>
	Grid<T, 1> get_subgrid(GridType<T, 1> const& src, len_type* extent, iter_type n)
	{
		len_type ds = extent[0] + extent[0] + 1;
		len_type s2 = extent[0];
		len_type cdims[] = { ds };

		Grid<T, 1> grid(cdims);
		T* dest = grid.values;

		if (n < s2)
		{
			iter_type offset = src.dims[0] - (s2 - n);
			for (iter_type i = 0; i < s2 - n; ++i)
			{
				dest[i] = src[offset + i];
			}

			for (iter_type i = 0; i < ds - (s2 - n); ++i)
			{
				dest[i + s2 - n] = src[i];
			}
		}
		else if (n > src.dims[0] - s2)
		{
			iter_type offset = n - s2;
			for (iter_type i = 0; i < src.dims[0] - (n - s2); ++i)
			{
				dest[i] = src[offset + i];
			}
			offset = src.dims[0] - n + s2;
			for (iter_type i = 0; i < ds - offset; ++i)
			{
				dest[offset + i] = src[i];
			}
		}
		else
		{
			iter_type index = n - s2;
			for (iter_type i = 0; i < ds; ++i)
			{
				dest[i] = src[index + i];
			}
		}

		return grid;
	}





	//! Specialization based on grid::get_subgrid().
	template<template<typename, size_t> typename GridType, typename T>
	Grid<T, 2> get_subgrid(GridType<T, 2> const& src, len_type* extent, iter_type n)
	{
		len_type dsx = extent[0] + extent[0] + 1;
		len_type dsy = extent[1] + extent[1] + 1;
		len_type s2x = extent[0];
		len_type s2y = extent[1];
		len_type cdims[] = { dsx, dsy };

		Grid<T, 2> grid(cdims);
		T* dest = grid.values;

		iter_type x = n % src.dims[0] - s2x;
		iter_type y = (n / src.dims[0]) - s2y;

		iter_type index = 0;
		for (iter_type j = 0; j < dsy; ++j)
		{
			y %= src.dims[1];
			for (iter_type i = 0; i < dsx; ++i)
			{
				x %= src.dims[0];
				dest[index++] = src[x + y * src.dims[0]];
				++x;
			}
			++y;
		}



		return grid;
	}




	//! Specialization based on grid::get_subgrid().
	template<template<typename, size_t> typename GridType, typename T>
	Grid<T, 3> get_subgrid(GridType<T, 3> const& src, len_type* extent, iter_type n)
	{

		len_type dsx = extent[0] + extent[0] + 1;
		len_type dsy = extent[1] + extent[1] + 1;
		len_type dsz = extent[2] + extent[2] + 1;
		len_type s2x = extent[0];
		len_type s2y = extent[1];
		len_type s2z = extent[2];
		len_type cdims[] = { dsx, dsy, dsz };

		Grid<T, 3> grid(cdims);
		T* dest = grid.values;

		iter_type x = n % src.dims[0] - s2x;
		iter_type y = (n / src.dims[0]) % src.dims[1] - s2y;
		iter_type z = n / (src.dims[0] * src.dims[1]) - s2z;
		iter_type index = 0;
		for (iter_type k = 0; k < dsz; ++k)
		{
			z %= src.dims[2];
			for (iter_type j = 0; j < dsy; ++j)
			{
				y %= src.dims[1];
				for (iter_type i = 0; i < dsx; ++i)
				{
					x %= src.dims[0];
					dest[index++] = src[x + y * src.dims[0] + z * src.dims[0] * src.dims[1]];
					++x;
				}
				++y;
			}
			++z;
		}


		return grid;
	}

}

//! @}


