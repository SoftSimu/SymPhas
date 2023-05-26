
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
 */


#include "readxdr.h"
#include <xdrfile.h>


XDRFILE* symphas::io::xdr::open_xdrgridf(const char* name)
{
	XDRFILE* f = xdrfile_open(name, "r");
	return f;
}

void symphas::io::xdr::xdrfile_close(XDRFILE* f)
{
	::xdrfile_close(f);
}

template<>
symphas::grid_info symphas::io::xdr::read_header(XDRFILE* f, int* index)
{
	symphas::grid_info ginfo{ nullptr, 0 };

	iter_type read_index = (index) ? *index : BAD_INDEX;
	if (read_index < 0)
	{
		/* test the first input to see if there are characters
			* if there nothing is read, then return index -1
			*/
		int dim;
		if (xdrfile_read_int(&dim, 1, f) != 1)
		{
			if (index)
			{
				*index = BAD_INDEX;
			}
			return { nullptr, dim };
		}

		int* dims = new int[dim];
		if (xdrfile_read_int(dims, dim, f) != dim)
		{
			if (index)
			{
				*index = BAD_INDEX;
			}
			return { nullptr, dim };
		}

		ginfo = symphas::grid_info{ dims, dim };

		/* the output of the intervals always goes x, y, z
			* keep checking size of dimension and print corresponding interval
			*/

		double v[2];
		for (iter_type i = 0; i < dim; ++i)
		{
			if (xdrfile_read_double(v, 2, f) != 2)
			{
				if (index)
				{
					*index = BAD_INDEX;
				}
				return { nullptr, dim };
			}
			ginfo.at(symphas::index_to_axis(i)).set_count(v[0], v[1], dims[i]);
		}

		delete[] dims;
	}

	/* if nothing was read, then -1 should be returned
		*/
	if (index)
	{
		if (xdrfile_read_int(index, 1, f) != 1)
		{
			*index = BAD_INDEX;
		}
	}
	return ginfo;
}



/*
 * functions to write data from the grid to a file
 */
template<>
void symphas::io::xdr::read_block(scalar_t* grid, symphas::grid_info ginfo, XDRFILE* f)
{
	xdrfile_read_double(grid, ginfo.num_points(), f);
}

template<>
void symphas::io::xdr::read_block(complex_t* grid, symphas::grid_info ginfo, XDRFILE* f)
{
	xdrfile_read_double(reinterpret_cast<double*>(grid), ginfo.num_points() * 2, f);
}

template<>
void symphas::io::xdr::read_block(double_arr2 *grid, symphas::grid_info ginfo, XDRFILE* f)
{
	xdrfile_read_double(reinterpret_cast<double*>(grid), ginfo.num_points() * 2, f);
}

template<>
void symphas::io::xdr::read_block(vector_t<3>* grid, symphas::grid_info ginfo, XDRFILE* f)
{
	for (iter_type k = 0, ii = 0; k < ginfo.at(Axis::Z).get_count(); ++k)
	{
		for (iter_type j = 0; j < ginfo.at(Axis::Y).get_count(); ++j)
		{
			for (iter_type i = 0; i < ginfo.at(Axis::X).get_count(); ++i, ++ii)
			{
				int d[3];
				xdrfile_read_int(d, 3, f);
				grid[ii] = vector_t<3>{
					d[0] / XDR_COORD_COMPRESSION,
					d[1] / XDR_COORD_COMPRESSION,
					d[2] / XDR_COORD_COMPRESSION };
			}
		}
	}
}

template<>
void symphas::io::xdr::read_block(vector_t<2>* grid, symphas::grid_info ginfo, XDRFILE* f)
{
	for (iter_type j = 0, ii = 0; j < ginfo.at(Axis::Y).get_count(); ++j)
	{
		for (iter_type i = 0; i < ginfo.at(Axis::X).get_count(); ++i, ++ii)
		{
			int d[2];
			xdrfile_read_int(d, 2, f);
			grid[ii] = vector_t<2>{
				d[0] / XDR_COORD_COMPRESSION,
				d[1] / XDR_COORD_COMPRESSION };
		}
	}
}

template<>
void symphas::io::xdr::read_block(vector_t<1>* grid, symphas::grid_info ginfo, XDRFILE* f)
{
	for (iter_type i = 0; i < ginfo.at(Axis::X).get_count(); ++i)
	{
		int d;
		xdrfile_write_int(&d, 1, f);
		grid[i] = vector_t<1>{
			d / XDR_COORD_COMPRESSION };
	}
}



template<>
void symphas::io::xdr::read_block(scalar_ptr_t(&grid)[3], symphas::grid_info ginfo, XDRFILE* f)
{
	for (iter_type k = 0, ii = 0; k < ginfo.at(Axis::Z).get_count(); ++k)
	{
		for (iter_type j = 0; j < ginfo.at(Axis::Y).get_count(); ++j)
		{
			for (iter_type i = 0; i < ginfo.at(Axis::X).get_count(); ++i, ++ii)
			{
				int d[3];
				xdrfile_read_int(d, 3, f);
				grid[0][ii] = d[0] / XDR_COORD_COMPRESSION;
				grid[1][ii] = d[1] / XDR_COORD_COMPRESSION;
				grid[2][ii] = d[2] / XDR_COORD_COMPRESSION;
			}
		}
	}
}

template<>
void symphas::io::xdr::read_block(scalar_ptr_t(&grid)[2], symphas::grid_info ginfo, XDRFILE* f)
{
	for (iter_type j = 0, ii = 0; j < ginfo.at(Axis::Y).get_count(); ++j)
	{
		for (iter_type i = 0; i < ginfo.at(Axis::X).get_count(); ++i, ++ii)
		{
			int d[2];
			xdrfile_read_int(d, 2, f);
			grid[0][ii] = d[0] / XDR_COORD_COMPRESSION;
			grid[1][ii] = d[1] / XDR_COORD_COMPRESSION;
		}
	}
}

template<>
void symphas::io::xdr::read_block(scalar_ptr_t(&grid)[1], symphas::grid_info ginfo, XDRFILE* f)
{
	for (iter_type i = 0; i < ginfo.at(Axis::X).get_count(); ++i)
	{
		int d;
		xdrfile_write_int(&d, 1, f);
		grid[0][i] = d / XDR_COORD_COMPRESSION;
	}
}
