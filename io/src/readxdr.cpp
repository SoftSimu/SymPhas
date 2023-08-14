
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

		for (iter_type i = 0; i < dim; ++i)
		{
			double v[2];
			if (xdrfile_read_double(v, 2, f) != 2)
			{
				if (index)
				{
					*index = BAD_INDEX;
				}
				return { nullptr, dim };
			}
			ginfo[symphas::index_to_axis(i)].set_count(v[0], v[1], dims[i]);
		}

		delete[] dims;
	}

	/* if nothing was read, then -1 should be returned
		*/

	int get = BAD_INDEX;
	bool index_scanned = (xdrfile_read_int(&get, 1, f) == 1);
	if (index)
	{
		*index = (index_scanned) ? get : BAD_INDEX;
	}

	for (auto const& [axis, interval] : ginfo)
	{
		double data[2]{};
		xdrfile_write_double(data, 2, f);
		ginfo.intervals[axis].set_interval(data[0], data[1]);
	}

	return ginfo;
}


void adjust_data_format(symphas::grid_info const& ginfo, scalar_t* from, scalar_t* to)
{
	auto helper = symphas::io::gp::new_helper(ginfo);
	iter_type n = 0;
	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				to[ii] = from[n++];
			}
		}
	}
}


/*
 * functions to write data from the grid to a file
 */
template<>
void symphas::io::xdr::read_block(scalar_t* grid, symphas::grid_info ginfo, XDRFILE* f)
{
	double* raw_data = new double[ginfo.num_interval_points()] {};
	xdrfile_read_double(raw_data, ginfo.num_interval_points(), f);
	adjust_data_format(ginfo, raw_data, grid);
	delete[] raw_data;
}

template<>
void symphas::io::xdr::read_block(complex_t* grid, symphas::grid_info ginfo, XDRFILE* f)
{
	xdrfile_read_double(reinterpret_cast<double*>(grid), ginfo.num_interval_points() * 2, f);
}

template<>
void symphas::io::xdr::read_block(double_arr2 *grid, symphas::grid_info ginfo, XDRFILE* f)
{
	xdrfile_read_double(reinterpret_cast<double*>(grid), ginfo.num_interval_points() * 2, f);
}

template<>
void symphas::io::xdr::read_block(vector_t<3>* grid, symphas::grid_info ginfo, XDRFILE* f)
{
	xdrfile_read_double(reinterpret_cast<double*>(grid), ginfo.num_interval_points() * 3, f);
}

template<>
void symphas::io::xdr::read_block(vector_t<2>* grid, symphas::grid_info ginfo, XDRFILE* f)
{
	xdrfile_read_double(reinterpret_cast<double*>(grid), ginfo.num_interval_points() * 2, f);
}

template<>
void symphas::io::xdr::read_block(vector_t<1>* grid, symphas::grid_info ginfo, XDRFILE* f)
{
	xdrfile_read_double(reinterpret_cast<double*>(grid), ginfo.num_interval_points() * 1, f);
}



template<>
void symphas::io::xdr::read_block(scalar_ptr_t(&grid)[3], symphas::grid_info ginfo, XDRFILE* f)
{
	xdrfile_read_double(reinterpret_cast<double*>(grid), ginfo.num_interval_points() * 3, f);
}

template<>
void symphas::io::xdr::read_block(scalar_ptr_t(&grid)[2], symphas::grid_info ginfo, XDRFILE* f)
{
	xdrfile_read_double(reinterpret_cast<double*>(grid), ginfo.num_interval_points() * 2, f);
}

template<>
void symphas::io::xdr::read_block(scalar_ptr_t(&grid)[1], symphas::grid_info ginfo, XDRFILE* f)
{
	xdrfile_read_double(reinterpret_cast<double*>(grid), ginfo.num_interval_points() * 1, f);
}
