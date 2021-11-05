
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


#include "writexdr.h"
#include <xdrfile.h>


XDRFILE* open_xdrgridf(symphas::io::write_info winfo, bool is_checkpoint)
{
	char name[BUFFER_LENGTH];
	symphas::io::copy_data_file_name(winfo.dir_str_ptr, winfo.index, winfo.id, (is_checkpoint ? DataFileType::CHECKPOINT_DATA : winfo.type), name);

	char mode[2];
	if (params::single_output_file)
	{
		mode[0] = 'a';
	}
	else
	{
		mode[0] = 'w';
	}
	mode[1] = '\0';
	XDRFILE* f = xdrfile_open(name, mode);
	return f;
}

void print_xdr_header(int index, size_t id, symphas::grid_info const& ginfo, symphas::io::write_info const& winfo, XDRFILE* f)
{
	static std::vector<std::tuple<std::string, size_t>> idlist;
	if (!params::single_output_file || (std::find(idlist.begin(), idlist.end(), std::make_tuple(winfo.dir_str_ptr, id)) == idlist.end()))
	{

		int dim = static_cast<int>(ginfo.dimension());
		xdrfile_write_int(&dim, 1, f);

		int* dims = new int[dim];
		for (iter_type i = 0; i < dim; ++i)
		{
			dims[i] = ginfo.at(symphas::index_to_axis(i)).count();
		}
		xdrfile_write_int(dims, dim, f);
		delete[] dims;

		/* the output of the intervals always goes x, y, z
		 * keep checking size of dimension and print corresponding interval
		 */

		for (iter_type i = 0; i < dim; ++i)
		{
			double v[2] = { 
				ginfo.intervals.at(symphas::index_to_axis(i)).left(),
				ginfo.intervals.at(symphas::index_to_axis(i)).right() };
			xdrfile_write_double(v, 2, f);
		}


		if (std::find(idlist.begin(), idlist.end(), std::make_tuple(winfo.dir_str_ptr, id)) == idlist.end())
		{
			idlist.emplace_back(winfo.dir_str_ptr, id);
		}
	}
	xdrfile_write_int(&index, 1, f);
}



/*
 * functions to write data from the grid to a file
 */
void save_xdr(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(grid), ginfo.num_points(), f);
	xdrfile_close(f);
}

void save_xdr(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	const double* data = reinterpret_cast<const double*>(grid);
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(data), ginfo.num_points() * 2, f);
	xdrfile_close(f);
}

void save_xdr(const double(*grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	const double* data = reinterpret_cast<const double*>(grid);
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(data), ginfo.num_points() * 2, f);
	xdrfile_close(f);
}

void save_xdr(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);

	for (iter_type k = 0, ii = 0; k < ginfo.at(Axis::Z).count(); ++k)
	{
		for (iter_type j = 0; j < ginfo.at(Axis::Y).count(); ++j)
		{
			for (iter_type i = 0; i < ginfo.at(Axis::X).count(); ++i, ++ii)
			{
				int d[3];
				d[0] = static_cast<int>((grid[ii] * XDR_COORD_COMPRESSION).v[0]);
				d[1] = static_cast<int>((grid[ii] * XDR_COORD_COMPRESSION).v[1]);
				d[2] = static_cast<int>((grid[ii] * XDR_COORD_COMPRESSION).v[2]);

				xdrfile_write_int(d, 3, f);
			}
		}
	}
	xdrfile_close(f);
}

void save_xdr(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);

	for (iter_type j = 0, ii = 0; j < ginfo.at(Axis::Y).count(); ++j)
	{
		for (iter_type i = 0; i < ginfo.at(Axis::X).count(); ++i, ++ii)
		{
			int d[2];
			d[0] = static_cast<int>((grid[ii] * XDR_COORD_COMPRESSION).v[0]);
			d[1] = static_cast<int>((grid[ii] * XDR_COORD_COMPRESSION).v[1]);

			xdrfile_write_int(d, 2, f);
		}
	}
	xdrfile_close(f);
}

void save_xdr(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);

	for (iter_type i = 0; i < ginfo.at(Axis::X).count(); ++i)
	{
		int d = static_cast<int>((grid[i] * XDR_COORD_COMPRESSION).v[0]);
		xdrfile_write_int(&d, 1, f);
	}
	xdrfile_close(f);
}




/*
 * functions to checkpoint the data to a file
 */
template<>
void symphas::io::xdr::save_grid_plotting<scalar_t>(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

template<>
void symphas::io::xdr::save_grid_plotting<complex_t>(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

template<>
void symphas::io::xdr::save_grid_plotting<double[2]>(const double(*grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

template<>
void symphas::io::xdr::save_grid_plotting<vector_t<3>>(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

template<>
void symphas::io::xdr::save_grid_plotting<vector_t<2>>(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

template<>
void symphas::io::xdr::save_grid_plotting<vector_t<1>>(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}





/*
 * functions to checkpoint the data to a file
 */
template<>
void symphas::io::xdr::save_grid<scalar_t>(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

template<>
void symphas::io::xdr::save_grid<complex_t>(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

template<>
void symphas::io::xdr::save_grid<double[2]>(const double(*grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

template<>
void symphas::io::xdr::save_grid<vector_t<3>>(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

template<>
void symphas::io::xdr::save_grid<vector_t<2>>(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

template<>
void symphas::io::xdr::save_grid<vector_t<1>>(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}




