
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
	symphas::io::copy_data_file_name(winfo.dir, winfo.index, winfo.id, (is_checkpoint ? DataFileType::CHECKPOINT_DATA : winfo.type), name);

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
	
	int dim = static_cast<int>(ginfo.dimension());
	if (!params::single_output_file || (std::find(idlist.begin(), idlist.end(), std::make_tuple(winfo.dir, id)) == idlist.end()))
	{
		xdrfile_write_int(&dim, 1, f);

		int* dims = new int[dim];
		for (iter_type i = 0; i < dim; ++i)
		{
			dims[i] = ginfo.at(symphas::index_to_axis(i)).get_count();
		}
		xdrfile_write_int(dims, dim, f);
		delete[] dims;


		/* the output of the intervals always goes x, y, z
		 * keep printing corresponding interval
		 */
		for (const auto& [axis, interval] : ginfo)
		{
			double domain[2]{ interval.domain_left(), interval.domain_right() };
			xdrfile_write_double(&domain[0], 2, f);
		}

		if (std::find(idlist.begin(), idlist.end(), std::make_tuple(winfo.dir, id)) == idlist.end())
		{
			idlist.emplace_back(winfo.dir, id);
		}
	}
	xdrfile_write_int(&index, 1, f);
	xdrfile_write_int(&dim, 1, f);

	for (auto const& [axis, interval] : ginfo)
	{
		double data[2]{ interval.left(), interval.right() };
		xdrfile_write_double(&data[0], 2, f);
	}
}



/*
 * functions to write data from the grid to a file
 */
void save_xdr(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(grid), ginfo.num_interval_points(), f);
	xdrfile_close(f);
}

void save_xdr(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	const double* data = reinterpret_cast<const double*>(grid);
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(data), ginfo.num_interval_points() * 2, f);
	xdrfile_close(f);
}

void save_xdr(const double(*grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	const double* data = reinterpret_cast<const double*>(grid);
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(data), ginfo.num_interval_points() * 2, f);
	xdrfile_close(f);
}

void save_xdr(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	const double* data = reinterpret_cast<const double*>(grid);
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(data), ginfo.num_interval_points() * 3, f);
	xdrfile_close(f);
}

void save_xdr(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	const double* data = reinterpret_cast<const double*>(grid);
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(data), ginfo.num_interval_points() * 2, f);
	xdrfile_close(f);
}

void save_xdr(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	const double* data = reinterpret_cast<const double*>(grid);
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(data), ginfo.num_interval_points() * 1, f);
	xdrfile_close(f);
}

void save_xdr(const scalar_ptr_t(&grid)[3], symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	const double* data = reinterpret_cast<const double*>(grid);
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(data), ginfo.num_interval_points() * 3, f);
	xdrfile_close(f);
}

void save_xdr(const scalar_ptr_t(&grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	const double* data = reinterpret_cast<const double*>(grid);
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(data), ginfo.num_interval_points() * 2, f);
	xdrfile_close(f);
}

void save_xdr(const scalar_ptr_t(&grid)[1], symphas::io::write_info winfo, symphas::grid_info ginfo, bool is_checkpoint)
{
	const double* data = reinterpret_cast<const double*>(grid);
	XDRFILE* f = open_xdrgridf(winfo, is_checkpoint);
	print_xdr_header(winfo.index, winfo.id, ginfo, winfo, f);
	xdrfile_write_double(const_cast<double*>(data), ginfo.num_interval_points() * 1, f);
	xdrfile_close(f);
}



/*
 * functions to checkpoint the data to a file
 */
void symphas::io::xdr::save_grid_plotting(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

void symphas::io::xdr::save_grid_plotting(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

void symphas::io::xdr::save_grid_plotting(const double_arr2* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

void symphas::io::xdr::save_grid_plotting(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

void symphas::io::xdr::save_grid_plotting(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

void symphas::io::xdr::save_grid_plotting(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}


void symphas::io::xdr::save_grid_plotting(const scalar_ptr_t(&grid)[3], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

void symphas::io::xdr::save_grid_plotting(const scalar_ptr_t(&grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}

void symphas::io::xdr::save_grid_plotting(const scalar_ptr_t(&grid)[1], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, false);
}




/*
 * functions to checkpoint the data to a file
 */
void symphas::io::xdr::save_grid(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

void symphas::io::xdr::save_grid(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

void symphas::io::xdr::save_grid(const double_arr2* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

void symphas::io::xdr::save_grid(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

void symphas::io::xdr::save_grid(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

void symphas::io::xdr::save_grid(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

void symphas::io::xdr::save_grid(const scalar_ptr_t(&grid)[3], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

void symphas::io::xdr::save_grid(const scalar_ptr_t(&grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}

void symphas::io::xdr::save_grid(const scalar_ptr_t(&grid)[1], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	save_xdr(grid, winfo, ginfo, true);
}






