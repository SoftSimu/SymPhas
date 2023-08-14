
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
 * MODULE:  io
 * PURPOSE: Implements the write methods, where the write algorithm is chosen
 * based on the writer chosen. Also includes overloads of the write methods.
 *
 * ***************************************************************************
 */

#pragma once

#include "writeincludes.h"

template<typename T, typename>
void symphas::io::save_grid_plotting(const T* values, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	SWITCH_IO_WRITE(save_grid_plotting(values, winfo, ginfo))
}

template<typename T, size_t N>
void symphas::io::save_grid_plotting(const T(*values)[N], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	using data_unit = T[N];
	SWITCH_IO_WRITE(save_grid_plotting(values, winfo, ginfo))
}

template<typename T, size_t N>
void symphas::io::save_grid_plotting(T* const(&values)[N], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	using data_unit = T*[N];
	SWITCH_IO_WRITE(save_grid_plotting(values, winfo, ginfo))
}


template<typename T, typename>
void symphas::io::save_grid(const T* values, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	SWITCH_IO_WRITE(save_grid(values, winfo, ginfo))
}

template<typename T, size_t N>
void symphas::io::save_grid(const T(*values)[N], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	SWITCH_IO_WRITE(save_grid(values, winfo, ginfo))
}

template<typename T, size_t N>
void symphas::io::save_grid(T* const(&values)[N], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	SWITCH_IO_WRITE(save_grid(values, winfo, ginfo))
}


//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, typename>
void symphas::io::save_grid(const T* values, symphas::grid_info ginfo)
{
	save_grid(values, symphas::io::write_info{}, ginfo);
}

//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, size_t N>
void symphas::io::save_grid(const T(*values)[N], symphas::grid_info ginfo)
{
	using data_unit = T[N];
	save_grid(values, symphas::io::write_info{}, ginfo);
}

//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, size_t N>
void symphas::io::save_grid(T* const(&values)[N], symphas::grid_info ginfo)
{
	using data_unit = T[N];
	save_grid(values, symphas::io::write_info{}, ginfo);
}



//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, typename>
void symphas::io::save_grid(const T* values, symphas::io::write_info winfo, const len_type* dims, size_t dimension)
{
	symphas::io::save_grid(values, winfo, symphas::grid_info{ dims, dimension });
}

//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, size_t N>
void symphas::io::save_grid(const T(*values)[N], symphas::io::write_info winfo, const len_type* dims, size_t dimension)
{
	symphas::io::save_grid(values, winfo, symphas::grid_info{ dims, dimension });
}

//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, size_t N>
void symphas::io::save_grid(T* const(&values)[N], symphas::io::write_info winfo, const len_type* dims, size_t dimension)
{
	symphas::io::save_grid(values, winfo, symphas::grid_info{ dims, dimension });
}


//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, typename>
void symphas::io::save_grid(const T* values, const len_type* dims, size_t dimension)
{
	symphas::io::save_grid(values, symphas::io::write_info{}, dims, dimension);
}

//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, size_t N>
void symphas::io::save_grid(const T(*values)[N], const len_type* dims, size_t dimension)
{
	symphas::io::save_grid(values, symphas::io::write_info{}, dims, dimension);
}

//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, size_t N>
void symphas::io::save_grid(T* const(&values)[N], const len_type* dims, size_t dimension)
{
	symphas::io::save_grid(values, symphas::io::write_info{}, dims, dimension);
}




template<typename M>
void symphas::io::write_plot_config(M const& model, const char* directory, const char* const* names, SaveParams const& save)
{
	SWITCH_IO_WRITE(write_plot_config(model, directory, names, save))
}




template<typename T, size_t D>
void symphas::io::print_grid(FILE* out, Grid<T, D> const& grid)
{
	iter_type last_pos[D]{};

	iter_type last_x = grid.dims[0] - 1;
	iter_type last_y = grid.dims[1] - 2;

	grid::get_grid_position(last_pos, grid.dims, 0);
	printf("/");
	for (iter_type n = 0; n < grid::length<D>(grid.dims); ++n)
	{
		iter_type pos[D]{};
		grid::get_grid_position(pos, grid.dims, n);
		if (pos[1] != last_pos[1] && last_pos[0] == last_x)
		{
			if (last_pos[1] == 0)
			{
				printf(" \\ \n| ");
			}
			else if (last_pos[1] == last_y)
			{
				printf(" | \n\\ ");
			}
			else
			{
				printf(" | \n| ");
			}
		}
		else
		{
			printf(" ");
		}
		write_data_entry(out, grid::value_at(grid, pos));
		for (iter_type i = 0; i < D; ++i)
		{
			last_pos[i] = pos[i];
		}
	}
	printf(" /\n");
}

template<typename T, size_t D>
void symphas::io::print_grid(Grid<T, D> const& grid)
{
	print_grid(stdout, grid);
}

template<typename T, size_t D>
void symphas::io::print_grid(FILE* out, RegionalGrid<T, D> const& grid)
{
	iter_type last_pos[D]{};
	iter_type dims[D]{};

	for (iter_type i = 0; i < D; ++i)
	{
		dims[i] = grid.dims[i] - grid.region.boundary_size * 2;
	}

	iter_type last_x = dims[0] - 1;
	iter_type last_y = dims[1] - 2;

	grid::get_grid_position(last_pos, dims, 0);
	fprintf(out, "/");
	for (iter_type n = 0; n < grid::length<D>(dims); ++n)
	{
		iter_type pos[D]{};
		grid::get_grid_position(pos, dims, n);
		if (pos[1] != last_pos[1] && last_pos[0] == last_x)
		{
			if (last_pos[1] == 0)
			{
				fprintf(out, " \\ \n| ");
			}
			else if (last_pos[1] == last_y)
			{
				fprintf(out, " | \n\\ ");
			}
			else
			{
				fprintf(out, " | \n| ");
			}
		}
		else
		{
			fprintf(out, " ");
		}
		write_data_entry(out, grid::value_at(grid, pos));
		for (iter_type i = 0; i < D; ++i)
		{
			last_pos[i] = pos[i];
		}
	}
	fprintf(out, " /\n");
}

template<typename T, size_t D>
void symphas::io::print_grid(RegionalGrid<T, D> const& grid)
{
	print_grid(stdout, grid);
}

template<typename T, size_t D>
void symphas::io::print_region(FILE* out, RegionalGrid<T, D> const& grid)
{
	bool skip = false;

	iter_type pos[D]{};
	iter_type last_pos[D]{};
	iter_type inner_dims[D]{};
	iter_type inner_domain_dims[D]{};
	for (iter_type i = 0; i < D; ++i)
	{
		inner_dims[i] = grid.region.dims[i] - grid.region.boundary_size * 2;
		inner_domain_dims[i] = grid.dims[i] - grid.region.boundary_size * 2;

		skip = (inner_dims[i] < 0) ? true : skip;
	}
	 

	if (skip)
	{
		fprintf(out, "< >\n");
		return;
	}

	iter_type last_x = (inner_dims[0] + grid.region.origin[0]) % inner_domain_dims[0] - 1;
	iter_type last_y = (inner_dims[1] + grid.region.origin[1]) % inner_domain_dims[1] - 2;

	grid::get_grid_position_offset(last_pos, inner_dims, grid.region.origin, 0);
	fprintf(out, "/");
	for (iter_type n = 0; n < grid::length<D>(inner_dims); ++n)
	{
		grid::get_grid_position_offset(pos, inner_dims, grid.region.origin, n);
		for (iter_type i = 0; i < D; ++i)
		{
			pos[i] = (pos[i] >= inner_domain_dims[i]) ? pos[i] - inner_domain_dims[i] : pos[i];
		}
		if (pos[1] != last_pos[1] && last_pos[0] == last_x)
		{
			if (last_pos[1] == grid.region.origin[1] && last_pos[1] == last_y)
			{
				fprintf(out, " \\ \n\\ ");
			}
			else if (last_pos[1] == grid.region.origin[1])
			{
				fprintf(out, " \\ \n| ");
			}
			else if (last_pos[1] == last_y)
			{
				fprintf(out, " | \n\\ ");
			}
			else
			{
				fprintf(out, " | \n| ");
			}
		}
		else
		{
			fprintf(out, " ");
		}
		write_data_entry(out, grid::value_at(grid, pos));
		for (iter_type i = 0; i < D; ++i)
		{
			last_pos[i] = pos[i];
		}
	}
	fprintf(out, " /\n");
}

template<typename T, size_t D>
void symphas::io::print_region(RegionalGrid<T, D> const& grid)
{
	print_region(stdout, grid);
}


template<size_t D>
void symphas::io::print_region(FILE* out, grid::region_interval<D> const& region, char inside, char outside)
{
	len_type pos[D]{};
	len_type last_pos[D]{};
	for (iter_type n = 0; n < grid::length<D>(region.dims); ++n)
	{
		grid::get_grid_position(pos, region.dims, n);
		if (pos[1] != last_pos[1])
		{
			fprintf(out, "\n");
		}

		fprintf(out, "%c ", (grid::is_in_region(pos, region.intervals)) ? inside : outside);
		
		for (iter_type i = 0; i < D; ++i)
		{
			last_pos[i] = pos[i];
		}
	}
	fprintf(out, "\n");
}

template<size_t D>
void symphas::io::print_region(grid::region_interval<D> const& region, char inside, char outside)
{
	print_region(stdout, region, inside, outside);
}

template<size_t D>
void symphas::io::print_region(FILE* out, grid::region_interval_multiple<D> const& regions, char inside, char outside)
{
	len_type pos[D]{};
	len_type last_pos[D]{};
	for (iter_type n = 0; n < grid::length<D>(regions.dims); ++n)
	{
		grid::get_grid_position(pos, regions.dims, n);
		if (pos[1] != last_pos[1])
		{
			fprintf(out, "\n");
		}

		bool is_inside = false;
		for (grid::region_interval<D> region : regions)
		{
			is_inside = (grid::is_in_region(pos, region.intervals)) ? true : is_inside;
		}
		fprintf(out, "%c ", (is_inside) ? inside : outside);

		for (iter_type i = 0; i < D; ++i)
		{
			last_pos[i] = pos[i];
		}
	}
	fprintf(out, "\n");
}

template<size_t D>
void symphas::io::print_region(grid::region_interval_multiple<D> const& regions, char inside, char outside)
{
	print_region(stdout, regions, inside, outside);
}


