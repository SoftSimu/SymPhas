
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


#include "writegpc.h"


/*
 * functions to write data from the grid to a file
 */
void symphas::io::gp::col::save_grid_plotting(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (auto const& [axis, _] : ginfo)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POS(axis, { i, j, k }));
				}
				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				fprintf(f, "%lE\n", grid[ii]);
			}
		}
	}
	symphas::io::free_helper(helper);
	fprintf(f, "\n");
	fclose(f);
}

void symphas::io::gp::col::save_grid_plotting(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (auto const& [axis, _] : ginfo)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POS(axis, { i, j, k }));
				}
				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				fprintf(f, "%lE+i%lE\n", grid[ii].real(), grid[ii].imag());
			}
		}
	}
	symphas::io::free_helper(helper);
	fprintf(f, "\n");
	fclose(f);
}

void symphas::io::gp::col::save_grid_plotting(const double_arr2 *grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (auto const& [axis, _] : ginfo)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POS(axis, { i, j, k }));
				}
				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				fprintf(f, "%lE %lE\n", grid[ii][0], grid[ii][1]);
			}
		}
	}
	symphas::io::free_helper(helper);
	fprintf(f, "\n");
	fclose(f);
}

void symphas::io::gp::col::save_grid_plotting(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (auto const& [axis, _] : ginfo)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POS(axis, { i, j, k }));
				}

				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				vector_t<3> const &v = grid[ii];

				double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2]),
					dx = v.v[0] / m,
					dy = v.v[1] / m,
					dz = v.v[2] / m;

				fprintf(f,
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f\n",
					dx, dy, dz, m);
			}
		}
	}
	symphas::io::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}

void symphas::io::gp::col::save_grid_plotting(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (auto const& [axis, _] : ginfo)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POS(axis, { i, j, k }));
				}

				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				vector_t<2> const &v = grid[ii];

				double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1]),
					dx = v.v[0] / m,
					dy = v.v[1] / m;

				fprintf(f,
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f\n",
					dx, dy, m);
			}
		}
	}
	symphas::io::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}


void symphas::io::gp::col::save_grid_plotting(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (auto const& [axis, _] : ginfo)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POS(axis, { i, j, k }));
				}
                
				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				const vector_t<1> &data = grid[ii];
				fprintf(f, 	"%" DATA_OUTPUT_ACCURACY_STR "f\n", *data.v);
			}
		}
	}
	symphas::io::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}


void symphas::io::gp::col::save_grid_plotting(const scalar_ptr_t(&grid)[3], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (auto const& [axis, _] : ginfo)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POS(axis, { i, j, k }));
				}

				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				double m = sqrt(grid[0][ii] * grid[0][ii] + grid[1][ii] * grid[1][ii] + grid[2][ii] * grid[2][ii]),
					dx = grid[0][ii] / m,
					dy = grid[1][ii] / m,
					dz = grid[2][ii] / m;

				fprintf(f,
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f\n",
					dx, dy, dz, m);

			}
		}
	}
	symphas::io::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}

void symphas::io::gp::col::save_grid_plotting(const scalar_ptr_t(&grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (auto const& [axis, _] : ginfo)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POS(axis, { i, j, k }));
				}

				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				double m = sqrt(grid[0][ii] * grid[0][ii] + grid[1][ii] * grid[1][ii]),
					dx = grid[0][ii] / m,
					dy = grid[1][ii] / m;

				fprintf(f,
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f\n",
					dx, dy, m);

			}
		}
	}
	symphas::io::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}


void symphas::io::gp::col::save_grid_plotting(const scalar_ptr_t(&grid)[1], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir, winfo.index, winfo.id, winfo.type);
	auto helper = symphas::io::new_helper(winfo, ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (auto const& [axis, _] : ginfo)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", GP_HELPER_POS(axis, { i, j, k }));
				}
				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				fprintf(f, "%" DATA_OUTPUT_ACCURACY_STR "f\n", grid[0][ii]);
			}
		}
	}
	symphas::io::free_helper(helper);
	fprintf(f, "\n\n");
	fclose(f);
}


