
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


#include "writecsv.h"
#include "writegp.h"




/*
 * functions to write data from the grid to a file
 */
void symphas::io::csv::save_grid_plotting(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	if (ginfo.dimension() > 1)
	{
		len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).get_count();
		len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).get_count();
		len_type L = ginfo.at(Axis::X).get_count();

		for (iter_type k = 0; k < N; k++)
		{

			for (iter_type j = 0; j < M; j++)
			{
				for (iter_type i = 0; i < L; i++)
				{
					iter_type ii = i + j * L + k * L * M;
					fprintf(f, "% 10lE", grid[ii]);
					if (i < L - 1)
					{
						fprintf(f, ",");
					}
				}
				fprintf(f, "\n");
			}
		}
	}
	else
	{
		len_type L = ginfo.at(Axis::X).get_count();
		for (iter_type i = 0; i < L; i++)
		{
			fprintf(f, "% 10lE", grid[i]);
			if (i < L - 1)
			{
				fprintf(f, ",");
			}
		}
		fprintf(f, "\n");
	}

	fclose(f);
}

void symphas::io::csv::save_grid_plotting(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	if (ginfo.dimension() > 1)
	{
		len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).get_count();
		len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).get_count();
		len_type L = ginfo.at(Axis::X).get_count();

		for (iter_type k = 0; k < N; k++)
		{

			for (iter_type j = 0; j < M; j++)
			{
				for (iter_type i = 0; i < L; i++)
				{
					iter_type ii = i + j * L + k * L * M;
					fprintf(f, "\"% 10lE,%10lE\"", grid[ii].real(), grid[ii].imag());
					if (i < L - 1)
					{
						fprintf(f, ",");
					}
				}
				fprintf(f, "\n");
			}
		}
	}
	else
	{
		len_type L = ginfo.at(Axis::X).get_count();
		for (iter_type i = 0; i < L; i++)
		{
			fprintf(f, "\"% 10lE,%10lE\"", grid[i].real(), grid[i].imag());
			if (i < L - 1)
			{
				fprintf(f, ",");
			}
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void symphas::io::csv::save_grid_plotting(const double_arr2 *grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	if (ginfo.dimension() > 1)
	{
		len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).get_count();
		len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).get_count();
		len_type L = ginfo.at(Axis::X).get_count();

		for (iter_type k = 0; k < N; k++)
		{

			for (iter_type j = 0; j < M; j++)
			{
				for (iter_type i = 0; i < L; i++)
				{
					iter_type ii = i + j * L + k * L * M;
					fprintf(f, "\"% 10lE,%10lE\"", grid[ii][0], grid[ii][1]);
					if (i < L - 1)
					{
						fprintf(f, ",");
					}
				}
				fprintf(f, "\n");
			}
		}
	}
	else
	{
		len_type L = ginfo.at(Axis::X).get_count();
		for (iter_type i = 0; i < L; i++)
		{
			fprintf(f, "\"% 10lE,%10lE\"", grid[i][0], grid[i][1]);
			if (i < L - 1)
			{
				fprintf(f, ",");
			}
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void symphas::io::csv::save_grid_plotting(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	len_type N = ginfo.at(Axis::Z).get_count();
	len_type M = ginfo.at(Axis::Y).get_count();
	len_type L = ginfo.at(Axis::X).get_count();

	for (iter_type k = 0; k < N; k++)
	{

		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;

				vector_t<3> v = grid[ii];
				double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2]),
					dx = v.v[0] / m,
					dy = v.v[1] / m,
					dz = v.v[2] / m;

				fprintf(f,
					"\"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f\"",
					dx, dy, dz, m);
					
				if (i < L - 1)
				{
					fprintf(f, ",");
				}
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

void symphas::io::csv::save_grid_plotting(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	len_type N = ginfo.at(Axis::Z).get_count();
	len_type M = ginfo.at(Axis::Y).get_count();
	len_type L = ginfo.at(Axis::X).get_count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;

				vector_t<2> v = grid[ii];
				double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1]),
					dx = v.v[0] / m,
					dy = v.v[1] / m;

				fprintf(f,
					"\"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f\"",
					dx, dy, m);

				if (i < L - 1)
				{
					fprintf(f, ",");
				}
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}


void symphas::io::csv::save_grid_plotting(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	len_type N = ginfo.at(Axis::Z).get_count();
	len_type M = ginfo.at(Axis::Y).get_count();
	len_type L = ginfo.at(Axis::X).get_count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;

				vector_t<1> v = grid[ii];

				fprintf(f,
					"\"%" DATA_OUTPUT_ACCURACY_STR "f\"",
					v.v[0]);

				if (i < L - 1)
				{
					fprintf(f, ",");
				}
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}


void symphas::io::csv::save_grid_plotting(const scalar_ptr_t(&grid)[3], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	len_type N = ginfo.at(Axis::Z).get_count();
	len_type M = ginfo.at(Axis::Y).get_count();
	len_type L = ginfo.at(Axis::X).get_count();

	for (iter_type k = 0; k < N; k++)
	{

		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;

				double m = sqrt(grid[0][ii] * grid[0][ii] + grid[1][ii] * grid[1][ii] + grid[2][ii] * grid[2][ii]),
					dx = grid[0][ii] / m,
					dy = grid[1][ii] / m,
					dz = grid[2][ii] / m;

				fprintf(f,
					"\"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f\"",
					dx, dy, dz, m);

				if (i < L - 1)
				{
					fprintf(f, ",");
				}
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

void symphas::io::csv::save_grid_plotting(const scalar_ptr_t(&grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	len_type N = ginfo.at(Axis::Z).get_count();
	len_type M = ginfo.at(Axis::Y).get_count();
	len_type L = ginfo.at(Axis::X).get_count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;

				double m = sqrt(grid[0][ii] * grid[0][ii] + grid[1][ii] * grid[1][ii]),
					dx = grid[0][ii] / m,
					dy = grid[1][ii] / m;

				fprintf(f,
					"\"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f "
					"%" DATA_OUTPUT_ACCURACY_STR "f\"",
					dx, dy, m);

				if (i < L - 1)
				{
					fprintf(f, ",");
				}
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}


void symphas::io::csv::save_grid_plotting(const scalar_ptr_t(&grid)[1], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	len_type N = ginfo.at(Axis::Z).get_count();
	len_type M = ginfo.at(Axis::Y).get_count();
	len_type L = ginfo.at(Axis::X).get_count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + j * L + k * L * M;


				fprintf(f,
					"\"%" DATA_OUTPUT_ACCURACY_STR "f\"",
					grid[0][ii]);

				if (i < L - 1)
				{
					fprintf(f, ",");
				}
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}







