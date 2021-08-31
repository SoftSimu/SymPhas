
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
template<>
void symphas::io::gpc::save_grid_plotting<scalar_t>(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);
	
	double dZ = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Zh;
	double dY = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Yh;
	double dX = ginfo.INTERVAL_Xh;
	double h[]{ dX, dY, dZ };

	double Z0 = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Z0;
	double Y0 = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Y0;
	double X0 = ginfo.INTERVAL_X0;
	double p0[]{ X0, Y0, Z0 };

	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	iter_type index = 0;
	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type p[]{ i, j, k };
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", p0[n] + h[n] * p[n]);
				}
				fprintf(f, "%lE\n", grid[index++]);
			}
		}
	}
	fprintf(f, "\n");
	fclose(f);
}

template<>
void symphas::io::gpc::save_grid_plotting<complex_t>(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	double dZ = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Zh;
	double dY = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Yh;
	double dX = ginfo.INTERVAL_Xh;
	double h[]{ dX, dY, dZ };

	double Z0 = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Z0;
	double Y0 = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Y0;
	double X0 = ginfo.INTERVAL_X0;
	double p0[]{ X0, Y0, Z0 };

	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	iter_type index = 0;
	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type p[]{ i, j, k };
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", p0[n] + h[n] * p[n]);
				}
				fprintf(f, "%lE+i%lE\n", grid[index].real(), grid[index].imag());
				++index;
			}
		}
	}
	fprintf(f, "\n");
	fclose(f);
}

template<>
void symphas::io::gpc::save_grid_plotting<double[2]>(const double (*grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	double dZ = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Zh;
	double dY = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Yh;
	double dX = ginfo.INTERVAL_Xh;
	double h[]{ dX, dY, dZ };

	double Z0 = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Z0;
	double Y0 = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Y0;
	double X0 = ginfo.INTERVAL_X0;
	double p0[]{ X0, Y0, Z0 };

	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	iter_type index = 0;
	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type p[]{ i, j, k };
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", p0[n] + h[n] * p[n]);
				}

				fprintf(f, "%lE %lE\n", grid[index][0], grid[index][1]);
				++index;
			}
		}
	}
	fprintf(f, "\n");
	fclose(f);
}

template<>
void symphas::io::gpc::save_grid_plotting<vector_t<3>>(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	double dZ = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Zh;
	double dY = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Yh;
	double dX = ginfo.INTERVAL_Xh;
	double h[]{ dX, dY, dZ };

	double Z0 = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Z0;
	double Y0 = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Y0;
	double X0 = ginfo.INTERVAL_X0;
	double p0[]{ X0, Y0, Z0 };

	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	iter_type index = 0;
	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type p[]{ i, j, k };
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", p0[n] + h[n] * p[n]);
				}

				vector_t<3> v = grid[index++];

				double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2]),
					dx = v.v[0] / m,
					dy = v.v[1] / m,
					dz = v.v[2] / m;

				fprintf(f,
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f\n",
					dx, dy, dz, m);
			}
		}
	}
	fprintf(f, "\n");
	fclose(f);
}

template<>
void symphas::io::gpc::save_grid_plotting<vector_t<2>>(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	double dZ = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Zh;
	double dY = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Yh;
	double dX = ginfo.INTERVAL_Xh;
	double h[]{ dX, dY, dZ };

	double Z0 = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Z0;
	double Y0 = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Y0;
	double X0 = ginfo.INTERVAL_X0;
	double p0[]{ X0, Y0, Z0 };

	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	iter_type index = 0;
	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type p[]{ i, j, k };
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", p0[n] + h[n] * p[n]);
				}

				vector_t<2> v = grid[index++];

				double m = sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1]),
					dx = v.v[0] / m,
					dy = v.v[1] / m;

				fprintf(f,
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f "
					"%." DATA_OUTPUT_ACCURACY_STR "f\n",
					dx, dy, m);
			}
		}
	}
	fprintf(f, "\n");
	fclose(f);
}


template<>
void symphas::io::gpc::save_grid_plotting<vector_t<1>>(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	FILE* f = symphas::io::open_data_file(winfo.dir_str_ptr, winfo.index, winfo.id, winfo.type);

	double dZ = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Zh;
	double dY = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Yh;
	double dX = ginfo.INTERVAL_Xh;
	double h[]{ dX, dY, dZ };

	double Z0 = (ginfo.dimension() < 3) ? 0 : ginfo.INTERVAL_Z0;
	double Y0 = (ginfo.dimension() < 2) ? 0 : ginfo.INTERVAL_Y0;
	double X0 = ginfo.INTERVAL_X0;
	double p0[]{ X0, Y0, Z0 };

	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	iter_type index = 0;
	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type p[]{ i, j, k };
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fprintf(f, "%." AXIS_OUTPUT_ACCURACY_STR "f ", p0[n] + h[n] * p[n]);
				}

				vector_t<1> v = grid[index++];
				fprintf(f, 	"%." DATA_OUTPUT_ACCURACY_STR "f\n", v.v);
			}
		}
	}
	fprintf(f, "\n");
	fclose(f);
}


