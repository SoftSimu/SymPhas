
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



#include "readgpc.h"




template<>
void symphas::io::gp::col::read_block<scalar_t>(scalar_t* grid, symphas::grid_info ginfo, FILE* f)
{
	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}
				iter_type ii = i + j * L + k * L * M;
				fscanf(f, "%lf", grid + ii);
			}
		}
	}
}

template<>
void symphas::io::gp::col::read_block<complex_t>(complex_t* grid, symphas::grid_info ginfo, FILE* f)
{
	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}
				double re, im;
				fscanf(f, "%lf+i%lf", &re, &im);

				iter_type ii = i + j * L + k * L * M;
				grid[ii] = complex_t{ re, im };
			}
		}
	}
}

template<>
void symphas::io::gp::col::read_block<double[2]>(double(*grid)[2], symphas::grid_info ginfo, FILE* f)
{
	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}
				iter_type ii = i + j * L + k * L * M;
				fscanf(f, "%lf %lf", grid[ii], grid[ii] + 1);
			}
		}
	}
}

template<>
void symphas::io::gp::col::read_block<vector_t<3>>(vector_t<3>* grid, symphas::grid_info ginfo, FILE* f)
{
	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}

				double m, dx, dy, dz;

				fscanf(f,
					"%lf %lf %lf %lf",
					&dx, &dy, &dz, &m);

				iter_type ii = i + j * L + k * L * M;
				grid[ii] = vector_t<3>{ dx * m, dy * m, dz * m };
			}
		}
	}
}

template<>
void symphas::io::gp::col::read_block<vector_t<2>>(vector_t<2>* grid, symphas::grid_info ginfo, FILE* f)
{
	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}

				double m, dx, dy;

				fscanf(f,
					"%lf %lf %lf",
					&dx, &dy, &m);

				iter_type ii = i + j * L + k * L * M;
				grid[ii] = vector_t<2>{ dx * m, dy * m };
			}
		}
	}
}


template<>
void symphas::io::gp::col::read_block<vector_t<1>>(vector_t<1>* grid, symphas::grid_info ginfo, FILE* f)
{
	len_type N = (ginfo.dimension() < 3) ? 1 : ginfo.at(Axis::Z).count();
	len_type M = (ginfo.dimension() < 2) ? 1 : ginfo.at(Axis::Y).count();
	len_type L = ginfo.at(Axis::X).count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}
				double m;
				fscanf(f, "%lf ", &m);

				iter_type ii = i + j * L + k * L * M;
				grid[ii] = vector_t<1>{ m };
			}
		}
	}
}





