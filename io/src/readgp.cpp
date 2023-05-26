
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



#include "readgp.h"


FILE* symphas::io::gp::open_gpgridf(const char* name)
{
	FILE* f;
	if ((f = fopen(name, "r")) == 0)
	{
		symphas::lib::make_directory_for_file(name);
		if ((f = fopen(name, "r")) == 0)
		{
			fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_OPEN, name);
			exit(ERR_CODE_FILE_OPEN);
		}
	}
	return f;
}

template<>
symphas::grid_info symphas::io::gp::read_header(FILE* f, int* index)
{
	symphas::grid_info ginfo{ nullptr, 0 };

	iter_type read_index = (index) ? *index : BAD_INDEX;
	if (read_index < 0)
	{
		int dim;
		if (fscanf(f, "%d", &dim) != 1)
		{
			if (index)
			{
				*index = BAD_INDEX;
			}
			ginfo = { nullptr, dim };
		}

		len_type* dims = new len_type[dim];
		for (iter_type i = 0; i < dim; ++i)
		{
			if (fscanf(f, "%d", dims + i) != 1)
			{
				if (index)
				{
					*index = BAD_INDEX;
				}
				return { nullptr, dim };
			}
		}

		ginfo = symphas::grid_info{ dims, dim };


		/* the output of the intervals always goes x, y, z
			* keep checking size of dimension and print corresponding interval
			*/

		double left, right;
		for (iter_type i = 0; i < dim; ++i)
		{
			if (fscanf(f, "%lf", &left) + fscanf(f, "%lf", &right) != 2)
			{
				if (index)
				{
					*index = BAD_INDEX;
				}
				return { nullptr, dim };
			}
			ginfo.at(symphas::index_to_axis(i)).set_count(left, right, dims[i]);
		}

		delete[] dims;
	}
		
	if (index)
	{
		if (fscanf(f, "%d", index) != 1)
		{
			*index = BAD_INDEX;
		}
	}

	return ginfo;
}




template<>
void symphas::io::gp::read_block(scalar_t* grid, symphas::grid_info ginfo, FILE* f)
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
				fscanf(f, "%lf", grid + ii);
			}
		}
	}
}

template<>
void symphas::io::gp::read_block(complex_t* grid, symphas::grid_info ginfo, FILE* f)
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
				double re, im;
				fscanf(f, "%lf %lf", &re, &im);
				grid[ii] = complex_t{ re, im };
			}
		}
	}
}

template<>
void symphas::io::gp::read_block(double_arr2*grid, symphas::grid_info ginfo, FILE* f)
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
				fscanf(f, "%lf %lf", grid[ii], grid[ii] + 1);
			}
		}
	}

}

template<>
void symphas::io::gp::read_block(vector_t<3>* grid, symphas::grid_info ginfo, FILE* f)
{
	len_type N = ginfo.at(Axis::Z).get_count();
	len_type M = ginfo.at(Axis::Y).get_count();
	len_type L = ginfo.at(Axis::X).get_count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + (j * L) + (k * L * M);
				double m, dx, dy, dz;

				fscanf(f,
					"%lf %lf %lf %lf",
					&dx, &dy, &dz, &m);

				grid[ii] = vector_t<3>{ dx * m, dy * m, dz * m };
			}
		}
	}

}

template<>
void symphas::io::gp::read_block(vector_t<2>* grid, symphas::grid_info ginfo, FILE* f)
{
	for (iter_type j = 0; j < ginfo.at(Axis::Y).get_count(); j++)
	{
		for (iter_type i = 0; i < ginfo.at(Axis::X).get_count(); i++)
		{
			iter_type ii = i + j * ginfo.at(Axis::X).get_count();
			double m, dx, dy;

			fscanf(f,
				"%lf %lf %lf",
				&dx, &dy, &m);

			grid[ii] = vector_t<2>{ dx * m, dy * m };
		}
	}

}


template<>
void symphas::io::gp::read_block(vector_t<1>* grid, symphas::grid_info ginfo, FILE* f)
{
	for (iter_type i = 0; i < ginfo.at(Axis::X).get_count(); i++)
	{
		double m;
		fscanf(f, "%lf", &m);
		grid[i] = vector_t<1>{ m };
	}

}

template<>
void symphas::io::gp::read_block(scalar_ptr_t(&grid)[3], symphas::grid_info ginfo, FILE* f)
{
	len_type N = ginfo.at(Axis::Z).get_count();
	len_type M = ginfo.at(Axis::Y).get_count();
	len_type L = ginfo.at(Axis::X).get_count();

	for (iter_type k = 0; k < N; k++)
	{
		for (iter_type j = 0; j < M; j++)
		{
			for (iter_type i = 0; i < L; i++)
			{
				iter_type ii = i + (j * L) + (k * L * M);
				double m, dx, dy, dz;

				fscanf(f,
					"%lf %lf %lf %lf",
					&dx, &dy, &dz, &m);

				grid[0][ii] = dx * m;
				grid[1][ii] = dy * m;
				grid[2][ii] = dz * m;
			}
		}
	}

}

template<>
void symphas::io::gp::read_block(scalar_ptr_t(&grid)[2], symphas::grid_info ginfo, FILE* f)
{
	for (iter_type j = 0; j < ginfo.at(Axis::Y).get_count(); j++)
	{
		for (iter_type i = 0; i < ginfo.at(Axis::X).get_count(); i++)
		{
			iter_type ii = i + j * ginfo.at(Axis::X).get_count();
			double m, dx, dy;

			fscanf(f,
				"%lf %lf %lf",
				&dx, &dy, &m);

			grid[0][ii] = dx * m;
			grid[1][ii] = dy * m;
		}
	}

}


template<>
void symphas::io::gp::read_block(scalar_ptr_t(&grid)[1], symphas::grid_info ginfo, FILE* f)
{
	for (iter_type i = 0; i < ginfo.at(Axis::X).get_count(); i++)
	{
		double m;
		fscanf(f, "%lf", &m);
		grid[0][i] = m;
	}

}




