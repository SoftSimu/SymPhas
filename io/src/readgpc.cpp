
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
void symphas::io::gp::col::read_block(scalar_t* grid, symphas::grid_info ginfo, FILE* f)
{
	auto helper = symphas::io::gp::new_helper(ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}
				scalar_t value;
				fscanf(f, "%lf", &value);
				if (grid != nullptr)
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[ii] = value;
				}
			}
		}
	}
	symphas::io::gp::free_helper(helper);
}

template<>
void symphas::io::gp::col::read_block(complex_t* grid, symphas::grid_info ginfo, FILE* f)
{
	auto helper = symphas::io::gp::new_helper(ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}
				double re, im;
				fscanf(f, "%lf+i%lf", &re, &im);

				if (grid != nullptr)
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[ii] = complex_t{ re, im };
				}
			}
		}
	}
	symphas::io::gp::free_helper(helper);
}


template<>
void symphas::io::gp::col::read_block(double_arr2*grid, symphas::grid_info ginfo, FILE* f)
{
	auto helper = symphas::io::gp::new_helper(ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}

				double a, b;
				fscanf(f, "%lf %lf", &a, &b);
				if (grid != nullptr)
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[ii][0] = a;
					grid[ii][1] = b;
				}
			}
		}
	}
	symphas::io::gp::free_helper(helper);
}

template<>
void symphas::io::gp::col::read_block(vector_t<3>* grid, symphas::grid_info ginfo, FILE* f)
{
	auto helper = symphas::io::gp::new_helper(ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}

				double m, dx, dy, dz;

				fscanf(f,
					"%lf %lf %lf %lf",
					&dx, &dy, &dz, &m);

				if (grid != nullptr)
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[ii] = vector_t<3>{ dx * m, dy * m, dz * m };
				}
			}
		}
	}
	symphas::io::gp::free_helper(helper);
}

template<>
void symphas::io::gp::col::read_block(vector_t<2>* grid, symphas::grid_info ginfo, FILE* f)
{
	auto helper = symphas::io::gp::new_helper(ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}

				double m, dx, dy;

				fscanf(f,
					"%lf %lf %lf",
					&dx, &dy, &m);

				if (grid != nullptr)
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[ii] = vector_t<2>{ dx * m, dy * m };
				}
			}
		}
	}
	symphas::io::gp::free_helper(helper);
}


template<>
void symphas::io::gp::col::read_block(vector_t<1>* grid, symphas::grid_info ginfo, FILE* f)
{
	auto helper = symphas::io::gp::new_helper(ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}
				double m;
				fscanf(f, "%lf ", &m);

				iter_type ii = GP_HELPER_INDEX({ i, j, k });
				grid[ii] = vector_t<1>{ m };
			}
		}
	}
	symphas::io::gp::free_helper(helper);
}


template<>
void symphas::io::gp::col::read_block(scalar_ptr_t(&grid)[3], symphas::grid_info ginfo, FILE* f)
{
	auto helper = symphas::io::gp::new_helper(ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}

				double m, dx, dy, dz;

				fscanf(f,
					"%lf %lf %lf %lf",
					&dx, &dy, &dz, &m);

				if (grid != nullptr)
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[0][ii] = dx * m;
					grid[1][ii] = dy * m;
					grid[2][ii] = dz * m;
				}
			}
		}
	}
	symphas::io::gp::free_helper(helper);
}

template<>
void symphas::io::gp::col::read_block(scalar_ptr_t(&grid)[2], symphas::grid_info ginfo, FILE* f)
{
	auto helper = symphas::io::gp::new_helper(ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}

				double m, dx, dy;

				fscanf(f,
					"%lf %lf %lf",
					&dx, &dy, &m);

				if (grid != nullptr)
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[0][ii] = dx * m;
					grid[1][ii] = dy * m;
				}
			}
		}
	}
	symphas::io::gp::free_helper(helper);
}


template<>
void symphas::io::gp::col::read_block(scalar_ptr_t(&grid)[1], symphas::grid_info ginfo, FILE* f)
{
	auto helper = symphas::io::gp::new_helper(ginfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				for (len_type n = 0; n < ginfo.dimension(); ++n)
				{
					fscanf(f, "%*f");
				}
				double m;
				fscanf(f, "%lf ", &m);

				if (grid != nullptr)
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[0][ii] = m;
				}
			}
		}
	}
	symphas::io::gp::free_helper(helper);
}






