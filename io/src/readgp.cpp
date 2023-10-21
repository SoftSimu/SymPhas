
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


FILE* symphas::io::gp::open_gpgridf_nofail(const char* name)
{
	FILE* f;
	if ((f = fopen(name, "r")) == 0)
	{
		symphas::lib::make_directory_for_file(name);
		return fopen(name, "r");
	}
	else
	{
		return f;
	}
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
			if (fscanf(f, "%lf %lf", &left, &right) != 2)
			{
				if (index)
				{
					*index = BAD_INDEX;
				}
				return { nullptr, dim };
			}
			ginfo.at(symphas::index_to_axis(i)).set_count(left, right, dims[i]);
			ginfo.at(symphas::index_to_axis(i)).interval_to_domain();
		}

		delete[] dims;
	}

	int get;
	bool index_scanned = (fscanf(f, "%d", &get) == 1);

	if (index)
	{
		*index = (index_scanned) ? get : BAD_INDEX;
	}

	char buffer[4]{};
	if (fscanf(f, " %s", buffer) == 1)
	{
		if (*buffer != CONFIG_OPTION_PREFIX_C)
		{
			size_t dimension;
			if (sscanf(buffer, "%zd", &dimension) != 1)
			{
				if (index)
				{
					*index = BAD_INDEX;
				}
				return { nullptr, ginfo.dimension() };
			}
			else
			{
				for (iter_type i = 0; i < dimension; ++i)
				{
					double value[2]{};
					if (fscanf(f, "%lf %lf", value, value + 1) == 2)
					{
						Axis axis(symphas::index_to_axis(i));

						if (ginfo.intervals.count(axis) == 0)
						{
							symphas::interval_element_type interval;
							interval.set_interval(value[0], value[1]);
							interval.domain_to_interval();
							ginfo[axis] = interval;
						}
						else
						{
							ginfo[axis].set_interval(value[0], value[1]);
						}
					}
					else
					{
						if (index)
						{
							*index = BAD_INDEX;
						}
						return { nullptr, ginfo.dimension() };
					}
				}
			}
		}
		return ginfo;
	}
	else
	{
		if (index)
		{
			*index = BAD_INDEX;
		}
		return { nullptr, ginfo.dimension() };
	}

}




template<>
void symphas::io::gp::read_block(scalar_t* grid, symphas::io::block_info binfo, FILE* f)
{
	auto helper = symphas::io::new_helper(binfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				scalar_t value;
				fscanf(f, "%lf", &value);

				if (grid != nullptr && helper->in_bounds({ i, j, k }))
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[ii] = value;
				}
			}
		}
	}
	symphas::io::free_helper(helper);
}

template<>
void symphas::io::gp::read_block(complex_t* grid, symphas::io::block_info binfo, FILE* f)
{
	auto helper = symphas::io::new_helper(binfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				double re, im;
				fscanf(f, "%lf %lf", &re, &im);

				if (grid != nullptr && helper->in_bounds({ i, j, k }))
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[ii] = complex_t{ re, im };
				}
			}
		}
	}
	symphas::io::free_helper(helper);
}

template<>
void symphas::io::gp::read_block(double_arr2*grid, symphas::io::block_info binfo, FILE* f)
{
	auto helper = symphas::io::new_helper(binfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				double a, b;
				fscanf(f, "%lf %lf", &a, &b);

				if (grid != nullptr && helper->in_bounds({ i, j, k }))
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[ii][0] = a;
					grid[ii][1] = b;
				}
			}
		}
	}
	symphas::io::free_helper(helper);
}

template<>
void symphas::io::gp::read_block(vector_t<3>* grid, symphas::io::block_info binfo, FILE* f)
{
	auto helper = symphas::io::new_helper(binfo);

	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				double m, dx, dy, dz;

				fscanf(f,
					"%lf %lf %lf %lf",
					&dx, &dy, &dz, &m);

				if (grid != nullptr && helper->in_bounds({ i, j, k }))
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[ii] = vector_t<3>{ dx * m, dy * m, dz * m };
				}
			}
		}
	}
	symphas::io::free_helper(helper);
}

template<>
void symphas::io::gp::read_block(vector_t<2>* grid, symphas::io::block_info binfo, FILE* f)
{
	auto helper = symphas::io::new_helper(binfo);
	for (iter_type j = 0; j < GP_HELPER_LENY; j++)
	{
		for (iter_type i = 0; i < GP_HELPER_LENX; i++)
		{
			double m, dx, dy;

			fscanf(f,
				"%lf %lf %lf",
				&dx, &dy, &m);

			if (grid != nullptr && helper->in_bounds({ i, j }))
			{
				iter_type ii = GP_HELPER_INDEX({ i, j });
				grid[ii] = vector_t<2>{ dx * m, dy * m };
			}
		}
	}
	symphas::io::free_helper(helper);
}


template<>
void symphas::io::gp::read_block(vector_t<1>* grid, symphas::io::block_info binfo, FILE* f)
{
	auto helper = symphas::io::new_helper(binfo);
	for (iter_type i = 0; i < GP_HELPER_LENX; i++)
	{
		double m;
		fscanf(f, "%lf", &m);

		if (grid != nullptr && helper->in_bounds({ i }))
		{
			iter_type ii = GP_HELPER_INDEX({ i });
			grid[ii] = vector_t<1>{ m };
		}
	}
	symphas::io::free_helper(helper);
}

template<>
void symphas::io::gp::read_block(scalar_ptr_t(&grid)[3], symphas::io::block_info binfo, FILE* f)
{
	auto helper = symphas::io::new_helper(binfo);
	for (iter_type k = 0; k < GP_HELPER_LENZ; k++)
	{
		for (iter_type j = 0; j < GP_HELPER_LENY; j++)
		{
			for (iter_type i = 0; i < GP_HELPER_LENX; i++)
			{
				double m, dx, dy, dz;

				fscanf(f,
					"%lf %lf %lf %lf",
					&dx, &dy, &dz, &m);
				
				if (*grid != nullptr && helper->in_bounds({ i, j, k }))
				{
					iter_type ii = GP_HELPER_INDEX({ i, j, k });
					grid[0][ii] = dx * m;
					grid[1][ii] = dy * m;
					grid[2][ii] = dz * m;
				}
			}
		}
	}
	symphas::io::free_helper(helper);
}

template<>
void symphas::io::gp::read_block(scalar_ptr_t(&grid)[2], symphas::io::block_info binfo, FILE* f)
{
	auto helper = symphas::io::new_helper(binfo);
	for (iter_type j = 0; j < GP_HELPER_LENY; j++)
	{
		for (iter_type i = 0; i < GP_HELPER_LENX; i++)
		{
			double m, dx, dy;

			fscanf(f,
				"%lf %lf %lf",
				&dx, &dy, &m);

			if (*grid != nullptr && helper->in_bounds({ i, j }))
			{
				iter_type ii = GP_HELPER_INDEX({ i, j });
				grid[0][ii] = dx * m;
				grid[1][ii] = dy * m;
			}
		}
	}
	symphas::io::free_helper(helper);
}


template<>
void symphas::io::gp::read_block(scalar_ptr_t(&grid)[1], symphas::io::block_info binfo, FILE* f)
{
	auto helper = symphas::io::new_helper(binfo);
	for (iter_type i = 0; i < GP_HELPER_LENX; i++)
	{
		double m;
		fscanf(f, "%lf", &m);
		if (*grid != nullptr && helper->in_bounds({ i }))
		{
			iter_type ii = GP_HELPER_INDEX({ i });
			grid[0][ii] = m;
		}
	}
	symphas::io::free_helper(helper);
}




