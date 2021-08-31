
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



#include "systemlib.h"

symphas::problem_parameters_type::problem_parameters_type(problem_parameters_type const& other) : problem_parameters_type(other.len)
{
	set_initial_data(other.tdata, len);
	set_interval_data(other.vdata, len);
	set_boundary_data(other.bdata, len);
	set_problem_time_step(other.dt);
}

void symphas::problem_parameters_type::set_initial_data(symphas::init_data_type* tdata_set, size_t n)
{
	if (n <= len)
	{
		for (size_t i = 0; i < n; ++i)
		{
			tdata[i] = tdata_set[i];
		}
		for (size_t i = n; i < len; ++i)
		{
			tdata[i] = tdata_set[0];
		}
	}
	else
	{
		fprintf(SYMPHAS_LOG, "the number of given initial data is larger than the number of systems\n");
	}
}

void symphas::problem_parameters_type::set_interval_data(symphas::interval_data_type* vdata_set, size_t n)
{
	if (n <= len)
	{
		for (size_t i = 0; i < n; ++i)
		{
			vdata[i] = vdata_set[i];
		}
		for (size_t i = n; i < len; ++i)
		{
			vdata[i] = vdata_set[0];
		}
	}
	else
	{
		fprintf(SYMPHAS_LOG, "the number of given interval data is larger than the number of systems\n");
	}
}

void symphas::problem_parameters_type::set_boundary_data(symphas::b_data_type* bdata_set, size_t n)
{
	if (n <= len)
	{
		for (size_t i = 0; i < n; ++i)
		{
			bdata[i] = bdata_set[i];
		}
		for (size_t i = n; i < len; ++i)
		{
			bdata[i] = bdata_set[0];
		}
	}
	else
	{
		fprintf(SYMPHAS_LOG, "the number of given boundary data is larger than the number of systems\n");
	}
}

void symphas::problem_parameters_type::equalize_discretization(const len_type* dims)
{
	for (iter_type i = 0; i < len; ++i)
	{
		size_t dim = vdata[i].size();
		for (iter_type n = 0; n < dim; ++n)
		{
			auto& interval = vdata[i].at(symphas::index_to_axis(n));

			if (interval.count() != dims[n])
			{
				interval.set_interval_count(interval.right(), interval.left(), dims[n]);
			}
		}
	}
}

void symphas::problem_parameters_type::equalize_discretization(iter_type i)
{
	size_t dim = vdata[i].size();
	len_type* dims = new len_type[dim];

	for (iter_type n = 0; n < dim; ++n)
	{
		dims[n] = vdata[i].at(symphas::index_to_axis(n)).count();
	}

	equalize_discretization(dims);
	delete[] dims;
}

symphas::problem_parameters_type::~problem_parameters_type()
{
	delete[] tdata;
	delete[] vdata;
	delete[] bdata;
}


void symphas::swap(symphas::problem_parameters_type& first, symphas::problem_parameters_type& second)
{
	using std::swap;
	swap(first.tdata, second.tdata);
	swap(first.vdata, second.vdata);
	swap(first.bdata, second.bdata);
	swap(first.dt, second.dt);
	swap(first.len, second.len);
}


void swap(symphas::problem_parameters_type& first, symphas::problem_parameters_type& second)
{
	symphas::swap(first, second);
}




