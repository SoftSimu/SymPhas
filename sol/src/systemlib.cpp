
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


double symphas::time_step_list::get_time_step(double time) const
{
	if (dts_len > 0)
	{
		if (time == t_dts[0])
		{
			return t_dts[1];
		}
		for (iter_type i = 1; i < dts_len; ++i)
		{
			if (time == t_dts[i * 2])
			{
				return t_dts[i * 2 + 1];
			}
			else if (time < t_dts[i * 2] && time >= t_dts[(i - 1) * 2])
			{
				return t_dts[(i - 1) * 2 + 1];
			}
		}
		return t_dts[(dts_len - 1) * 2 + 1];
	}
	else
	{
		return 0;
	}
}


 //! Provide the time step used in the solution of the problem.
void symphas::time_step_list::set_time_step(double dt, double time, bool insert)
{
	if (dts_len > 0)
	{
		if (t_dts[0] == time)
		{
			t_dts[0] = time;
			t_dts[1] = dt;
		}
		else
		{
			for (iter_type i = 1; i < dts_len + 1; ++i)
			{
				if (i == dts_len)
				{
					if (insert)
					{
						double* t_dts_updated = new double[(dts_len + 1) * 2] {};

						std::copy(t_dts, t_dts + dts_len * 2, t_dts_updated);

						t_dts_updated[dts_len * 2] = time;
						t_dts_updated[dts_len * 2 + 1] = dt;

						std::swap(t_dts, t_dts_updated);
						delete[] t_dts_updated;
						++dts_len;
					}
					else
					{
						t_dts[(i - 1) * 2] = time;
						t_dts[(i - 1) * 2 + 1] = dt;
					}
					break;
				}
				else
				{
					if (t_dts[i * 2] == time)
					{
						t_dts[i * 2] = time;
						t_dts[i * 2 + 1] = dt;
						break;
					}
					else if (time < t_dts[i * 2] && time > t_dts[(i - 1) * 2])
					{
						if (insert)
						{
							double* t_dts_updated = new double[(dts_len + 1) * 2] {};

							std::copy(t_dts, t_dts + i * 2, t_dts_updated);
							std::copy(t_dts + i * 2, t_dts + dts_len, t_dts_updated + (i + 1) * 2);

							t_dts_updated[i * 2] = time;
							t_dts_updated[i * 2 + 1] = dt;

							std::swap(t_dts, t_dts_updated);
							delete[] t_dts_updated;
							++dts_len;
						}
						else
						{
							t_dts[i * 2] = time;
							t_dts[i * 2 + 1] = dt;
						}
						break;
					}
				}
			}
		}
	}
}

void symphas::time_step_list::clear_time_steps(double default_dt)
{
	delete[] t_dts;
	dts_len = 1;

	t_dts = new double[dts_len * 2] { TIME_INIT, 1.0 };
}

//! Provide the time step used in the solution of the problem.
void symphas::time_step_list::set_time_steps(const double* dts, const double* t_dts, size_t dts_len)
{
	if (dts_len > 0)
	{
		delete[] t_dts;
		this->dts_len = dts_len;

		this->t_dts = new double[this->dts_len * 2] {};

		for (iter_type i = 0; i < dts_len; ++i)
		{
			this->t_dts[i * 2] = t_dts[i];
			this->t_dts[i * 2 + 1] = dts[i];
		}
	}
}



void symphas::problem_parameters_type::set_modifiers(const char* const* str, size_t n)
{
	n = (n > 0) ? std::min(n, num_fields_len) : num_fields_len;
	for (iter_type i = 0; i < n; ++i)
	{
		set_modifier(str[i], i);
	}
}

void symphas::problem_parameters_type::set_modifier(const char* str, iter_type i)
{
	if (str)
	{
		size_t next = 0;
		size_t end = std::strlen(str);
		char buffer[BUFFER_LENGTH_R2];


		std::map<std::string, ModelModifiers> map = {
			{ "plot.default" , ModelModifiers::PLOT_DEFAULT },
			{ "plot.max" , ModelModifiers::PLOT_MAX },
			{ "plot.min" , ModelModifiers::PLOT_MIN },
			{ "plot.contours" , ModelModifiers::PLOT_CONTOURS } };

		while (next < end)
		{
			sscanf(str, "%s %zn", buffer, &next);
			for (auto [k, v] : map)
			{
				if (k == buffer)
				{
					set_modifier(v, i);
				}
			}
		}
	}
}

symphas::problem_parameters_type::problem_parameters_type(problem_parameters_type const& other) : problem_parameters_type(other.len)
{
	set_initial_data(other.tdata, len);
	set_interval_data(other.vdata, len);
	set_boundary_data(other.bdata, len);
	set_time_steps(other.dt_list);
	time = other.time;
	index = other.index;
}

void symphas::problem_parameters_type::set_initial_data(const symphas::init_data_type* tdata_set, size_t n)
{
	if (n <= len && n > 0)
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
		fprintf(SYMPHAS_ERR, "the number of given initial data is larger than the number of systems\n");
	}
}

void symphas::problem_parameters_type::set_interval_data(const symphas::interval_data_type* vdata_set, size_t n)
{
	if (n <= len && n > 0)
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
		fprintf(SYMPHAS_ERR, "the number of given interval data is larger than the number of systems\n");
	}
}

void symphas::problem_parameters_type::set_boundary_data(const symphas::b_data_type* bdata_set, size_t n)
{
	if (n <= len && n > 0)
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
		fprintf(SYMPHAS_ERR, "the number of given boundary data is larger than the number of systems\n");
	}
}

void symphas::problem_parameters_type::set_num_fields(len_type* num_fields_set, size_t n)
{
	if (n <= num_fields_len && n > 0)
	{
		for (size_t i = 0; i < n; ++i)
		{
			num_fields[i] = num_fields_set[i];
		}
		for (size_t i = n; i < num_fields_len; ++i)
		{
			num_fields[i] = num_fields_set[0];
		}
		len = std::reduce(num_fields, num_fields + num_fields_len);
	}
	else
	{
		fprintf(SYMPHAS_ERR, "the number of given field lengths is larger than the number of lengths\n");
	}
}


//! Provide the time step used in the solution of the problem.
void symphas::problem_parameters_type::set_time_step(double dt, double time)
{
	dt_list.set_time_step(dt, time);
}


//! Provide the time step used in the solution of the problem.
void symphas::problem_parameters_type::set_time_steps(const double* dts, const double* t_dts, size_t dts_len)
{
	dt_list.set_time_steps(dts, t_dts, dts_len);
}

void symphas::problem_parameters_type::clear_time_steps(double default_dt)
{
	dt_list.clear_time_steps(default_dt);
}

void symphas::problem_parameters_type::equalize_discretization(const len_type* dims)
{
	for (iter_type i = 0; i < len; ++i)
	{
		size_t dim = vdata[i].size();
		for (iter_type n = 0; n < dim; ++n)
		{
			auto& interval = vdata[i].at(symphas::index_to_axis(n));

			if (interval.get_count() != dims[n])
			{
				interval.set_count(interval.right(), interval.left(), dims[n]);
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
		dims[n] = vdata[i].at(symphas::index_to_axis(n)).get_count();
	}

	equalize_discretization(dims);
	delete[] dims;
}

symphas::problem_parameters_type::~problem_parameters_type()
{
	delete[] tdata;
	delete[] vdata;
	delete[] bdata;
	delete[] num_fields;
	delete[] modifiers;
}


void symphas::swap(symphas::problem_parameters_type& first, symphas::problem_parameters_type& second)
{
	using std::swap;
	swap(first.tdata, second.tdata);
	swap(first.vdata, second.vdata);
	swap(first.bdata, second.bdata);
	swap(first.dt_list, second.dt_list);
	swap(first.num_fields, second.num_fields);
	swap(first.modifiers, second.modifiers);
	swap(first.len, second.len);
	swap(first.num_fields_len, second.num_fields_len);
	swap(first.time, second.time);
	swap(first.index, second.index);
}


void swap(symphas::problem_parameters_type& first, symphas::problem_parameters_type& second)
{
    symphas::swap(first, second);
}



