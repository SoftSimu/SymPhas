
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
 * MODULE:  sol
 * PURPOSE: Defines a virtual model which does no work, but is only used
 * to step through already precomputed data and mimic the work of a real
 * model for postprocessing purposes.
 *
 * ***************************************************************************
 */

#pragma once


#ifdef USING_IO

#include "io.h"
#include "model.h"

 //! \cond
#define SN sizeof...(S)
//


//! Mimics functionality of the solver to 
/*!
 * An object which mimics the functionality of a solver, but data is sourced
 * from data files. The data files are taken from the initial data parameters,
 * which specify that initial conditions are loaded from a file. 
 * 
 * When this is used in a virtual model, the save parameters must be
 * initialized. If the save parameters are not initialized, the index 
 */
struct DataStepper
{
	using id_type = void;
	using init_data_read = symphas::init_data_read;

	SaveParams save;
	init_data_read *data;
	size_t len;

	DataStepper(SaveParams const& save, init_data_read* data, size_t len = 1) : 
		save{ save }, data{ (len > 0) ? new init_data_read[len] : nullptr }, len{ len }
	{
		if (this->data)
		{
			std::copy(data, data + len, this->data);
		}
	}
	DataStepper(size_t len) :
		save{ {} }, data{ (len > 0) ? new init_data_read[len] : nullptr }, len{ len } {}
	DataStepper() : DataStepper({}, nullptr, 0) {}

	DataStepper(DataStepper const& other) : DataStepper(other.save, other.data, other.len) {}
	DataStepper(DataStepper&& other) : DataStepper() { swap(*this, other); }
	DataStepper& operator=(DataStepper other) { swap(*this, other); return *this; }

	friend void swap(DataStepper& first, DataStepper& second)
	{
		using std::swap;
		swap(first.save, second.save);
		swap(first.data, second.data);
		swap(first.len, second.len);
	}



	template<typename S>
	void step(S& system, double)
	{
		iter_type next_index = save.next_save(data[system.id].index);
		iter_type read_index = symphas::io::read_grid(
			system.values,
			symphas::io::read_info{ next_index, system.id, data[system.id].name });
		data[system.id].index = read_index;
		
	}

	void set_save_object(SaveParams const& save)
	{
		this->save = save;
	}




	iter_type index() const
	{
		return data->index;
	}

	static DataStepper make_solver(symphas::problem_parameters_type const& pp)
	{
		DataStepper ds{ pp.length() };
		for (iter_type n = 0; n < pp.length(); ++n)
		{
			if (pp.get_initial_data()[n].in == Inside::FILE)
			{
				ds.data[n] = pp.get_initial_data()[n].file;
			}
		}
		return ds;
	}
};


 //! A representation of the phase field crystal model.
 /*!
  * A representation of the phase field crystal model. Implements the phase
  * field crystal model equation, and allows
  *
  * \tparam PFC The phase field crystal parameters specialized object.
  * \tparam D The dimension of the phase field crystal problem.
  * \tparam Sp The solver type for numerically solving the phase field crystal
  * problem.
  * \tparam S... The types of the order parameters.
  */
template<size_t D, typename... S>
struct ModelVirtual : Model<D, DataStepper, S...>
{
	using parent_type = Model<D, DataStepper, S...>;
	using parent_type::solver;
	using parent_type::lastindex;

	using parent_type::parent_type;

	ModelVirtual(symphas::problem_parameters_type const& parameters, SaveParams const& save) :
		parent_type(nullptr, 0, parameters) 
	{
		solver.set_save_object(save);
	}

	//! Advances to the next solution iteration.
	/*!
	 * Calls the `step` function of the solver, stepping forward the solution
	 * iteration. This also increments Moodel::lastindex.
	 *
	 * \param dt The time increment.
	 */
	void step(double dt)
	{
		parent_type::step(std::make_index_sequence<SN>{}, dt);
		lastindex = solver.index();
	}

	void update(double time)
	{
		parent_type::update_systems(time);
	}

	void equation() {}
};


#undef SN

#endif




