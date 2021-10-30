
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
 * model for postprocessing purposes or saving as a new format.
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




//! Mimics solver functionality and reads existing datafiles instead.
/*!
 * An object which mimics the functionality of a solver, but data is sourced
 * from existing solution data. The data files are taken from the initial data 
 * parameters, which specify that initial conditions are loaded from a file. 
 * 
 * When this is used in a virtual model, the save parameters must be
 * initialized.
 */
template<size_t D>
struct DataStepper
{
	using id_type = void;

	symphas::init_data_type data[D];
	SaveParams save;

	DataStepper(SaveParams const& save = {}) :
		save{ save }, data{ {} } {}


	friend void swap(DataStepper& first, DataStepper& second)
	{
		using std::swap;
		swap(first.save, second.save);
		swap(first.data, second.data);
	}


	template<typename S>
	void step(S& system, double)
	{
		iter_type next_index = save.next_save(data[system.id].file.get_index());

		iter_type read_index = symphas::io::read_grid(
			system.values,
			symphas::io::read_info{
				next_index, system.id, 
				data[system.id].file.get_name(), data[system.id].in == Inside::CHECKPOINT 
			});
		data[system.id].file.set_index(read_index);
		
	}

	void set_save_object(SaveParams const& save)
	{
		this->save = save;
	}



	iter_type index() const
	{
		return data->file.get_index();
	}

	static DataStepper<D> make_solver(symphas::problem_parameters_type const& pp)
	{
		DataStepper<D> ds;
		for (iter_type n = 0; n < pp.length(); ++n)
		{
			if (pp.get_initial_data()[n].in == Inside::CHECKPOINT
				|| pp.get_initial_data()[n].in == Inside::FILE)
			{
				ds.data[n] = pp.get_initial_data()[n];
			}
			else
			{
				fprintf(SYMPHAS_ERR, "a virtual model can only be initialized "
					"with FILE or CHECKPOINT initial conditions type\n");
			}
		}
		return ds;
	}
};


 //! A model which will read from computed solution data.
 /*!
  * This model takes a precomputed solution and reads all the input data. The
  * purpose is to rewrite all the data to a new simulation directory, following 
  * the usual model workflow, typically in order to write in a different
  * format or run postprocessing on existing data.
  *
  * \tparam D The dimension of the phase field crystal problem.
  * \tparam S... The types of the order parameters.
  */
template<size_t D, typename... S>
struct ModelVirtual : Model<D, DataStepper<D>, S...>
{
	using parent_type = Model<D, DataStepper<D>, S...>;
	using parent_type::solver;
	using parent_type::lastindex;

	using parent_type::parent_type;

	symphas::problem_parameters_type build_parameters(const char* dir, iter_type start_index, bool checkpoint = true)
	{
		symphas::problem_parameters_type pp{ SN };

		symphas::init_data_type load;
		load.in = (checkpoint) ? Inside::CHECKPOINT : Inside::FILE;
		load.file = { dir, start_index };
		pp.set_initial_data(&load, 1);

		for (size_t n = 0; n < SN; ++n)
		{
			symphas::io::read_info rinfo{ start_index, n, dir, checkpoint };
			symphas::grid_info ginfo = symphas::io::read_header(rinfo);
			pp.set_interval_data(ginfo.intervals, static_cast<iter_type>(n));
		}

		symphas::b_data_type bdata;
		for (iter_type n = 0; n < D * 2; ++n)
		{
			bdata[symphas::index_to_side(n)] = BoundaryType::NONE;
		}
		pp.set_boundary_data(&bdata, 1);


		return pp;
	}

	symphas::problem_parameters_type build_parameters(const char* dir, bool checkpoint = true)
	{
		return build_parameters(dir, params::start_index, checkpoint);
	}


	//! Initialize system parameters using the checkpoint file. 
	/*!
	 * Reads the solution data from the given input and initializes the
	 * model system parameters.
	 * 
	 * \param solution_dir The input for the solution files.
	 * \param save The save parameters object for stepping through the solution 
	 * files. Importantly, specifies the frequency of the save intervals for the 
	 * solution files which cannot be inferred directly from the files in
	 * general.
	 * \param checkpoint Indicates the input is a checkpoint of a previous
	 * solution, rather than a plain file.
	 */
	ModelVirtual(const char* solution_input, SaveParams save, bool checkpoint = true) :
		parent_type(nullptr, 0, build_parameters(solution_input, save.get_start(), checkpoint))
	{
		solver.set_save_object(save);
		lastindex = save.start;
	}

	//! Initialize system parameters using the checkpoint file. 
	/*!
	 * Reads the solution data from the given input and initializes the
	 * model system parameters.
	 *
	 * \param solution_dir The input for the solution files.
	 * \param step_index The iteration steps between each saved solution data
	 * frame.
	 * \param stop_index The final index of the solution data.
	 * \param start_index The beginning index of the solution data.
	 * \param checkpoint Indicates the input is a checkpoint of a previous
	 * solution, rather than a plain file.
	 */
	ModelVirtual(const char* solution_input, iter_type step_index, iter_type start_index, iter_type stop_index, bool checkpoint = true) :
		ModelVirtual(solution_input, SaveParams{ step_index, start_index, stop_index }, checkpoint) {}

	//! Initialize system parameters using the checkpoint file. 
	/*!
	 * Reads the solution data from the given input and initializes the
	 * model system parameters.
	 *
	 * \param solution_dir The input for the solution files.
	 * \param step_index The iteration steps between each saved solution data
	 * frame.
	 * \param stop_index The final index of the solution data.
	 * \param checkpoint Indicates the input is a checkpoint of a previous
	 * solution, rather than a plain file.
	 */
	ModelVirtual(const char* solution_input, iter_type step_index, iter_type stop_index, bool checkpoint = true) :
		ModelVirtual(solution_input, step_index, stop_index, params::start_index, checkpoint) {}


	//! Advances to the next solution iteration.
	/*!
	 * Calls the `step` function of the solver, stepping forward the solution
	 * iteration. This also increments Model::lastindex.
	 *
	 * \param dt The time increment. Has no effect for this model.
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




