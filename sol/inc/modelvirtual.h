
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


namespace symphas::internal
{
	template<template<typename, size_t> typename G, typename T, size_t D>
	void update_for_regional(PhaseFieldSystem<G, T, D>& system, symphas::grid_info& info)
	{
		for (auto& [axis, interval] : info)
		{
			interval.domain_to_interval();
		}
	}

	template<typename T, size_t D>
	void update_for_regional(PhaseFieldSystem<RegionalGrid, T, D>& system, symphas::grid_info& info)
	{
		grid::resize_adjust_region(system, info);

		for (auto& [axis, interval] : info)
		{
			double offset = system.region.boundary_size * interval.width();
			system.info[axis].set_interval(interval.left() + offset, interval.right() - offset);

			interval.domain_to_interval();
		}
		info.update_strides();
	}
}


//! Mimics solver functionality and reads existing datafiles instead.
/*!
 * An object which mimics the functionality of a solver, but data is sourced
 * from existing solution data. The data files are taken from the initial data 
 * parameters, which specify that initial conditions are loaded from a file. 
 * 
 * When this is used in a virtual model, the save parameters must be
 * initialized.
 */
template<size_t D, typename Sp = void>
struct DataStepper
{
	template<typename T>
	using SolverSystemApplied = typename symphas::solver_system_type<Sp>::template type<T, D>;

	template<typename T>
	struct solver_system_impl : SolverSystemApplied<T>
	{
		solver_system_impl(SolverSystemApplied<T> const& s) : SolverSystemApplied<T>(s) {}
	};


	template<typename T>
	T _get_solver_system_type(SolverSystemApplied<T>) { return {}; }
	template<typename S>
	auto get_solver_system_type(S s) { return _get_solver_system_type(s); }

	template<typename S>
	struct solver_implicit_type
	{
		using type = std::invoke_result_t<decltype(&DataStepper<D, Sp>::get_solver_system_type<S>), DataStepper<D, Sp>, S>;
	};

	template<typename S>
	using solver_implicit_t = typename solver_implicit_type<S>::type;

	using id_type = Solver<Sp>;

	symphas::init_entry_type* data;
	len_type len;
	SaveParams save;
	double dt;

	DataStepper(len_type len, SaveParams const& save = {}) :
		data{ (len > 0) ? new symphas::init_entry_type[len]{} : nullptr }, len{ len }, save{ save }, dt{} {}


	friend void swap(DataStepper& first, DataStepper& second)
	{
		using std::swap;
		swap(first.save, second.save);
		swap(first.data, second.data);
	}


	template<typename S>
	void step(S& system)
	{
		iter_type next_index = save.next_save(data[system.id].file.get_index());

		symphas::io::read_info rinfo{
			next_index, system.id,
				data[system.id].file.get_name(), data[system.id].in == Inside::CHECKPOINT };
		rinfo.set_offset(false);

		symphas::grid_info ginfo(system.info);
		iter_type read_index = symphas::io::read_grid<solver_implicit_t<S>>(rinfo, &ginfo);
		symphas::internal::update_for_regional(system, ginfo);

		symphas::io::read_grid(system.values, rinfo, &ginfo);
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

	static DataStepper<D, Sp> make_solver(symphas::problem_parameters_type const& pp)
	{
		DataStepper<D, Sp> ds(pp.length());
		for (iter_type n = 0; n < ds.len; ++n)
		{
			if (pp.get_initial_data()[n].at(Axis::NONE).in == Inside::CHECKPOINT
				|| pp.get_initial_data()[n].at(Axis::NONE).in == Inside::FILE)
			{
				ds.data[n] = pp.get_initial_data()[n].at(Axis::NONE);
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


//
//template<typename T>
//struct symphas::solver_system_type<DataStepper<D>>
//{
//	template<typename Ty, size_t D>
//	using type = typename symphas::internal::solver_system_type_match<
//		typename T::id_type,
//		internal::solver_system_type_index<T>
//	>::template type<Ty, D>;
//};

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
template<size_t D, typename Sp = void, typename... S>
struct ModelVirtual : Model<D, DataStepper<D, Sp>, S...>
{
	using parent_type = Model<D, DataStepper<D, Sp>, S...>;
	using parent_type::solver;
	using parent_type::index;

	using parent_type::parent_type;

	auto modify_tdata(symphas::init_data_type& load, const char* dir, iter_type start_index, bool checkpoint) const
	{
		load[Axis::NONE].in = (checkpoint) ? Inside::CHECKPOINT : Inside::FILE;
		load[Axis::NONE].file = { dir, start_index };
	}

	symphas::problem_parameters_type build_parameters(const char* dir, iter_type start_index, bool checkpoint = true)
	{
		symphas::io::read_info rinfo{ start_index, 0, dir, checkpoint };
		symphas::grid_info ginfo = symphas::io::read_header(rinfo);

		if (ginfo.dimension() > 0)
		{
			symphas::problem_parameters_type pp(1);

			symphas::init_data_type load;
			modify_tdata(load, dir, start_index, checkpoint);
			pp.set_initial_data(&load, 1);

			symphas::b_data_type bdata;
			for (iter_type n = 0; n < D * 2; ++n)
			{
				bdata[symphas::index_to_side(n)] = BoundaryType::NONE;
			}
			pp.set_boundary_data(&bdata, 1);

			while (ginfo.dimension() > 0)
			{
				pp.extend(rinfo.get_id() + 1);
				pp.set_interval_data(ginfo.intervals, iter_type(rinfo.get_id()));
				++rinfo.get_id();

				ginfo = symphas::io::read_header(rinfo);
			}

			return pp;
		}
		else
		{
			return symphas::problem_parameters_type(0);
		}
	}

	symphas::problem_parameters_type build_parameters(const char* dir, symphas::problem_parameters_type const& parameters, bool checkpoint = true)
	{
		auto pp(parameters);

		for (auto& [tdata, vdata, bdata] : pp)
		{
			modify_tdata(tdata, dir, params::start_index, checkpoint);
		}
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
	ModelVirtual(const char* solution_dir, SaveParams save, bool checkpoint = true) :
		parent_type(nullptr, 0, build_parameters(solution_dir, save.get_start(), checkpoint))
	{
		solver.set_save_object(save);
		index = save.get_start();
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
	ModelVirtual(const char* solution_dir, iter_type step_index, iter_type start_index, iter_type stop_index, bool checkpoint = true) :
		ModelVirtual(solution_dir, SaveParams{ step_index, start_index, stop_index }, checkpoint) {}

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
	ModelVirtual(const char* solution_dir, iter_type step_index, iter_type stop_index, bool checkpoint = true) :
		ModelVirtual(solution_dir, step_index, stop_index, params::start_index, checkpoint) {}

	ModelVirtual(const char* solution_dir, symphas::problem_parameters_type const& parameters, bool checkpoint = true) :
		parent_type(nullptr, 0, build_parameters(solution_dir, parameters, checkpoint))
	{
		solver.set_save_object(SaveParams(index));
	}

	ModelVirtual(const char* solution_dir, bool checkpoint = true) :
		parent_type(nullptr, 0, build_parameters(solution_dir, checkpoint)) {}

	//! Advances to the next solution iteration.
	/*!
	 * Calls the `step` function of the solver, stepping forward the solution
	 * iteration. This also increments Model::lastindex.
	 *
	 * \param dt The time increment. Has no effect for this model.
	 */
	void step(double dt)
	{
		parent_type::step(dt);
		index = solver.index();
	}

	void update(double time)
	{
		parent_type::update_systems(time);
	}

	void equation() {}
};


template<>
struct ModelVirtual<0> 
{
	iter_type start;
	iter_type end;
	char* dir;

	ModelVirtual(const char* dir, iter_type start, iter_type end) :
		start{ start }, end{ end }, dir{ (dir && std::strlen(dir) > 0) ? new char[std::strlen(dir) + 1] {} : nullptr }
	{
		if (this->dir)
		{
			std::strcpy(this->dir, dir);
		}
	}

	ModelVirtual(const char* dir, iter_type start) : ModelVirtual(dir, start, start) {}
	ModelVirtual(const char* dir) : ModelVirtual(dir, 0, 0) {}
	ModelVirtual() : ModelVirtual(nullptr, 0, 0) {}

	ModelVirtual(ModelVirtual const& other) : ModelVirtual(other.dir, other.start, other.end) {}
	ModelVirtual(ModelVirtual&& other) : ModelVirtual() 
	{
		swap(*this, other);
	}

	ModelVirtual& operator=(ModelVirtual other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(ModelVirtual& first, ModelVirtual& second)
	{
		using std::swap;
		swap(first.start, second.start);
		swap(first.end, second.end);
		swap(first.dir, second.dir);
	}
};


template<size_t D, typename... S>
using ModelVirtualBasic = ModelVirtual<D, void, S...>;

template<size_t D, typename Sp, typename... S>
using ArrayModelVirtual = ModelVirtual<D, Sp, symphas::internal::field_array_t<void>, S...>;

template<size_t D, typename... S>
using ArrayModelVirtualBasic = ArrayModelVirtual<D, void, S...>;


#endif




