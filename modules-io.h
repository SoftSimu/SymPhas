
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
 * This file is part of the SymPhas API. It is used when the io module is 
 * enabled.
 *
 * ***************************************************************************
 */

// header used when io module is used

#pragma once

#include "prereq.h"
#include "write.h"
#include "read.h"
#include "modelvirtual.h"


#ifdef USING_PROC
#include "proc.h"
#endif

namespace symphas
{


	/* opens the record file for the simulation statistics and logging
	 * this is the backup configuration; data is appended to the end with a
	 * specified prefix depending on the data added
	 */


	 //! Obtain the file name in which simulation updates are published.
	 /*!
	  * Writes the record file name for the simulation statistics and logging.
	  *
	  * \param The directory of the record file.
	  */
	const char* get_record_name(const char* dir);

	//! Open the file in which simulation updates are published.
	/*!
	 * Opens the record file for the simulation progression.
	 * This is the backup configuration; data is appended to the end with a
	 * specified prefix depending on the type of data added.
	 *
	 * \param name The name of the record file.
	 */
	FILE* open_record(const char* name);

	//! Write the save index to the record file.
	/*!
	 * Opens the record file and appends the idnex to a new line. If the index
	 * matches the final save index, then additional information is added.
	 *
	 * \param dir The directory of the record file.
	 * \param save The save parameters of the simulation.
	 * \param index The current index of the simulation.
	 */
	void record_index(const char* dir, SaveParams const& save, int index);


	//! Write the checkpoint for the given system.
	/*!
	 * Create the checkpoint of the current solution data found in the given
	 * system.
	 *
	 * \param sys The system containing the solution data.
	 * \param dir The directory to which to save the checkpoint. This directory
	 * must exist.
	 * \param index The index of the solution progression.
	 */
	template<template<typename, size_t> typename G, typename T, size_t D>
	void checkpoint_system(PhaseFieldSystem<G, T, D> const& sys, const char* dir, iter_type index)
	{
		symphas::io::write_info w{ dir, index, sys.get_id(), DataFileType::CHECKPOINT_DATA };

		symphas::interval_data_type intervals;
		for (iter_type i = 0; i < D; ++i)
		{
			symphas::interval_element_type interval;
			interval.set_count(
				sys.get_info().intervals.at(symphas::index_to_axis(i)).left(),
				sys.get_info().intervals.at(symphas::index_to_axis(i)).right(),
				sys.dims[i]);
			intervals[symphas::index_to_axis(i)] = interval;
		}

		symphas::grid_info g{ intervals };
		symphas::io::save_grid(sys.values, w, g);
	}

	template<typename... Ts, size_t... Is>
	void checkpoint_systems(std::tuple<Ts...> const& sys, const char* dir, iter_type index, std::index_sequence<Is...>)
	{
		((checkpoint_system(std::get<Is>(sys), dir, index), ...));
	}

	//! Write the checkpoint for the list of systems.
	/*!
	 * Create the checkpoint of the current solution data found in the given
	 * systems. Each system will be individually saved.
	 *
	 * \param sys The system containing the solution data.
	 * \param dir The directory to which to save the checkpoint. This directory
	 * must exist.
	 * \param index The index of the solution progression.
	 */
	template<typename... Ts>
	void checkpoint_systems(std::tuple<Ts...> const& sys, const char* dir, iter_type index)
	{
		checkpoint_systems(sys, dir, index, std::make_index_sequence<sizeof...(Ts)>{});
	}

	//! Write the checkpoint for the list of systems.
	/*!
	 * Create the checkpoint of the current solution data found in the given
	 * systems. Each system will be individually saved.
	 *
	 * \param sys The system containing the solution data.
	 * \param dir The directory to which to save the checkpoint. This directory
	 * must exist.
	 * \param index The index of the solution progression.
	 */
	template<typename S>
	void checkpoint_systems(std::pair<S*, len_type> const& sys, const char* dir, iter_type index)
	{
		auto [systems, len] = sys;
		for (iter_type i = 0; i < len; ++i)
		{
			checkpoint_system(systems[i], dir, index);
		}
	}

	//! Save the checkpoint if the checkpoint index is reached.
	/*!
	 * If the solution index of the model corresponds to a checkpoint index,
	 * then the model phase fields are saved as a checkpoint.
	 *
	 * \param model The phase field problem data that is saved as a
	 * checkpoint.
	 * \param model_save The save parameters.
	 * \param dir The directory that data is persisted to.
	 */
	template<typename M>
	void save_if_checkpoint(M const& model, SaveParams const& model_save, const char* dir)
	{
		if (symphas::io::at_checkpoint(model.get_index(), model_save))
		{
			checkpoint_systems(model.systems_tuple(), dir, model.get_index());
			record_index(dir, model_save, model.get_index());
		}
	}

	//! Determine the solution of a phase field problem.
	/*!
	 * The phase field problem in the given model is run through the solver
	 * to compute the solution up to the index given by the save parameters.
	 * The function will return true when the index of the model after
	 * running the solution has not reached the stop index. This indicates
	 * that the model should still be solved.
	 *
	 * \param model The phase field problem data.
	 * \param save The save parameters.
	 * \param starttime The begin time of this solution.
	 * \param dt The time step between solution iterations.
	 */
	template<typename M>
	bool run_model(M& model, SaveParams const& save, double dt, double starttime = 0)
	{
		iter_type n = save.next_save(model.get_index()) - model.get_index();
		bool done = run_model(model, n, dt, starttime);
		return (done && !save.is_last_save(model.get_index()));
	}

	template<size_t D, typename... S>
	bool run_model(ModelVirtual<D, S...>& model)
	{
		return run_model(model, model.solver.save, 0);
	}



	//! Get the solution of the model using a standardized workflow.
	/*!
	 * Makes checkpoints and performs
	 * any postprocessing which is necessary (and if the module is enabled).
	 *
	 * \param model The model which is solved.
	 * \param dt The time stepping distance.
	 * \param save Save parameters dictating when the model is paused for
	 * persisting data.
	 * \param save_points Additional save points, which should correspond
	 * to postprocessing saves.
	 * \param dir The save directory.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M* models, len_type num_models, double dt, SaveParams const& save, std::vector<iter_type> save_points, const char* dir,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		M& model = *models;
		model.print_info(SYMPHAS_LOG);

		if (checkpoint)
		{
			fprintf(SYMPHAS_LOG, "Checkpoints are on.\n");
		}
		if (plotting_output)
		{
			fprintf(SYMPHAS_LOG, "Plotting data will be generated.\n");
		}
		if (!checkpoint && !plotting_output)
		{
			fprintf(SYMPHAS_LOG, "No output is generated.\n");
		}
		if (num_models > 1)
		{
			fprintf(SYMPHAS_LOG, "Running concurrent simulations of %d models.\n", num_models);
		}
		fprintf(SYMPHAS_LOG, OUTPUT_BANNER);


		double run_time = 0;
		double iteration_display_time = 0;
		auto print_progress = [&]()
		{
			double progress = model.get_index() / (0.01 * save.get_stop());
			fprintf(SYMPHAS_LOG, "[%5.1lf%%]", progress);
			fprintf(SYMPHAS_LOG, "%11d", model.get_index());
			fprintf(SYMPHAS_LOG, "%8s+%10.3lf", "", iteration_display_time);
			fprintf(SYMPHAS_LOG, "\n");

		};


		auto data_persistence = [&]()
		{
			if (checkpoint)
			{
				fprintf(SYMPHAS_LOG, "%25s", "checkpoint...");
				save_if_checkpoint(model, save, dir);
				fprintf(SYMPHAS_LOG, "done\n");
			}

			if (plotting_output)
			{
				if (save.current_save(model.get_index()))
				{
					fprintf(SYMPHAS_LOG, "%25s", "phase field output...");
					model.save_systems(dir);
					fprintf(SYMPHAS_LOG, "done\n");
				}

#ifdef USING_PROC
				if (model.get_index() > params::start_index
					|| (model.get_index() == params::start_index && save.get_init_flag()))
				{
					symphas::Time t;
					fprintf(SYMPHAS_LOG, "%25s", "postprocessing...");

					bool ran_procs = false;
					if (num_models == 1)
					{
						ran_procs = symphas::proc::collect_data(model, dir);
					}
					else
					{
						ran_procs = symphas::proc::collect_data_avg(models, dir, num_models);
					}

					if (ran_procs)
					{
						fprintf(SYMPHAS_LOG, "done (%.2lfs)\n", t.current_duration());
					}
					else
					{
						fprintf(SYMPHAS_LOG, "nothing to do\n");
					}
				}
#endif
			}

		};

		if (!symphas::io::is_last_save(model.get_index(), save_points))
		{
			data_persistence();
			print_progress();

			fprintf(SYMPHAS_LOG, "%7s", "Progress");
			fprintf(SYMPHAS_LOG, "%11s ", "Index");
			fprintf(SYMPHAS_LOG, "%19s\n", "Runtime (seconds)");

			bool running = true;
			for (iter_type i = 0; i < num_models; ++i)
			{
				models[i].update(TIME_INIT);
			}
			while (running)
			{
				{
					symphas::Time t;

					int next_model_save = save.next_save(model.get_index());
					int next_list_save = symphas::io::next_save_from_list(model.get_index(), save_points);
					int n = std::min(next_model_save, next_list_save) - model.get_index();

					for (iter_type i = 0; i < num_models; ++i)
					{
						running = run_model(models[i], n, dt, models[i].get_time() + dt) && running;
					}
					running = !save.is_last_save(model.get_index()) && running;

					run_time += t.current_duration();
					iteration_display_time = t.current_duration();
				
				}

				print_progress();
				data_persistence();
			}
		}
		else
		{
			data_persistence();
		}

		fprintf(SYMPHAS_LOG, "Completed %d iterations in %.2lf seconds.\n", model.get_index(), run_time);
		fprintf(SYMPHAS_LOG, OUTPUT_BANNER);
	}


	template<typename M>
	void find_solution(M& model, double dt, SaveParams const& save, std::vector<iter_type> save_points, const char* dir,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		find_solution(&model, 1, dt, save, save_points, dir, plotting_output, checkpoint);
	}

	//! See symphas::find_solution().
	/*!
	 * Overload which passes an empty list of save points. Thus the save will
	 * be dictated only by the given save parameter.
	 * 
	 * \param models The phase field problems which are solved.
	 * \param num_models The number of models in the pointer array.
	 * \param dt The time step increment of the solution.
	 * \param save The save parameters for persisting the solution.
	 * \param dir The directory to put the persistent data.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M* models, len_type num_models, double dt, SaveParams const& save, const char* dir,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		find_solution(models, num_models, dt, save, { save.get_stop() }, dir, plotting_output, checkpoint);
	}

	//! See symphas::find_solution().
	/*!
	 * Overload which passes an empty list of save points. Thus the save will
	 * be dictated only by the given save parameter.
	 * 
	 * \param model The phase field problem which is solved.
	 * \param dt The time step increment of the solution.
	 * \param save The save parameters for persisting the solution.
	 * \param dir The directory to put the persistent data.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M& model, double dt, SaveParams const& save, const char* dir,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		find_solution(model, dt, save, { save.get_stop() }, dir, plotting_output, checkpoint);
	}


	//! See symphas::find_solution().
	/*!
	 * Overload without directory specification. The directory is given to be
	 * the current directory.
	 * 
	 * \param models The phase field problems which are solved.
	 * \param num_models The number of models in the pointer array.
	 * \param dt The time step increment of the solution.
	 * \param save The save parameters for persisting the solution.
	 * \param save_points The list of points at which the solution will be
	 * persisted.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M* models, len_type num_models, double dt, SaveParams const& save, std::vector<iter_type> save_points,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		find_solution(models, num_models, dt, save, save_points, ".", plotting_output, checkpoint);
	}

	//! See symphas::find_solution().
	/*!
	 * Overload without directory specification. The directory is given to be
	 * the current directory.
	 * 
	 * \param model The phase field problem which is solved.
	 * \param dt The time step increment of the solution.
	 * \param save The save parameters for persisting the solution.
	 * \param save_points The list of points at which the solution will be
	 * persisted.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M& model, double dt, SaveParams const& save, std::vector<iter_type> save_points,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		find_solution(model, dt, save, save_points, ".", plotting_output, checkpoint);
	}


	//! See symphas::find_solution().
	/*!
	 * Overload which passes an empty list of save points and no directory. The
	 * directory is given to be the current directory.
	 * 
	 * \param models The phase field problems which are solved.
	 * \param num_models The number of models in the pointer array.
	 * \param dt The time step increment of the solution.
	 * \param save The save parameters for persisting the solution.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M* models, len_type num_models, double dt, SaveParams const& save,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		find_solution(models, num_models, dt, save, ".", plotting_output, checkpoint);
	}

	//! See symphas::find_solution().
	/*!
	 * Overload which passes an empty list of save points and no directory. The
	 * directory is given to be the current directory.
	 * 
	 * \param model The phase field problem which is solved.
	 * \param dt The time step increment of the solution.
	 * \param save The save parameters for persisting the solution.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M& model, double dt, SaveParams const& save,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		find_solution(model, dt, save, ".", plotting_output, checkpoint);
	}



	//! See symphas::find_solution().
	/*!
	 * Takes the number of iterations to apply the solver.
	 * The directory is not specified and is set to the default value of
	 * current directory. Solves multiple models simultaneously. Only one model
	 * is saved.
	 *
	 * \param models The phase field problems which are solved.
	 * \param num_models The number of models in the pointer array.
	 * \param dt The time step increment of the solution.
	 * \param stop The terminating index.
	 * \param dir The directory to put the persistent data.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M* models, len_type num_models, double dt, iter_type stop, const char* dir,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		SaveParams save{ SaveType::DEFAULT, stop, 0, stop };
		find_solution(models, num_models, dt, save, { save.get_stop() }, dir, plotting_output, checkpoint);
	}

	//! See symphas::find_solution().
	/*!
	 * Takes the number of iterations to apply the solver.
	 * The directory is not specified and is set to the default value of
	 * current directory. Solves multiple models simultaneously. Only one model
	 * is saved.
	 *
	 * \param model The phase field problem which is solved.
	 * \param dt The time step increment of the solution.
	 * \param stop The terminating index.
	 * \param dir The directory to put the persistent data.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M& model, double dt, iter_type stop, const char* dir,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		SaveParams save{ SaveType::DEFAULT, static_cast<double>(stop), 0, stop };
		find_solution(model, dt, save, { save.get_stop() }, dir, plotting_output, checkpoint);
	}



	//! See symphas::find_solution().
	/*!
	 * Takes the number of iterations to apply the solver.
	 * The directory is not specified and is set to the default value of
	 * current directory. Solves multiple models simultaneously. Only one model
	 * is saved.
	 * 
	 * \param models The phase field problems which are solved.
	 * \param num_models The number of models in the pointer array.
	 * \param dt The time step increment of the solution.
	 * \param stop The terminating index.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M* models, len_type num_models, double dt, iter_type stop,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		find_solution(models, num_models, dt, stop, ".", plotting_output, checkpoint);
	}

	//! See symphas::find_solution().
	/*!
	 * Takes the number of iterations to apply the solver.
	 * The directory is not specified and is set to the default value of 
	 * current directory.
	 * 
	 * \param model The phase field problem which is solved.
	 * \param dt The time step increment of the solution.
	 * \param stop The terminating index.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M& model, double dt, iter_type stop,
		bool plotting_output = true, bool checkpoint = symphas::io::is_checkpoint_set())
	{
		find_solution(model, dt, stop, ".", plotting_output, checkpoint);
	}
}