
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
 * This file is part of the SymPhas API. It is used when the configuration
 * module is enabled.
 *
 * ***************************************************************************
 */

// header used when configuration headers are included

#pragma once

#include "modules-io.h"
#include "conf.h"
#include "model.h"


namespace symphas
{

	//! Execute a solution to the phase field problem from its given state.
	/*!
	 * The solution will continue to run until the next save index
	 * is reached. The function will return true until the model is complete.
	 *
	 * \param model The phase field problem data.
	 * \param config The configuration dictating parameters of the solution
	 * and from where save points are derived.
	 */
	template<typename M>
	bool run_model(M& model, Conf const& config)
	{
		iter_type n = symphas::conf::next_save(model.index(), config) - model.index();
		run_model(model, n, config.get_time_step_list());

		return !(symphas::conf::is_last_save(model.index(), config));
	}

	//! Execute a solution to the phase field problem from its given state.
	/*!
	 * The solution will continue to run until the next save index
	 * is reached. A truth value representing whether the solution has been
	 * entirely completed is returned. The global configuration is used
	 * to determine parameters of the solution.
	 *
	 * \param model The phase field problem data.
	 */
	template<typename M>
	bool run_model(M& model)
	{
		return run_model(model, symphas::conf::config());
	}




	//! Get the solution of the model using a standardized workflow.
	/*!
	 * Uses the global configuration, makes checkpoints and
	 * saves output.
	 * 
	 * \param models The phase field problems which are solved.
	 * \param num_models The number of models in the pointer array.
	 * \param config The configuration parameters of the solution.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M* models, len_type num_models, Conf const& config, bool plotting_output = true, bool checkpoint = true)
	{
		auto save_points = symphas::conf::get_save_indices(config);
		find_solution(models, num_models, config.get_time_step_list(), config.save, save_points, config.get_result_dir(), plotting_output, checkpoint);
	}

	//! Get the solution of the model using a standardized workflow.
	/*!
	 * Uses the global configuration, makes checkpoints and
	 * saves output.
	 * 
	 * \param model The phase field problem which is solved.
	 * \param config The configuration parameters of the solution.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M& model, Conf const& config, bool plotting_output = true, bool checkpoint = true)
	{
		auto save_points = symphas::conf::get_save_indices(config);
		find_solution(model, config.get_time_step_list(), config.save, save_points, config.get_result_dir(), plotting_output, checkpoint);
	}


	//! Get the solution of the model using a standardized workflow.
	/*!
	 * Uses the global configuration, and if selected, makes checkpoints and
	 * saves output. Solves multiple models simultaneously. Only one model
	 * is saved.
	 *
	 * \param models The phase field problems which are solved.
	 * \param num_models The number of models in the pointer array.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M* models, len_type num_models, bool plotting_output = true, bool checkpoint = true)
	{
		find_solution(models, num_models, symphas::conf::config(), plotting_output, checkpoint);
	}

	//! Get the solution of the model using a standardized workflow.
	/*!
	 * Uses the global configuration, and if selected, makes checkpoints and
	 * saves output.
	 *
	 * \param model The phase field problem which is solved.
	 * \param plotting_output If true, output will be generated to be used
	 * by a plotting utility.
	 * \param checkpoint If true, checkpoints will be saved, which are raw
	 * data outputs of the phase field.
	 */
	template<typename M>
	void find_solution(M& model, bool plotting_output = true, bool checkpoint = true)
	{
		find_solution(model, symphas::conf::config(), plotting_output, checkpoint);
	}

	namespace io
	{

		//! Uses configuration data to write the plot configuration.
		/*!
		 * Calls the existing function to write the plotting configuration for the given model,
		 * passing parameters from the configuration. See: write_plot_config().
		 * Uses the given names of the phase fields.
		 *
		 * \param model The model for which to write the plotting configuration for.
		 */
		template<typename M>
		void write_plot_config(M const& model, const char* directory, const char* const* names)
		{
			write_plot_config(model, directory, names, symphas::conf::config().save);
		}

		//! Uses configuration data to write the plot configuration.
		/*!
		 * Calls the existing function to write the plotting configuration for the given model,
		 * passing parameters from the configuration. See: write_plot_config().
		 * Uses the given names of the phase fields.
		 *
		 * \param model The model for which to write the plotting configuration for.
		 */
		template<typename M>
		void write_plot_config(M const& model, const char* const* names)
		{
			write_plot_config(model, symphas::conf::config().get_result_dir(), names, symphas::conf::config().save);
		}


		//! Uses configuration data to write the plot configuration.
		/*!
		 * Calls the existing function to write the plotting configuration for the given model,
		 * passing parameters from the configuration. See: write_plot_config().
		 * For the phase field names, substitutes the configuration defined
		 * names.
		 *
		 * \param model The model for which to write the plotting configuration for.
		 */
		template<typename M>
		void write_plot_config(M const& model, const char* directory, SaveParams const& save)
		{
			len_type len = symphas::model_num_fields(model);
			char** names = new char* [len] {};

			iter_type i = 0;
			for (auto name : symphas::conf::config().get_names(len))
			{
				names[i] = nullptr;
				std::swap(names[i], name.data);
				++i;
			}

			symphas::io::write_plot_config(model, directory, names, symphas::conf::config().save);
			for (iter_type i = 0; i < len; ++i)
			{
				delete[] names[i];
			}
			delete[] names;
		}


		//! Uses configuration data to write the plot configuration.
		/*!
		 * Calls the existing function to write the plotting configuration for the given model,
		 * passing parameters from the configuration. See: write_plot_config().
		 * Uses the given names of the phase fields.
		 *
		 * \param model The model for which to write the plotting configuration for.
		 */
		template<typename M>
		void write_plot_config(M const& model, const char* directory)
		{
			write_plot_config(model, directory, symphas::conf::config().save);
		}

		//! Uses configuration data to write the plot configuration.
		/*!
		 * Calls the existing function to write the plotting configuration for the given model,
		 * passing parameters from the configuration. See: write_plot_config().
		 * Uses the given names of the phase fields.
		 *
		 * \param model The model for which to write the plotting configuration for.
		 */
		template<typename M>
		void write_plot_config(M const& model, SaveParams const& save)
		{
			write_plot_config(model, symphas::conf::config().get_result_dir(), save);
		}

		//! Uses configuration data to write the plot configuration.
		/*!
		 * Calls the existing function to write the plotting configuration for the given model,
		 * passing parameters from the configuration. See: write_plot_config().
		 * Uses the given names of the phase fields.
		 *
		 * \param model The model for which to write the plotting configuration for.
		 */
		template<typename M>
		void write_plot_config(M const& model)
		{
			write_plot_config(model, symphas::conf::config().get_result_dir());
		}

	}
}
