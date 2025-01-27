
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
 * This file is part of the SymPhas API. This is the main header file that is
 * configured by the CMake build process.
 *
 * ***************************************************************************
 */

#pragma once

#include "definitions.h"

#ifdef USING_CONF
#include "modules-conf.h"
#elif defined(USING_IO)
#include "modules-io.h"
#else
#include "modules-none.h"
#endif

#include "modules-models.h"


#ifdef MODULES_EXPORT
#define DLLMOD DLLEXPORT
#else
#define DLLMOD DLLIMPORT
#endif


 // include all the models
#ifdef MODEL_TESTS
#include "testdefs.h"
#endif


#include "timer.h"


namespace symphas
{

#ifdef PRINT_TIMINGS
	DLLMOD extern double init_time;
	DLLMOD extern double model_update_time;
	DLLMOD extern double model_equation_time;
	DLLMOD extern double model_step_time;
	DLLMOD extern double iteration_time;

	DLLMOD extern int iteration_count;

	void print_timings(FILE* out);
#endif


	//! One iteration of the solution loop.
	/*!
	 * The phase field problem represented by the given model is iterated
	 * through one solution step.
	 *
	 * \param model The phase field problem.
	 * \param dt The time step of this iteration.
	 * \param time The current solution time.
	 */
	template<typename M>
	void model_iteration(M& model, double dt, double time)
	{

#ifdef PRINT_TIMINGS
		symphas::Time t;

		{
			symphas::Time tt;
			model.equation();
			model_equation_time += tt.current_duration();
		}

		{
			symphas::Time tt;
			model.step(dt);
			model_step_time += tt.current_duration();
		}

		{
			symphas::Time tt;
			model.update(time + dt);
			model_update_time += tt.current_duration();
		}


		iteration_time += t.current_duration();
		iteration_count += 1;
#else

		model.update(time);
		model.equation();
		model.step(dt);

#endif

	}

	//! Determine the solution of a phase field problem.
	/*!
	 * The phase field problem in the given model is run through the solver
	 * to compute the solution after the specified number of iterations.
	 * The function will always return true.
	 *
	 * \param model The phase field problem data.
	 * \param n The number of iterations of the solution to compute.
	 * \param starttime The begin time of this solution.
	 * \param dt The time step between solution iterations.
	 */
	template<typename M>
	bool run_model(M& model, iter_type n, symphas::time_step_list const& dts, double starttime = 0)
	{
		double time = starttime;
		iter_type end = model.get_index() + n;

		for (iter_type i = model.get_index(); i < end; i = model.get_index())
		{
			double dt = dts.get_time_step(time);
			model_iteration(model, dt, time);
			time += dt * (model.get_index() - i);
		}
		return true;
	}

	//! Initialize the program parameters.
	/*!
	 * Initialize the program parameters and the command line parameters.
	 *
	 * \param config The name of the configuration file, or name of the title
	 * if there no configuration.
	 * \param param_list The list of strings containing the key value pairs
	 * in the format: "key=value" from which command line parameters
	 * are extracted.
	 * \param num_params The number of command line arguments in the list.
	 */
	void init(const char* config, const char* const* param_list, int num_params);
	inline void init(const char* param)
	{
		init(param, nullptr, 0);
	}
	inline void init()
	{
		init("--" ARGUMENT_HELP_STRING);
	}

	void finalize();

}


namespace params
{

	enum program_params_value
	{
		PARAMS
	};

	inline void operator+=(program_params_value, const char* arg)
	{
		params::parse_params(symphas::build_params(), arg);
	}

	inline void operator+=(program_params_value, std::string const& str)
	{
		params::parse_params(symphas::build_params(), str.c_str());
	}

	template<typename T>
	void operator+=(program_params_value, T arg)
	{
		params::set_param(arg);
	}

	template<typename T>
	void operator,(program_params_value, T arg)
	{
		params::set_param(arg);
	}

}

#undef BUFFER_LENGTH_R4
#undef BUFFER_LENGTH_R2
#undef BUFFER_LENGTH_R1
#undef BUFFER_LENGTH
#undef BUFFER_LENGTH_L2
#undef BUFFER_LENGTH_L3
#undef BUFFER_LENGTH_L4
#undef BUFFER_LENGTH_L5
#undef BUFFER_LENGTH_R
#undef BUFFER_LENGTH_L
#undef LINE_READ_BUFFER
#undef FILECHUNK_READ_BUFFER

#undef ERR_CODE_FILE_OPEN
#undef ERR_MSG_FILE_OPEN




