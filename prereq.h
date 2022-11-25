
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
 * This file is part of the SymPhas API. This file includes objects
 * and files required for compilation.
 *
 * ***************************************************************************
 */

#pragma once

#include "prereq-defs.h"
#include "system.h"
#include "solver.h"
#include "stencilincludes.h"

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
	//! Construct a map of all the parameter keys and initialization strategies.
	/*!
	 * The parameter map which specifies the command line key string and the
	 * corresponding initialization strategy for the variable it is supposed
	 * to initialize is constructed and returned.
	 * 
	 * If the _io_ module is enabled, it will also add those parameters.
	 */
	inline param_map_type build_param_map();


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
			model.update(time);
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
	bool run_model(M& model, iter_type n, double dt, double starttime = 0)
	{
		double time = starttime;
		iter_type end = model.get_index() + n;

		for (iter_type i = model.get_index(); i < end; i = model.get_index())
		{
			time += dt * (model.get_index() - i);
			model_iteration(model, dt, time);
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
	void init(const char* config, const char* const* param_list, size_t num_params);
}



// \cond

template<>
struct Solver<GeneralizedStencil<0>> : Solver<Solver<GeneralizedStencil<0>>>
{
	using id_type = void;

	static auto make_solver(...)
	{
		return Solver<GeneralizedStencil<0>>{};
	}

	template<size_t En, typename SS, typename S, typename E>
	auto form_expr_one(SS&&, std::pair<S, E>&& e)
	{
		return e;
	}

	/*
	 * forward euler method for actually updating the grid
	 */
	template<typename S>
	void step(S& sys, double dt) {}


	/*
	 * the parameter is a tuple with the first parameter being a Variable object with ref
	 * base type, and the ref is a reference to the system; like this, the variable index is
	 * packaged and the equality operator is deferred to the oplvariable
	 */
	template<typename S, typename E>
	inline void equation(std::pair<S, E>& r)
	{
		auto& [sys, equation] = r;
		expr::prune::update(equation);
		expr::result(equation, sys.get().frame);
	}
};

// \endcond

//! Placeholder solver.
/*!
 * A placeholder solver which will avoid errors when no concrete solver is implemented
 * or provided in the driver implementation.
 */
template<>
struct Solver<void> : Solver<Solver<GeneralizedStencil<0>>>
{
	using id_type = void;
};


template<>
struct symphas::internal::solver_supported_type_match<void, void>
{
	static const bool value = true;
};


