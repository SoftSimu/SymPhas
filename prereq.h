
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


#ifdef MODULES_EXPORT
#define DLLMOD DLLEXPORT
#else
#define DLLMOD DLLIMPORT
#endif

template<typename T, typename Ty>
struct symphas::internal::solver_supported_type_match
{
	static const bool value = false;
};


/* implementation of suported system, after the specialization of individual supported types
 */
template<typename T, typename S>
struct symphas::internal::solver_supported_system
{

protected:

	template<typename Ty>
	struct supported_wrap
	{
		static const bool value = solver_supported_type<T, void>::value || solver_supported_type<T, Ty>::value;
	};

	template<typename Ty>
	static constexpr auto is_supported(Block<Ty>)
	{
		return supported_wrap<Ty>{};
	}

	static constexpr auto get_value(S s)
	{
		return is_supported(s);
	}

	using wrap_type = typename std::invoke_result_t<decltype(&solver_supported_system<T, S>::get_value), S>;

public:

	static const size_t value = wrap_type::value;
};



template<typename T, typename... S>
struct symphas::internal::solver_supported_systems
{
	static const size_t value = ((solver_supported_system<T, S>::value && ...));
};


template<typename T>
struct symphas::internal::solver_system_type_match
{
	template<typename Ty, size_t D>
	using type = SolverSystem<Ty, D>;
};

template<typename T>
struct symphas::internal::provisional_system_type_match
{
	template<typename Ty, size_t D>
	using type = ProvisionalSystem<Ty, D>;
};



#ifdef USING_PROC
#include "proc.h"
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
		iter_type end = model.index() + n;

		for (iter_type i = model.index(); i < end; i = model.index())
		{
			model_iteration(model, dt, time);
			time += dt * (model.index() - i);
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



