
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
 * This file is part of the SymPhas API. It is used when model definitions
 * are included from a separate file.
 *
 * ***************************************************************************
 */

 // header used when model headers are included

#pragma once

#ifdef USING_CONF
#include "conf.h"
#endif

#include "modelmacros.h"
#ifdef MODEL_INCLUDE_HEADER

#include MODEL_INCLUDE_HEADER

#ifdef USING_MODEL_SELECTION

namespace symphas::internal
{

#ifdef USING_CONF

	template<size_t D, typename Sp, typename... S>
	auto run_virtual_model(ModelVirtual<D, Sp, S...>& virt, Conf const& conf)
	{
		len_type len = symphas::model_num_fields(virt);
		char** names = new char* [len] {};

		iter_type i = 0;
		for (auto name : conf.get_names(len))
		{
			names[i] = nullptr;
			std::swap(names[i], name.data);
			++i;
		}

		symphas::io::write_plot_config(virt, conf.get_result_dir(), names, conf.save);
		for (iter_type i = 0; i < len; ++i)
		{
			delete[] names[i];
		}
		delete[] names;


		find_solution(virt, conf);
		return 1;
	}

	template<size_t D, typename Sp, typename... S>
	auto initialize_virtual_model(Model<D, Sp, S...>* model, ModelVirtual<0> const& setup, Conf const& conf)
	{
		return ModelVirtual<D, Sp, S...>(setup.dir, conf.save);
	}

	template<size_t D, typename Sp, typename... S>
	auto initialize_virtual_model(ArrayModel<D, Sp, S...>* model, ModelVirtual<0> const& setup, Conf const& conf)
	{
		auto save(conf.save);
		save.set_start(setup.start);
		return ArrayModelVirtual<D, Sp, S...>(setup.dir, save);
	}

	template<typename M, typename... Ts>
	int handle_simulation_creation(ModelVirtual<0> const& setup, Conf const& conf, Ts&&... args)
	{
		if (params::plots_only)
		{
			M model(conf.get_coeff_list(), conf.get_coeff_len(), conf.get_problem_parameters());
			symphas::io::write_plot_config(model);
			return 1;
		}
		else
		{
			auto virt = initialize_virtual_model((M*)0, setup, conf);
			return run_virtual_model(virt, conf);
		}
	}

#else

	template<size_t D, typename Sp, typename... S>
	auto run_virtual_model(ModelVirtual<D, Sp, S...>& virt, iter_type start, iter_type end)
	{
		symphas::io::write_plot_config(virt, '.', SaveParams(start, end));

		symphas::io::write_plot_config(virt);
		find_solution(virt, 1.0, end - start);
		return 1;
	}

	template<size_t D, typename Sp, typename... S>
	auto initialize_virtual_model(Model<D, Sp, S...>* model, ModelVirtual<0> const& setup)
	{
		return ModelVirtual<D, Sp, S...>(setup.dir);
	}

	template<size_t D, typename Sp, typename... S>
	auto initialize_virtual_model(ArrayModel<D, Sp, S...>* model, ModelVirtual<0> const& setup)
	{
		return ArrayModelVirtual<D, Sp, S...>(setup.dir);
	}

	template<typename M, typename... Ts, size_t D = model_dimension<M>::value>
	int handle_simulation_creation(ModelVirtual<0> const& setup, int, Ts&&... args)
	{
		if (params::plots_only)
		{
			M model(nullptr, 0, {});
			symphas::io::write_plot_config(model);
		}
		else
		{
			auto virt = initialize_virtual_model<D>(model_parameter_types_t<M>{}, setup);
			return run_virtual_model(virt, setup.start, setup.end);
		}
	}

#endif

	template<typename M, typename T0, typename T1, typename... Ts>
	auto handle_simulation_creation(T0 const& arg0, T1 const& arg1, Ts&&... args)
	{
		return MODEL_APPLY_CALL<M>(arg0, arg1, std::forward<Ts>(args)...);
	}

	template<typename M, typename T0>
	auto handle_simulation_creation(T0 const& arg0)
	{
		return MODEL_APPLY_CALL<M>(arg0);
	}

	template<typename M>
	auto handle_simulation_creation()
	{
		return MODEL_APPLY_CALL<M>();
	}

	template<size_t D>
	void print_stencil_message(const char* (&deriv_names)[D], size_t(&stencil_values)[D])
	{
		size_t max_name_len = 0;
		for (iter_type i = 0; i < D; ++i)
		{
			size_t name_len = std::strlen(deriv_names[i]);
			max_name_len = (name_len > max_name_len) ? name_len : max_name_len;
		}

		fprintf(SYMPHAS_LOG, "the finite difference stencils of the solver use the "
			"following point values:\n");

		for (iter_type i = 0; i < D; ++i)
		{
			fprintf(SYMPHAS_LOG, "\t%-*s : %zd\n", static_cast<int>(max_name_len), deriv_names[i], stencil_values[i]);
		}
	}


	template<size_t D>
	void print_stencil_message(const char* (&deriv_names)[D], size_t(&stencil_values)[0])
	{
		fprintf(SYMPHAS_LOG, "no stencil point values are defined\n");
	}

	template<
		template<size_t, typename> typename Model,
		template<typename, size_t = 0> typename Solver,
		size_t N, size_t D, size_t O, size_t... Ps,
		typename... Ts
	>
	auto run_model_call(std::index_sequence<N, D, O, Ps...>, Ts&& ...args)
	{
		const char* names[]{ "laplacian", "gradlaplacian", "bilaplacian" };
		size_t values[]{ Ps... };
		print_stencil_message(names, values);

		using type = typename SelfSelectingStencil<D, O>::template Points<Ps...>;
		return handle_simulation_creation<Model<D, Solver<type, N>>>(std::forward<Ts>(args)...);
	}

	template<
		template<size_t, typename> typename Model,
		template<size_t> typename Solver, size_t N, size_t D, typename... Ts
	>
	auto run_model_call(std::index_sequence<N, D>, Ts&& ...args)
	{
		return handle_simulation_creation<Model<D, Solver<N>>>(std::forward<Ts>(args)...);
	}
}

struct model_select
{
	model_select(size_t dimension, StencilParams stp) :
		dimension{ dimension }, stp{ stp } {}

	size_t dimension;
	StencilParams stp;


#ifdef USING_CONF

	auto load_virtual_conf(const char* spec, iter_type& start, iter_type& end)
	{
		char* dir = new char[std::strlen(spec) + 1] {};
		std::strcpy(dir, spec);
		auto dir_it = std::strchr(dir, VIRTUAL_MODEL_SEP_KEY_C);

		start = 0;
		end = start;
		if (dir_it != NULL)
		{
			*dir_it = '\0';
			sscanf(dir_it + 1, "%d", &start);
			if ((dir_it = std::strchr(dir_it + 1, VIRTUAL_MODEL_SEP_KEY_C)) != NULL)
			{
				sscanf(dir_it + 1, "%d", &end);
			}
		}


		param_map_type param_map = symphas::build_params();
		auto conf = symphas::conf::restore_checkpoint(param_map, dir, start);

#ifdef USING_PROC
		symphas::internal::update_data_stop(conf.save.get_stop());
#endif

		delete[] dir;
		return conf;
	}

	auto load_virtual_conf(const char* name)
	{
		iter_type start, end;
		return load_virtual_conf(name, start, end);
	}


#else

	int load_virtual_conf(const char* name, iter_type& start, iter_type& end)
	{
		start = 0;
		end = 0;
		return {};
	}

	auto load_virtual_conf(const char* name)
	{
		iter_type start, end;
		return load_virtual_conf(name, start, end);
	}

#endif

	inline bool check_virtual_model(const char* name)
	{
		auto sep_it = std::strchr(name, VIRTUAL_MODEL_SEP_KEY_C);
		if (sep_it != NULL)
		{
			char buffer[STR_ARR_LEN(STR(VIRTUAL_MODEL_KEYWORD))]{};
			std::copy(name, sep_it, buffer);
			return (std::strcmp(buffer, STR(VIRTUAL_MODEL_KEYWORD)) == 0);
		}
		else
		{
			return false;
		}
	}

	template<template<typename, size_t> typename AppliedSolver, typename... Ts>
	inline int call_virtual_model(const char* name, Ts&&... args)
	{
		auto sep_it = std::strchr(name, VIRTUAL_MODEL_SEP_KEY_C) + 1;
		auto dir_it = std::strchr(sep_it, VIRTUAL_MODEL_SEP_KEY_C);
		auto end_it = (dir_it == NULL) ? name + std::strlen(name) + 1 : dir_it;

		char* dir = new char[end_it - sep_it + 1] {};
		std::copy(sep_it, end_it, dir);
		dir[end_it - sep_it] = '\0';

		iter_type start, end;
		auto load_conf = load_virtual_conf(sep_it, start, end);
		load_conf.set_title(STR(VIRTUAL_MODEL_KEYWORD));

		constexpr int last_index = decltype(symphas::internal::model_counter(symphas::internal::model_count_index<MAX_DEFINED_MODELS>{}))::value;
		auto result = model_call_wrapper<last_index - 1>::call<AppliedSolver>(load_conf.dimension, load_conf.stp, load_conf.get_model_name(), ModelVirtual<0>(dir, start, end), load_conf);

		delete[] dir;
		return result;
	}

	/* returns false if there is no model with the given name
	 */
	template<template<typename, size_t> typename AppliedSolver, typename... Ts>
	auto call(const char* name, Ts&& ...args)
	{
		if (check_virtual_model(name))
		{
			return call_virtual_model<AppliedSolver>(name, std::forward<Ts>(args)...);
		}
		else
		{
			constexpr int last_index = decltype(symphas::internal::model_counter(symphas::internal::model_count_index<MAX_DEFINED_MODELS>{}))::value;
			return model_call_wrapper<last_index - 1>::call<AppliedSolver>(dimension, stp, name, std::forward<Ts>(args)...);
		}
	}

	template<template<size_t> typename AppliedSolver, typename... Ts>
	auto call(const char* name, Ts&& ...args)
	{
		constexpr int last_index = decltype(symphas::internal::model_counter(symphas::internal::model_count_index<MAX_DEFINED_MODELS>{}))::value;
		return model_call_wrapper<last_index - 1>::call<AppliedSolver>(dimension, name, std::forward<Ts>(args)...);
	}

};

namespace symphas
{
	template<size_t D, typename Sp, typename... S>
	bool run_model(ModelVirtual<D, Sp, S...>& model)
	{
		return run_model(model, model.solver.save, 0);
	}
}

// free energy parameters

//#undef ii_
//#undef ii
//#undef jj_
//#undef jj
//#undef kk_
//#undef kk
//
//#undef op_ii_
//#undef op_ii
//#undef op_jj_
//#undef op_jj

// functions

#undef Gaussian
#undef smoothing
#undef cross


//#undef CONSERVED
//#undef NONCONSERVED


#endif

#endif


template<size_t D>
struct init_expr_select
{
	/* returns false if there is no model with the given name
	 */
	template<typename T>
	auto call(
		const char* name, Axis ax, T* values,
		grid::region_interval<D> const& interval, symphas::interval_data_type const& vdata, double const* coeff, size_t num_coeff)
	{
		using namespace symphas::internal;
		constexpr int last_index = decltype(init_expr_counter(init_expr_count_index<255>{}))::value;
		return init_expr_call_wrapper<D, last_index - 1>::template call(name, ax, values, interval, vdata, coeff, num_coeff);
	}
};

template<size_t D, typename T>
bool match_init_expr(
	const char* initname, Axis ax, T* values,
	grid::region_interval<D> const& interval, symphas::interval_data_type const& vdata, double const* coeff, size_t num_coeff)
{
	init_expr_select<D> ie;

	if (ie.call(initname, ax, values, interval, vdata, coeff, num_coeff) < 0)
	{
		fprintf(SYMPHAS_ERR, "unknown initialization equation provided, '%s'\n", initname);
		return false;
	}
	return true;
}



// undefine all the parameters

// order parameter names

//#undef op
//#undef dop
//#undef var

//#undef lap
//#undef bilap
//#undef gradlap
//#undef grad
//#undef param


#undef x
#undef y
#undef z
#undef t





