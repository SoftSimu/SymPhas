
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
#include "modules-conf.h"
#endif

#include "prereq-defs.h"

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

namespace symphas {
//! Construct a map of all the parameter keys and initialization strategies.
/*!
 * The parameter map which specifies the command line key string and the
 * corresponding initialization strategy for the variable it is supposed
 * to initialize is constructed and returned.
 *
 * If the _io_ module is enabled, it will also add those parameters.
 */
param_map_type build_params();
}  // namespace symphas

namespace symphas::internal {
template <typename>
struct infer_model_from_apply;

template <template <typename> typename model_apply_specialized,
          typename model_type>
struct infer_model_from_apply<model_apply_specialized<model_type>> {
  using type = model_type;
};

#ifdef USING_CONF

template <size_t D, typename Sp, typename... S>
auto run_virtual_model(ModelVirtual<D, Sp, S...>& virt, Conf const& conf) {
  len_type len = len_type(symphas::model_num_fields(virt));
  char** names = new char* [len] {};

  iter_type i = 0;
  for (auto name : conf.name_settings.get_names(len)) {
    names[i] = nullptr;
    std::swap(names[i], name.data);
    ++i;
  }

  symphas::io::write_plot_config(virt, conf.get_result_dir(), names,
                                 conf.simulation_settings.save);
  for (iter_type i = 0; i < len; ++i) {
    delete[] names[i];
  }
  delete[] names;

  find_solution(virt, conf);
  return 1;
}

template <size_t D, typename Sp, typename... S>
auto initialize_virtual_model(Model<D, Sp, S...>* model,
                              ModelVirtual<0> const& setup, Conf const& conf) {
  return ModelVirtual<D, Sp, S...>(setup.dir, conf.simulation_settings.save);
}

template <size_t D, typename Sp, typename... S>
auto initialize_virtual_model(ArrayModel<D, Sp, S...>* model,
                              ModelVirtual<0> const& setup, Conf const& conf) {
  auto save(conf.simulation_settings.save);
  save.set_start(setup.start);
  return ArrayModelVirtual<D, Sp, S...>(setup.dir, save);
}

template <typename model_apply_specialized, typename... Ts>
int handle_simulation_creation(ModelVirtual<0> const& setup, Conf const& conf,
                               Ts&&... args) {
  using model_type =
      typename infer_model_from_apply<model_apply_specialized>::type;
  if (params::plots_only) {
    model_type model(conf.model_settings.coeff, conf.model_settings.coeff_len,
                     conf.get_problem_parameters());
    symphas::io::write_plot_config(model);
    return 1;
  } else {
    auto virt = initialize_virtual_model((model_type*)0, setup, conf);
    return run_virtual_model(virt, conf);
  }
}

#else

template <size_t D, typename Sp, typename... S>
auto run_virtual_model(ModelVirtual<D, Sp, S...>& virt, iter_type start,
                       iter_type end) {
  symphas::io::write_plot_config(virt, '.', SaveParams(start, end));

  symphas::io::write_plot_config(virt);
  find_solution(virt, 1.0, end - start);
  return 1;
}

template <size_t D, typename Sp, typename... S>
auto initialize_virtual_model(Model<D, Sp, S...>* model,
                              ModelVirtual<0> const& setup) {
  return ModelVirtual<D, Sp, S...>(setup.dir);
}

template <size_t D, typename Sp, typename... S>
auto initialize_virtual_model(ArrayModel<D, Sp, S...>* model,
                              ModelVirtual<0> const& setup) {
  return ArrayModelVirtual<D, Sp, S...>(setup.dir);
}

template <typename model_apply_specialized, typename... Ts,
          size_t D = model_dimension<M>::value>
int handle_simulation_creation(ModelVirtual<0> const& setup, int,
                               Ts&&... args) {
  if (params::plots_only) {
    M model(nullptr, 0, {});
    symphas::io::write_plot_config(model);
  } else {
    auto virt =
        initialize_virtual_model<D>(model_parameter_types_t<M>{}, setup);
    return run_virtual_model(virt, setup.start, setup.end);
  }
}

#endif

template <typename model_apply_specialized, typename T0, typename T1,
          typename... Ts>
auto handle_simulation_creation(T0 const& arg0, T1 const& arg1, Ts&&... args) {
  return model_apply_specialized{}(arg0, arg1, std::forward<Ts>(args)...);
}

template <typename model_apply_specialized, typename T0>
auto handle_simulation_creation(T0 const& arg0) {
  return model_apply_specialized{}(arg0);
}

template <typename model_apply_specialized>
auto handle_simulation_creation() {
  return model_apply_specialized{}();
}

template <size_t D>
void print_stencil_message(const char* (&deriv_names)[D],
                           size_t (&stencil_values)[D]) {
  size_t max_name_len = 0;
  for (iter_type i = 0; i < D; ++i) {
    size_t name_len = std::strlen(deriv_names[i]);
    max_name_len = (name_len > max_name_len) ? name_len : max_name_len;
  }

  fprintf(SYMPHAS_LOG,
          "The finite difference stencils of the solver use the "
          "following point values:\n");

  for (iter_type i = 0; i < D; ++i) {
    fprintf(SYMPHAS_LOG, "\t%-*s : %zd\n", static_cast<int>(max_name_len),
            deriv_names[i], stencil_values[i]);
  }
}

template <size_t D, size_t D0, std::enable_if_t<(D != D0), int> = 0>
void print_stencil_message(const char* (&deriv_names)[D],
                           size_t (&stencil_values)[D0]) {
  fprintf(SYMPHAS_LOG, "no stencil point values are defined\n");
}

template <template <typename> typename model_apply_type,
          template <size_t, typename> typename Model,
          template <typename, size_t = 0> typename Solver, size_t N, size_t D,
          size_t O, size_t... Ps, typename... Ts>
auto run_model_call(std::index_sequence<N, D, O, Ps...>, Ts&&... args) {
  fprintf(SYMPHAS_LOG, "The simulation is using solver variation %ld\n", N);

  const char* names[]{"laplacian", "gradlaplacian", "bilaplacian"};
  size_t values[]{Ps...};
  print_stencil_message(names, values);

  using type = typename SelfSelectingStencil<D, O>::template Points<Ps...>;
  using model_apply_specialized = model_apply_type<Model<D, Solver<type, N>>>;
  return handle_simulation_creation<model_apply_specialized>(
      std::forward<Ts>(args)...);
}

template <template <typename> typename model_apply_type,
          template <size_t, typename> typename Model,
          template <size_t> typename Solver, size_t N, size_t D, typename... Ts>
auto run_model_call(std::index_sequence<N, D>, Ts&&... args) {
  using model_apply_specialized = model_apply_type<Model<D, Solver<N>>>;
  return handle_simulation_creation<model_apply_specialized>(
      std::forward<Ts>(args)...);
}
}  // namespace symphas::internal

//! Used for runtime-based model selection.
/*!
 * Used for runtime-based model selection and requires selectable models to have
 * #LINK_WITH_NAME defined. Using the dimension and stencil parameters, a model
 * with the given name will be selected from the list of compiled models. The
 * name must match exactly.
 */
template <template <typename> typename model_apply_type>
struct model_select {
  model_select(size_t dimension, StencilParams stp)
      : dimension{dimension}, stp{stp} {}

  size_t dimension;
  StencilParams stp;

#ifdef USING_CONF

  auto load_virtual_conf(const char* spec, iter_type& start, iter_type& end) {
    char* dir = new char[std::strlen(spec) + 1]{};
    std::strcpy(dir, spec);
    auto dir_it = std::strchr(dir, VIRTUAL_MODEL_SEP_KEY_C);

    start = 0;
    end = start;
    if (dir_it != NULL) {
      *dir_it = '\0';
      sscanf(dir_it + 1, "%d", &start);
      if ((dir_it = std::strchr(dir_it + 1, VIRTUAL_MODEL_SEP_KEY_C)) != NULL) {
        sscanf(dir_it + 1, "%d", &end);
      }
    }

    param_map_type param_map = symphas::build_params();
    auto conf = symphas::conf::restore_checkpoint(param_map, dir, start);

#ifdef USING_PROC
    symphas::internal::update_data_stop(
        conf.simulation_settings.save.get_stop());
#endif

    delete[] dir;
    return conf;
  }

  auto load_virtual_conf(const char* name) {
    iter_type start, end;
    return load_virtual_conf(name, start, end);
  }

#else

  int load_virtual_conf(const char* name, iter_type& start, iter_type& end) {
    start = 0;
    end = 0;
    return {};
  }

  auto load_virtual_conf(const char* name) {
    iter_type start, end;
    return load_virtual_conf(name, start, end);
  }

#endif

  inline bool check_virtual_model(const char* name) {
    if (name != nullptr) {
      auto sep_it = std::strchr(name, VIRTUAL_MODEL_SEP_KEY_C);
      if (sep_it != NULL) {
        char buffer[STR_ARR_LEN(STR(VIRTUAL_MODEL_KEYWORD))]{};
        std::copy(name, sep_it, buffer);
        return (std::strcmp(buffer, STR(VIRTUAL_MODEL_KEYWORD)) == 0);
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  template <template <typename, size_t> typename AppliedSolver, typename... Ts>
  inline int call_virtual_model(const char* name, Ts&&... args) {
    auto sep_it = std::strchr(name, VIRTUAL_MODEL_SEP_KEY_C) + 1;
    auto dir_it = std::strchr(sep_it, VIRTUAL_MODEL_SEP_KEY_C);
    auto end_it = (dir_it == NULL) ? name + std::strlen(name) + 1 : dir_it;

    char* dir = new char[end_it - sep_it + 1]{};
    std::copy(sep_it, end_it, dir);
    dir[end_it - sep_it] = '\0';

    iter_type start, end;
    auto load_conf = load_virtual_conf(sep_it, start, end);
    load_conf.name_settings.set_title(STR(VIRTUAL_MODEL_KEYWORD));

    constexpr int last_index = decltype(symphas::internal::model_counter(
        symphas::internal::model_count_index<MAX_DEFINED_MODELS>{}))::value;
    auto result =
        model_call_wrapper<model_apply_type, last_index - 1>::template call<
            AppliedSolver>(load_conf.simulation_settings.dimension,
                           load_conf.simulation_settings.stp,
                           load_conf.model_settings.model,
                           ModelVirtual<0>(dir, start, end), load_conf);

    delete[] dir;
    return result;
  }

  /* returns false if there is no model with the given name
   */
  template <template <typename, size_t> typename AppliedSolver, typename... Ts>
  auto call(const char* name, Ts&&... args) {
    if (name != nullptr) {
      if (check_virtual_model(name)) {
        return call_virtual_model<AppliedSolver>(name,
                                                 std::forward<Ts>(args)...);
      } else {
        constexpr int last_index = decltype(symphas::internal::model_counter(
            symphas::internal::model_count_index<MAX_DEFINED_MODELS>{}))::value;
        return model_call_wrapper<model_apply_type, last_index - 1>::
            template call<AppliedSolver>(dimension, stp, name,
                                         std::forward<Ts>(args)...);
      }
    } else {
      return 0;
    }
  }

  template <template <size_t> typename AppliedSolver, typename... Ts>
  auto call(const char* name, Ts&&... args) {
    if (name != nullptr) {
      constexpr int last_index = decltype(symphas::internal::model_counter(
          symphas::internal::model_count_index<MAX_DEFINED_MODELS>{}))::value;
      return model_call_wrapper<model_apply_type, last_index - 1>::
          template call<AppliedSolver>(dimension, name,
                                       std::forward<Ts>(args)...);
    } else {
      return 0;
    }
  }
};

namespace symphas {
template <size_t D, typename Sp, typename... S>
bool run_model(ModelVirtual<D, Sp, S...>& model) {
  return run_model(model, model.solver.simulation_settings.save, 0);
}
}  // namespace symphas

template <size_t D>
struct init_expr_select {
  /* returns false if there is no model with the given name
   */
  template <typename T>
  auto call(const char* name, Axis ax, T* values,
            grid::region_interval<D> const& interval,
            symphas::interval_data_type const& vdata, double const* coeff,
            size_t num_coeff) {
    using namespace symphas::internal;
    constexpr int last_index =
        decltype(init_expr_counter(init_expr_count_index<128>{}))::value;
    return init_expr_call_wrapper<D, last_index - 1>::call(
        name, ax, values, interval, vdata, coeff, num_coeff);
  }
};

template <size_t D, typename T>
bool match_init_expr(const char* initname, Axis ax, T* values,
                     grid::region_interval<D> const& interval,
                     symphas::interval_data_type const& vdata,
                     double const* coeff, size_t num_coeff) {
  init_expr_select<D> ie;

  if (ie.call(initname, ax, values, interval, vdata, coeff, num_coeff) < 0) {
    fprintf(SYMPHAS_ERR, "unknown initialization equation provided, '%s'\n",
            initname);
    return false;
  }
  return true;
}

namespace symphas {

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
template <typename M>
void model_iteration(M& model, double dt, double time) {
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
template <typename M>
bool run_model(M& model, iter_type n, symphas::time_step_list const& dts,
               double starttime = 0) {
  double time = starttime;
  iter_type end = model.get_index() + n;

  for (iter_type i = model.get_index(); i < end; i = model.get_index()) {
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
inline void init(const char* param) { init(param, nullptr, 0); }
inline void init() { init("--" ARGUMENT_HELP_STRING); }

void finalize();

}  // namespace symphas

namespace params {

enum program_params_value { PARAMS };

inline void operator+=(program_params_value, const char* arg) {
  params::parse_params(symphas::build_params(), arg);
}

inline void operator+=(program_params_value, std::string const& str) {
  params::parse_params(symphas::build_params(), str.c_str());
}

template <typename T>
void operator+=(program_params_value, T arg) {
  params::set_param(arg);
}

template <typename T>
void operator,(program_params_value, T arg) {
  params::set_param(arg);
}

}  // namespace params
