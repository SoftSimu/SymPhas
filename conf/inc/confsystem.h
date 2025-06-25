
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
 * MODULE:  conf
 * PURPOSE: Defines the configuration functionality for defining parameters
 * of a phase field system.
 *
 * ***************************************************************************
 */

#pragma once

#include <functional>
#include <stack>
#include <string>

#include "io.h"
#include "stencilincludes.h"
#include "systemlib.h"

//! The extension of the configuration file for SymPhas.
#define CONFIG_EXTENSION "in"

//! Key used to prefix comments in the configuration file.
#define CONFIG_COMMENT_PREFIX "#"
#define CONFIG_SUBSTITUTION_PREFIX "$"
#define CONFIG_SUBSTITUTION_PREFIX_C ((CONFIG_SUBSTITUTION_PREFIX)[0])

//! The character usage of the comment key.
#define CONFIG_COMMENT_PREFIX_C ((CONFIG_COMMENT_PREFIX)[0])

//! Key used to prefix the title in the configuration file.
/*!
 * The title is typically specified at the top of the configuration file,
 * but can appear anywhere.
 */
#define CONFIG_TITLE_PREFIX "!"

//! The character usage of the title key.
#define CONFIG_TITLE_PREFIX_C ((CONFIG_TITLE_PREFIX)[0])

//! The key which separates keys and values in the configuration.
#define CONFIG_SEP_KEY ":"

//! The character usage of the separation key.
#define CONFIG_SEP_KEY_C ((CONFIG_SEP_KEY)[0])

//! The way the keys are formatted in the configuration.
#define CONFIG_NAME_FMT "%-12s" CONFIG_SEP_KEY

//! Parameters are added to the configuration using this format.
/*!
 * Parameters are added to the configuration using this format. This is
 * typically only when a backup configuration is written that needs to persist
 * the command line parameters passed to the program.
 */
#define CONFIG_PARAM_PREFIX CONFIG_COMMENT_PREFIX "*"

//! The index at which backups are written are saved with this format.
/*!
 * When backup data is written, the backup configuration file contents are
 * appended to add the current index. The index that is written needs to be
 * preceded by a specific key.
 */
#define CONFIG_INDEX_PREFIX CONFIG_COMMENT_PREFIX "$"

//! Extension given to the file specifying the coefficients.
#define COEFF_SPEC_EXTENSION "constants"

#define COEFF_STR_PREPEND_NAME "c"

//! The format string for a default phase field name.
/*!
 * The format for the default phase field name string. It only takes
 * an index format specification.
 */
#define DEFAULT_PF_NAME_FMT "Phase Field %zd"

namespace symphas {
//! Specifies elements used in the management of SymPhas configuration.
/*!
 * Specifies elements used in the management of SymPhas configuration.
 */
namespace conf {}
}  // namespace symphas

namespace symphas::internal {
char option_open_bracket();
char option_close_bracket();
}  // namespace symphas::internal

namespace symphas::conf {

std::vector<std::string> parse_options(const char* options,
                                       bool spaces_are_delimiters = false,
                                       const char* extra_delimiters = "");
inline std::vector<std::string> parse_options(const char* options,
                                              const char* extra_delimiters) {
  return parse_options(options, false, extra_delimiters);
}

//! Append default phase field names starting at the given pointer.
/*!
 * The iterator to the names list is provided, which is expected to point to the
 * start of an array of at least length `num`. The array is then filled with
 * that many default names, with starting index which may be chosen to be
 * greater than 0.
 *
 * \param names_it The iterator to the position in the names list to begin
 * inserting the default names.
 * \param num The number of default names to insert, and the expected length of
 * the array.
 * \param start_id The starting index for the default names, which is 0 by
 * default.
 */
inline void insert_default_pf_names(char** names_it, size_t num,
                                    size_t start_id = 0) {
  for (size_t id = start_id; id < start_id + num; ++id) {
    char name_buffer[BUFFER_LENGTH_R2];
    sprintf(name_buffer, DEFAULT_PF_NAME_FMT, id);
    names_it[id] = new char[std::strlen(name_buffer) + 1];
    std::strcpy(names_it[id], name_buffer);
  }
}

//! Expand memory and append default names to the given names list.
/*
 * Given a list of names of phase fields, append an additional number of names
 * given by num, in the default phase field name format. The default format
 * is given by #DEFAULT_PF_NAME_FMT. The list is expanded to fit the additional
 * names.
 *
 * \param names The list of strings with the names of the phase fields, which
 * will be appended with num default names. It will be reallocated to the
 * expanded length to accommodate additional names.
 * \param names_len The current length of the names list.
 * \param num The number of default names to append.
 */
inline void append_default_pf_names(char**& names, size_t names_len,
                                    size_t num) {
  char** extend = new char*[names_len + num];
  for (iter_type i = 0; i < names_len; ++i) {
    extend[i] = names[i];
  }

  insert_default_pf_names(extend + names_len, num);

  delete[] names;
  names = extend;
}
}  // namespace symphas::conf

// ***************************************************************************

//! Object storing configurable properties.
/*!
 * The configuration specifies a number of properties which could be used
 * throughout the workflow of a simulation. The properties are entirely managed
 * by this object, so  for example, copies of a configuration can be made in
 * which case all the properties will be correctly constructed and assigned.
 * Also allows the configuration to be persisted.
 *
 * Most of the members of the configuration are publicly accessible, and these
 * members can safely be modified. Members which are which are accessed through
 * function interfaces include the boundaries, initial conditions and
 * system intervals (these are arrays, such that each element represents
 * the description of one field, initialized in the order of the array).
 *
 * The configuration has a default initialization, which does not represent
 * a valid set of properties of a phase field problem.
 *
 * The parameters that are not detailed in the function docs are mentioned here:
 *
 * |Property            | Configuration Key             | Format |
 * |--------------------|-------------------------------|--------|
 * | Dimension          | `DIM` | Either `1`, `2` or `3`. |
 * | FD Stencil Order   | `ORDER` | Either `2` or `4`. |
 * | `N`-point Laplacian Stencil | `PTL` | Refer to available stencils. |
 * | `N`-point Biaplacian Stencil | `PTB` | Refer to available stencils. |
 * | `N`-point Gradlaplacian Stencil | `PTG` | Refer to available stencils. |
 * | Time Step (\f$\Delta t\f$)  | `DELTA` | Any positive floating point number.
 * | | Number of Solution Iterations  | `FRAMES` | Any positive integer. | |
 * Number of Concurrent Simulations  | `FRAMES` | Any positive integer. | |
 * Directory of Results | `DIR` | Any valid system path (relative or absolute).
 * | | Save Interval | `SAVE` | See symphas::io::parse_save_str(). | | Save
 * Index 0? | `SAVEINIT` | `YES` or `NO`. |
 *
 * The format and key of each possible configurable property are described
 * for each function. The general configuration formatting rules are the
 * following:
 * - Each configuration property line consists of the key, the character `:`
 * and then the value of the property. The `:` character separates the key
 * from the value, and is called the _option delimiter_.
 * - Each item in a specification is separated by any number of spaces or tabs.
 * - Numeric values can be specified using scientific notation.
 * - There is a character used for special configuration options, `@`. This
 * character is called the _detail specifier_. It changes the behaviour
 * of some properties, which are described in their respective format
 * explanations.
 * - The title of the simulation is specified on its own line, preceded by
 * the character `!`, the _title specifier_.
 *
 * Some configurable properties support multiple options, where
 * options are simply bracket surrounded elements. Each element of the options
 * follows the format of the property it is written for. Typically, options are
 * given to specify properties about multiple phase fields. This is
 * indicated in the format specification of the respective property.
 *
 * Options may be surrounded in either brackets, `(` and `)`, or braces,
 * `{` and `}`. Options may also be surrounded in quotes, `"` or `'`. When
 * options are surrounded in brackets, nested options can be used.
 *
 * Additionally, the configuration property key is case insensitive. Any
 * number of spaces and tabs may separate the key and the option delimiter, and
 * any number of spaces may separate the option delimiter and the value.
 * Moreover, as long as the key is followed by at least one tab or space,
 * then anything may be written up to the option delimiter (for example,
 * to provide a hint as to the meaning of that key).
 *
 */

using settings_spec_type = std::vector<std::pair<std::string, std::string>>;

struct DirectorySettings {
  char* root_dir;    //!< The directory to save all the simulation information.
  char* result_dir;  //!< The directory to save all the simulation information.

  DirectorySettings() : root_dir{nullptr}, result_dir{nullptr} {}

  DirectorySettings(const char* directory) : DirectorySettings() {
    set_directory(directory);
  }

  DirectorySettings(DirectorySettings const& other)
      : root_dir{new char[std::strlen(other.root_dir) + 1]{}},
        result_dir{new char[std::strlen(other.result_dir) + 1]{}} {
    std::strcpy(root_dir, other.root_dir);
    std::strcpy(result_dir, other.result_dir);
  }

  DirectorySettings(DirectorySettings&& other) : DirectorySettings() {
    swap(*this, other);
  }

  friend void swap(DirectorySettings& first, DirectorySettings& second) {
    using std::swap;
    swap(first.root_dir, second.root_dir);
    swap(first.result_dir, second.result_dir);
  }

  void set_directory(const char* directory, const char* title = nullptr);

  ~DirectorySettings() {
    delete[] root_dir;
    delete[] result_dir;
  }
};

struct SimulationSettings {
  symphas::time_step_list dt_list;  //!< The list of time steps.
  SaveParams save;                  //!< Encapsulation of save parameters.
  StencilParams stp;                //!< Parameters characterizing the stencils.
  size_t dimension;  //!< The dimension of the system, can be 1 or 2.

  SimulationSettings() : dt_list{}, stp{}, dimension{} {}

  void parse_dimension(const char* value);
  void parse_stencil(size_t order, const char* value);
  void parse_stencil_accuracy(const char* str);
  void parse_solver_type(const char* str);
  void parse_dt(const char* str);

  void parse_save(const char* save_spec) {
    symphas::io::parse_save_str(save_spec, &save);
  }

  void parse_save_init(const char* value) {
    save.set_init_flag(std::strcmp(value, "NO") != 0);
  }

  void parse_frame_count(const char* value) { save.set_stop(atoi(value)); }

  //! Provide the time step used in the solution of the problem.
  void set_time_step(double dt, double time) {
    dt_list.set_time_step(dt, time);
  }

  //! Provide the time step used in the solution of the problem.
  void set_time_steps(const double* dts, const double* t_dts, size_t dts_len) {
    dt_list.set_time_steps(dts, t_dts, dts_len);
  }

  void clear_time_steps(double default_dt) {
    dt_list.clear_time_steps(default_dt);
  }

  //! Provide the time step used in the solution of the problem.
  void set_time_steps(symphas::time_step_list const& list) { dt_list = list; }

  inline auto get_time_steps() const { return dt_list.get_time_steps(); }

  inline auto get_times_of_steps() const {
    return dt_list.get_times_of_steps();
  }

  inline size_t get_num_time_steps() const {
    return dt_list.get_num_time_steps();
  }

  inline auto& get_time_step_list() { return dt_list; }

  inline const auto& get_time_step_list() const { return dt_list; }
};

struct NameSettings {
  char* title;
  char** names;      //!< A list of names of each of the phase field.
  size_t names_len;  //!< The number of phase field names provided.

  NameSettings() : title{nullptr}, names{nullptr}, names_len{0} {}

  NameSettings(const char* title) : NameSettings() {
    this->title = new char[std::strlen(title) + 1]{};
    std::strcpy(this->title, title);
  }

  NameSettings(NameSettings const& other)
      : title{new char[std::strlen(other.title) + 1]{}}, names{}, names_len{0} {
    std::strcpy(title, other.title);
    symphas::lib::resize_array(other.names_len, names, names_len);
    for (int i = 0; i < names_len; ++i) {
      names[i] = new char[std::strlen(other.names[i]) + 1]{};
      std::strcpy(names[i], other.names[i]);
    }
  }

  NameSettings(NameSettings&& other) : NameSettings() { swap(*this, other); }

  friend void swap(NameSettings& first, NameSettings& second) {
    using std::swap;
    swap(first.title, second.title);
    swap(first.names, second.names);
    swap(first.names_len, second.names_len);
  }

  void parse_names(const char* value);

  void set_title(const char* title) {
    delete[] this->title;
    this->title = new char[std::strlen(title) + 1]{};
    std::strcpy(this->title, title);
  }

  void set_name(const char* name, size_t n) {
    symphas::lib::set_array_value(name, n, names, names_len);
  }

  symphas::lib::array_container<symphas::lib::string> get_names(
      len_type desired_len = -1) const;
  symphas::lib::string get_name(int n) const;

  const char* get_title() const { return title; }

  ~NameSettings() {
    delete[] title;
    symphas::lib::clear_array(names, names_len);
  }
};

struct ModelSettings {
  char* model;       //!< The string representing the model which is simulated.
  double* coeff;     //!< The coefficients in the equations of motion.
  char** modifiers;  //!< The coefficients in the equations of motion.
  len_type*
      num_fields;  //!< The number of fields to generate in the parameters.
  size_t num_fields_len;  //!< The number of groups of fields, when array model
                          //!< type is used.
  size_t coeff_len;       //!< The number of coefficients provided.
  bool init_coeff_copied;

  ModelSettings()
      : model{nullptr},
        coeff{nullptr},
        modifiers{nullptr},
        num_fields{nullptr},
        num_fields_len{0},
        coeff_len{0},
        init_coeff_copied{false} {}

  ModelSettings(ModelSettings const& other)
      : model{new char[std::strlen(other.model) + 1]{}},
        coeff{new double[other.coeff_len]{}},
        modifiers{new char* [other.num_fields_len] {}},
        num_fields{new len_type[other.num_fields_len]{}},
        num_fields_len{other.num_fields_len},
        coeff_len{other.coeff_len},
        init_coeff_copied{other.init_coeff_copied} {
    std::strcpy(model, other.model);
    for (int i = 0; i < num_fields_len; ++i) {
      num_fields[i] = other.num_fields[i];
      modifiers[i] = new char[std::strlen(other.modifiers[i]) + 1]{};
      std::strcpy(modifiers[i], other.modifiers[i]);
    }
    for (int i = 0; i < coeff_len; ++i) {
      coeff[i] = other.coeff[i];
    }
  }

  ModelSettings(ModelSettings&& other) : ModelSettings() { swap(*this, other); }

  friend void swap(ModelSettings& first, ModelSettings& second) {
    using std::swap;
    swap(first.model, second.model);
    swap(first.coeff, second.coeff);
    swap(first.modifiers, second.modifiers);
    swap(first.num_fields, second.num_fields);
    swap(first.num_fields_len, second.num_fields_len);
    swap(first.coeff_len, second.coeff_len);
    swap(first.init_coeff_copied, second.init_coeff_copied);
  }

  void parse_model_spec(const char* value, const char* dir);
  void set_model_name(const char* str);

  double get_coeff(size_t n) const {
    if (n < coeff_len) {
      return coeff[n];
    } else {
      return DEFAULT_COEFF_VALUE;
    }
  }

  ~ModelSettings() {
    delete[] model;
    delete[] num_fields;
    for (int i = 0; i < num_fields_len; ++i) {
      delete[] modifiers[i];
    }
    delete[] modifiers;
    delete[] coeff;
  }
};

struct DomainSettings {
  symphas::interval_data_type*
      intervals;  //!< A list of the intervals of all the axes of the system.
  symphas::b_data_type* bdata;  //!< A list of the boundary data corresponding
                                //!< to all sides of the system.
  symphas::init_data_type* tdata;  //!< A list of the parameters required for
                                   //!< generating the initial data.
  len_type** dims;  //!< An array of the full dimensions of the system, instead
                    //!< of single variables.

  size_t intervals_len;  //!< The number of interval data elements provided.
  size_t bdata_len;      //!< The number of boundary data elements provided.
  size_t tdata_len;      //!< The number of initial condition elements provided.

  DomainSettings()
      : intervals{nullptr},
        bdata{nullptr},
        tdata{nullptr},
        dims{nullptr},
        intervals_len{0},
        bdata_len{0},
        tdata_len{0} {}

  DomainSettings(DomainSettings const& other)
      : intervals{nullptr},
        bdata{nullptr},
        tdata{nullptr},
        dims{nullptr},
        intervals_len{0},
        bdata_len{0},
        tdata_len{0} {
    symphas::lib::resize_array(other.intervals_len, intervals, intervals_len);
    symphas::lib::resize_array(other.bdata_len, bdata, bdata_len);
    symphas::lib::resize_array(other.tdata_len, tdata, tdata_len);
    symphas::lib::resize_array(other.intervals_len, dims, 0,
                               other.get_dimension());

    for (int i = 0; i < intervals_len; ++i) {
      intervals[i] = other.intervals[i];

      len_type dimension = len_type(other.get_dimension());
      for (int n = 0; n < dimension; ++n) {
        dims[i][n] = other.dims[i][n];
      }
    }
    for (int i = 0; i < bdata_len; ++i) {
      bdata[i] = other.bdata[i];
    }
    for (int i = 0; i < tdata_len; ++i) {
      tdata[i] = other.tdata[i];
    }
  }

  DomainSettings(DomainSettings&& other) : DomainSettings() {
    swap(*this, other);
  }

  friend void swap(DomainSettings& first, DomainSettings& second) {
    using std::swap;
    swap(first.intervals, second.intervals);
    swap(first.bdata, second.bdata);
    swap(first.tdata, second.tdata);
    swap(first.dims, second.dims);
    swap(first.intervals_len, second.intervals_len);
    swap(first.bdata_len, second.bdata_len);
    swap(first.tdata_len, second.tdata_len);
  }

  void parse_boundaries_array(const char* const* b_specs, size_t dimension);
  void parse_boundaries_array(Side side, const char* b_specs, size_t dimension);
  void parse_boundaries(const char* const* b_spec, size_t n, size_t dimension);
  void parse_boundaries(Side side, const char* b_spec, size_t n);
  void parse_interval_array(const char* const* v_specs, size_t dimension);
  void parse_interval_array(Axis ax, const char* r_spec, size_t dimension);
  void parse_interval(const char* value, Axis ax, size_t n);
  void parse_initial_condition_array(const char* t_specs,
                                     ModelSettings* model_settings);
  void parse_initial_condition(const char* value, size_t n,
                               ModelSettings* model_settings);
  void parse_width(const char* str, size_t dimension);

  void set_interval(symphas::interval_data_type const& value, size_t n) {
    if (n >= intervals_len) {
      for (iter_type i = 0; i < intervals_len; ++i) {
        delete[] dims[i];
      }
      delete[] dims;
      dims = new len_type* [n + 1] {};
      for (iter_type i = 0; i < intervals_len; ++i) {
        set_dimensions(i);
      }
    }
    symphas::lib::set_array_value(value, n, intervals, intervals_len);
    set_dimensions(n);
  }
  void set_boundary(symphas::b_data_type const& value, size_t n) {
    symphas::lib::set_array_value(value, n, bdata, bdata_len);
  }
  void set_initial_condition(symphas::init_data_type const& value, size_t n) {
    symphas::lib::set_array_value(value, n, tdata, tdata_len);
  }

  void set_num_domains(size_t len) {
    intervals_len = len;
    bdata_len = len;
    tdata_len = len;
  }

  void set_interval_array(symphas::interval_data_type const* arr,
                          size_t len = 0) {
    for (iter_type i = 0; i < intervals_len; ++i) {
      delete[] dims[i];
    }
    delete[] dims;
    if (len == 0) {
      len = intervals_len;
    }
    std::copy(arr, arr + len, intervals);
    intervals_len = len;
    dims = new len_type* [len] {};
    for (iter_type i = 0; i < len; ++i) {
      set_dimensions(i);
    }
  }

  void set_boundary_array(symphas::b_data_type const* arr, size_t len = 0) {
    if (len == 0) {
      len = bdata_len;
    }
    std::copy(arr, arr + len, bdata);
    bdata_len = len;
  }

  void set_initial_condition_array(symphas::init_data_type const* arr,
                                   size_t len = 0) {
    if (len == 0) {
      len = tdata_len;
    }
    std::copy(arr, arr + len, tdata);
    tdata_len = len;
  }

  size_t get_dimension() const {
    if (intervals_len > 0) {
      auto* max_interval = std::max_element(
          intervals, intervals + intervals_len,
          [](auto& a, auto& b) { return a.size() < b.size(); });
      return max_interval->size();
    } else {
      return 0;
    }
  }

  ~DomainSettings() {
    size_t dimension = get_dimension();
    symphas::lib::clear_array(intervals, intervals_len);
    symphas::lib::clear_array(bdata, bdata_len);
    symphas::lib::clear_array(tdata, tdata_len);
    symphas::lib::clear_array(dims, intervals_len, dimension);
  }

 protected:
  void set_dimensions(size_t n);
  void parse_initial_condition_expression(char* input, size_t n,
                                          ModelSettings* model_settings);
  void parse_initial_condition_file(char* input, size_t n);
  void parse_initial_condition_params(char* input, size_t n);
};

struct SymPhasSettings {
  DirectorySettings directory_settings;
  NameSettings name_settings;
  SimulationSettings simulation_settings;
  ModelSettings model_settings;
  DomainSettings domain_settings;

  size_t runs;
  symphas::lib::string conf_file;

  SymPhasSettings(settings_spec_type const& settings_list, const char* title,
                  const char* dir = "");
  SymPhasSettings()
      : directory_settings{},
        name_settings{},
        simulation_settings{},
        model_settings{},
        domain_settings{},
        runs{1},
        conf_file{} {}

  symphas::problem_parameters_type get_problem_parameters() const;

  //! Return the problem parameters object.
  /*!
   * Constructs the problem parameters object in order to be used in model
   * initialization. The problem parameters are defined so that the maximum
   * amount of information is included.
   */
  void write(const char* savedir, const char* name = BACKUP_CONFIG_NAME) const;

 protected:
  void make_directories() const {
#ifdef USING_MPI
    if (symphas::parallel::is_host_node()) {
      init_work_dirs();
      MPI_Barrier(MPI_COMM_WORLD);
    } else {
      MPI_Barrier(MPI_COMM_WORLD);
      init_work_dirs();
    }
#else

    init_work_dirs();
#endif
  }
  void init_work_dirs() const;
};

struct SystemConf : SymPhasSettings {
  using SymPhasSettings::SymPhasSettings;

  SystemConf(SymPhasSettings const& settings) : SymPhasSettings(settings) {}

  //! Set the dimensions using the nth interval set.
  /*!
   * Set the SystemConf::dims member to be an array of length prescribed by
   * the member SystemConf::dimension. Then fill the array by computing the
   * length (number of grid points) along each of the intervals from the data
   * given by the nth element from the interval list SystemConf::intervals.
   *
   * \param n The element index from SystemConf::intervals to choose the
   * interval data from.
   */
  const char* get_result_dir() const { return directory_settings.result_dir; }
  size_t system_count() const {
    return std::max({domain_settings.bdata_len, domain_settings.intervals_len,
                     domain_settings.tdata_len});
  }
};
