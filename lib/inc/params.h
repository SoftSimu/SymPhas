
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
 * MODULE:  lib
 * PURPOSE: Manages functionality about command line parameters and
 * program global parameters used in changing the behaviour of the program.
 *
 * ***************************************************************************
 */

#pragma once

#include <algorithm>
#include <map>
#include <utility>

#include "spslib.h"

// \cond

#define ARGUMENT_HELP_STRING "help"

namespace symphas {
constexpr inline char symphas_usage_message[] = R"~(
symphas_impl [CONFIG_NAME] [OPTIONS]...
)~";

constexpr inline char symphas_description_message[] = R"~(
SymPhas is a symbolic algebra framework that formulates expressions at compile time for high-speed 
numerical solutions of phase-field problems. It is equipped with two base solvers: a finite 
difference solver with auto-generating stencils of any order of accuracy (up to compilation
resource limits), and a spectral solver. The finite difference solver additionally includes a
variation to find a minimal domain of a field to reduce the iterable region when evaluating an
expression.

If this SymPhas executable reports that a model doesn't exist, verify the CMake parameters are
correctly set. 
)~"
#ifdef DEBUG
                                                      R"~(
This build currently has DEBUG mode enabled, which may indicate optimizations have been turned off.
Parallelization is also off.
)~"
#endif
    ;

}  // namespace symphas

#define SYMPHAS_USAGE_MESSAGE symphas::symphas_usage_message
#define SYMPHAS_DESCRIPTION_MESSAGE symphas::symphas_description_message

#define VIRTUAL_MODEL_KEYWORD VIRTUAL

#define VIRTUAL_MODEL_SEP_KEY *
#define VIRTUAL_MODEL_SEP_KEY_C (STR(VIRTUAL_MODEL_SEP_KEY)[0])

// \endcond

/* **************************************************************************
 * Parameter assignment functionality.
 * **************************************************************************/

using config_key_value_type = std::vector<symphas::lib::string>;

namespace params {

//! Functionality object used to assign parameters.
/*!
 * Interface which is inherited by specialized classes that will assign
 * a parameter using the member function param_assign_base::assign(void*, const
 * char*).
 */
struct param_assign_base {
  //! Assigns the parameter to the given value.
  /*!
   * Assigns the value from the given string to parameter. Virtual function is
   * implemented within param_assign and specialized by parameter type.
   */
  virtual void assign(void* param, const char* value) = 0;
  virtual void set_default(void* param) = 0;

  //! Print the formatted name when writing to the help output.
  virtual size_t print_with_name(FILE* out, void* param,
                                 const char* name) const = 0;

  //! Print the formatted name when writing to the help output.
  virtual size_t print_value(int len, char* out, void* param) const = 0;

  virtual param_assign_base* clone() const = 0;

  virtual ~param_assign_base() {}
};

//! Implementation of parameter assignment functionality.
/*!
 * Child class of param_assign_base used for assigning a parameter of a
 * specified type to the correct value.
 *
 * Implements the virtual function param_assign_base::assign(void*, const
 * char*). Each specialization has a different implementation of how the
 * parameter is assigned depending on template parameter type `p_type`.
 *
 * \tparam p_type The parameter type name which will change how the parameter
 * string is parsed and assigned.
 */
template <typename p_type>
struct param_assign;

//! Implementation of multiple parameter assignment.
/*!
 * Uses the params::param_assign functionality to assign the same type to
 * multiple parameters. This is useful for options which typically
 * should control more than one parameter at once, such as simultaneously
 * setting read and write parameters.
 *
 * \tparam p_type The parameter type name which will change how the parameter
 * string is parsed and assigned.
 */
template <typename p_type, size_t N>
struct param_assign_multiple : params::param_assign_base {
  param_assign<p_type> default_assigns[N];

  //! Implementation to assign the same value to multiple parameters.
  /*!
   * A new character string is created with the same length as value. The
   * provided value string is copied into the new one, and the new one is
   * assigned to the parameter.
   *
   * \param param The parameter to assign.
   * \param value The value as a string to assign to `param`.
   */
  void assign(void* params, const char* value) override {
    using param_type = void* [N];
    void*(&param_arr)[N] = *static_cast<param_type*>(params);

    for (iter_type i = 0; i < N; ++i) {
      default_assigns[i].assign(param_arr[i], value);
    }
  }

  void set_default(void* params) override {
    using param_type = void* [N];
    void*(&param_arr)[N] = *static_cast<param_type*>(params);

    for (iter_type i = 0; i < N; ++i) {
      default_assigns[i].set_default(param_arr[i]);
    }
  }

  size_t print_with_name(FILE* out, void* param,
                         const char* name) const override {
    return fprintf(out, "%s", name);
  }

  size_t print_value(int len, char* out, void* params) const override {
    using param_type = void* [N];
    void*(&param_arr)[N] = *static_cast<param_type*>(params);

    size_t n = 0;
    for (iter_type i = 0; i < N; ++i) {
      if (len - n > 0) {
        n += default_assigns[i].print_value(int(len - n), out, param_arr[i]);
      }
    }
    return n;
  }

  param_assign_base* clone() const override {
    return new param_assign_multiple<p_type, N>{*this};
  }
};

template <typename... p_types>
struct param_assign_separate : params::param_assign_base {
  std::tuple<param_assign<p_types>...> default_assigns;

 protected:
  template <size_t... Is, size_t N = sizeof...(Is)>
  void assign(void* params, const char* value, std::index_sequence<Is...>) {
    using param_type = void* [N];
    void*(&param_arr)[N] = *static_cast<param_type*>(params);

    (std::get<Is>(default_assigns).assign(param_arr[Is], value), ...);
  }

  template <size_t... Is, size_t N = sizeof...(Is)>
  void set_default(void* params, std::index_sequence<Is...>) {
    using param_type = void* [N];
    void*(&param_arr)[N] = *static_cast<param_type*>(params);

    (std::get<Is>(default_assigns).set_default(param_arr[Is]), ...);
  }

  template <size_t... Is, size_t N = sizeof...(Is)>
  size_t print_value(int len, char* out, void* params,
                     std::index_sequence<Is...>) const {
    using param_type = void* [N];
    void*(&param_arr)[N] = *static_cast<param_type*>(params);

    size_t n = 0;
    ((n += (len - n > 0) ? param_assign<p_types>{}.print_value(
                               len - int(n), out, param_arr[Is])
                         : 0),
     ...);
    return n;
  }

 public:
  //! Implementation to assign the same value to multiple parameters.
  /*!
   * A new character string is created with the same length as value. The
   * provided value string is copied into the new one, and the new one is
   * assigned to the parameter.
   *
   * \param param The parameter to assign.
   * \param value The value as a string to assign to `param`.
   */
  void assign(void* params, const char* value) override {
    assign(params, value, std::make_index_sequence<sizeof...(p_types)>{});
  }

  void set_default(void* params) override {
    set_default(params, std::make_index_sequence<sizeof...(p_types)>{});
  }

  size_t print_with_name(FILE* out, void* param,
                         const char* name) const override {
    return fprintf(out, "%s", name);
  }

  size_t print_value(int len, char* out, void* params) const override {
    return print_value(len, out, params,
                       std::make_index_sequence<sizeof...(p_types)>{});
  }

  param_assign_base* clone() const override {
    return new param_assign_separate<p_types...>{*this};
  }
};

}  // namespace params

//! Specialization of parameter assignment functionality for `bool` type.
/*!
 * A parameter will be set to `true` if the string value matches either of
 * "true" or "yes", and set to `false` if it matches either of "false" or "no".
 * The value is case insensitive.
 */

template <>
struct params::param_assign<bool> : params::param_assign_base {
  bool default_value;

  //! Implementation to assign a bool from a string to a parameter.
  /*!
   * A truth value in a string is parsed and assigned to the parameter. If
   * the truth value could not be parsed, the program exits and reports the
   * error. The truth value can either be "true" or "yes" to indicate `true`, or
   * it can be "false" or "no" to indicate `false`. The input is not case
   * sensitive.
   *
   * \param param The parameter to assign.
   * \param value The value as a string to assign to `param`.
   */
  void assign(void* param, const char* value) override {
    *static_cast<bool*>(param) = extract_bool(value, default_value);
  }

  void set_default(void* param) override {
    default_value = *static_cast<bool*>(param);
  }

  size_t print_with_name(FILE* out, void* param,
                         const char* name) const override {
    return fprintf(out, "%s[=yes|no](default=%s)", name,
                   (*static_cast<bool*>(param)) ? "yes" : "no");
  }

  size_t print_value(int len, char* out, void* param) const override {
    return sprintf(out, "%.*s", len,
                   (*static_cast<bool*>(param)) ? "yes" : "no");
  }

  param_assign_base* clone() const override {
    return new param_assign<bool>{*this};
  }

 protected:
  bool extract_bool(const char* value, bool default_value) {
    size_t len = std::strlen(value) + 1;

    if (len > 1) {
      char* cpy = new char[len];
      std::strncpy(cpy, value, len);
      symphas::lib::to_lower(cpy);

      if (std::strcmp(cpy, "true") == 0 || std::strcmp(cpy, "yes") == 0) {
        return true;
      } else if (std::strcmp(cpy, "false") == 0 ||
                 std::strcmp(cpy, "no") == 0) {
        return false;
      } else if (atoi(cpy) > 0) {
        return true;
      } else {
        return default_value;
      }
    } else {
      return !default_value;
    }
  }
};

//! Specialization of parameter assignment functionality for `double` type.
template <>
struct params::param_assign<double> : params::param_assign_base {
  double default_value;

  //! Implementation to assign a double value from a string to a parameter.
  /*!
   * An double value in a string is parsed and assigned to the parameter. If the
   * string could not be parsed for a double, the parameter is simply assigned
   * `0.0`.
   *
   * \param param The parameter to assign.
   * \param value The value as a string to assign to `param`.
   */
  void assign(void* param, const char* value) override {
    *static_cast<double*>(param) = extract_double(value, default_value);
  }

  void set_default(void* param) override {
    default_value = *static_cast<double*>(param);
  }

  size_t print_with_name(FILE* out, void* param,
                         const char* name) const override {
    return fprintf(out, "%s=float(default=%.2lf)", name,
                   *static_cast<double*>(param));
  }

  size_t print_value(int len, char* out, void* param) const override {
    return sprintf(out, "%.2lf", *static_cast<double*>(param));
  }

  param_assign_base* clone() const override {
    return new param_assign<double>{*this};
  }

 protected:
  double extract_double(const char* str, double default_value) {
    char* end;
    double value = std::strtod(str, &end);
    if (end == str || *end != '\0') {
      return default_value;
    } else {
      return value;
    }
  }
};

//! Specialization of parameter assignment functionality for `int` type.
template <>
struct params::param_assign<int> : params::param_assign_base {
  int default_value;

  //! Implementation to assign an integer value from a string to a parameter.
  /*!
   * An integer value in a string is parsed and assigned to the parameter. If
   * the string could not be parsed for an integer, the parameter is simply
   * assigned `0`.
   *
   * \param param The parameter to assign.
   * \param value The value as a string to assign to `param`.
   */
  void assign(void* param, const char* value) override {
    *static_cast<int*>(param) = extract_int(value, default_value);
  }

  void set_default(void* param) override {
    default_value = *static_cast<int*>(param);
  }

  size_t print_with_name(FILE* out, void* param,
                         const char* name) const override {
    return fprintf(out, "%s=N(default=%d)", name, *static_cast<int*>(param));
  }

  size_t print_value(int len, char* out, void* param) const override {
    return sprintf(out, "%d", *static_cast<int*>(param));
  }

  param_assign_base* clone() const override {
    return new param_assign<int>{*this};
  }

 protected:
  int extract_int(const char* value, int default_value) {
    int result = atoi(value);
    if (result == 0 && *value != '0') {
      return default_value;
    } else {
      return result;
    }
  }
};

//! Specialization of parameter assignment functionality for `char*` type.
template <>
struct params::param_assign<char*> : params::param_assign_base {
  char* default_value;

  param_assign() : default_value{nullptr} {}
  param_assign(param_assign const& other)
      : default_value{(other.default_value)
                          ? new char[std::strlen(other.default_value) + 1]{}
                          : nullptr} {
    if (other.default_value) {
      std::strcpy(default_value, other.default_value);
    }
  }
  param_assign(param_assign&& other) noexcept : param_assign() {
    swap(*this, other);
  }

  friend void swap(param_assign& first, param_assign& second) {
    std::swap(first.default_value, second.default_value);
  }

  //! Implementation to assign a new char array to a parameter.
  /*!
   * A new character string is created with the same length as value. The
   * provided value string is copied into the new one, and the new one is
   * assigned to the parameter.
   *
   * \param param The parameter to assign.
   * \param value The value as a string to assign to `param`.
   */
  void assign(void* param, const char* value) override {
    char*& param_str = *static_cast<char**>(param);
    delete[] param_str;
    param_str = extract_string(value, default_value);
  }

  void set_default(void* param) override {
    char* ptr = static_cast<char*>(param);
    default_value = new char[std::strlen(ptr) + 1]{};
    std::strcpy(default_value, ptr);
  }

  size_t print_with_name(FILE* out, void* param,
                         const char* name) const override {
    return fprintf(out, "%s=NAME", name);
  }

  size_t print_value(int len, char* out, void* param) const override {
    char* str = *static_cast<char**>(param);
    if (str && std::strlen(str) > 0) {
      return sprintf(out, "%.*s", len, str);
    } else {
      return sprintf(out, "%.*s", len, "<empty>");
    }
  }

  ~param_assign() override { delete[] default_value; }
  param_assign_base* clone() const override {
    return new param_assign<char*>{*this};
  }

 protected:
  char* extract_string(const char* value, const char* default_value) {
    size_t len = std::strlen(value) + 1;
    char* save_value = new char[len];
    std::strncpy(save_value, value, len);
    return save_value;
  }
};

//! Specialization of parameter assignment functionality for lists of string
//! type.
template <>
struct params::param_assign<config_key_value_type> : params::param_assign_base {
  config_key_value_type default_value;

  //! Implementation to assign a substitution key-value list from a matrix
  //! parameter file
  /*!
   * A list of parameters is copied from a file containing a matrix of
   * parameters. A specific line is chosen from the list, which is by default
   * the first line.
   *
   * \param param The parameter to assign.
   * \param value The value as a string to assign to `param`.
   */
  void assign(void* param, const char* value) override {
    *static_cast<config_key_value_type*>(param) =
        extract_values(value, default_value);
  }

  void set_default(void* param) override {
    default_value = *static_cast<config_key_value_type*>(param);
  }

  size_t print_with_name(FILE* out, void* param,
                         const char* name) const override {
    return fprintf(out, "%s=FILE[,N]", name);
  }

  size_t print_value(int len, char* out, void* param) const override {
    auto& values = *static_cast<config_key_value_type*>(param);
    if (values.size() > 0) {
      size_t n = 0;
      n += sprintf(out, "(");
      len -= int(n);

      auto it = std::begin(values);
      for (; it < std::end(values) - 1; ++it) {
        char* str = it->data;
        if (len > 2) {
          size_t nn = sprintf(out + n, "%.*s, ", len - 2, str);
          len -= int(nn);
          n += nn;
        }
      }
      if (len > 1) {
        char* str = it->data;
        n += sprintf(out + n, "%.*s)", len, str);
      }
      return n;
    } else {
      return sprintf(out, "%.*s", len, "<empty>");
    }
  }

  param_assign_base* clone() const override {
    return new param_assign<config_key_value_type>{*this};
  }

 protected:
  config_key_value_type extract_values(
      const char* value, config_key_value_type const& default_value) {
    config_key_value_type key_values;

    iter_type index = 0;
    const char* pos_index = strrchr(value, ',');
    char* fname;
    if (pos_index != NULL) {
      index = atoi(pos_index + 1);
      fname = new char[int(pos_index - value) + 1];
      std::copy(value, pos_index, fname);
      fname[int(pos_index - value)] = '\0';
    } else {
      fname = new char[std::strlen(value) + 1];
      std::strcpy(fname, value);
    }

    FILE* f = fopen(fname, "r");
    if (f == NULL) {
      fprintf(SYMPHAS_ERR,
              "substitution parameter matrix file '%s' was not found\n", fname);
    } else {
      char* line = new char[BUFFER_LENGTH]{};
      do {
        std::fgets(line, BUFFER_LENGTH, f);
        symphas::lib::str_trim(line);
      } while (index-- > 0);

      char* token = std::strtok(line, " ");
      while (token != NULL) {
        // double value = atof(token);
        key_values.emplace_back(
            symphas::lib::string{token, std::strlen(token) + 1});
        token = std::strtok(NULL, " ");
      }
      delete[] line;
    }

    delete[] fname;
    return key_values;
  }
};

#ifdef EXECUTION_HEADER_AVAILABLE
#include <execution>

namespace symphas {
enum ParallelizationType { SEQ, PAR };
}

template <>
struct params::param_assign<symphas::ParallelizationType>
    : params::param_assign_base {
  symphas::ParallelizationType default_value;

  //! Implementation to assign a bool from a string to a parameter.
  /*!
   * A truth value in a string is parsed and assigned to the parameter. If
   * the truth value could not be parsed, the program exits and reports the
   * error. The truth value can either be "true" or "yes" to indicate `true`, or
   * it can be "false" or "no" to indicate `false`. The input is not case
   * sensitive.
   *
   * \param param The parameter to assign.
   * \param value The value as a string to assign to `param`.
   */
  void assign(void* param, const char* value) override {
    *static_cast<symphas::ParallelizationType*>(param) = extract_policy(
        value, *static_cast<symphas::ParallelizationType*>(param));
  }

  void set_default(void* param) override {
    default_value = *static_cast<symphas::ParallelizationType*>(param);
  }

  size_t print_with_name(FILE* out, void* param,
                         const char* name) const override {
    return fprintf(out, "%s[=seq|par](default=%s)", name,
                   (*static_cast<symphas::ParallelizationType*>(param) ==
                    symphas::ParallelizationType::SEQ)
                       ? "seq"
                       : "par");
  }

  size_t print_value(int len, char* out, void* param) const override {
    return sprintf(out, "%.*s", len,
                   (*static_cast<symphas::ParallelizationType*>(param) ==
                    symphas::ParallelizationType::SEQ)
                       ? "seq"
                       : "par");
  }

  param_assign_base* clone() const override {
    return new param_assign<symphas::ParallelizationType>{*this};
  }

 protected:
  symphas::ParallelizationType extract_policy(
      const char* value, symphas::ParallelizationType default_value) {
    size_t len = std::strlen(value) + 1;

    if (len > 1) {
      char* cpy = new char[len];
      std::strncpy(cpy, value, len);
      symphas::lib::to_lower(cpy);

      symphas::ParallelizationType policy =
          (std::strcmp(cpy, "seq") == 0)   ? symphas::ParallelizationType::SEQ
          : (std::strcmp(cpy, "par") == 0) ? symphas::ParallelizationType::PAR
                                           : default_value;

      delete[] cpy;
      return policy;
    } else {
      return (default_value == symphas::ParallelizationType::PAR)
                 ? symphas::ParallelizationType::SEQ
                 : symphas::ParallelizationType::PAR;
    }
  }
};
#endif

/* **************************************************************************
 * Parameter object typedefs
 * **************************************************************************/

struct param_map_element {
  param_map_element()
      : parameter{nullptr}, assign_method{nullptr}, alias{}, description{} {}

  template <typename specialized_assign_type>
  param_map_element(void* parameter,
                    specialized_assign_type const& assign_method)
      : param_map_element(parameter, assign_method, '\0', "") {}

  template <typename specialized_assign_type>
  param_map_element(void* parameter,
                    specialized_assign_type const& assign_method, char alias)
      : param_map_element(parameter, assign_method, alias, "") {}

  template <typename specialized_assign_type>
  param_map_element(void* parameter,
                    specialized_assign_type const& assign_method,
                    std::string description)
      : param_map_element(parameter, assign_method, '\0', description) {}

  template <typename specialized_assign_type>
  param_map_element(void* parameter,
                    specialized_assign_type const& assign_method, char alias,
                    std::string description)
      : parameter{parameter},
        assign_method{assign_method.clone()},
        alias{alias},
        description{description} {}

  param_map_element(param_map_element const& other)
      : parameter{other.parameter},
        assign_method{other.assign_method->clone()},
        alias{other.alias},
        description{other.description} {}

  param_map_element(param_map_element&& other) noexcept : param_map_element() {
    swap(*this, other);
  }
  param_map_element operator=(param_map_element other) noexcept {
    swap(*this, other);
    return *this;
  }

  friend void swap(param_map_element& first, param_map_element& second) {
    using std::swap;
    swap(first.parameter, second.parameter);
    swap(first.assign_method, second.assign_method);
    swap(first.alias, second.alias);
    swap(first.description, second.description);
  }

  template <typename T>
  void assign(T&& value) {
    assign_method->assign(parameter, std::forward<T>(value));
  }

  void print_help(FILE* out, const char* name) const;
  void print_help(const char* name) const { print_help(SYMPHAS_INFO, name); }
  void print_help(FILE* out, std::string name) const {
    print_help(out, name.c_str());
  }
  void print_help(std::string name) const { print_help(SYMPHAS_INFO, name); }

  symphas::lib::string value_str() const;

  void* parameter;
  params::param_assign_base* assign_method;
  char alias;
  std::string description;

  ~param_map_element() { delete assign_method; }
};

//! The type specifying the parameter map.
/*!
 * The parameter map is constructed with the command line parameter string key
 * as the map key, used to relate it to the program level variable that is
 * modified.
 *
 * The value elements of the parameter map are pair objects containing the
 * pointer to variable and variable interpretation strategy, respectively. The
 * variable interpretation strategy is simply the pointer to the concrete
 * class params::param_assign, specialized on the type of the variable that
 * it must assign.
 *
 * Once this map is populated, it is given to the function
 * params::parse_params() which initializes all the variables from the
 * command line.
 */
using param_map_type =
    std::map<std::string, param_map_element, symphas::lib::any_case_comparator>;

/* **************************************************************************
 * Parameters and parameter management functionality
 * **************************************************************************/

//! Adds parameters to the given parameter map.
/*!
 * Adds all the program parameters to the given parameter map. Upon program
 * initialization, the global parameter map at params::param_map is empty. In
 * order for command line parameters to be parsed and set, a complete map
 * representing the translation between user provided string and variable
 * needs to be constructed.
 *
 * \param param_map The map into which to add key/value pairs corresponding to
 * parameter string names and pointers to parameter variables.
 */
bool add_base_params(param_map_type& param_map);

//! Contains functionality for command line arguments.
/*!
 * API for dealing with the command line arguments and assigning program
 * level parameters.
 * The parameters controlled by this logic modify behaviour across the whole
 * SymPhas library for the required modules. Optional modules may introduce
 * additional parameters.
 */
namespace params {

//! Name of the simulation.
/*!
 * Global name of the simulation.
 */
DLLLIB extern char* title;

//! Flag indicating position of boundaries in the discrete grid.
/*!
 * Whether or not to extend the boundary when creating the boundary grid.
 * choosing this to be true means that the grid sides are increased
 * by BOUNDARY_DEPTH units (logically and in memory) rather than having the
 * interior be BOUNDARY_DEPTH units less on each side from the given dimensions.
 *
 * Briefly, this value will either:
 * - shrink the dimensions of the simulated
 * phase-field system and recompute the spatial discretization, thereby
 * maintaining the total grid size equal to what has been specified in
 * configuration; or
 * - extend the dimensions beyond the dimensions specified in the
 * configuration, thereby maintaining the value of the spatial discretization
 * the same as given in the configuration.
 *
 * **In more detail: (this only applies to the boundary system)**
 *
 * This will affect the way the width is computed. When the configuration
 * is specified and the left and right endpoints are provided along with
 * the width, this information is used to compute the system dimensions.
 * When params::extend_boundary is false, then the dimensions of the phase
 * field will be computed by reducing the computed system dimensions. Thus,
 * the system dimensions of the physical system will be reduced, and in
 * order to maintain the values of the left and right endpoints that were
 * provided, the spatial discretization will change.
 *
 * However, when this value is true, then the boundary will be extended
 * from the computed system dimensions, so the simulated system spatial
 * discretization will precisely match what is provided in the
 * configuration.
 *
 * This value is false by default. However,
 * the latter case may be desired, so this value should be changed on
 * the command line.
 * It is false by default since it is desirable for the system
 * dimensions to be exactly as provided in the configuration, for example,
 * so that dimensions that are a power of two are used to
 * maximize caching and other possible computations.
 */
DLLLIB extern bool extend_boundary;

//! Enable VTK output.
/*!
 * This program allows the user to visualize the phase field evolution with
 * vtk.
 */
DLLLIB extern bool viz_enabled;

//! Interval of updating VTK output.
/*!
 * This program allows the user to visualize the phase field evolution with
 * vtk, and by specifying this parameter for a value > 0, the plot will
 * update at the given value (in milliseconds).
 */
DLLLIB extern int viz_interval;

//! Index of the field to get VTK output.
/*!
 * This program allows the user to visualize the phase field evolution with
 * vtk, and will try to find a field matching that index, or the next
 * closest field, which it will display. By default, this is always
 * the first field.
 */
DLLLIB extern int viz_index;

DLLLIB extern void* viz_interval_enable[2];
DLLLIB extern void* viz_index_enable[2];

//! The values used as interior values to initial conditions.
/*!
 * These parameters let the user control the interior
 * values of the structured initial conditions.
 */
DLLLIB extern double init_inside_val;

//! The values used as exterior values to initial conditions.
/*!
 * These parameters let the user control the exterior
 * values of the structured initial conditions.
 */
DLLLIB extern double init_outside_val;

//! Tunes the magnitude of randomness in initial conditions.
/*!
 * Tunes the magnitude of randomness used in some algorithms that
 * generate the initial data.
 */
DLLLIB extern double init_rand_val;

//! A string for the name of a data file.
/*!
 * The name of a data file that can be used later in the program. The
 * configuration can use this parameter to read an input file as initial data.
 * The name can contain commas (with no spaces), in which case multiple files
 * can be given.
 */
DLLLIB extern char* input_data_file;

//! Determines the start index of the simulation.
/*!
 * The value of this parameter determines the simulation-wide starting
 * index. If the start index should be changed separately from the parsed
 * command line arguments, this value needs to be modified before the
 * simulation is initiated.
 *
 * This is #INDEX_INIT by default.
 */
DLLLIB extern int start_index;

#ifdef EXECUTION_HEADER_AVAILABLE
DLLLIB extern symphas::ParallelizationType parallelization;
#endif

DLLLIB extern config_key_value_type config_key_values;

//! List of all the key/value parameter strings.
/*!
 * Each of the parameters that are parsed and constructed are first
 * stored in the form of the initially provided string, which is first
 * split into a key (variable link name) and value (variable value) pair.
 *
 * This list can be used to refer to all the raw values processed by the
 * command line utility.
 */
DLLLIB extern std::vector<std::pair<std::string, std::string>> rawparams;

//! This enumeration is used as `ParamNames << "desired value"`.
/*!
 * This enumeration is used as `ParamNames << "desired value"`, where
 * desired value is either a boolean, int or floating point number depending
 * on the parameter that is being set.
 */
enum ParamNames {
  EXTEND_BOUNDARY,   //!< See params::extend_boundary.
  VIZ_ENABLED,       //!< See params::viz_enabled.
  VIZ_INTERVAL,      //!< See params::viz_interval.
  VIZ_INDEX,         //!< See params::viz_index.
  INIT_INSIDE_VAL,   //!< See params::init_inside_val.
  INIT_OUTSIDE_VAL,  //!< See params::init_outside_val.
  INIT_RAND_VAL,     //!< See params::init_rand_val.
  PARALLELIZATION    //!< See params::parallelization. Either true or false.
};

struct param_value_flag {
  ParamNames param;
  bool value;
};

struct param_value_count {
  ParamNames param;
  int value;
};

struct param_value_float {
  ParamNames param;
  double value;
};

inline param_value_flag operator<<(ParamNames param, bool value) {
  return {param, value};
}

inline param_value_count operator<<(ParamNames param, int value) {
  return {param, value};
}

inline param_value_float operator<<(ParamNames param, double value) {
  return {param, value};
}

inline void set_param(ParamNames param, bool value) {
  switch (param) {
    case ParamNames::EXTEND_BOUNDARY:
      extend_boundary = value;
      break;
    case ParamNames::VIZ_ENABLED:
      viz_enabled = value;
      break;
#ifdef EXECUTION_HEADER_AVAILABLE
    case ParamNames::PARALLELIZATION:
      parallelization = (value) ? symphas::ParallelizationType::PAR
                                : symphas::ParallelizationType::SEQ;
      break;
#endif
    default:
      break;
  }
}

inline void set_param(ParamNames param, int value) {
  switch (param) {
    case ParamNames::VIZ_INTERVAL:
      viz_enabled = true;
      viz_interval = value;
      break;
    case ParamNames::VIZ_INDEX:
      viz_enabled = true;
      viz_index = value;
      break;
    case ParamNames::INIT_INSIDE_VAL:
      init_inside_val = value;
      break;
    case ParamNames::INIT_OUTSIDE_VAL:
      init_outside_val = value;
      break;
    case ParamNames::INIT_RAND_VAL:
      init_rand_val = value;
      break;
    default:
      break;
  }
}

inline void set_param(ParamNames param, double value) {
  switch (param) {
    case ParamNames::INIT_INSIDE_VAL:
      init_inside_val = value;
      break;
    case ParamNames::INIT_OUTSIDE_VAL:
      init_outside_val = value;
      break;
    case ParamNames::INIT_RAND_VAL:
      init_rand_val = value;
      break;
    default:
      break;
  }
}

inline void set_param(param_value_flag const& input) {
  auto [param, value] = input;
  set_param(param, value);
}

inline void set_param(param_value_count const& input) {
  auto [param, value] = input;
  set_param(param, value);
}

inline void set_param(param_value_float const& input) {
  auto [param, value] = input;
  set_param(param, value);
}

//! Parse the arguments and assign them to the globally linked parameters.
/*!
 * Parsing arguments requires the number of arguments to parse and the
 * parameter map indicating the link between provided string input and the
 * parameter values.
 *
 * The list of arguments first requires the dictionary of keys in the form
 * of the ::param_map_type typedef. This is a map which specifies the
 * parameter name as a string, the pointer to the program variable (the
 * parameter) and the parameter type which is a specialized template class
 * of with parent param_assign_base. The parameter string name is the map
 * key, and the map value associated with the key consists of the pointer to
 * the program variable (parameter) represented by the key and a
 * param_assign template class specialized to the parameter type. The map
 * key selects the map element in order to be assigned the value that has
 * been provided on the command line.
 *
 * The command line arguments are parsed into key/value pairs, and the
 * provided keys are compared against the string keys from the parameter
 * map. Once the key is used to select the pointer to the parameter variable
 * from the map and its type, the param_assign specialization uses the
 * implementation of a method to assign a variable of the correct type by
 * processing the input and accordingly assign the parameter.
 *
 * \param param_map The dictionary between a key string and the variable
 * which the key refers to.
 * \param args The list of argument strings in the form `key=value`.
 * \param n The number of provided arguments, which is the length of the
 * list of arguments.
 */
void parse_params(param_map_type param_map, const char* const* args, int n);

//! See params::parse_params().
void parse_params(param_map_type param_map, const char* args, int n);

inline void parse_params(param_map_type param_map, const char* args) {
  parse_params(param_map, &args, 1);
}

template <typename... Ts>
void parse_params(param_map_type param_map, const char* arg0, const char* arg1,
                  Ts&&... args) {
  const char* list[]{arg0, arg1, std::forward<Ts>(args)...};
  parse_params(param_map, list, 2 + sizeof...(Ts));
}

//! Print the list of all arguments.
void print_params(FILE* out, param_map_type param_map);

//! Print the list of all arguments.
inline void print_params(param_map_type param_map) {
  print_params(SYMPHAS_INFO, param_map);
}

//! Print the list of all arguments.
void print_param_help(FILE* out, param_map_type param_map);

//! Print the list of all arguments.
void print_param_help(param_map_type param_map);

//! Assign the value of a parameter.
/*!
 * Assign a single parameter to its globally linked variable. Used outside
 * of the context of parsing arguments, as it assigns a particular
 * parameter to the given value. The value needs to correspond to the
 * correct type.
 *
 * \param param The parameter which is assigned, chosen from the parameter
 * map.
 * \param value The value to which to assign the parameter.
 */
template <typename T>
auto assign(param_map_element& param, T&& value)
    -> decltype(param.assign(std::forward<T>(value))) {
  param.assign(std::forward<T>(value));
}

}  // namespace params
