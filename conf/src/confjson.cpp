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

#include "confjson.h"

#include <algorithm>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "boundary.h"
#include "coeffparser.h"
#include "definitions.h"
#include "gridinfo.h"
#include "initialconditionslib.h"
#include "io.h"
#include "macros.h"
#include "params.h"
#include "stencildefs.h"

using json = nlohmann::json;

// Helper functions to validate stencil point configurations
namespace symphas {

std::vector<size_t> parse_comma_separated_list(const std::string& str) {
  std::vector<size_t> result;
  std::stringstream ss(str);
  std::string item;

  while (std::getline(ss, item, ',')) {
    // Trim whitespace
    item.erase(0, item.find_first_not_of(" \t"));
    item.erase(item.find_last_not_of(" \t") + 1);
    if (!item.empty()) {
      result.push_back(std::stoul(item));
    }
  }
  return result;
}

// Dimension-aware validation functions
template <size_t dimension>
bool is_valid_laplacian_points(size_t points) {
  if constexpr (dimension == 1) {
#ifdef AVAILABLE_LAPLACIAN_POINTS_1D
    std::string available_str = STR(AVAILABLE_LAPLACIAN_POINTS_1D);
    auto available = parse_comma_separated_list(available_str);
    return std::find(available.begin(), available.end(), points) !=
           available.end();
#else
    return true;  // Default behavior if not configured
#endif
  } else if constexpr (dimension == 2) {
#ifdef AVAILABLE_LAPLACIAN_POINTS_2D
    std::string available_str = STR(AVAILABLE_LAPLACIAN_POINTS_2D);
    auto available = parse_comma_separated_list(available_str);
    return std::find(available.begin(), available.end(), points) !=
           available.end();
#else
    return true;
#endif
  } else if constexpr (dimension == 3) {
#ifdef AVAILABLE_LAPLACIAN_POINTS_3D
    std::string available_str = STR(AVAILABLE_LAPLACIAN_POINTS_3D);
    auto available = parse_comma_separated_list(available_str);
    return std::find(available.begin(), available.end(), points) !=
           available.end();
#else
    return true;
#endif
  }
  return true;
}

template <size_t dimension>
bool is_valid_gradlaplacian_points(size_t points) {
  if constexpr (dimension == 1) {
#ifdef AVAILABLE_GRADLAPLACIAN_POINTS_1D
    std::string available_str = STR(AVAILABLE_GRADLAPLACIAN_POINTS_1D);
    auto available = parse_comma_separated_list(available_str);
    return std::find(available.begin(), available.end(), points) !=
           available.end();
#else
    return true;
#endif
  } else if constexpr (dimension == 2) {
#ifdef AVAILABLE_GRADLAPLACIAN_POINTS_2D
    std::string available_str = STR(AVAILABLE_GRADLAPLACIAN_POINTS_2D);
    auto available = parse_comma_separated_list(available_str);
    return std::find(available.begin(), available.end(), points) !=
           available.end();
#else
    return true;
#endif
  } else if constexpr (dimension == 3) {
#ifdef AVAILABLE_GRADLAPLACIAN_POINTS_3D
    std::string available_str = STR(AVAILABLE_GRADLAPLACIAN_POINTS_3D);
    auto available = parse_comma_separated_list(available_str);
    return std::find(available.begin(), available.end(), points) !=
           available.end();
#else
    return true;
#endif
  }
  return true;
}

template <size_t dimension>
bool is_valid_bilaplacian_points(size_t points) {
  if constexpr (dimension == 1) {
#ifdef AVAILABLE_BILAPLACIAN_POINTS_1D
    std::string available_str = STR(AVAILABLE_BILAPLACIAN_POINTS_1D);
    auto available = parse_comma_separated_list(available_str);
    return std::find(available.begin(), available.end(), points) !=
           available.end();
#else
    return true;
#endif
  } else if constexpr (dimension == 2) {
#ifdef AVAILABLE_BILAPLACIAN_POINTS_2D
    std::string available_str = STR(AVAILABLE_BILAPLACIAN_POINTS_2D);
    auto available = parse_comma_separated_list(available_str);
    return std::find(available.begin(), available.end(), points) !=
           available.end();
#else
    return true;
#endif
  } else if constexpr (dimension == 3) {
#ifdef AVAILABLE_BILAPLACIAN_POINTS_3D
    std::string available_str = STR(AVAILABLE_BILAPLACIAN_POINTS_3D);
    auto available = parse_comma_separated_list(available_str);
    return std::find(available.begin(), available.end(), points) !=
           available.end();
#else
    return true;
#endif
  }
  return true;
}

template <size_t dimension>
std::vector<size_t> get_available_laplacian_points() {
  if constexpr (dimension == 1) {
#ifdef AVAILABLE_LAPLACIAN_POINTS_1D
    std::string available_str = STR(AVAILABLE_LAPLACIAN_POINTS_1D);
    return parse_comma_separated_list(available_str);
#else
    return {3};  // Default 1D
#endif
  } else if constexpr (dimension == 2) {
#ifdef AVAILABLE_LAPLACIAN_POINTS_2D
    std::string available_str = STR(AVAILABLE_LAPLACIAN_POINTS_2D);
    return parse_comma_separated_list(available_str);
#else
    return {5, 9};  // Default 2D
#endif
  } else if constexpr (dimension == 3) {
#ifdef AVAILABLE_LAPLACIAN_POINTS_3D
    std::string available_str = STR(AVAILABLE_LAPLACIAN_POINTS_3D);
    return parse_comma_separated_list(available_str);
#else
    return {7, 15, 19, 21, 27};  // Default 3D
#endif
  }
  return {};
}

template <size_t dimension>
std::vector<size_t> get_available_gradlaplacian_points() {
  if constexpr (dimension == 1) {
#ifdef AVAILABLE_GRADLAPLACIAN_POINTS_1D
    std::string available_str = STR(AVAILABLE_GRADLAPLACIAN_POINTS_1D);
    return parse_comma_separated_list(available_str);
#else
    return {4};
#endif
  } else if constexpr (dimension == 2) {
#ifdef AVAILABLE_GRADLAPLACIAN_POINTS_2D
    std::string available_str = STR(AVAILABLE_GRADLAPLACIAN_POINTS_2D);
    return parse_comma_separated_list(available_str);
#else
    return {6, 8, 12, 16};
#endif
  } else if constexpr (dimension == 3) {
#ifdef AVAILABLE_GRADLAPLACIAN_POINTS_3D
    std::string available_str = STR(AVAILABLE_GRADLAPLACIAN_POINTS_3D);
    return parse_comma_separated_list(available_str);
#else
    return {10, 12, 28, 36, 40};
#endif
  }
  return {};
}

template <size_t dimension>
std::vector<size_t> get_available_bilaplacian_points() {
  if constexpr (dimension == 1) {
#ifdef AVAILABLE_BILAPLACIAN_POINTS_1D
    std::string available_str = STR(AVAILABLE_BILAPLACIAN_POINTS_1D);
    return parse_comma_separated_list(available_str);
#else
    return {5};
#endif
  } else if constexpr (dimension == 2) {
#ifdef AVAILABLE_BILAPLACIAN_POINTS_2D
    std::string available_str = STR(AVAILABLE_BILAPLACIAN_POINTS_2D);
    return parse_comma_separated_list(available_str);
#else
    return {13, 17, 21};
#endif
  } else if constexpr (dimension == 3) {
#ifdef AVAILABLE_BILAPLACIAN_POINTS_3D
    std::string available_str = STR(AVAILABLE_BILAPLACIAN_POINTS_3D);
    return parse_comma_separated_list(available_str);
#else
    return {21, 25, 41, 52, 57};
#endif
  }
  return {};
}

// Non-template wrapper functions that infer dimension from context
// For now, assume 2D as default (most common case)
bool is_valid_laplacian_points(size_t points, size_t dimension) {
  if (dimension == 1) {
    return is_valid_laplacian_points<1>(points);
  } else if (dimension == 2) {
    return is_valid_laplacian_points<2>(points);
  } else if (dimension == 3) {
    return is_valid_laplacian_points<3>(points);
  }
  return false;
}

bool is_valid_gradlaplacian_points(size_t points, size_t dimension) {
  if (dimension == 1) {
    return is_valid_gradlaplacian_points<1>(points);
  } else if (dimension == 2) {
    return is_valid_gradlaplacian_points<2>(points);
  } else if (dimension == 3) {
    return is_valid_gradlaplacian_points<3>(points);
  }
  return false;
}

bool is_valid_bilaplacian_points(size_t points, size_t dimension) {
  if (dimension == 1) {
    return is_valid_bilaplacian_points<1>(points);
  } else if (dimension == 2) {
    return is_valid_bilaplacian_points<2>(points);
  } else if (dimension == 3) {
    return is_valid_bilaplacian_points<3>(points);
  }
  return false;
}

std::vector<size_t> get_available_laplacian_points(size_t dimension) {
  if (dimension == 1) {
    return get_available_laplacian_points<1>();
  } else if (dimension == 2) {
    return get_available_laplacian_points<2>();
  } else if (dimension == 3) {
    return get_available_laplacian_points<3>();
  }
  return {};
}

std::vector<size_t> get_available_gradlaplacian_points(size_t dimension) {
  if (dimension == 1) {
    return get_available_gradlaplacian_points<1>();
  } else if (dimension == 2) {
    return get_available_gradlaplacian_points<2>();
  } else if (dimension == 3) {
    return get_available_gradlaplacian_points<3>();
  }
  return {};
}

std::vector<size_t> get_available_bilaplacian_points(size_t dimension) {
  if (dimension == 1) {
    return get_available_bilaplacian_points<1>();
  } else if (dimension == 2) {
    return get_available_bilaplacian_points<2>();
  } else if (dimension == 3) {
    return get_available_bilaplacian_points<3>();
  }
  return {};
}

}  // namespace symphas

template <typename T>
T parse_definition(json const& entry,
                   std::map<std::string, json> const& definitions) {
  if (entry.is_string()) {
    std::string str = entry.get<std::string>();
    if (str.size() >= 4 && str[0] == '$' && str[1] == '{' &&
        str.back() == '}') {
      std::string def_key = str.substr(2, str.size() - 3);
      if (definitions.count(def_key) > 0) {
        return definitions.at(def_key).get<T>();
      } else {
        throw std::runtime_error("Key not found in definitions: " + def_key);
      }
    }
  }
  return entry.get<T>();
}

template <typename T>
T parse_definition(json const& entry, std::string key, T const& fallback,
                   std::map<std::string, json> const& definitions) {
  if (!entry.contains(key)) {
    return fallback;
  }
  auto sub_entry = entry[key];
  if (sub_entry.is_string()) {
    std::string str = sub_entry.get<std::string>();
    if (str.size() >= 4 && str[0] == '$' && str[1] == '{' &&
        str.back() == '}') {
      std::string def_key = str.substr(2, str.size() - 3);
      if (definitions.count(def_key) > 0) {
        return definitions.at(def_key).get<T>();
      } else {
        throw std::runtime_error("Key not found in definitions: " + def_key);
      }
    }
  }
  return sub_entry.get<T>();
}

symphas::interval_element_type parse_interval_data(
    json const& entry, std::map<std::string, json> const& definitions) {
  if (entry.contains("start") && entry.contains("end") &&
      entry.contains("points") && entry.contains("width")) {
    throw std::runtime_error(
        "Interval is overdefined. Please specify only three of the four "
        "parameters: start, end, points, width.");
  }

  // Get parameters from JSON, or use default values
  double start = parse_definition(entry, "start", 0.0, definitions);
  double end = parse_definition(entry, "end", 1.0, definitions);
  int points = parse_definition(entry, "points", 0, definitions);
  double width = parse_definition(entry, "width", 0.0, definitions);

  symphas::interval_element_type interval;
  if (!points && !width) {
    interval.set_domain(start, end);
  } else if (!points) {
    interval.set_domain(start, end, width);
  } else if (!width) {
    interval.set_domain_count(start, end, points);
  } else {
    interval.set_domain(start, end, width);
    interval.set_domain_count_from_r(points, (!entry.contains("end")) ? 0 : 1);
  }
  interval.interval_to_domain();
  return interval;
}

symphas::b_element_type parse_boundary_data(
    json const& entry, std::map<std::string, json> const& definitions) {
  BoundaryType type = BoundaryType::DEFAULT;
  symphas::boundary_tag_list tags(BoundaryTag::NONE, BoundaryTag::NONE);

  // Parse boundary type - handle both explicit type and implicit DEFAULT
  if (entry.contains("type")) {
    std::string type_str = entry["type"].get<std::string>();

    // Parse boundary type and tags from space-separated string (legacy format)
    std::istringstream iss(type_str);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token) {
      tokens.push_back(token);
    }

    if (!tokens.empty()) {
      // First token is the boundary type (or first tag if DEFAULT is implicit)
      type = symphas::boundary_from_str(tokens[0].c_str());

      size_t tag_start = 0;
      if (type != BoundaryType::DEFAULT) {
        // Explicit boundary type, tags start from index 1
        tag_start = 1;
      } else {
        // Check if first token is actually a tag (implicit DEFAULT)
        if (symphas::boundary_tag_from_str(tokens[0].c_str()) !=
            BoundaryTag::NONE) {
          tag_start = 0;  // All tokens are tags
        } else {
          tag_start = 1;  // First token was explicit DEFAULT
        }
      }

      // Parse up to 2 boundary tags from type string
      for (size_t i = tag_start; i < tokens.size() && i - tag_start < 2; ++i) {
        BoundaryTag tag = symphas::boundary_tag_from_str(tokens[i].c_str());
        if (tag != BoundaryTag::NONE) {
          tags[i - tag_start] = tag;
        }
      }
    }
  }

  // Parse separate tags array (new format) - this takes precedence over type
  // string tags
  if (entry.contains("tags")) {
    // Reset tags if separate tags array is provided
    tags = symphas::boundary_tag_list{};

    auto const& tag_array = entry["tags"];
    size_t tag_count = 0;

    for (auto const& tag_entry : tag_array) {
      if (tag_count >= 2) break;  // Maximum 2 tags supported

      std::string tag_str = tag_entry.get<std::string>();
      BoundaryTag tag = symphas::boundary_tag_from_str(tag_str.c_str());

      if (tag != BoundaryTag::NONE) {
        tags[tag_count] = tag;
        tag_count++;
      }
    }

    // If no explicit type was provided and tags are present, default to DEFAULT
    if (!entry.contains("type")) {
      type = BoundaryType::DEFAULT;
    }
  }

  // Handle case where neither type nor tags are provided
  if (!entry.contains("type") && !entry.contains("tags")) {
    return {BoundaryType::NONE, nullptr, 0};
  }

  // Parse parameters
  std::vector<double> parameters;
  if (entry.contains("parameters")) {
    for (auto const& param : entry["parameters"]) {
      parameters.push_back(parse_definition<double>(param, definitions));
    }
  }
  int argc = int(parameters.size());

  // Check if any tags were specified
  if (tags[0] != BoundaryTag::NONE || tags[1] != BoundaryTag::NONE) {
    return {type, tags, parameters.data(), argc};
  }

  return {type, parameters.data(), argc};
}

symphas::init_data_type parse_initial_conditions(
    json const& entry, std::map<std::string, json> const& definitions) {
  symphas::init_data_type init_data;

  for (auto const& [axis_str, init_config] : entry.items()) {
    Axis axis = symphas::axis_from_str(axis_str.c_str());

    Inside init_type =
        symphas::in_from_str(init_config["type"].get<std::string>().c_str());

    std::vector<double> parameters;
    if (init_config.contains("parameters")) {
      for (auto const& param : init_config["parameters"]) {
        parameters.push_back(parse_definition<double>(param, definitions));
      }
    }

    symphas::init_data_parameters params(parameters.data(), parameters.size());
    symphas::init_entry_type init_entry(init_type, params);

    init_data[axis] = init_entry;
  }

  return init_data;
}

void JsonConfManager::parse_single_initial_condition(
    const json& init_config, std::map<std::string, json> const& definitions,
    int field_idx) {
  // Determine axis (default to NONE for all axes)
  Axis axis = Axis::NONE;
  if (init_config.contains("axis")) {
    std::string axis_str = parse_value_or_definition<std::string>(
        init_config["axis"], definitions);
    if (axis_str == "X")
      axis = Axis::X;
    else if (axis_str == "Y")
      axis = Axis::Y;
    else if (axis_str == "Z")
      axis = Axis::Z;
    else if (axis_str == "R")
      axis = Axis::R;
    else if (axis_str == "T")
      axis = Axis::T;
    else if (axis_str == "S")
      axis = Axis::S;
    else if (axis_str == "C")
      axis = Axis::C;
    else if (axis_str == "I")
      axis = Axis::I;
    else if (axis_str == "NONE")
      axis = Axis::NONE;
    else {
      fprintf(SYMPHAS_ERR,
              "Invalid axis specification for initial condition: '%s'\n",
              axis_str.c_str());
      return;
    }
  }

  // Parse the initial condition type
  if (!init_config.contains("type")) {
    fprintf(SYMPHAS_ERR, "Initial condition must specify a 'type'\n");
    return;
  }

  std::string type_str =
      parse_value_or_definition<std::string>(init_config["type"], definitions);
  Inside init_type = symphas::in_from_str(type_str.c_str());

  // Check for NONE type (invalid)
  if (init_type == Inside::NONE) {
    fprintf(SYMPHAS_ERR,
            "Invalid initial condition type: '%s'. A valid initial condition "
            "must be provided.\n",
            type_str.c_str());
    return;
  }

  symphas::init_entry_type init_entry;
  init_entry.in = init_type;

  // Handle different types of initial conditions
  if (init_type == Inside::EXPRESSION) {
    // Expression type - requires expression name
    if (!init_config.contains("expression")) {
      fprintf(
          SYMPHAS_ERR,
          "EXPRESSION type initial condition requires 'expression' field\n");
      return;
    }
    std::string expr_name = parse_value_or_definition<std::string>(
        init_config["expression"], definitions);
    init_entry.expr_data = expr_name.c_str();

    // Handle coefficients if provided
    if (init_config.contains("coefficients")) {
      std::vector<double> coeffs;
      for (auto const& coeff : init_config["coefficients"]) {
        coeffs.push_back(parse_definition<double>(coeff, definitions));
      }
      // Note: This would need integration with model_settings coefficients
    }

  } else if (init_type == Inside::FILE) {
    // File type - requires filename
    if (!init_config.contains("file")) {
      fprintf(SYMPHAS_ERR,
              "FILE type initial condition requires 'file' field\n");
      return;
    }
    std::string filename = parse_value_or_definition<std::string>(
        init_config["file"], definitions);
    int index = 0;
    if (init_config.contains("index")) {
      index = parse_value_or_definition<int>(init_config["index"], definitions);
    }
    init_entry.file = {filename.c_str(), index};

  } else {
    // Standard types with parameters
    // Set default parameters based on type
    switch (init_type) {
      case Inside::GAUSSIAN:
        init_entry.data.gp[0] = 0;     // mean
        init_entry.data.gp[1] = 0.25;  // std deviation
        break;
      case Inside::UNIFORM:
        init_entry.data.gp[0] = -1;  // minimum value
        init_entry.data.gp[1] = 1;   // maximum value
        break;
      case Inside::CAPPED:
        init_entry.data.gp[0] = -1;  // minimum value
        init_entry.data.gp[1] = 1;   // maximum value
        break;
      case Inside::CONSTANT:
        init_entry.data.gp[0] = 0;  // the constant value
        break;
      case Inside::CIRCLE:
        init_entry.data.gp[0] = 2;  // ratio of the x width
        init_entry.data.gp[1] = 2;  // ratio of the y width
        break;
      case Inside::VORONOI:
        init_entry.data.gp[0] = 10;  // number of crystals
        init_entry.data.gp[1] = -1;  // lower range
        init_entry.data.gp[2] = 1;   // upper range
        break;
      default:
        init_entry.data = symphas::init_data_parameters::one();
        break;
    }

    // Override defaults with provided parameters
    if (init_config.contains("parameters")) {
      int param_idx = 0;
      for (auto const& param : init_config["parameters"]) {
        init_entry.data.gp[param_idx++] =
            parse_value_or_definition<double>(param, definitions);
      }
    }
  }
  // Set the initial condition for the specified field and axis
  // Create or update the axis mapping for this field
  symphas::init_data_type field_init_data;
  if (domain_settings.tdata_len > field_idx) {
    field_init_data = domain_settings.tdata[field_idx];
  }
  field_init_data[axis] = init_entry;
  domain_settings.set_initial_condition(field_init_data, field_idx);
}

template <typename T>
T JsonConfManager::parse_value_or_definition(
    const json& entry, std::map<std::string, json> const& definitions) {
  return parse_definition<T>(entry, definitions);
}

template <typename T>
T JsonConfManager::parse_value_or_definition(
    const json& entry, T const& fallback,
    std::map<std::string, json> const& definitions) {
  return parse_definition<T>(entry, "", fallback, definitions);
}

std::string JsonConfManager::build_save_spec_from_json(
    const json& save, std::map<std::string, json> const& definitions) {
  std::string spec;

  // Handle LIST type - comma-separated indices
  if (save.contains("indices")) {
    bool first = true;
    for (auto const& index : save["indices"]) {
      if (!first) spec += ",";
      spec += std::to_string(parse_definition<int>(index, definitions));
      first = false;
    }
    return spec;
  }

  // Handle other types with interval value (renamed from base)
  double interval = 1.0;
  if (save.contains("interval")) {
    interval = parse_definition<double>(save["interval"], definitions);
    spec += std::to_string(interval);
  }

  // Add save type if specified
  if (save.contains("type")) {
    std::string type = parse_definition<std::string>(save["type"], definitions);
    if (type != "DEFAULT") {
      spec += " " + type;
    }
  }

  return spec;
}

void JsonConfManager::parse_simulation_settings(
    const json& config_json, std::map<std::string, json> const& definitions) {
  if (!config_json.contains("simulation")) return;

  auto const& sim = config_json["simulation"];

  // Parse dimension
  if (sim.contains("dimension")) {
    simulation_settings.dimension =
        parse_value_or_definition<size_t>(sim["dimension"], definitions);
  }

  // Parse time steps
  if (sim.contains("time_step")) {
    double dt =
        parse_value_or_definition<double>(sim["time_step"], definitions);
    double time = parse_value_or_definition<double>(
        sim.value("time", json(1.0)), definitions);
    simulation_settings.set_time_step(dt, time);
  }
  if (sim.contains("time_steps")) {
    auto const& time_steps = sim["time_steps"];
    std::vector<double> dts, times;
    iter_type stop = 0;
    for (auto const& step : time_steps) {
      double dt = parse_value_or_definition<double>(step["dt"], definitions);
      dts.push_back(dt);

      // Support both "time" and "iterations" specifications
      if (step.contains("time")) {
        double next_time =
            parse_value_or_definition<double>(step["time"], definitions);
        times.push_back(next_time);
        stop += next_time / dt;
      } else if (step.contains("iterations")) {
        double iterations =
            parse_value_or_definition<double>(step["iterations"], definitions);
        times.push_back(iterations * dt);
        stop += iterations;
      } else {
        throw std::runtime_error(
            "time_steps entry must contain either 'time' or 'iterations'");
      }
    }
    simulation_settings.set_time_steps(dts.data(), times.data(), dts.size());
    simulation_settings.save.set_stop(stop);
    simulation_settings.save.set_params(stop);
  }  // Parse stencil parameters
  if (sim.contains("stencil")) {
    auto const& stencil = sim["stencil"];

    // Parse order (accuracy)
    if (stencil.contains("order")) {
      size_t order =
          parse_value_or_definition<size_t>(stencil["order"], definitions);
      simulation_settings.parse_stencil_accuracy(std::to_string(order).c_str());
    } else {
      // Use default order
      simulation_settings.parse_stencil_accuracy("@");
    }

    // Parse PTL (laplacian stencil points)
    if (stencil.contains("ptl")) {
      size_t ptl =
          parse_value_or_definition<size_t>(stencil["ptl"], definitions);

      // Validate that the requested PTL is available
      if (!symphas::is_valid_laplacian_points(ptl,
                                              simulation_settings.dimension)) {
        auto available = symphas::get_available_laplacian_points(
            simulation_settings.dimension);
        std::string available_str = "";
        for (size_t i = 0; i < available.size(); ++i) {
          if (i > 0) available_str += ", ";
          available_str += std::to_string(available[i]);
        }
        throw std::runtime_error("Invalid laplacian stencil points (" +
                                 std::to_string(ptl) +
                                 "). Available options: " + available_str);
      }

      simulation_settings.parse_stencil(2, std::to_string(ptl).c_str());
    } else {
      // Use default PTL
      simulation_settings.parse_stencil(2, "@");
    }

    // Parse PTG (gradlaplacian stencil points)
    if (stencil.contains("ptg")) {
      size_t ptg =
          parse_value_or_definition<size_t>(stencil["ptg"], definitions);

      // Validate that the requested PTG is available
      if (!symphas::is_valid_gradlaplacian_points(
              ptg, simulation_settings.dimension)) {
        auto available = symphas::get_available_gradlaplacian_points(
            simulation_settings.dimension);
        std::string available_str = "";
        for (size_t i = 0; i < available.size(); ++i) {
          if (i > 0) available_str += ", ";
          available_str += std::to_string(available[i]);
        }
        throw std::runtime_error("Invalid gradlaplacian stencil points (" +
                                 std::to_string(ptg) +
                                 "). Available options: " + available_str);
      }

      simulation_settings.parse_stencil(3, std::to_string(ptg).c_str());
    } else {
      // Use default PTG
      simulation_settings.parse_stencil(3, "@");
    }

    // Parse PTB (bilaplacian stencil points)
    if (stencil.contains("ptb")) {
      size_t ptb =
          parse_value_or_definition<size_t>(stencil["ptb"], definitions);

      // Validate that the requested PTB is available
      if (!symphas::is_valid_bilaplacian_points(
              ptb, simulation_settings.dimension)) {
        auto available = symphas::get_available_bilaplacian_points(
            simulation_settings.dimension);
        std::string available_str = "";
        for (size_t i = 0; i < available.size(); ++i) {
          if (i > 0) available_str += ", ";
          available_str += std::to_string(available[i]);
        }
        throw std::runtime_error("Invalid bilaplacian stencil points (" +
                                 std::to_string(ptb) +
                                 "). Available options: " + available_str);
      }

      simulation_settings.parse_stencil(4, std::to_string(ptb).c_str());
    } else {
      // Use default PTB
      simulation_settings.parse_stencil(4, "@");
    }
  }  // Parse save parameters
  if (sim.contains("save")) {
    auto const& save = sim["save"];

    // Parse structured save format
    std::string save_spec = build_save_spec_from_json(save, definitions);
    simulation_settings.parse_save(save_spec.c_str());

    if (save.contains("save_initial")) {
      bool save_init =
          parse_value_or_definition<bool>(save["save_initial"], definitions);
      simulation_settings.parse_save_init(save_init ? "YES" : "NO");
    }
  }

  // Parse solver variation
  if (sim.contains("solver_variation")) {
    int solver_variation =
        parse_value_or_definition<int>(sim["solver_variation"], definitions);
    simulation_settings.parse_solver_type(
        std::to_string(solver_variation).c_str());
  }
}

void JsonConfManager::parse_model_settings(
    const json& config_json, std::map<std::string, json> const& definitions) {
  if (!config_json.contains("model")) return;

  auto const& model = config_json["model"];

  // Parse model name
  if (model.contains("name")) {
    std::string model_name =
        parse_value_or_definition<std::string>(model["name"], definitions);
    model_settings.set_model_name(model_name.c_str());
  }

  // Parse coefficients: support 'coefficients_file' or inline 'coefficients'
  std::vector<double> coeffs;
  if (model.contains("coefficients_file")) {
    std::string coeff_file = parse_value_or_definition<std::string>(
        model["coefficients_file"], definitions);
    // If not absolute, resolve relative to config_dir
    if (!(coeff_file.size() > 0 &&
          (coeff_file[0] == '/' || coeff_file[0] == '\\' ||
           (coeff_file.size() > 1 && coeff_file[1] == ':')))) {
      coeff_file = this->config_dir + "/" + coeff_file;
    }
    coeffs = parse_pfc_coeff_file(coeff_file);
  } else if (model.contains("coefficients")) {
    if (model["coefficients"].is_array()) {
      for (auto const& coeff : model["coefficients"]) {
        coeffs.push_back(parse_value_or_definition<double>(coeff, definitions));
      }
    } else if (model["coefficients"].is_object()) {
      // Handle named coefficients with support for "c1", "c2", etc. syntax
      std::map<int, double> coefficient_map;
      int max_index = -1;
      for (auto const& [name, value] : model["coefficients"].items()) {
        double coeff_value =
            parse_value_or_definition<double>(value, definitions);
        std::regex c_pattern("^c(\\d+)$");
        std::smatch match;
        if (std::regex_match(name, match, c_pattern)) {
          int index = std::stoi(match.str(1)) - 1;
          if (index >= 0) {
            coefficient_map[index] = coeff_value;
            max_index = std::max(max_index, index);
          }
        } else {
          coeffs.push_back(coeff_value);
        }
      }
      if (!coefficient_map.empty()) {
        int array_size =
            std::max(max_index + 1, static_cast<int>(coeffs.size()));
        std::vector<double> indexed_coeffs(array_size, DEFAULT_COEFF_VALUE);
        for (size_t i = 0; i < coeffs.size(); ++i) {
          indexed_coeffs[i] = coeffs[i];
        }
        for (auto const& [index, value] : coefficient_map) {
          if (index < array_size) {
            indexed_coeffs[index] = value;
          }
        }
        coeffs = std::move(indexed_coeffs);
      }
    }
  }
  if (!coeffs.empty()) {
    model_settings.coeff = new double[coeffs.size()];
    std::copy(coeffs.begin(), coeffs.end(), model_settings.coeff);
    model_settings.coeff_len = coeffs.size();
  }

  // Parse field specifications
  if (model.contains("fields")) {
    auto const& fields = model["fields"];
    if (fields.is_array()) {
      model_settings.num_fields_len = fields.size();
      model_settings.num_fields = new len_type[model_settings.num_fields_len];
      model_settings.modifiers = new char*[model_settings.num_fields_len];

      for (size_t i = 0; i < fields.size(); ++i) {
        if (fields[i].is_number()) {
          model_settings.num_fields[i] =
              parse_value_or_definition<len_type>(fields[i], definitions);
          model_settings.modifiers[i] = new char[1];
          model_settings.modifiers[i][0] = '\0';
        } else if (fields[i].is_object()) {
          model_settings.num_fields[i] = parse_value_or_definition<len_type>(
              fields[i]["count"], definitions);
          std::string modifier = parse_value_or_definition<std::string>(
              fields[i].value("modifier", json("")), definitions);
          model_settings.modifiers[i] = new char[modifier.length() + 1];
          std::strcpy(model_settings.modifiers[i], modifier.c_str());
        }
      }
    }
  }
}

void JsonConfManager::parse_domain_settings(
    const json& config_json, std::map<std::string, json> const& definitions) {
  if (!config_json.contains("domains")) return;

  auto const& domains = config_json["domains"];
  // Get domain mappings from model
  std::map<std::string, std::string> domain_mappings;
  if (config_json.contains("model") &&
      config_json["model"].contains("domains")) {
    auto const& model = config_json["model"];
    auto const& mappings = model["domains"];
    for (auto const& [field_index, domain_config] : mappings.items()) {
      auto domain_name = domain_config["domain"];
      domain_mappings[field_index] =
          parse_value_or_definition<std::string>(domain_name, definitions);

      if (domain_config.contains("initial_condition")) {
        auto const& init_condition = domain_config["initial_condition"];
        int field_idx = std::stoi(field_index);

        // Handle array of initial conditions for different axes
        if (init_condition.is_array()) {
          for (auto const& init_config : init_condition) {
            parse_single_initial_condition(init_config, definitions, field_idx);
          }
        } else {
          // Single initial condition (applies to all axes)
          parse_single_initial_condition(init_condition, definitions,
                                         field_idx);
        }
      }  // Parse boundaries (new axis-based structure)
    }
  }

  // Parse each domain
  for (auto const& [field_index, domain_name] : domain_mappings) {
    if (!domains.contains(domain_name)) {
      continue;  // Skip if domain not found
    }

    auto const& domain =
        domains[domain_name];  // Parse initial conditions (enhanced with axis
                               // support, expressions, and files)
    if (domain.contains("boundaries")) {
      auto const& boundaries = domain["boundaries"];

      if (boundaries.is_object()) {
        symphas::b_data_type bdata;

        // Parse boundaries by axis (x, y, z)
        for (auto const& [axis_str, axis_boundaries] : boundaries.items()) {
          if (axis_str == "all") {
            // Apply to all sides
            auto boundary = parse_boundary_data(axis_boundaries, definitions);
            for (size_t i = 0; i < simulation_settings.dimension * 2; ++i) {
              Side side = symphas::index_to_side(i);
              bdata[side] = boundary;
            }
          } else {
            // Parse left and right boundaries for this axis
            Axis axis = symphas::axis_from_str(axis_str.c_str());

            if (axis_boundaries.contains("left")) {
              Side left_side;
              if (axis == Axis::X)
                left_side = Side::LEFT;
              else if (axis == Axis::Y)
                left_side = Side::BOTTOM;
              else if (axis == Axis::Z)
                left_side = Side::BACK;

              auto boundary =
                  parse_boundary_data(axis_boundaries["left"], definitions);
              bdata[left_side] = boundary;
            }

            if (axis_boundaries.contains("right")) {
              Side right_side;
              if (axis == Axis::X)
                right_side = Side::RIGHT;
              else if (axis == Axis::Y)
                right_side = Side::TOP;
              else if (axis == Axis::Z)
                right_side = Side::FRONT;

              auto boundary =
                  parse_boundary_data(axis_boundaries["right"], definitions);
              bdata[right_side] = boundary;
            }
          }
        }

        // Set boundaries for this field
        int field_idx = std::stoi(field_index);
        domain_settings.set_boundary(bdata, field_idx);
      }
    }
    if (domain.contains("intervals")) {
      auto const& intervals = domain["intervals"];
      if (intervals.is_object()) {
        symphas::interval_data_type domain_intervals;
        for (auto const& [axis_str, interval_config] : intervals.items()) {
          Axis axis = symphas::axis_from_str(axis_str.c_str());
          if (interval_config.is_string()) {
            auto interval_name = interval_config.get<std::string>();
            auto const& all_intervals = config_json["intervals"];
            if (all_intervals.contains(interval_name)) {
              auto interval_data = parse_interval_data(
                  all_intervals[interval_name], definitions);
              domain_intervals[axis] = interval_data;
            } else {
              throw std::runtime_error("Interval definition not found: " +
                                       interval_name);
            }
          } else {
            auto interval_data =
                parse_interval_data(interval_config, definitions);
            domain_intervals[axis] = interval_data;
          }
        }
        // Set intervals for the first domain (index 0)
        int field_idx = std::stoi(field_index);
        domain_settings.set_interval(domain_intervals, field_idx);
      }
    }
  }
}

void JsonConfManager::parse_name_settings(
    const json& config_json, std::map<std::string, json> const& definitions) {
  if (config_json.contains("title")) {
    std::string title = parse_value_or_definition<std::string>(
        config_json["title"], definitions);
    name_settings.set_title(title.c_str());
  }
  if (config_json.contains("names")) {
    auto const& names = config_json["names"];
    if (names.is_array()) {
      for (size_t i = 0; i < names.size(); ++i) {
        std::string name =
            parse_value_or_definition<std::string>(names[i], definitions);
        name_settings.set_name(name.c_str(), i);
      }
    }
  }
}

void JsonConfManager::parse_directory_settings(
    const json& config_json, std::map<std::string, json> const& definitions) {
  if (config_json.contains("output")) {
    auto const& output = config_json["output"];
    if (output.contains("title")) {
      std::string title = parse_value_or_definition<std::string>(
          output.value("title", config_json.value("title", json(""))),
          definitions);
      name_settings.set_title(title.c_str());
    } else {
      name_settings.set_title(model_settings.model);
    }
    if (output.contains("directory")) {
      std::string dir = parse_value_or_definition<std::string>(
          output["directory"], definitions);
      directory_settings.set_directory(dir.c_str(), name_settings.title);
    } else {
      directory_settings.set_directory("symphas_output", name_settings.title);
    }
  }
}

JsonConfManager::JsonConfManager(const char* config_file) : SymPhasSettings{} {
  std::ifstream file(config_file);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open configuration file: " +
                             std::string(config_file));
  }
  std::string config_path(config_file);
  std::string config_dir = ".";
  size_t last_slash = config_path.find_last_of("/\\");
  if (last_slash != std::string::npos) {
    config_dir = config_path.substr(0, last_slash);
  }
  this->config_dir = config_dir;

  json config;
  file >> config;

  std::ifstream schema_file("schema.json");
  json schema;
  schema_file >> schema;
  try {
    nlohmann::json_schema::json_validator validator;
    validator.set_root_schema(schema);
    validator.validate(config);
    fprintf(SYMPHAS_LOG, "Configuration validated successfully against schema");
  } catch (const std::exception& e) {
    fprintf(SYMPHAS_ERR, "Configuration validation failed: %s\n", e.what());
    throw std::runtime_error("Invalid configuration file");
  }

  std::map<std::string, json> definitions;
  if (config.contains("definitions")) {
    definitions = config["definitions"].get<std::map<std::string, json>>();
  }

  // Parse all configuration sections
  parse_simulation_settings(config, definitions);
  parse_model_settings(config, definitions);
  parse_domain_settings(config, definitions);
  parse_name_settings(config, definitions);
  parse_directory_settings(config, definitions);
  SymPhasSettings::conf_file = config_file;

  make_directories();
}

void JsonConfManager::write(const char* savedir, const char* name) const {
  json config;

  // Add simulation settings
  json simulation;
  simulation["dimension"] = simulation_settings.dimension;

  // Add time steps - handle multiple time steps correctly
  json time_steps = json::array();

  if (simulation_settings.dt_list.get_num_time_steps() > 0) {
    auto dt_values = simulation_settings.dt_list.get_time_steps();
    auto time_values = simulation_settings.dt_list.get_times_of_steps();

    for (size_t i = 0; i < simulation_settings.dt_list.get_num_time_steps();
         ++i) {
      json time_step;
      time_step["dt"] = dt_values[i];

      if (i == 0) {
        // First time step uses the time value directly
        time_step["time"] = time_values[i];
      } else {
        // Subsequent time steps: calculate iterations from time difference
        double time_diff = time_values[i] - time_values[i - 1];
        iter_type iterations =
            static_cast<iter_type>(time_diff / dt_values[i - 1]);
        time_step["iterations"] = iterations;
      }

      time_steps.push_back(time_step);
    }
  } else {
    // Fallback: create a single time step with default values
    json time_step;
    time_step["dt"] = 1.0;
    time_step["iterations"] = simulation_settings.save.get_stop();
    time_steps.push_back(time_step);
  }

  simulation["time_steps"] = time_steps;

  // Add stencil settings
  json stencil;
  stencil["order"] = simulation_settings.stp.ord;
  stencil["ptl"] = simulation_settings.stp.ptl;
  stencil["ptb"] = simulation_settings.stp.ptb;
  stencil["ptg"] = simulation_settings.stp.ptg;
  simulation["stencil"] = stencil;

  simulation["solver_variation"] = simulation_settings.stp.type;
  config["simulation"] = simulation;

  // Add model settings
  json model;
  if (model_settings.model) {
    model["name"] = model_settings.model;
  }

  // Add intervals
  json intervals;
  for (size_t d = 0; d < simulation_settings.dimension; ++d) {
    Axis axis = symphas::index_to_axis(d);
    const char* axis_char = symphas::str_from_axis(axis);

    if (axis_char) {
      char axis_key[]{*axis_char, '\0'};
      symphas::lib::to_lower(axis_key);
      json interval;

      // Get interval data from the first domain (index 0)
      if (domain_settings.intervals &&
          domain_settings.intervals[0].find(axis) !=
              domain_settings.intervals[0].end()) {
        auto& interval_elem = domain_settings.intervals[0].at(axis);
        interval["start"] = interval_elem.domain_left();
        interval["points"] = interval_elem.get_count();
        interval["width"] = interval_elem.width();
      } else {
        // Fallback to simple array access if interval data structure is
        // different
        interval["start"] = 0.0;
        interval["points"] = 64;
        interval["width"] = 1.0;
      }

      intervals[axis_key] = interval;
    }
  }
  model["intervals"] = intervals;

  // Add coefficients if present
  if (model_settings.coeff_len > 0 && model_settings.coeff) {
    json coefficients;
    for (size_t i = 0; i < model_settings.coeff_len; ++i) {
      std::string coeff_name = "c" + std::to_string(i + 1);
      coefficients[coeff_name] = model_settings.coeff[i];
    }
    model["coefficients"] = coefficients;
  }

  model["domains"] = json::object();

  // Add domains (simplified - creates a single domain for each field)
  json domains;
  for (size_t field_idx = 0; field_idx < domain_settings.tdata_len;
       ++field_idx) {
    std::string domain_name = "field" + std::to_string(field_idx) + "_domain";
    model["domains"][std::to_string(field_idx)] = domain_name;

    json domain;

    // Add initial conditions
    if (domain_settings.tdata && field_idx < domain_settings.tdata_len) {
      auto& tdata = domain_settings.tdata[field_idx];

      // Find initial condition for this field
      if (!tdata.empty()) {
        // Get the first initial condition (usually for Axis::NONE)
        auto it = tdata.begin();
        json ic;

        if (it->second.in == Inside::FILE) {
          ic["type"] = "FILE";
          ic["file"] = it->second.file.get_name();
        } else if (it->second.in == Inside::EXPRESSION) {
          ic["type"] = "EXPRESSION";
          ic["expression"] = it->second.expr_data.get_name();
          if (it->second.expr_data.get_num_coeff() > 0) {
            json coefficients = json::array();
            for (iter_type i = 0; i < it->second.expr_data.get_num_coeff();
                 ++i) {
              coefficients.push_back(it->second.expr_data.get_coeff()[i]);
            }
            ic["coefficients"] = coefficients;
          }
        } else {
          // Convert Inside enum to string
          const char* init_type_str = symphas::str_from_in(it->second.in);
          if (init_type_str) {
            ic["type"] = init_type_str;

            // Add parameters if available
            if (it->second.data.gp && it->second.data.N > 0) {
              json parameters = json::array();
              for (size_t j = 0; j < it->second.data.N; ++j) {
                parameters.push_back(it->second.data.gp[j]);
              }
              ic["parameters"] = parameters;
            }
          }
        }
        domain["initial_conditions"] = ic;
      }
    }

    // Add boundaries
    json boundaries;
    if (domain_settings.bdata && field_idx < domain_settings.bdata_len) {
      auto& bdata = domain_settings.bdata[field_idx];

      for (size_t d = 0; d < simulation_settings.dimension; ++d) {
        Axis axis = symphas::index_to_axis(d);
        const char* axis_char = symphas::str_from_axis(axis);
        if (axis_char) {
          json axis_boundaries;

          // For 2D: X axis uses LEFT/RIGHT, Y axis uses BOTTOM/TOP
          // For 3D: X axis uses LEFT/RIGHT, Y axis uses BOTTOM/TOP, Z axis uses
          // BACK/FRONT
          Side left_side, right_side;
          if (d == 0) {  // X axis
            left_side = Side::LEFT;
            right_side = Side::RIGHT;
          } else if (d == 1) {  // Y axis
            left_side = Side::BOTTOM;
            right_side = Side::TOP;
          } else {  // Z axis
            left_side = Side::BACK;
            right_side = Side::FRONT;
          }

          // Left boundary
          if (bdata.find(left_side) != bdata.end()) {
            json left_boundary;
            const char* left_type =
                symphas::str_from_boundary(bdata.at(left_side).type);
            if (left_type) {
              left_boundary["type"] = left_type;
            }
            axis_boundaries["left"] = left_boundary;
          }

          // Right boundary
          if (bdata.find(right_side) != bdata.end()) {
            json right_boundary;
            const char* right_type =
                symphas::str_from_boundary(bdata.at(right_side).type);
            if (right_type) {
              right_boundary["type"] = right_type;
            }
            axis_boundaries["right"] = right_boundary;
          }

          char axis_key[]{*axis_char, '\0'};
          symphas::lib::to_lower(axis_key);
          boundaries[axis_key] = axis_boundaries;
        }
      }
    }
    domain["boundaries"] = boundaries;
    domains[domain_name] = domain;
  }
  config["model"] = model;
  config["domains"] = domains;

  // Add save settings if present
  json save;
  save["interval"] = simulation_settings.save.get_base();
  save["save_initial"] = simulation_settings.save.get_init_flag();
  if (simulation_settings.save.get_type() != SaveType::DEFAULT) {
    const char* save_type_str =
        symphas::io::get_save_str(simulation_settings.save.get_type());
    if (save_type_str) {
      save["type"] = save_type_str;
    }
  }
  config["simulation"]["save"] = save;

  // Add raw parameters as definitions if present
  if (!params::rawparams.empty()) {
    json parameters;
    for (const auto& param : params::rawparams) {
      if (!param.second.empty()) {
        // Try to parse as number first, otherwise keep as string
        try {
          double num_val = std::stod(param.second);
          parameters[param.first] = num_val;
        } catch (...) {
          parameters[param.first] = param.second;
        }
      } else {
        parameters[param.first] = nullptr;
      }
    }
    config["parameters"] = parameters;
  }

  json output;
  output["directory"] = directory_settings.root_dir;
  if (name_settings.title) {
    output["title"] = name_settings.title;
  }
  config["output"] = output;

  json names;
  if (name_settings.names) {
    for (size_t i = 0; i < name_settings.names_len; ++i) {
      names.push_back(name_settings.names[i]);
    }
  }
  if (!names.empty()) {
    config["names"] = names;
  }

  // Write to file
  char bname[BUFFER_LENGTH]{};
  snprintf(bname, BUFFER_LENGTH, "%s/%s.json", savedir, name);

  std::ofstream file(bname);
  if (!file.is_open()) {
    fprintf(SYMPHAS_ERR, "error opening write JSON configuration file '%s'\n",
            bname);
    exit(ERR_CODE_FILE_OPEN);
  }

  file << config.dump(2);  // Pretty print with 2-space indentation
  file.close();
}
