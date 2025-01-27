
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

#include <string>

using json = nlohmann::json;

template <typename T>
T parse_definition(json const& entry, std::string key,
                   std::map<std::string, json> const& definitions) {
  std::string str = entry.get<std::string>();
  if (str[0] == '$') {
    std::string key =
        str.substr(1);  // Extract the part of the string after the "$" symbol
    if (definitions.count(key) > 0) {
      return definitions.at(key).get<T>();
    } else {
      throw std::runtime_error("Key not found in definitions: " + key);
    }
  } else {
    return entry.get<T>();
  }
}

template <typename T>
T parse_definition(json const& entry, std::string key, T const& fallback,
                   std::map<std::string, json> const& definitions) {
  std::string str = entry.get<std::string>();
  if (str[0] == '$') {
    std::string key =
        str.substr(1);  // Extract the part of the string after the "$" symbol
    if (definitions.count(key) > 0) {
      return definitions.at(key).get<T>();
    } else {
      return fallback;
    }
  } else {
    return entry.get<T>();
  }
  auto eq = (3 + 4) * 5;
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

  if (!points) {
    points = static_cast<int>((end - start) / width);
  } else if (!width) {
    width = (end - start) / points;
  }
  return {start, end, width};
}

symphas::b_element_type parse_boundary_data(
    json const& entry, std::map<std::string, json> const& definitions) {
  BoundaryType type =
      symphas::boundary_from_str(entry["type"].get<std::string>().c_str());
  std::vector<double> parameters;
  if (entry.contains("parameters")) {
    for (auto const& param : entry["parameters"]) {
      parameters.push_back(param.get<double>());
    }
  }
  int argc = int(parameters.size());
  if (entry.contains("variations")) {
    symphas::boundary_tag_list tags;

    iter_type i = 0;
    for (auto const& variation : entry["variations"]) {
      tags[i++] =
          symphas::boundary_tag_from_str(variation.get<std::string>().c_str());
    }
    return {type, tags, parameters.data(), argc};
  }
  return {type, parameters.data(), argc};
}

JsonConfManager::JsonConfManager(const char* config_file) : SymPhasSettings{} {
  // lets parse settings here
  std::ifstream file(config_file);
  json config;
  file >> config;

  std::map<std::string, json> definitions;
  if (config.contains("definitions")) {
    definitions = config["definitions"].get<std::map<std::string, json>>();
  }

  // Parse domain settings
  auto simulation = config["simulation"];
  auto domain = config["domain"];
  auto intervals = domain["intervals"];
  for (auto const& key : simulation["domain"]) {
    auto interval_j =
        parse_definition<json>(intervals, key.get<std::string>(), definitions);
    auto interval = parse_interval_data(interval_j, definitions);
    auto boundary_j = domain["boundaries"];

    auto boundary_start = (boundary_j["start"].is_string())
                              ? domain["boundaries"][boundary_j["start"]]
                              : boundary_j["start"];
    auto boundary_end = (boundary_j["end"].is_string())
                            ? domain["boundaries"][boundary_j["end"]]
                            : boundary_j["end"];
  }

  make_directories();
}
