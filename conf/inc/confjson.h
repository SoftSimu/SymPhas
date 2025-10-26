
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

#include <fstream>
#include <nlohmann/json.hpp>
#include <nlohmann/json-schema.hpp>


#include "confsystem.h"

using json = nlohmann::json;

namespace symphas {
//! Specifies elements used in the management of SymPhas configuration.
/*!
 * Specifies elements used in the management of SymPhas configuration.
 */
namespace conf {
struct SymPhasJsonKeys {
  static constexpr char model[] = "model";
  static constexpr char initial_conditions[] = "initial_conditions";
  static constexpr char domain[] = "domain";
  static constexpr char boundaries[] = "boundaries";
  static constexpr char intervals[] = "intervals";
  static constexpr char parameters[] = "parameters";
};

}  // namespace conf
}  // namespace symphas

struct JsonConfManager : SymPhasSettings {
  JsonConfManager() : SymPhasSettings{} {}

  JsonConfManager(const char* config_file);

  //! Write a JSON backup of the current configuration.
  /*!
   * Write a JSON backup file of the current configuration to the specified
   * directory with the given name. This reconstructs the JSON from the current
   * settings and saves it as a .json file.
   *
   * \param savedir The directory where the backup should be saved.
   * \param name The base name for the backup file (without extension).
   */
  void write(const char* savedir, const char* name = BACKUP_CONFIG_NAME) const;

  std::string config_dir;
 private:
  void parse_simulation_settings(
      const json& config_json, std::map<std::string, json> const& definitions);
  void parse_model_settings(const json& config_json,
                            std::map<std::string, json> const& definitions);
  void parse_domain_settings(const json& config_json,
                             std::map<std::string, json> const& definitions);
  void parse_name_settings(const json& config_json,
                           std::map<std::string, json> const& definitions);
  void parse_directory_settings(const json& config_json,
                                std::map<std::string, json> const& definitions);

  void parse_single_initial_condition(
      const json& init_config, std::map<std::string, json> const& definitions,
      int field_idx);

  std::string build_save_spec_from_json(
      const json& save, std::map<std::string, json> const& definitions);

  template <typename T>
  T parse_value_or_definition(const json& entry,
                              std::map<std::string, json> const& definitions);
  template <typename T>
  T parse_value_or_definition(const json& entry, T const& fallback,
                              std::map<std::string, json> const& definitions);
};
