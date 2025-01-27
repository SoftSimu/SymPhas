
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

  //! Generate the configuration from the given options.
  /*!
   * Generate the configuration from the given options. The name of the
   * configuration is also provided.
   *
   * \param options The list of key-value pairs for all the options used
   * in the configuration.
   * \param title The title of the configuration.
   * \param dir The absolute parent directory of the configuration, which is
   * optional.
   */
  JsonConfManager(std::vector<std::pair<std::string, std::string>> const& options,
                  const char* title, const char* dir = "") {}
};
