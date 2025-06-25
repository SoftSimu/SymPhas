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
 */

#include "conf.h"

#include "confjson.h"

namespace symphas::internal {
inline Conf *global_config;
}

const Conf &symphas::conf::config() {
  return *symphas::internal::global_config;
}

/*!
 * Read the configuration.
 */
std::tuple<std::vector<std::pair<std::string, std::string>>, std::string>
parse_conf_option_list(const char *file, param_map_type const &param_map) {
  FILE *f;
  if ((f = fopen(file, "r")) == nullptr) {
    fprintf(SYMPHAS_ERR, "error opening configuration file '%s'\n", file);
    exit(444);
  }

  char line[BUFFER_LENGTH_L3]{}, title[BUFFER_LENGTH]{};

  std::vector<std::pair<std::string, std::string>> options;

  while (fgets(line, BUFFER_LENGTH_L3, f) != nullptr) {
    symphas::lib::str_trim(line);
    int pos = 0;
    if ((pos = symphas::lib::pos_after_token(line, CONFIG_TITLE_PREFIX)) > 0) {
      if (*title) {
        fprintf(SYMPHAS_ERR, "redefinition of title!\n");
        exit(445);
      } else if (sscanf(line, CONFIG_TITLE_PREFIX "%[^\n]", title) < 1) {
        fprintf(SYMPHAS_ERR, "given title is in incorrect format\n");
        exit(446);
      }
    } else if (!param_map.empty() && (pos = symphas::lib::pos_after_token(
                                          line, CONFIG_PARAM_PREFIX)) > 0) {
      int start = 0, end = 0;
      sscanf(line, CONFIG_PARAM_PREFIX "%n%*[^\n]%n", &start, &end);

      if (end > start) {
        char param[BUFFER_LENGTH];
        std::copy(line + start, line + end, param);
        param[end - start] = '\0';

        params::parse_params(param_map, param, 1);
      }
    } else if ((pos = symphas::lib::pos_after_token(
                    line, CONFIG_COMMENT_PREFIX)) > 0) {
      // comment line
    } else {
      char *sep_pos = std::strchr(line, CONFIG_SEP_KEY_C);
      size_t len = std::strlen(line);
      if (sep_pos) {
        auto name_len = static_cast<size_t>(sep_pos - line);

        /* copy the key (before the separation character) into k
         */
        char *k = new char[name_len + 1];
        std::copy(line, sep_pos, k);
        k[name_len] = '\0';
        symphas::lib::str_trim(k);

        /* copies the value which is everything after the separation char
         */
        char *v = new char[len - name_len];
        std::copy(sep_pos + 1, line + len, v);
        v[len - name_len - 1] = '\0';
        symphas::lib::str_trim(v);

        /* get only the first word in k, find first character that isn't in
         * the capital alphabetical range, since those are the only valid
         * characters for a key
         */
        char *it = k;
        for (; *it >= 'A' && *it <= 'Z'; ++it);
        *it = '\0';

        options.emplace_back(k, v);

        delete[] v;
        delete[] k;
      } else {
        if (*line) {
          fprintf(SYMPHAS_WARN,
                  "unidentified configuration entry found: '%s'\n", line);
        }
      }
    }
  }
  fclose(f);

  return {options, title};
}

std::pair<std::vector<iter_type>, bool> get_completed_indices(
    const char *parent_dir) {
  std::vector<iter_type> indices;
  bool sim_done = false;

  const char file_name[] = "progress.txt ";
  size_t out_len = std::strlen(parent_dir) + std::strlen(file_name) + 2;
  char *indices_file_name = new char[out_len];
  snprintf(indices_file_name, out_len, "%s/%s", parent_dir, file_name);

  FILE *indices_file;
  if ((indices_file = fopen(indices_file_name, "r")) != nullptr) {
    char indices_line[BUFFER_LENGTH_L3]{};

    while (fgets(indices_line, BUFFER_LENGTH_L3, indices_file) != nullptr) {
      symphas::lib::str_trim(indices_line);
      int pos = 0;

      if ((pos = symphas::lib::pos_after_token(indices_line,
                                               CONFIG_INDEX_PREFIX)) > 0) {
        char *endptr;
        int index = strtol(indices_line + pos, &endptr, 10);

        if (!*endptr) {
          indices.emplace_back(index);
        } else if (std::strcmp(indices_line + pos, SIMULATION_DONE_KEY) == 0) {
          sim_done = true;
        } else {
          fprintf(SYMPHAS_WARN,
                  "the written completion index '%s' is not "
                  "in the correct format\n",
                  indices_line + pos);
        }
      }
    }
    fclose(indices_file);
  }

  delete[] indices_file_name;

  return {indices, sim_done};
}

Conf symphas::conf::make_config(const char *file,
                                param_map_type const &param_map) {
  if (file == nullptr) {
    throw std::invalid_argument("config file name was not specified");
  } else {
    const char *ext_pos = std::strrchr(file, '.');
    if (ext_pos == nullptr) {
      fprintf(SYMPHAS_ERR, "no file extension found in '%s'\n", file);
      exit(443);
    }

    char *ext = new char[std::strlen(ext_pos) + 1]{};
    std::strcpy(ext, ext_pos);
    for (int i = 0; ext[i]; ++i) ext[i] = std::tolower(ext[i]);

    bool is_json = (std::strcmp(ext, ".json") == 0);
    delete[] ext;

    symphas::lib::string parent_dir(symphas::lib::get_parent_directory(file));
    auto [indices, sim_done] = get_completed_indices(parent_dir);

    if (is_json) {
      JsonConfManager json_conf(file);
      char savedir[BUFFER_LENGTH];
      snprintf(savedir, BUFFER_LENGTH, "%s/" CHECKPOINT_DIR,
               json_conf.directory_settings.result_dir);

      json_conf.write(savedir);
      Conf c{json_conf};

      c.append_computed_index(indices);
      if (sim_done) {
        c.set_done();
      }

      fprintf(SYMPHAS_LOG, "configuration loaded from '%s'\n", file);
      fprintf(SYMPHAS_LOG, OUTPUT_BANNER);

      return c;
    } else {
      auto [options, title] = parse_conf_option_list(file, param_map);

      Conf c{options, title.c_str(), parent_dir};
      c.conf_file = file;

      c.append_computed_index(indices);
      if (sim_done) {
        c.set_done();
      }

      char savedir[BUFFER_LENGTH];
      snprintf(savedir, BUFFER_LENGTH, "%s/" CHECKPOINT_DIR,
               c.get_result_dir());
      c.write(savedir);

      fprintf(SYMPHAS_LOG, "configuration loaded from '%s'\n", file);
      fprintf(SYMPHAS_LOG, OUTPUT_BANNER);
      return c;
    }
  }
}

Conf symphas::conf::make_config(const char *file) {
  return symphas::conf::make_config(file, {});
}

Conf::Conf(std::vector<std::pair<std::string, std::string>> const &options,
           const char *title, const char *dir)
    : SystemConf(options, title, dir)
#ifdef WITH_PROC_LIB
      ,
      DataConf(options)
#endif
      ,
      computed_indices{{}},
      sim_done{false} {
}

Conf::Conf(const char *file) : Conf{symphas::conf::make_config(file)} {}

void Conf::append_computed_index(std::vector<iter_type> indices) {
  computed_indices.insert(computed_indices.end(), indices.begin(),
                          indices.end());
  std::sort(computed_indices.begin(), computed_indices.end());
}

void Conf::append_computed_index(iter_type index) {
  computed_indices.emplace_back(index);
  std::sort(computed_indices.begin(), computed_indices.end());
}

auto Conf::get_computed_indices() const { return computed_indices; }

//! Return the last index in the list.
/*!
 * Since the list is sorted, this is always the
 * largest index.
 */
iter_type Conf::get_last_computed_index() const {
  if (computed_indices.empty()) {
    return 0;
  } else {
    return computed_indices.back();
  }
}

//! Set the simulation for this configuration as complete.
void Conf::set_done() { sim_done = true; }

//! Check if the configuration is complete.
bool Conf::is_done() const { return sim_done; }

void Conf::write(const char *savedir) const {
  SystemConf::write(savedir);
#ifdef WITH_PROC_LIB
  DataConf::write(savedir);
#endif
}

void symphas::conf::setup_global_config(Conf const &configuration) {
  symphas::internal::global_config = new Conf{configuration};

#ifdef WITH_PROC_LIB
  symphas::internal::update_data_stop(
      symphas::internal::global_config->simulation_settings.save.get_stop());
#endif
}

void symphas::conf::setup_global_config(const char *file,
                                        param_map_type const &param_map) {
  setup_global_config(make_config(file, param_map));
}

std::vector<iter_type> symphas::conf::get_save_indices(Conf const &config) {
#ifdef WITH_PROC_LIB
  std::vector<SaveParams> data_saves;
  for (auto [_, datas] : symphas::internal::get_data_map()) {
    for (auto &data : *datas) {
      data_saves.push_back(data.save);
    }
  }
  return symphas::io::get_save_indices(config.simulation_settings.save,
                                       data_saves);
#else
  return symphas::io::get_save_indices(config.simulation_settings.save);
#endif
}

iter_type symphas::conf::next_save(iter_type index, Conf const &config) {
  auto saves = get_save_indices(config);
  auto lower = std::lower_bound(saves.begin(), saves.end(), index + 1);
  if (lower != saves.end()) {
    return *lower;
  } else {
    return index;
  }
}

bool symphas::conf::is_last_save(iter_type index, Conf const &config) {
  std::vector<iter_type> saves = symphas::conf::get_save_indices(config);
  return symphas::io::is_last_save(index, saves);
}

iter_type symphas::conf::num_saves(Conf const &config) {
  return static_cast<iter_type>(
      symphas::io::get_save_indices(config.simulation_settings.save).size());
}

Conf symphas::conf::restore_checkpoint(param_map_type param_map) {
  char buffer[BUFFER_LENGTH_L4];
  std::strncpy(buffer, params::checkpoint, sizeof(buffer) / sizeof(char) - 1);

  char *tok = std::strtok(buffer, ",");
  char dir[BUFFER_LENGTH_L4];
  std::strncpy(dir, tok, sizeof(dir) / sizeof(char) - 1);

  if ((tok = std::strtok(NULL, ",")) != 0) {
    char *endptr;
    int index = strtol(tok, &endptr, 10);
    if (*endptr) {
      fprintf(SYMPHAS_WARN, "error reading the given checkpoint index '%s'\n",
              tok);
      return restore_checkpoint(param_map, dir);
    } else {
      return restore_checkpoint(param_map, dir, index);
    }
  }

  return restore_checkpoint(param_map, dir);
}

Conf symphas::conf::restore_checkpoint(param_map_type const &param_map,
                                       const char *dir, int index) {
  Conf configuration = symphas::conf::make_config_auto(dir, param_map);
  if (configuration.is_done()) {
    fprintf(SYMPHAS_LOG,
            "the provided backup configuration refers to a solution "
            "which has been completed\n");
  }

  if (index < 0) {
    params::start_index = configuration.get_last_computed_index();
  } else {
    auto indices = configuration.get_computed_indices();
    auto lower = std::lower_bound(indices.begin(), indices.end(), index);
    if (lower != indices.end() && *lower > index) {
      fprintf(SYMPHAS_LOG,
              "the provided index does not match any candidates "
              "in the list of completed indices, so using index '%d' instead\n",
              *lower);
      params::start_index = *lower;
    } else {
      params::start_index = index;
    }
  }
  fprintf(SYMPHAS_LOG, "loading from index '%d'\n", params::start_index);

  char *load_name = new char[std::strlen(dir) + 1];
  std::strcpy(load_name, dir);

  symphas::init_data_type load;
  load[Axis::NONE].in = Inside::CHECKPOINT;
  load[Axis::NONE].file = {load_name, params::start_index};

  /* change the initial data so it is loaded from the saved name
   */
  for (iter_type i = 0; i < configuration.system_count(); ++i) {
    configuration.domain_settings.set_initial_condition(load, i);
  }

  return configuration;
}

Conf symphas::conf::restore_checkpoint(const char *dir, int index) {
  return restore_checkpoint({}, dir, index);
}

/*!
 * Find and return the path to either configuration.json or configuration.in
 * file. Returns the first file found, preferring .json over .in.
 *
 * \param base_dir The directory to search in (optional, defaults to current
 * directory).
 * \return A string containing the path to the found configuration file, or
 * empty string if none found.
 */
symphas::lib::string symphas::conf::find_configuration_file(
    const char *base_dir) {
  char json_path[BUFFER_LENGTH];
  char in_path[BUFFER_LENGTH];

  // Construct paths for both configuration files
  snprintf(json_path, BUFFER_LENGTH, "%s/%s/%s.json", base_dir,
           CHECKPOINT_DIR, BACKUP_CONFIG_NAME);
  snprintf(in_path, BUFFER_LENGTH, "%s/%s/%s.in", base_dir,
           CHECKPOINT_DIR, BACKUP_CONFIG_NAME);

  // Check for configuration.json first (preferred)
  FILE *json_file = fopen(json_path, "r");
  if (json_file != nullptr) {
    fclose(json_file);
    return symphas::lib::string(json_path);
  }

  // Check for configuration.in as fallback
  FILE *in_file = fopen(in_path, "r");
  if (in_file != nullptr) {
    fclose(in_file);
    return symphas::lib::string(in_path);
  }

  // Neither file found
  return symphas::lib::string();
}

/*!
 * Create configuration using auto-detected configuration file.
 * Searches for either configuration.json or configuration.in in the specified
 * directory.
 *
 * \param base_dir The directory to search for configuration files (optional,
 * defaults to current directory).
 * \param param_map The parameter keys which may be set by this configuration.
 * \return A Conf object created from the found configuration file.
 */
Conf symphas::conf::make_config_auto(const char *base_dir,
                                     param_map_type const &param_map) {
  symphas::lib::string config_file = find_configuration_file(base_dir);

  if (config_file.empty()) {
    fprintf(SYMPHAS_ERR,
            "No configuration file found. Looking for '%s.json' or '%s.in' in "
            "'%s'\n",
            BACKUP_CONFIG_NAME, BACKUP_CONFIG_NAME, base_dir);
    exit(447);
  }

  fprintf(SYMPHAS_LOG, "Found configuration file: '%s'\n", config_file.data);
  return make_config(config_file, param_map);
}
