
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

namespace symphas::internal
{
	inline Conf* global_config;
}

const Conf& symphas::conf::config()
{
	return *symphas::internal::global_config;
}



Conf symphas::conf::make_config(const char* file, param_map_type const& param_map)
{

	char line[BUFFER_LENGTH_L3]{},
		title[BUFFER_LENGTH]{};
	std::vector<std::pair<std::string, std::string>> options;


	std::vector<iter_type> indices;
	bool sim_done = false;

	FILE* f;
	if ((f = fopen(file, "r")) == 0)
	{
		fprintf(SYMPHAS_ERR, "error opening configuration file '%s'\n", file);
		exit(444);
	}

	char* parent_dir = new char[std::strlen(file) + 1];
	std::strcpy(parent_dir, file);
	symphas::lib::get_parent_directory(parent_dir, parent_dir);


	while (fgets(line, BUFFER_LENGTH_L3, f) != NULL)
	{
		symphas::lib::str_trim(line);
		int pos = 0;
		if ((pos = symphas::lib::pos_after_token(line, CONFIG_TITLE_PREFIX)) > 0)
		{
			if (*title)
			{
				fprintf(SYMPHAS_ERR, "redefinition of title!\n");
				exit(445);
			}
			else if (sscanf(line, CONFIG_TITLE_PREFIX "%[^\n]", title) < 1)
			{
				fprintf(SYMPHAS_ERR, "given title is in incorrect format\n");
				exit(446);
			}
		}
		else if (
			!param_map.empty() &&
			(pos = symphas::lib::pos_after_token(line, CONFIG_PARAM_PREFIX)) > 0)
		{
			int 
				start = 0, 
				end = 0;
			sscanf(line, CONFIG_PARAM_PREFIX "%n%*[^\n]%n", &start, &end);


			if (end > start)
			{
				char param[BUFFER_LENGTH];
				std::copy(line + start, line + end, param);
				param[end - start] = '\0';

				params::parse_arguments(param_map, param, 1);

			}
		}
		else if ((pos = symphas::lib::pos_after_token(line, CONFIG_INDEX_PREFIX)) > 0)
		{
			char* endptr;
			int index = strtol(line + pos, &endptr, 10);
			
			if (!*endptr)
			{
				indices.emplace_back(index);
			}
			else if (std::strcmp(line + pos, SIMULATION_DONE_KEY) == 0)
			{
				sim_done = true;
			}
			else
			{
				fprintf(SYMPHAS_WARN, "the written completion index '%s' is not "
					"in the correct format\n", line);
			}

		}
		else if ((pos = symphas::lib::pos_after_token(line, CONFIG_COMMENT_PREFIX)) > 0)
		{
			// comment line
		}
		else
		{
			char* sep_pos = std::strchr(line, CONFIG_SEP_KEY_C);
			size_t len = std::strlen(line);
			if (sep_pos)
			{
				size_t name_len = static_cast<size_t>(sep_pos - line);

				/* copy the key (before the separation character) into k
				 */
				char* k = new char[name_len + 1];
				std::copy(line, sep_pos, k);
				k[name_len] = '\0';
				symphas::lib::str_trim(k);


				/* copies the value which is everything after the separation char
				 */
				char* v = new char[len - name_len];
				std::copy(sep_pos + 1, line + len, v);
				v[len - name_len - 1] = '\0';
				symphas::lib::str_trim(v);


				/* get only the first word in k, find first character that isn't in
				 * the capital alphabetical range, since those are the only valid
				 * characters for a key
				 */
				char* it = k;
				for (; *it >= 'A' && *it <= 'Z'; ++it);
				*it = '\0';

				options.emplace_back(k, v);

				delete[] v;
				delete[] k;
			}
			else
			{
				if (*line)
				{
					fprintf(SYMPHAS_WARN, "unidentified configuration entry found: '%s'\n", line);
				}
			}
		}
	}
	fclose(f);

	Conf c{ options, title, parent_dir };
	c.append_computed_index(indices);
	if (sim_done)
	{
		c.set_done();
	}

	char savedir[BUFFER_LENGTH];
	snprintf(savedir, BUFFER_LENGTH, "%s/" CHECKPOINT_DIR, c.get_result_dir());
	c.write(savedir);

	fprintf(SYMPHAS_LOG, "configuration loaded from '%s'\n", file);
	fprintf(SYMPHAS_LOG, OUTPUT_BANNER);
	return c;
}


Conf symphas::conf::make_config(const char* file)
{
	return symphas::conf::make_config(file, {});
}






Conf::Conf(std::vector<std::pair<std::string, std::string>> options, const char* title, const char* dir) :
	SystemConf(options, title, dir)
#ifdef USING_PROC
	, DataConf(options)
#endif
	, sim_done{ false }
	, computed_indices{ {} }
{}

Conf::Conf(const char* file) : Conf{ symphas::conf::make_config(file) } {}


Conf::Conf(symphas::problem_parameters_type const& parameters, const char* title, const char* dir)
	: SystemConf(parameters, title, dir)
#ifdef USING_PROC
	, DataConf()
#endif
	, sim_done{ false }
	, computed_indices{ {} }
{}


void Conf::append_computed_index(std::vector<iter_type> indices)
{
	computed_indices.insert(computed_indices.end(), indices.begin(), indices.end());
	std::sort(computed_indices.begin(), computed_indices.end());
}

void Conf::append_computed_index(iter_type index)
{
	computed_indices.emplace_back(index);
	std::sort(computed_indices.begin(), computed_indices.end());
}

auto Conf::get_computed_indices() const
{
	return computed_indices;
}

//! Return the last index in the list.
/*!
	* Since the list is sorted, this is always the
	* largest index.
	*/
iter_type Conf::get_last_computed_index() const
{
	if (computed_indices.empty())
	{
		return 0;
	}
	else
	{
		return computed_indices.back();
	}
}


//! Set the simulation for this configuration as complete.
void Conf::set_done()
{
	sim_done = true;
}

//! Check if the configuration is complete.
bool Conf::is_done() const
{
	return sim_done;
}


void Conf::write(const char* savedir) const
{

	SystemConf::write(savedir);
#ifdef USING_PROC
	DataConf::write(savedir);
#endif
}



void symphas::conf::setup_global_config(Conf const& configuration)
{
	symphas::internal::global_config = new Conf{ configuration };


#ifdef USING_PROC
	symphas::internal::update_data_stop(symphas::internal::global_config->save.get_stop());
#endif
}

void symphas::conf::setup_global_config(const char* file, param_map_type param_map)
{
	setup_global_config(make_config(file, param_map));
}




std::vector<iter_type> symphas::conf::get_save_indices(Conf const& config)
{
#ifdef USING_PROC
	std::vector<SaveParams> data_saves;
	for (auto [_, datas] : symphas::internal::get_data_map())
	{
		for (auto& data : *datas)
		{
			data_saves.push_back(data.save);
		}
	}
	return symphas::io::get_save_indices(config.save, data_saves);
#else
	return symphas::io::get_save_indices(config.save);
#endif
}

iter_type symphas::conf::next_save(iter_type index, Conf const& config)
{
	auto saves = get_save_indices(config);
	auto lower = std::lower_bound(saves.begin(), saves.end(), index + 1);
	if (lower != saves.end())
	{
		return *lower;
	}
	else
	{
		return index;
	}
}


bool symphas::conf::is_last_save(iter_type index, Conf const& config)
{
	std::vector<iter_type> saves = symphas::conf::get_save_indices(config);
	return symphas::io::is_last_save(index, saves);
}



iter_type symphas::conf::num_saves(Conf const& config)
{
	return static_cast<iter_type>(symphas::io::get_save_indices(config.save).size());
}


Conf symphas::conf::restore_checkpoint(param_map_type param_map)
{
	char buffer[BUFFER_LENGTH_L4];
	std::strncpy(buffer, params::checkpoint, sizeof(buffer) / sizeof(char) - 1);

	char* tok = std::strtok(buffer, ",");
	char dir[BUFFER_LENGTH_L4];
	std::strncpy(dir, tok, sizeof(dir) / sizeof(char) - 1);

	if ((tok = std::strtok(NULL, ",")) != 0)
	{
		char* endptr;
		int index = strtol(tok, &endptr, 10);
		if (*endptr)
		{
			fprintf(SYMPHAS_WARN, "error reading the given checkpoint index '%s'\n", tok);
			return restore_checkpoint(param_map, dir);
		}
		else
		{
			return restore_checkpoint(param_map, dir, index);
		}
	}

	return restore_checkpoint(param_map, dir);
}


Conf symphas::conf::restore_checkpoint(param_map_type param_map, const char* dir, int index)
{
	char name[BUFFER_LENGTH_L4];
	sprintf(name, BACKUP_CONFIG_LOC_FMT, dir);


	Conf configuration = symphas::conf::make_config(name, param_map);
	if (configuration.is_done())
	{
		fprintf(SYMPHAS_LOG, "the provided backup configuration refers to a solution "
			"which has been completed\n");
	}

	if (index < 0)
	{
		params::start_index = configuration.get_last_computed_index();
	}
	else
	{
		auto indices = configuration.get_computed_indices();
		auto lower = std::lower_bound(indices.begin(), indices.end(), index);
		if (*lower > index)
		{
			fprintf(SYMPHAS_LOG, "the provided index does not match any candidates "
				"in the list of completed indices, so loading from '%d' instead\n", *lower);
			params::start_index = *lower;
		}
		else
		{
			params::start_index = index;
		}
	}
	fprintf(SYMPHAS_LOG, "loading from index '%d'\n", params::start_index);

	char* load_name = new char[std::strlen(dir) + 1];
	std::strcpy(load_name, dir);

	symphas::init_data_type load;
	load.in = Inside::CHECKPOINT;
	load.file.index = params::start_index;
	load.file.name = load_name;

	/* change the initial data so it is loaded from the saved name
	 */
	for (iter_type i = 0; i < configuration.system_count(); ++i)
	{
		configuration.set_initial_condition(load, i);
	}

	return configuration;
}

