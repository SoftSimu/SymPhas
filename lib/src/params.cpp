
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


#include "params.h"





DLLLIB std::vector<std::pair<std::string, std::string>> params::rawparams;

bool add_base_params(param_map_type& param_map)
{
	using namespace params;

	param_map["title"] = std::make_pair(&title, new param_assign<char*>);
	param_map["extend-boundary"] = std::make_pair(&extend_boundary, new param_assign<bool>);
	param_map["extend-boundaryend"] = std::make_pair(&extend_boundary, new param_assign<bool>);
	param_map["extend-boundary"] = std::make_pair(&extend_boundary, new param_assign<bool>);
	param_map["extend-boundaries"] = std::make_pair(&extend_boundary, new param_assign<bool>);
	param_map["visualization"] = std::make_pair(&viz_interval, new param_assign<int>);
	param_map["visualize"] = std::make_pair(&viz_interval, new param_assign<int>);
	param_map["viz-interval"] = std::make_pair(&viz_interval, new param_assign<int>);
	param_map["init-inside_val"] = std::make_pair(&init_inside_val, new param_assign<double>);
	param_map["init-inside"] = std::make_pair(&init_inside_val, new param_assign<double>);
	param_map["init-outside_val"] = std::make_pair(&init_outside_val, new param_assign<double>);
	param_map["init-outside"] = std::make_pair(&init_outside_val, new param_assign<double>);
	param_map["init-rand_val"] = std::make_pair(&init_rand_val, new param_assign<double>);
	param_map["init-rand"] = std::make_pair(&init_rand_val, new param_assign<double>);
	param_map["data-file"] = std::make_pair(&input_data_file, new param_assign<char*>);
	param_map["data"] = std::make_pair(&input_data_file, new param_assign<char*>);

	return true;
}


DLLLIB char* params::title = nullptr;
DLLLIB bool params::extend_boundary = false;
DLLLIB int params::viz_interval = 0;
DLLLIB double params::init_inside_val = 1;
DLLLIB double params::init_outside_val = -1;
DLLLIB double params::init_rand_val = 1;
DLLLIB char* params::input_data_file = nullptr;
DLLLIB int params::start_index = INDEX_INIT;



void params::parse_arguments(param_map_type param_map, const char* args, size_t n)
{
	parse_arguments(param_map, &args, n);
}


void params::parse_arguments(param_map_type param_map, const char* const* args, size_t n)
{
	const char err_msg[] = "command line parameter '%s' cannot be interpreted, skipping...\n";

	for (iter_type c = 0; c < n; ++c)
	{
		char* tok;
		char argument[BUFFER_LENGTH_L4];
		std::strncpy(argument, args[c], sizeof(argument) / sizeof(char) - 1);
		char key[BUFFER_LENGTH_R2], value[BUFFER_LENGTH_L4];

		if ((tok = std::strtok(argument, "=")) == 0)
		{
			fprintf(SYMPHAS_WARN, err_msg, args[c]);
			continue;
		}
		std::strncpy(key, tok, sizeof(key) / sizeof(char) - 1);

		if ((tok = std::strtok(NULL, "=")) == 0)
		{
			fprintf(SYMPHAS_WARN, err_msg, args[c]);
			continue;
		}
		std::strncpy(value, tok, sizeof(value) / sizeof(char) - 1);

		if (std::find_if(rawparams.begin(), rawparams.end(),
			[&](auto e) { return std::strcmp(e.first.c_str(), key) == 0; })
			== rawparams.end())
		{
			rawparams.emplace_back(key, value);
		}

		if (param_map.count(key) > 0)
		{
			auto [param, type] = param_map[key];
			type->assign(param, value);
		}
		else
		{
			fprintf(SYMPHAS_LOG, "unrecognized command argument '%s'\n", key);
		}
	}
}

