
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


#define OPTION_DESCRIPTION_LEN 55
#define OPTION_DESCRIPTION_TAB_BUFFER 26
#define OPTION_DESCRIPTION_TAB 2
#define OPTION_DESCRIPTION_ALIAS_BUFFER 6
#define OPTION_DESCRIPTION_NAME_LEN 19


DLLLIB std::vector<std::pair<std::string, std::string>> params::rawparams;

bool add_base_params(param_map_type& param_map)
{
	using namespace params;

	param_map["title"] = { &title, new param_assign<char*>, 'T' };
	param_map["extend-boundary"] = { &extend_boundary, new param_assign<bool>, 'b' };
	param_map["extend-boundaries"] = { &extend_boundary, new param_assign<bool>, 'b' };
	param_map["visualization"] = { &viz_enabled, new param_assign<bool>,
		"enables visualization with vtk, and requires the interval of visualization"
		" to be set" };
	param_map["visualize"] = { &viz_enabled, new param_assign<bool> };
	param_map["viz-interval-set"] = { &viz_interval, new param_assign<int>,
		"specifies directly the visualization interval" };
	param_map["viz-index-set"] = { &viz_index, new param_assign<int>,
		"specifies directly the visualization index" };
	param_map["viz-interval"] = { &viz_interval_enable, new param_assign_separate<int, bool>, 'v',
		"enables visualization with vtk if not already enabled, and specifies the interval of"
		" visualization to be set" };
	param_map["viz-index"] = { &viz_index_enable, new param_assign_separate<int, bool>, 'V',
		"enables visualization with vtk if not already enabled, and specifies the index of the"
		" field to visualize" };
	param_map["init-inside-val"] = { &init_inside_val, new param_assign<double>, '1' };
	param_map["init-inside"] = { &init_inside_val, new param_assign<double>, '1' };
	param_map["init-outside-val"] = { &init_outside_val, new param_assign<double>, '0' };
	param_map["init-outside"] = { &init_outside_val, new param_assign<double>, '0' };
	param_map["init-rand-val"] = { &init_rand_val, new param_assign<double>, 's' };
	param_map["init-rand"] = { &init_rand_val, new param_assign<double>, 's' };
	param_map["data-file"] = { &input_data_file, new param_assign<char*>, 'd', 
		"specifies an input data file for initial conditions, read by the program when '!' is"
		" provided in the configuration at the initial condition parameter" };
	param_map["data"] = { &input_data_file, new param_assign<char*>, 'd' };


#ifdef EXECUTION_HEADER_AVAILABLE
	param_map["parallelization"] = { &parallelization, new param_assign<symphas::ParallelizationType>, 'P' };
#endif

	return true;
}


DLLLIB char* params::title = nullptr;
DLLLIB bool params::extend_boundary = false;
DLLLIB bool params::viz_enabled = false;
DLLLIB int params::viz_interval = 5;
DLLLIB int params::viz_index = 0;
DLLLIB double params::init_inside_val = 1;
DLLLIB double params::init_outside_val = -1;
DLLLIB double params::init_rand_val = 1;
DLLLIB char* params::input_data_file = nullptr;
DLLLIB int params::start_index = INDEX_INIT;

#ifdef EXECUTION_HEADER_AVAILABLE
DLLLIB symphas::ParallelizationType params::parallelization = symphas::ParallelizationType::PAR;
#endif

DLLLIB void* params::viz_interval_enable[2]{ (void*)&params::viz_interval, (void*)&params::viz_enabled };
DLLLIB void* params::viz_index_enable[2]{ (void*)&params::viz_index, (void*)&params::viz_enabled };


void param_map_element::print_help(FILE* out, const char* name) const
{
	len_type tab_buffer = 26;
	len_type tab = 2;
	len_type len = len_type(description.size());

	if (alias)
	{
		fprintf(out, "  -%c, ", alias);
	}
	else
	{
		fprintf(out, "%6s", "");
	}

	fprintf(out, "--");
	len_type size = len_type(assign_method->print_with_name(out, parameter, name));
	if (size < OPTION_DESCRIPTION_TAB_BUFFER - OPTION_DESCRIPTION_TAB - OPTION_DESCRIPTION_ALIAS_BUFFER)
	{
		if (size < OPTION_DESCRIPTION_NAME_LEN)
		{
			fprintf(out, "%*s", OPTION_DESCRIPTION_NAME_LEN - size, "");
		}
	}
	else
	{
		if (len > 0)
		{
			fprintf(out, "\n%*s", OPTION_DESCRIPTION_TAB_BUFFER + OPTION_DESCRIPTION_TAB, "");
		}
	}
	



	char buffer[OPTION_DESCRIPTION_LEN + 1]{};
	int end = std::min(len, OPTION_DESCRIPTION_LEN);
	std::copy(description.c_str(), description.c_str() + end, buffer);

	char* space = (end == len) ? buffer + end : std::strrchr(buffer, ' ');
	*space = '\0';

	fprintf(out, "%s\n", buffer);
	int start = space - buffer;
	while (end < len)
	{
		char buffer0[OPTION_DESCRIPTION_LEN + 1]{};
		end = std::min(len, start + OPTION_DESCRIPTION_LEN);
		std::copy(description.c_str() + start, description.c_str() + end, buffer0);
		
		char* space = (end == len) ? buffer0 + end - start : std::strrchr(buffer0, ' ');
		*space = '\0';
		start += space - buffer0;

		fprintf(out, "%*s%s\n", OPTION_DESCRIPTION_TAB_BUFFER + OPTION_DESCRIPTION_TAB, "", buffer0);
	}

}


void params::parse_arguments(param_map_type param_map, const char* args, size_t n)
{
	parse_arguments(param_map, &args, n);
}


void params::parse_arguments(param_map_type param_map, const char* const* args, size_t n)
{
	const char err_msg[] = "command line parameter '%s' cannot be interpreted, skipping...\n";

	if (n == 1)
	{
		if (std::strcmp(args[0], "--" ARGUMENT_HELP_STRING) == 0 || std::strcmp(args[0], ARGUMENT_HELP_STRING) == 0)
		{
			print_argument_help(SYMPHAS_INFO, param_map);
			exit(1);
		}
	}

	for (iter_type c = 0; c < n; ++c)
	{
		char key[BUFFER_LENGTH_R2]{}, value[BUFFER_LENGTH_L4]{};


		const char* it;
		const char* start = args[c];
		const char* end = std::strchr(args[c], '\0');
		
		if (std::strlen(args[c]) > 1 && args[c][0] == '-')
		{
			char short_key = args[c][1];
			if (args[c][2] == '=')
			{
				it = args[c] + 3;
				std::copy(it, end, value);
				value[(end - it)] = '\0';
			}

			for (auto& [name, element] : param_map)
			{
				if (element.alias == short_key)
				{
					element.assign(value);
					std::strcpy(key, name.c_str());
					break;
				}
			}
			if (*key == '\0')
			{
				fprintf(SYMPHAS_WARN, err_msg, args[c]);
				continue;
			}
		}
		else if ((it = std::strchr(args[c], '=')) != NULL
			|| (std::strlen(args[c]) > 2 && args[c][0] == '-' && args[c][1] == '-'))
		{
			if (std::strlen(args[c]) > 2 && args[c][0] == '-' && args[c][1] == '-')
			{
				start += 2;
			}

			if (it != NULL)
			{

				std::copy(start, it, key);
				key[(it - start)] = '\0';

				std::copy(it + 1, end, value);
				value[(end - it - 1)] = '\0';
			}
			else
			{
				std::copy(start, end, key);
				key[(end - start)] = '\0';
			}
		}
		else
		{
			fprintf(SYMPHAS_WARN, err_msg, args[c]);
			continue;
		}

		if (std::find_if(rawparams.begin(), rawparams.end(),
			[&](auto e) { return std::strcmp(e.first.c_str(), key) == 0; })
			== rawparams.end())
		{
			rawparams.emplace_back(key, value);
		}

		if (param_map.count(key) > 0)
		{
			param_map[key].assign(value);
		}
		else
		{
			fprintf(SYMPHAS_LOG, "unrecognized command argument '%s'\n", key);
		}
	}
}


void params::print_argument_help(FILE* out, param_map_type param_map)
{
	fprintf(out, SYMPHAS_USAGE_MESSAGE);
	fprintf(out, SYMPHAS_DESCRIPTION_MESSAGE);
	fprintf(out, "\n");

	for (const auto& [param, element] : param_map)
	{
		element.print_help(out, param);
	}
}

