
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

#include "prereq.h"

#ifdef PRINT_TIMINGS
DLLMOD double symphas::iteration_time = 0;
DLLMOD double symphas::init_time = 0;
DLLMOD double symphas::model_update_time = 0;
DLLMOD double symphas::model_equation_time = 0;
DLLMOD double symphas::model_step_time = 0;

DLLMOD int symphas::iteration_count = 0;

void symphas::print_timings(FILE* out)
{
	if (iteration_count > 0)
	{
		fprintf(out, OUTPUT_BANNER);
		fprintf(out, "Simulation timings:\n\n");

		fprintf(out, "\tConfiguration and parameter initialization time:\n");
		fprintf(out, "\t\t... %lf\n\n", init_time);

		fprintf(out, "\tCombined / average time for model updating:\n");
		fprintf(out, "\t\t... %lf / %lf\n\n", model_update_time, model_update_time / iteration_count);

		fprintf(out, "\tCombined / average time for model equation processing:\n");
		fprintf(out, "\t\t... %lf / %lf\n\n", model_equation_time, model_equation_time / iteration_count);

		fprintf(out, "\tCombined / average time for model stepping:\n");
		fprintf(out, "\t\t... %lf / %lf\n\n", model_step_time, model_step_time / iteration_count);

		fprintf(out, "\tCombined / average time for whole solution iteration:\n");
		fprintf(out, "\t\t... %lf / %lf\n\n", iteration_time, iteration_time / iteration_count);

		fprintf(out, OUTPUT_BANNER);
	}
}

#endif




param_map_type symphas::build_params()
{
	param_map_type param_map;

	add_base_params(param_map);
#ifdef USING_IO
	add_io_params(param_map);
#endif

	add_solution_params(param_map);

	return param_map;
}



#ifdef USING_CONF

#include "conf.h"

//! Initialize the program parameters.
/*!
 * Given the name of the configuration, the program level (global)
 * configuration will be loaded and initialized. The command line arguments
 * will also be parsed and the variables corresponding to the command
 * line arguments will be initialized to the given values. The command
 * line parameter list, which gives the mapping between command line keys
 * and program level variables, is initialized using
 * symphas::build_params().
 *
 * \param config The name of the configuration file.
 * \param param_list The list of strings containing the key value pairs
 * in the format: "key=value" from which command line parameters
 * are extracted.
 * \param num_params The number of command line arguments in the list.
 */
void symphas::init(const char* config, const char* const* param_list, int num_params)
{

#ifdef USING_MPI
	MPI_Init(NULL, NULL);
#endif

#ifdef PRINT_TIMINGS
	symphas::Time t;
#endif

	param_map_type param_map = build_params();
	params::parse_params(param_map, param_list, num_params);

	if (config != NULL)
	{
		if (strchr(config, '=') != NULL || (std::strlen(config) > 2 && config[0] == '-'))
		{
			params::parse_params(param_map, &config, 1);
		}
	}

#ifdef PRINT_TIMINGS
	init_time += t.current_duration();
#endif

	if (symphas::io::is_checkpoint_set())
	{
		Conf backup_conf = symphas::conf::restore_checkpoint(param_map);
		symphas::conf::setup_global_config(backup_conf);
	}
	else
	{
		symphas::conf::setup_global_config(config, param_map);
	}

#ifdef DEBUG
	params::print_params(param_map);
#endif
	params::assign(param_map["title"], symphas::conf::config().get_title());

}

#else

//! Initialize the program parameters.
/*!
 * Initialize the name of the solution. The command line parameters
 * will also be parsed and the variables corresponding to the command
 * line arguments will be initialized to the given values. The command
 * line parameter list, which gives the mapping between command line keys
 * and program level variables, is initialized using
 * symphas::build_params().
 *
 * \param title The name of the solution.
 * \param param_list The list of strings containing the key value pairs
 * in the format: "key=value" from which command line parameters
 * are extracted.
 * \param num_params The number of command line arguments in the list.
 */
void symphas::init(const char* title, const char* const* param_list, size_t num_params)
{

#ifdef USING_MPI
	MPI_Init(NULL, NULL);
#endif

#ifdef PRINT_TIMINGS
	symphas::Time t;
#endif

	param_map_type param_map = build_params();
	params::parse_params(param_map, param_list, num_params);

	if (strchr(config, '=') != NULL || (std::strlen(config) > 2 && config[0] == '-'))
	{
		params::parse_params(param_map, &config, 1);
	}

#ifdef PRINT_TIMINGS
	init_time += t.current_duration();
#endif

	params::assign(param_map["title"], title);
}

#endif



void symphas::finalize()
{
#ifdef USING_MPI
	MPI_Finalize();
#endif
}