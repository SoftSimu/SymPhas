
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


#include "processtypes.h"

void symphas::internal::parse_proc_job(const char* v, DataParams* d, ProcessType default_proc, char split_char)
{
	char str[BUFFER_LENGTH]{};
	std::strncpy(str, v, sizeof(str) / sizeof(char) - 1);


	/* set defaults
	 */
	d->type = default_proc;
	d->tag = ProcessTag::NONE;

	d->save = SaveParams{};
	d->save.set_start(params::start_index);

	char delim[]{ split_char, '\0' };

	char* tok;
	if (str[0] == split_char)
	{
		tok = std::strtok(str, delim);
		symphas::io::parse_save_str(tok, &d->save);
	}
	else
	{
		tok = std::strtok(str, delim);
		if (tok != NULL)
		{
			symphas::internal::parse_proc_spec(tok, d, default_proc);

			if ((tok = std::strtok(NULL, delim)) != NULL)
			{
				symphas::io::parse_save_str(tok, &d->save);
			}
		}
	}
}

void symphas::internal::parse_proc_spec(const char* str, DataParams* d, ProcessType default_proc)
{
	char const err_msg[] = "incorrect format for postprocessing configuration, token '%s'\n";
	char type[BUFFER_LENGTH_R4]{}, tag[BUFFER_LENGTH_R4]{};

	size_t n = sscanf(str, "%s %s", type, tag);

	if (n == 0)
	{
		d->type = default_proc;
	}
	else
	{
		if (std::strcmp(type, "YES") == 0)
		{
			d->type = default_proc;
		}
		else if (std::strcmp(type, "NO") == 0)
		{
			d->type = ProcessType::NO;
		}
		else if (std::strcmp(type, "ABS") == 0)
		{
			d->type = ProcessType::SCALAR;
		}
		else if (std::strcmp(type, "VEC") == 0)
		{
			d->type = ProcessType::VECTOR;
		}
		else if (std::strcmp(type, "PT") == 0)
		{
			d->type = ProcessType::POINT;
		}
		else
		{
			fprintf(SYMPHAS_WARN, err_msg, type);
		}

		if (n > 1)
		{
			if (std::strcmp(tag, "DFT") == 0)
			{
				d->tag = ProcessTag::DFT;
			}
			else
			{
				fprintf(SYMPHAS_WARN, err_msg, tag);
			}
		}
	}
}

void symphas::internal::parse_proc(const char* proc_name, char split_char)
{
	char name[BUFFER_LENGTH_R2]{}, args[BUFFER_LENGTH_R2]{};
	int pos = 0;
	size_t n = sscanf(proc_name, "%s %n", name, &pos);
	if (n == 0)
	{
		fprintf(SYMPHAS_WARN, "incorrect format for postprocessing configuration, "
			"read '%s'\n", proc_name);
		return;
	}

	if (pos < std::strlen(proc_name))
	{
		std::strncpy(args, proc_name + pos, sizeof(args) / sizeof(char) - 1);
	}
	else
	{
		std::strcpy(args, "");
	}


	if (symphas::internal::data_map.count(name) != 0)
	{
		data_map[name]->emplace_back(DataParams{});
		auto d = &symphas::internal::data_map[name]->back();
		auto p = symphas::internal::default_procs_map[name];
		parse_proc_job(args, d, p, split_char);

		if (d->type != ProcessType::NO)
		{
			fprintf(SYMPHAS_LOG, "added postprocessing job: '%s'\n", name);
		}
	}
	else
	{
		fprintf(SYMPHAS_WARN, "name not identified as a postprocess '%s'\n", name);
	}
}

void symphas::internal::update_data_stop(iter_type stop)
{
	for (auto& [_, datas] : data_map)
	{
		for (auto& data : *datas)
		{
			data.save.set_stop(stop);
		}
	}
}


std::map<std::string, std::vector<DataParams>*, symphas::internal::any_case_comparator> symphas::internal::get_data_map()
{
	return data_map;
}

std::map<std::string, ProcessType, symphas::internal::any_case_comparator> symphas::internal::get_default_procs_map()
{
	return default_procs_map;
}