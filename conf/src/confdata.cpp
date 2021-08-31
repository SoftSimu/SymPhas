
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



#include "confdata.h"


char C_PROC[] = "POSTPROC";


DataConf::DataConf(std::vector<std::pair<std::string, std::string>> params)
{
	for (auto&& [k, v] : params)
	{
		if (strcasecmp(C_PROC, k.c_str()) == 0)
		{
			char str[BUFFER_LENGTH_L2]{};
			std::strncpy(str, v.c_str(), sizeof(str) / sizeof(char) - 1);

			symphas::lib::str_trim(str);
			symphas::lib::to_upper(str);

			// check to make sure the NO parameter is not given
			// if NO is given, then nothing should be added
			if (std::strlen(str) > 0 && std::strcmp(str, "NO") != 0)
			{
				auto procs = symphas::conf::parse_options(str);
				for (auto&& proc : procs)
				{
					symphas::internal::parse_proc(proc.c_str(), CONFIG_OPTION_PREFIX_C);
				}
			}
		}
	}
}


void DataConf::write(const char* savedir, const char* name) const
{
	FILE* f;
	char bname[BUFFER_LENGTH];
	sprintf(bname, "%s/%s"  "." CONFIG_EXTENSION, savedir, name);
	if ((f = fopen(bname, "a")) == 0)
	{
		fprintf(SYMPHAS_ERR, "error opening write configuration file '%s'\n", bname);
		exit(ERR_CODE_FILE_OPEN);
	}


	char open = symphas::internal::option_open_bracket();
	char close = symphas::internal::option_close_bracket();

	fprintf(f, CONFIG_NAME_FMT, C_PROC);
	for (auto& [k, v] : symphas::internal::get_data_map())
	{
		if (!v->empty())
		{
			fprintf(f, "%c ", open);
			for (auto& proc : *v)
			{
				fprintf(f, "%s ", k.c_str());
				if (proc.type == ProcessType::SCALAR)
				{
					fprintf(f, "%s ", "ABS");
				}
				else if (proc.type == ProcessType::VEC)
				{
					fprintf(f, "%s ", "VEC");
				}
				else if (proc.type == ProcessType::POINT)
				{
					fprintf(f, "%s ", "PT");
				}

				if (proc.tag == ProcessTag::DFT)
				{
					fprintf(f, "%s ", "DFT");
				}

				char save_spec[BUFFER_LENGTH];
				symphas::io::save_as_str(&proc.save, save_spec, BUFFER_LENGTH);
				fprintf(f, STR(CONFIG_OPTION_PREFIX) " %s ", save_spec);
			}
			fprintf(f, "%c", close);
		}
	}
	fprintf(f, "\n");
	fclose(f);
}






