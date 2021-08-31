
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

#ifdef USING_CONF
#include "modules-conf.h"
#else
#include "modules-io.h"
#endif

const char* symphas::get_record_name(const char* dir)
{

#ifdef USING_CONF
	size_t out_len = std::strlen(dir) + sizeof(BACKUP_CONFIG_LOC_FMT) / sizeof(char) - 1;
	char* out = new char[out_len];
	snprintf(out, out_len, BACKUP_CONFIG_LOC_FMT, dir);
	return out;

#else
	char* out = new char[BUFFER_LENGTH];
	snprintf(out, BUFFER_LENGTH, "%s/computed_indices.txt", dir);
	return out;

#endif
}

FILE* symphas::open_record(const char* name)
{
	FILE* f;
	if ((f = fopen(name, "a")) == 0)
	{
		fprintf(SYMPHAS_ERR, "error opening backup configuration file '%s'\n", name);
		exit(1001);
	}
	return f;
}

void symphas::record_index(const char* dir, SaveParams const& save, int index)
{
	const char* name = get_record_name(dir);
	FILE* f = open_record(name);
#ifdef USING_CONF
	fprintf(f, CONFIG_INDEX_PREFIX "%d\n", index);
	if (index == save.get_stop())
	{
		fprintf(f, CONFIG_INDEX_PREFIX SIMULATION_DONE_KEY "\n");
	}
#else
	fprintf(f, "%d\n", index);
#endif
	fclose(f);

	delete[] name;
}

