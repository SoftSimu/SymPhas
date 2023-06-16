
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


#include "writedefines.h"
#include <mutex>


namespace symphas::internal
{
	struct nameidstore;
}

DLLIO std::vector<symphas::internal::nameidstore> nameid_list{};

void symphas::io::copy_data_file_name(const char* dir, const char* data_name, int index, size_t id, DataFileType type, char* out)
{
	if (type == DataFileType::NAMED_DATA)
	{
		snprintf(out, BUFFER_LENGTH, "%s", dir);
	}
	else if (type == DataFileType::CHECKPOINT_DATA)
	{
		if (params::single_output_file)
		{
			snprintf(out, BUFFER_LENGTH, OUTPUT_CHECKPOINT_FMT, dir, id);
		}
		else
		{
			snprintf(out, BUFFER_LENGTH, OUTPUT_CHECKPOINT_INDEX_FMT, dir, id, index);
		}
	}
	else
	{
		char data_dir[BUFFER_LENGTH],
			data_file[BUFFER_LENGTH],
			name_buffer[BUFFER_LENGTH_R1],
			fid[BUFFER_LENGTH_R2];
		
		snprintf(fid, BUFFER_LENGTH_R2, POSTFIX_ID_FMT, id);

		/* generic data means the file name is specified entirely by dataname
		 * the generic data is placed in a subdirectory 
		 */
		if (type == DataFileType::GENERIC_DATA)
		{


			if (params::single_output_file)
			{
				snprintf(data_dir, sizeof(data_dir) / sizeof(char), OUTPUT_DATA_DIR "/%zd/%d", dir, id, index);
			}
			else
			{
				snprintf(data_dir, sizeof(data_dir) / sizeof(char), OUTPUT_DATA_DIR "%zd", dir, id);
			}

			if (*data_name)
			{
				std::strncpy(data_file, data_name, sizeof(data_file) / sizeof(char) - 1);
			}
			else
			{
				std::strcpy(data_file, "out");
			}
		}
		else
		{
			snprintf(data_dir, sizeof(data_dir) / sizeof(char), OUTPUT_DATA_DIR, dir);


			/* in the case of latex plot output, the name is preceded by the title, which will
			 * be printed to the parameter
			 */

			switch (type)
			{
			case DataFileType::SOLUTION_DATA:
			{
				if (*data_name)
				{
					snprintf(name_buffer, sizeof(name_buffer) / sizeof(char), "%s_%s", PHASEFIELD_DATA_NAME, data_name);
				}
				else
				{
					std::strcpy(name_buffer, PHASEFIELD_DATA_NAME);
				}
				break;
			}
			case DataFileType::POSTPROC_DATA:
			{
				std::strncpy(name_buffer, data_name, sizeof(name_buffer) / sizeof(char) - 1);
				break;
			}
			default:
				break;
			}

#ifdef LATEX_PLOT
			char title_name[BUFFER_LENGTH_R2];
			symphas::lib::to_file_name(params::title, title_name, BUFFER_LENGTH_R2);
			if (params::single_output_file)
			{
				snprintf(data_file, sizeof(data_file) / sizeof(char), OUTPUT_DATA_FILE_LATEX_FMT, title_name, name_buffer, fid);
			}
			else
			{
				snprintf(data_file, sizeof(data_file) / sizeof(char), OUTPUT_DATA_FILE_LATEX_INDEX_FMT, title_name, name_buffer, fid, index);
			}
#else

			if (params::single_output_file)
			{
				snprintf(data_file, sizeof(data_file) / sizeof(char), OUTPUT_DATA_FILE_FMT, name_buffer, fid);
			}
			else
			{
				snprintf(data_file, sizeof(data_file) / sizeof(char), OUTPUT_DATA_FILE_INDEX_FMT, name_buffer, fid, index);
			}

#endif

		}

		sprintf(out, "%s/%s", data_dir, data_file);

	}


}

void symphas::io::copy_data_file_name(const char* dir, int index, size_t id, DataFileType type, char* out)
{
	copy_data_file_name(dir, "", index, id, type, out);
}

/*
 * returns the opened file pointer (if it is successful) with the given name at the index and id
 */

FILE* symphas::io::open_data_file(const char* dir, const char* data_name, int index, size_t id, DataFileType type)
{
	static std::vector<std::tuple<std::string, size_t, DataFileType>> idlist;
	static std::mutex idlist_mtx;

	bool first_open;
	if (idlist.empty())
	{
		first_open = true;
	}
	else
	{
		first_open = (std::find(idlist.begin(), idlist.end(), std::make_tuple(dir, id, type)) == idlist.end());
	}
	const char* mode = (params::single_output_file) ? "a" : "w";

	if (first_open)
	{
		idlist_mtx.lock();
		idlist.emplace_back( dir, id, type );
		idlist_mtx.unlock();
	}

	const char* sw_mode = (first_open)
		? "w"
		: mode;

	
	char name[BUFFER_LENGTH];
	copy_data_file_name(dir, data_name, index, id, type, name);

	FILE* f;
	if ((f = fopen(name, sw_mode)) == 0)
	{
		symphas::lib::make_directory_for_file(name);
		if ((f = fopen(name, sw_mode)) == 0)
		{
			fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_OPEN, name);
			exit(ERR_CODE_FILE_OPEN);
		}
	}
	return f;
}

FILE* symphas::io::open_data_file(const char* dir, int index, size_t id, DataFileType type)
{
	return open_data_file(dir, "", index, id, type);
}



FILE* symphas::io::open_postproc_plot_file(const char* dir, const char* name, size_t id, const char* mode)
{
	FILE* plot;
	char fname[BUFFER_LENGTH];

	sprintf(fname, POSTPROC_PLOT_LOC_FMT, dir, name, id);

	if ((plot = fopen(fname, mode)) == 0)
	{
		symphas::lib::make_directory_for_file(fname);
		if ((plot = fopen(fname, mode)) == 0)
		{
			fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_OPEN, fname);
			exit(ERR_CODE_FILE_OPEN);
		}
	}

	return plot;
}


namespace symphas::internal
{

	/*
	 * represents the name/id combination for a given data
	 */
	struct nameidstore
	{
		nameidstore(const char* str, size_t id) : id{ id }, name{ (str && std::strlen(str) > 0) ? new char[std::strlen(str) + 1] : nullptr }, entry_index{ 0 }
		{
			if (name)
			{
				std::strcpy(name, str);
			}
		}
		nameidstore() : nameidstore(nullptr, 0) {}

		nameidstore(const nameidstore& other) : nameidstore(other.name, other.id)
		{
			entry_index = other.entry_index;
		}

		nameidstore(nameidstore&& other) : nameidstore()
		{
			swap(*this, other);
		}

		nameidstore& operator=(nameidstore other)
		{
			swap(*this, other);
			return *this;
		}

		friend void swap(nameidstore& first, nameidstore& second)
		{
			using std::swap;
			swap(first.name, second.name);
			swap(first.id, second.id);
			swap(first.entry_index, second.entry_index);
		}

		~nameidstore()
		{
			delete[] name;
		}

		/* comparison between nameid stores only considers the name itself, rather than
		 * both name and id matching, since only name uniquely differentiates files
		 * if the format later changes to include id, then this can compare on the id as well
		 */
		bool operator==(const nameidstore& other) const
		{
			return (std::strcmp(name, other.name) == 0) && (id == other.id);
		}

		size_t id;
		char* name;// [BUFFER_LENGTH_R4] ;
		int entry_index;
	};

}

/* given a new name, a new file will be created
 * after the plot file is written for the first time, each consecutive call will prepend the sets[n-1] format to the plot file
 */
void symphas::io::write_postproc_plot_file(const char* (*sets), len_type n, const char* dir, const char* name, int index, size_t id, const char* prepend)
{
	using symphas::internal::nameidstore;

	/* contains the list of all the name/id combinations which have been written by the function
	 */


	/* if there is no unique name (since comparison only checks name), then a new file
	 * needs to be created
	 */
	nameidstore nid(name, id);
	auto el = std::find(nameid_list.begin(), nameid_list.end(), nid);
	if (el == nameid_list.end())
	{
		nameid_list.push_back(nid);
		overwrite_postproc_plot_file(sets, n, dir, name, index, id, nameid_list.back().entry_index++, prepend);
	}
	else
	{
		append_postproc_plot_file(sets[n - 1], dir, name, index, id, el->entry_index++, prepend);
	}
}

/*! 
 * The sets in the plot file will always follow the given format:
 * all the sets except the last one will have no formatting to provide
 * the last set will have a single format into which will be printed the location of the datafile to load
 *
 * in the case of the #LATEX_PLOT functionality, the final set format will include, before the datafile name:
 * the title
 * latex output filename, which is formatted by the title, id and index
 */
void symphas::io::overwrite_postproc_plot_file(
	const char* (*sets), len_type n, const char* dir, const char* name, 
	int index, size_t id, int entry_index, const char* prepend)
{
	FILE* plot = open_postproc_plot_file(dir, name, id, "w");
	char fid[BUFFER_LENGTH_R2];
	sprintf(fid, POSTFIX_ID_FMT, id);

	for (iter_type i = 0; i < n - 1; ++i)
	{
		fprintf(plot, "%s", sets[i]);
	}
	fclose(plot);

	append_postproc_plot_file(sets[n - 1], dir, name, index, id, entry_index, prepend);
}

void symphas::io::append_postproc_plot_file(
	const char* const set, const char* dir, const char* name, 
	int index, size_t id, int entry_index, const char* prepend)
{
	FILE* plot = open_postproc_plot_file(dir, name, id, "a");
	fprintf(plot, "\n%s", prepend);


	char fid[BUFFER_LENGTH_R2], data_loc[BUFFER_LENGTH];
	sprintf(fid, POSTFIX_ID_FMT, id);
	copy_data_file_name(DATA_DIR_RELATIVE_PLOT, name, index, id, DataFileType::POSTPROC_DATA, data_loc);
	fprintf(plot, "\n");

#ifdef LATEX_PLOT
	char title_name[BUFFER_LENGTH_R2];
	char latex_name[BUFFER_LENGTH_R2];
	symphas::lib::to_file_name(params::title, title_name);

	sprintf(latex_name, OUTPUT_LATEX_FILE_FMT, title_name, name, id, index);
	fprintf(plot, set, title_name, index, latex_name, data_loc, entry_index);
#else

	fprintf(plot, set, data_loc, entry_index);

#endif

	fclose(plot);
}





