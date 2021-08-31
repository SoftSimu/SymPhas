
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
 * MODULE:  io
 * PURPOSE: Defines basic parts of the input (reading) functionality.
 *
 * ***************************************************************************
 */

#pragma once


#include "io.h"





/* definitions for output error messages from reading data
 */
#define SYMPHAS_MSG_BAD_DIM_READ \
"the datafile being loaded has data which is an inconsistent " \
"dimension with the one specified in the configuration: found dimension " \
"%d, but expected dimension %d\n"
#define SYMPHAS_MSG_READ_NO_DATA \
"reading from the datafile failed after EOF was reached\n"
#define SYMPHAS_MSG_BAD_INDEX_READ \
"the specified index '%d' was not found while parsing the datafile '%s'\n"
#define SYMPHAS_MSG_READ_DIFFERENT_INDEX \
SYMPHAS_MSG_BAD_INDEX_READ "the closest index found was '%d'\n"

#define BAD_INDEX -1





namespace symphas::io
{
	//! Data used to read information from a datafile.
	/*!
	 * Groups together information that is used to read data from a file, typically
	 * grid data, which is indexed within a file with the given name, where the
	 * name also has an id associated with it depending on the number of grids.
	 */
	struct read_info
	{
		iter_type index;		//!< Index to retrieve from the file.
		size_t id;				//!< Id of the system that is being loaded.
		char* name;				//!< Name of the file to retrieve.
	};

	
	//! A data source will be accessed and the given data array initialized.
	/*!
	 * The values from the data source given by the read information parameter 
	 * and accessed and copied into the array. Typically this is a data file.
	 * 
	 * A data file will be opened to be parsed by the correct read utility,
	 * and the values of the given array will be initialized based on the data.
	 * The length of the array has to be correctly sized beforehand.
	 * 
	 * Reads data from the grid to a file. It will ignore the header and other
	 * file information until it reaches the index, and then it will read the 
	 * data corresponding to that index and return the index which was read
	 * from.
	 *
	 * if the index is not found, because either data could not be read or the 
	 * index found does not end up matching the requested index, that index 
	 * will instead be returned.
	 * 
	 * \param values The arrays into which the values are read into.
	 * \param rinfo Information about how to access the persistent information.
	 */
	template<typename T>
	int read_grid(T* values, symphas::io::read_info rinfo);
	

	//! A data source will be accessed and the given data array initialized.
	/*!
	 * The values from the data source given by the read information parameter
	 * and accessed and copied into the array. Typically this is a data file.
	 *
	 * A data file will be opened to be parsed by the correct read utility,
	 * and the values of the given array will be initialized based on the data.
	 * The length of the array has to be correctly sized beforehand.
	 *
	 * Reads data from the grid to a file. It will ignore the header and other
	 * file information until it reaches the index, and then it will read the
	 * data corresponding to that index and return the index which was read
	 * from.
	 *
	 * if the index is not found, because either data could not be read or the
	 * index found does not end up matching the requested index, that index
	 * will instead be returned.
	 *
	 * \param values The arrays into which the values are read into.
	 * \param rinfo Information about how to access the persistent information.
	 */
	template<typename T, size_t N>
	int read_grid(T(*values)[N], symphas::io::read_info rinfo);
	



	//! Read a grid from a file in a standardized header-block way.
	/*!
	 * Read a grid from a file in a standardized header-block way. The
	 * methods to read the header, open the file and the block are provided
	 * as parameters.
	 * 
	 * \param grid The array to where the data from file is written.
	 * \param rinfo Information about the file to be read.
	 * \param open_file_f The function to open the file to be read.
	 * \param read_header_f The function to read the header from the file.
	 * \param read_block_f The function to read the block data, containing
	 * the array data.
	 */
	template<typename T, typename Fo, typename Fc, typename Fh, typename Fb>
	int read_grid_standardized(T* grid, symphas::io::read_info rinfo, Fo open_file_f, Fc close_file_f, Fh read_header_f, Fb read_block_f)
	{
		auto f = open_file_f(rinfo.name);

		int index;
		auto bginfo = read_header_f(&index, f, true);
		int prev;

		if (!params::single_output_file)
		{
			read_block_f(grid, bginfo, f);
			prev = index;
		}
		else
		{
			do
			{
				read_block_f(grid, bginfo, f);
				prev = index;
				read_header_f(&index, f, false);
			} while (index <= rinfo.index && index > BAD_INDEX);
		}
		close_file_f(f);

		if (prev != rinfo.index)
		{
			if (index == BAD_INDEX)
			{
				fprintf(SYMPHAS_ERR, SYMPHAS_MSG_BAD_INDEX_READ, rinfo.index, rinfo.name);
				return BAD_INDEX;
			}
			else
			{
				int pmax = std::abs(rinfo.index - prev);
				int cmax = std::abs(rinfo.index - index);
				int close_index = (pmax > cmax) ? index : prev;

				fprintf(SYMPHAS_LOG, SYMPHAS_MSG_READ_DIFFERENT_INDEX, rinfo.index, rinfo.name, close_index);
				return close_index;
			}
		}
		else
		{
			return prev;
		}
	}
}


