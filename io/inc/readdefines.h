
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



// \cond

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

#define SPECIALIZE_READ_BLOCK_FILE(NAMESPACE, F) \
template<> void symphas::io::NAMESPACE::read_block(scalar_t* grid, symphas::grid_info ginfo, F* f); \
template<> void symphas::io::NAMESPACE::read_block(complex_t* grid, symphas::grid_info ginfo, F* f); \
template<> void symphas::io::NAMESPACE::read_block(double_arr2* grid, symphas::grid_info ginfo, F* f); \
template<> void symphas::io::NAMESPACE::read_block(vector_t<3>* grid, symphas::grid_info ginfo, F* f); \
template<> void symphas::io::NAMESPACE::read_block(vector_t<2>* grid, symphas::grid_info ginfo, F* f); \
template<> void symphas::io::NAMESPACE::read_block(vector_t<1>* grid, symphas::grid_info ginfo, F* f); \
template<> void symphas::io::NAMESPACE::read_block(scalar_ptr_t (&grid)[3], symphas::grid_info ginfo, F* f); \
template<> void symphas::io::NAMESPACE::read_block(scalar_ptr_t (&grid)[2], symphas::grid_info ginfo, F* f); \
template<> void symphas::io::NAMESPACE::read_block(scalar_ptr_t (&grid)[1], symphas::grid_info ginfo, F* f);

#define SPECIALIZE_READ_BLOCK(NAMESPACE) SPECIALIZE_READ_BLOCK_FILE(NAMESPACE, FILE)


// \endcond


namespace symphas::io
{
	struct read_info;
}

void swap(symphas::io::read_info& first, symphas::io::read_info& second);

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

		read_info(iter_type index, size_t id, const char* name)
			: index{ index }, id{ id }, name{ (name) ? new char[std::strlen(name) + 1] : nullptr }
		{
			if (name)
			{
				std::strcpy(this->name, name);
			}
		}

		read_info(iter_type index, size_t id, const char* input, bool checkpoint)
			: index{ index }, id{ id }, name{ nullptr }
		{
			char buffer[BUFFER_LENGTH];
			if (checkpoint)
			{
				if (params::single_input_file)
				{
					snprintf(buffer, STR_ARR_LEN(buffer), OUTPUT_CHECKPOINT_FMT,
						input, id);
				}
				else
				{
					snprintf(buffer, STR_ARR_LEN(buffer), OUTPUT_CHECKPOINT_INDEX_FMT,
						input, id, index);
				}
			}
			else
			{
				snprintf(buffer, STR_ARR_LEN(buffer), "%s", input);
			}

			name = new char[std::strlen(buffer) + 1];
			std::strcpy(name, buffer);

		}

		read_info(read_info const& other) : read_info{ other.index, other.id, other.name } {}
		read_info(read_info&& other) noexcept : read_info()
		{
			::swap(*this, other);
		}
		read_info& operator=(read_info other)
		{
			::swap(*this, other);
			return *this;
		}

		friend void ::swap(symphas::io::read_info& first, symphas::io::read_info& second);

		~read_info()
		{
			delete[] name;
		}


	protected:

		read_info() : read_info(0, 0, nullptr) {}
	};

}

inline void swap(symphas::io::read_info& first, symphas::io::read_info& second)
{
	using std::swap;
	swap(first.index, second.index);
	swap(first.id, second.id);
	swap(first.name, second.name);
}



namespace symphas::io
{
	
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
	int read_grid(T* values, symphas::io::read_info const& rinfo);
	

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
	int read_grid(T(*values)[N], symphas::io::read_info const& rinfo);

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
	 * \param values The arrays into which the values are read into. The values
	 * have the type vector of arrays.
	 * \param rinfo Information about how to access the persistent information.
	 */
	template<typename T, size_t N>
	int read_grid(T* (&values)[N], symphas::io::read_info const& rinfo);

	//! Read the header from the given data source. 
	/*!
	 * The header will be read from the given data source, initializing
	 * the grid info object and returning it, if the information is available.
	 * This is typically used if the grid parameters contained by a file are 
	 * not known before reading the file.
	 * 
	 * The full header is only read when the given index is `nullptr` or
	 * initialized to less than 0. Otherwise, only the index will be parsed.
	 * This is to support the params::single_input_file and params::single_output_file
	 * configurations.
	 * 
	 * \param rinfo Information about how to access the persistent information.
	 * \param index Returns the index from the header that is parsed. When this
	 * is `nullptr` or initialized to less than 0, the full header format is
	 * parsed.
	 */
	symphas::grid_info read_header(symphas::io::read_info const& rinfo, iter_type* index = nullptr);
	

	//! Read the header from the given data source. 
	/*!
	 * The header will be read from the given data source, initializing
	 * the grid info object and returning it, if the information is available.
	 * This is typically used if the grid parameters contained by a file are 
	 * not known before reading the file.
	 *
	 * The full header is only read when the given index is `nullptr` or
	 * initialized to less than 0. Otherwise, only the index will be parsed.
	 * This is to support the params::single_input_file and params::single_output_file
	 * configurations.
	 *
	 * \param f The file pointer from where the header is read.
	 * \param index Returns the index from the header that is parsed. When this
	 * is `nullptr` or initialized to less than 0, the full header format is
	 * parsed.
	 * 
	 * \tparam F The type of the file pointer from where the header is read.
	 */
	template<typename F>
	symphas::grid_info read_header(F* f, iter_type* index = nullptr);



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
	template<typename value_type, typename Fo, typename Fc, typename Fb>
	int read_grid_standardized(value_type grid, symphas::io::read_info const& rinfo, Fo open_file_f, Fc close_file_f, Fb read_block_f)
	{
		auto f = open_file_f(rinfo.name);

		int index = -1;
		symphas::grid_info bginfo = read_header(f, &index);
		int prev;

		if (!params::single_input_file)
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
				read_header(f, &index);
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



