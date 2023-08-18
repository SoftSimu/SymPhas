
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
 * PURPOSE: Defines basic parts of the input (reading) functionality for
 * xdr (binary) files.
 *
 * ***************************************************************************
 */

#pragma once


#include "readdefines.h"

struct XDRFILE;

namespace symphas::io::xdr
{

	//! Open the file which was written in the xdr (binary) format.
	XDRFILE* open_xdrgridf(const char* name);

	//! Close the xdr (binary) file.
	void xdrfile_close(XDRFILE*);

	//! Read one segment of data from the input file.
	/*!
	 * Read one data block from the saved output, into the grid parameter.
	 * This is used when there are multiple data indices written to the same
	 * file, such as when params::single_input_file is chosen.
	 *
	 * \param grid The values into which to read the data.
	 * \param ginfo The grid parameter specification.
	 * \param f The pointer to the file that is accessed.
	 */
	template<typename value_type>
	void read_block(value_type grid, symphas::grid_info ginfo, XDRFILE* f);

	//! Read the xdr data into the given array.
	/*!
	 * This is an xdr implementation of symphas::io::read_grid().
	 *
	 * It will read the header and other file information until it reaches the
	 * index, and then it will read the data corresponding to that index and
	 * return the index which was read from.
	 * If the index was not found, because either data could not be read or the
	 * index found does not end up matching the requested index, that index
	 * will instead be returned.
	 *
	 * \param grid The array into which values are read.
	 * \param rinfo Information about the data file that is read from.
	 */
	template<typename value_type>
	int read_grid(value_type grid, symphas::io::read_info const& rinfo)
	{
		return read_grid_standardized(grid, rinfo, open_xdrgridf, xdrfile_close, read_block<value_type>);
	}


	//! Read the header of the xdr datafile.
	/*!
	 * When the given index is less than 0, then the full file
	 * header is parsed.
	 */
	template<typename F>
	symphas::grid_info read_header(F* f, int* index) { return { nullptr, 0 }; }
	template<>
	symphas::grid_info read_header(XDRFILE* f, int* index);

	inline symphas::grid_info read_header(symphas::io::read_info const& rinfo, iter_type* index = nullptr)
	{
		XDRFILE* x = open_xdrgridf(rinfo.get_name());
		symphas::grid_info ginfo = read_header(x, index);
		xdrfile_close(x);

		return ginfo;
	}

}

//! \cond
SPECIALIZE_READ_BLOCK_FILE(xdr, XDRFILE)
//! \endcond



