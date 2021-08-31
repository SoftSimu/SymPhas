
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
	//! Read the header of the xdr datafile.
	symphas::grid_info read_xdr_header(int* ginfo, XDRFILE* f, bool read_header = true);

	//! Open the file which was written in the xdr (binary) format.
	XDRFILE* open_xdrgridf(const char* name);

	//! Close the xdr (binary) file.
	void xdrfile_close(XDRFILE*);

	//! Read one segment of data from the input file.
	/*!
	 * Read one data block from the saved output, into the grid parameter.
	 * This is used when there are multiple data indices written to the same
	 * file, such as when params::single_output_file is chosen.
	 *
	 * \param grid The values into which to read the data.
	 * \param ginfo The grid parameter specification.
	 * \param f The pointer to the file that is accessed.
	 */
	template<typename T>
	void read_block(T* grid, symphas::grid_info ginfo, XDRFILE* f);

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
	template<typename T>
	int read_grid(T* grid, symphas::io::read_info rinfo)
	{
		return read_grid_standardized(grid, rinfo, open_xdrgridf, xdrfile_close, read_xdr_header, read_block<T>);
	}
}

//! \cond

//! Specialization based on symphas::io::xdr::read_block.
template<>
void symphas::io::xdr::read_block<scalar_t>(scalar_t* grid, symphas::grid_info ginfo, XDRFILE* f);
//! Specialization based on symphas::io::xdr::read_block.
template<>
void symphas::io::xdr::read_block<complex_t>(complex_t* grid, symphas::grid_info ginfo, XDRFILE* f);
//! Specialization based on symphas::io::xdr::read_block.
template<>
void symphas::io::xdr::read_block<double[2]>(double(*grid)[2], symphas::grid_info ginfo, XDRFILE* f);
//! Specialization based on symphas::io::xdr::read_block.
template<>
void symphas::io::xdr::read_block<vector_t<3>>(vector_t<3>* grid, symphas::grid_info ginfo, XDRFILE* f);
//! Specialization based on symphas::io::xdr::read_block.
template<>
void symphas::io::xdr::read_block<vector_t<2>>(vector_t<2>* grid, symphas::grid_info ginfo, XDRFILE* f);
//! Specialization based on symphas::io::xdr::read_block.
template<>
void symphas::io::xdr::read_block<vector_t<1>>(vector_t<1>* grid, symphas::grid_info ginfo, XDRFILE* f);

//! \endcond



