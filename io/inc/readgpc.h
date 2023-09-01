
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
 * PURPOSE: Defines the text input functionality for plain text files.
 *
 * ***************************************************************************
 */

#pragma once

#include "readdefines.h"
#include "readgp.h"


namespace symphas::io::gp::col
{

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
	void read_block(value_type values, symphas::grid_info ginfo, FILE* f);

	//! Plain text implementation of reading data.
	template<typename value_type>
	int read_grid(value_type values, symphas::io::read_info const& rinfo, symphas::grid_info* ginfo = nullptr)
	{
		return read_grid_standardized(values, rinfo, ginfo, symphas::io::gp::open_gpgridf, fclose, read_block<value_type>);
	}

	DECLARE_GP_HEADER_FUNCTIONS
}

//! \cond

SPECIALIZE_READ_BLOCK(gp::col)

//! \endcond

