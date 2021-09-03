
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
	 * file, such as when params::single_output_file is chosen.
	 * 
	 * \param grid The values into which to read the data.
	 * \param ginfo The grid parameter specification.
	 * \param f The pointer to the file that is accessed.
	 */
	template<typename T>
	void read_block(T* grid, symphas::grid_info ginfo, FILE* f);

	//! Plain text implementation of reading data.
	template<typename T>
	int read_grid(T* grid, symphas::io::read_info rinfo)
	{
		return read_grid_standardized(grid, rinfo, symphas::io::gp::open_gpgridf, fclose, 
			symphas::io::gp::read_gp_header, read_block<T>);
	}

}

//! \cond

//! Specialization of symphas::io::gp::col::read_block().
template<>
void symphas::io::gp::col::read_block<scalar_t>(scalar_t* grid, symphas::grid_info ginfo, FILE* f);
//! Specialization of symphas::io::gp::col::read_block().
template<>
void symphas::io::gp::col::read_block<complex_t>(complex_t* grid, symphas::grid_info ginfo, FILE* f);
//! Specialization of symphas::io::gp::col::read_block().
template<>
void symphas::io::gp::col::read_block<double[2]>(double(*grid)[2], symphas::grid_info ginfo, FILE* f);
//! Specialization of symphas::io::gp::col::read_block().
template<>
void symphas::io::gp::col::read_block<vector_t<3>>(vector_t<3>* grid, symphas::grid_info ginfo, FILE* f);
//! Specialization of symphas::io::gp::col::read_block().
template<>
void symphas::io::gp::col::read_block<vector_t<2>>(vector_t<2>* grid, symphas::grid_info ginfo, FILE* f);
//! Specialization of symphas::io::gp::col::read_block().
template<>
void symphas::io::gp::col::read_block<vector_t<1>>(vector_t<1>* grid, symphas::grid_info ginfo, FILE* f);

//! \endcond

