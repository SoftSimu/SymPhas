
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
 * PURPOSE: Defines the text output functionality for binary text files 
 * in xdr format.
 *
 * ***************************************************************************
 */

#pragma once

#include "writedefines.h"

namespace symphas::io
{
	//! Defines elements used in input and output for the xdr format.
	/*!
	 * The xdr format is a binary format, useful for large data files.
	 */
	namespace xdr {}
}

namespace symphas::io::xdr
{


	//! Writes a data array to a file for plotting.
	/*!
	 * The given data is written to a file set by \p winfo. The information
	 * about the grid and the grid type is used to select a writing method. The
	 * output format and amount of data written is chosen to allow the resulting
	 * file to be ingested by a plotting program.
	 *
	 * \param grid The data which is written to the file.
	 * \param winfo Information about the file that is written.
	 * \param ginfo Information about the grid.
	 */
	template<typename T>
	void save_grid_plotting(const T* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);


	//! Writes a data array to a file.
	/*!
	 * The given data is written to a file set by \p winfo. The information
	 * about the grid and the grid type is used to select a writing method.
	 *
	 * \param grid The data which is written to the file.
	 * \param winfo Information about the file that is written.
	 * \param ginfo Information about the grid.
	 */
	template<typename T>
	void save_grid(const T* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);

	/*
	 * prints the plotting data file
	 * takes as parameters the configuration and the number of systems that
	 * will be plotted
	 */

	template<typename M>
	void write_plot_config(M const&, const char*, const char* const*, SaveParams const&) {}

}


//! \cond

//! Specialization of symphas::io::xdr::save_grid_plotting().
template<>
void symphas::io::xdr::save_grid_plotting<scalar_t>(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);
//! Specialization of symphas::io::xdr::save_grid_plotting().
template<>
void symphas::io::xdr::save_grid_plotting<complex_t>(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);
//! Specialization of symphas::io::xdr::save_grid_plotting().
template<>
void symphas::io::xdr::save_grid_plotting<double[2]>(const double(*grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo);
template<>
void symphas::io::xdr::save_grid_plotting<vector_t<3>>(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);
//! Specialization of symphas::io::xdr::save_grid_plotting().
template<>
void symphas::io::xdr::save_grid_plotting<vector_t<2>>(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);
//! Specialization of symphas::io::xdr::save_grid_plotting().
template<>
void symphas::io::xdr::save_grid_plotting<vector_t<1>>(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);

//! Specialization of symphas::io::xdr::save_grid().
template<>
void symphas::io::xdr::save_grid<scalar_t>(const scalar_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);
//! Specialization of symphas::io::xdr::save_grid().
template<>
void symphas::io::xdr::save_grid<complex_t>(const complex_t* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);
//! Specialization of symphas::io::xdr::save_grid().
template<>
void symphas::io::xdr::save_grid<double[2]>(const double(*grid)[2], symphas::io::write_info winfo, symphas::grid_info ginfo);
//! Specialization of symphas::io::xdr::save_grid().
template<>
void symphas::io::xdr::save_grid<vector_t<3>>(const vector_t<3>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);
//! Specialization of symphas::io::xdr::save_grid().
template<>
void symphas::io::xdr::save_grid<vector_t<2>>(const vector_t<2>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);
//! Specialization of symphas::io::xdr::save_grid().
template<>
void symphas::io::xdr::save_grid<vector_t<1>>(const vector_t<1>* grid, symphas::io::write_info winfo, symphas::grid_info ginfo);

//! \endcond

