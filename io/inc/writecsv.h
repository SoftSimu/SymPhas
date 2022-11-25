
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
 * PURPOSE: Defines the text output functionality for writing data csv.
 *
 * ***************************************************************************
 */


#pragma once


#include "writedefines.h"
#include "writegp.h"

//! \cond


namespace symphas::io
{
	//! Defines elements used in input and output for the text format.
	/*!
	 * The format of input/output is based on using it with the gnuplot
	 * utility.
	 */
	namespace csv {}
}

namespace symphas::io::csv
{

	template<size_t D, typename Sp, typename... S>
	void write_plot_config(Model<D, Sp, S...> const& model, const char* directory, const char* const* names, SaveParams const& save) {}

	//! \cond
	DECLARE_SAVE_GRID_PLOTTING_FUNCTIONS
	//! \endcond
	
	template<typename value_type>
	void save_grid(value_type&& grid, symphas::io::write_info winfo, symphas::grid_info ginfo)
	{
		symphas::io::gp::save_grid(std::forward<value_type>(grid), winfo, ginfo);
	}

}


