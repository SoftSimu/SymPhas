
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


	//! \cond
	DECLARE_SAVE_GRID_ALL_FUNCTIONS
	//! \endcond

	/*
	 * prints the plotting data file
	 * takes as parameters the configuration and the number of systems that
	 * will be plotted
	 */

	template<typename M>
	void write_plot_config(M const&, const char*, const char* const*, SaveParams const&) {}

}

