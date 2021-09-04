
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


namespace symphas::io::gp::mov
{
	//! Plain text implementation of reading data.
	template<typename T>
	int read_grid(T* grid, symphas::io::read_info const& rinfo)
	{
		return symphas::io::gp::read_grid(grid, rinfo);
	}
	
	DECLARE_GP_HEADER_FUNCTIONS
}

