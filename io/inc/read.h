
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
 * PURPOSE: Implements the read methods, where the read algorithm is chosen
 * based on the writer chosen.
 *
 * ***************************************************************************
 */

#pragma once

#include "readincludes.h"


template<typename T>
int symphas::io::read_grid(T* values, symphas::io::read_info rinfo)
{
	SWITCH_IO_READ(read_grid<T>(values, rinfo))
}
template<typename T, size_t N>
int symphas::io::read_grid(T(*values)[N], symphas::io::read_info rinfo)
{
	using data_unit = T[N];
	SWITCH_IO_READ(read_grid<data_unit>(values, rinfo))
}


