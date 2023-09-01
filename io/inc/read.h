
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
struct nullptr_for_read
{
	auto operator()(symphas::io::read_info const& rinfo, symphas::grid_info* ginfo) const
	{
		SWITCH_IO_READ(read_grid<T*>(nullptr, rinfo, ginfo), BAD_INDEX)
	}
};

template<typename T, size_t D>
struct nullptr_for_read<T[D]>
{
	auto operator()(symphas::io::read_info const& rinfo, symphas::grid_info* ginfo) const
	{
		T* data[D]{};
		for (iter_type i = 0; i < D; ++i)
		{
			data[i] = nullptr;
		}
		SWITCH_IO_READ(read_grid<T*(&)[D]>(data, rinfo, ginfo), BAD_INDEX)
	}
};

template<typename T>
int symphas::io::read_grid(symphas::io::read_info const& rinfo, symphas::grid_info* ginfo)
{
	return nullptr_for_read<T>{}(rinfo, ginfo);
}

template<typename T>
int symphas::io::read_grid(T* &values, symphas::io::read_info const& rinfo, symphas::grid_info* ginfo)
{
	using data_unit = T*;
	SWITCH_IO_READ(read_grid<data_unit>(values, rinfo, ginfo), BAD_INDEX)
}

template<typename T, size_t N>
int symphas::io::read_grid(T(*values)[N], symphas::io::read_info const& rinfo, symphas::grid_info* ginfo)
{
	using data_unit = T(*)[N];
	SWITCH_IO_READ(read_grid<data_unit>(values, rinfo, ginfo), BAD_INDEX)
}

template<typename T, size_t N>
int symphas::io::read_grid(T* (&values)[N], symphas::io::read_info const& rinfo, symphas::grid_info* ginfo)
{
	using data_unit = T * (&)[N];
	SWITCH_IO_READ(read_grid<data_unit>(values, rinfo, ginfo), BAD_INDEX)
}

inline symphas::grid_info symphas::io::read_header(symphas::io::read_info const& rinfo, iter_type* index)
{
	SWITCH_IO_READ(read_header(rinfo, index), grid::dim_list())
}

template<typename F>
inline symphas::grid_info symphas::io::read_header(F* f, iter_type* index)
{
	SWITCH_IO_READ(read_header(f, index), grid::dim_list())
}


