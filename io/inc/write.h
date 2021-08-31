
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
 * PURPOSE: Implements the write methods, where the write algorithm is chosen
 * based on the writer chosen. Also includes overloads of the write methods.
 *
 * ***************************************************************************
 */

#pragma once

#include "writeincludes.h"


template<typename T>
void symphas::io::save_grid_plotting(const T* values, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	SWITCH_IO_WRITE(save_grid_plotting<T>(values, winfo, ginfo))
}

template<typename T, size_t N>
void symphas::io::save_grid_plotting(T(*values)[N], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	using data_unit = T[N];
	SWITCH_IO_WRITE(save_grid_plotting<data_unit>(values, winfo, ginfo))
}

template<typename T>
void symphas::io::save_grid(const T* values, symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	SWITCH_IO_WRITE(save_grid<T>(values, winfo, ginfo))
}

template<typename T, size_t N>
void symphas::io::save_grid(T(*values)[N], symphas::io::write_info winfo, symphas::grid_info ginfo)
{
	using data_unit = T[N];
	SWITCH_IO_WRITE(save_grid<data_unit>(values, winfo, ginfo))
}


//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T>
void symphas::io::save_grid(const T(*values), symphas::grid_info ginfo)
{
	save_grid<T>(values, symphas::io::write_info{}, ginfo);
}

//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, size_t N>
void symphas::io::save_grid(T(*values)[N], symphas::grid_info ginfo)
{
	using data_unit = T[N];
	save_grid<data_unit>(values, symphas::io::write_info{}, ginfo);
}




//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T>
void symphas::io::save_grid(const T(*values), symphas::io::write_info winfo, const len_type* dims, size_t dimension)
{
	symphas::io::save_grid(values, winfo, symphas::grid_info{ dims, dimension });
}

//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, size_t N>
void symphas::io::save_grid(T(*values)[N], symphas::io::write_info winfo, const len_type* dims, size_t dimension)
{
	using data_unit = T[N];
	symphas::io::save_grid(values, winfo, symphas::grid_info{ dims, dimension });
}




//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T>
void symphas::io::save_grid(const T(*values), const len_type* dims, size_t dimension)
{
	symphas::io::save_grid(values, symphas::io::write_info{}, dims, dimension);
}

//! Default argument overload of symphas::io::save_grid(const T*, symphas::io::write_info, symphas::grid_info).
template<typename T, size_t N>
void symphas::io::save_grid(T(*values)[N], const len_type* dims, size_t dimension)
{
	using data_unit = T[N];
	symphas::io::save_grid(values, symphas::io::write_info{}, dims, dimension);
}


template<typename M>
void symphas::io::write_plot_config(M const& model, const char* directory, const char* const* names, SaveParams const& save)
{
	SWITCH_IO_WRITE(write_plot_config(model, directory, names, save))
}









