
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
 * MODULE:  sol
 * PURPOSE: Used to generate a VTK visualization of the simulated phase-field
 * model.
 *
 * ***************************************************************************
 */

#pragma once

#include <utility>

#include "gridfunctions.h"

struct ColourPlotUpdater
{
	virtual void update() {}
	virtual ~ColourPlotUpdater() {}
};


struct ColourPlot2d
{
	void init(scalar_t* (&values), len_type* dims, iter_type& index, ColourPlotUpdater* (&updater));

	template<size_t D>
	void init(scalar_t* (&values), grid::region_interval<D> const& interval, iter_type& index, ColourPlotUpdater* (&updater))
	{
		len_type dims[D]{};
		for (iter_type i = 0; i < D; ++i)
		{
			dims[i] = interval[i][1] - interval[i][0];
		}
		init(values, &dims[0], index, updater);
	}
};
