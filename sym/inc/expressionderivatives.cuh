
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
 * MODULE:  expr
 * PURPOSE: Implements derivatives of the symbolic algebra library.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressionderivatives.h"
#include "expressionlogic.cuh"
#include "gridfunctions.cuh"

namespace symphas::internal {

template <typename T, size_t D, typename E>
void update_temporary_grid(GridCUDA<T, D>& grid, OpEvaluable<E> const& e) {}

template <typename T, size_t D>
void update_temporary_grid(RegionalGridCUDA<T, D>& grid, ...) {
  grid::region_interval<D> interval;
  grid::resize_adjust_region(grid, interval);
}

template <typename T, size_t D>
void update_temporary_grid(RegionalGridCUDA<T, D>& grid,
                           grid::region_interval<D> interval) {
  for (iter_type i = 0; i < D; ++i) {
    interval[i][0] -= grid.region.boundary_size;
    interval[i][1] += grid.region.boundary_size;
  }
  grid::resize_adjust_region(grid, interval);
}

template <typename T, size_t D>
void update_temporary_grid(RegionalGridCUDA<T, D>& grid,
                           grid::region_interval_multiple<D> const& regions) {
  update_temporary_grid(grid, +regions);
}

template <typename T, size_t D, typename E>
void update_temporary_grid(RegionalGridCUDA<T, D>& grid,
                           OpEvaluable<E> const& e) {
  update_temporary_grid(grid,
                        expr::iterable_domain(*static_cast<E const*>(&e)));
}
}  // namespace symphas::internal
