
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
 * MODULE:  datatypes
 * PURPOSE: Defines how to update the boundaries, depending on the type
 * and dimension.
 *
 * ***************************************************************************
 */

#pragma once

#include "boundaryupdatecuda.cuh"

#ifdef USING_CUDA

#include <cuda_runtime.h>

// *********************************************************************
/* DEFAULT BOUNDARY ALGORITHMS
 */

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::FRONT, 2>::
operator()(const grid::Boundary<T, 2>* b, RegionalGridCUDA<T, 3>& grid,
           double time) {
  regional_update_boundary(symphas::lib::side_list<Side::FRONT, Side::FRONT>{},
                           b, grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::BACK, 2>::
operator()(const grid::Boundary<T, 2>* b, RegionalGridCUDA<T, 3>& grid,
           double time) {
  regional_update_boundary(symphas::lib::side_list<Side::BACK, Side::BACK>{}, b,
                           grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 2>::
operator()(const grid::Boundary<T, 2>* b, RegionalGridCUDA<T, 3>& grid,
           double time) {
  regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, b, grid,
                           time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 2>::
operator()(const grid::Boundary<T, 2>* b, RegionalGridCUDA<T, 3>& grid,
           double time) {
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, b, grid,
                           time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::TOP, 2>::
operator()(const grid::Boundary<T, 2>* b, RegionalGridCUDA<T, 3>& grid,
           double time) {
  regional_update_boundary(
      symphas::lib::side_list<Side::TOP, Side::TOP, Side::TOP>{}, b, grid,
      time);
  regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::BACK>{}, b,
                           grid, time);
  regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::FRONT>{}, b,
                           grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::DEFAULT, Side::BOTTOM,
    2>::operator()(const grid::Boundary<T, 2>* b, RegionalGridCUDA<T, 3>& grid,
                   double time) {
  regional_update_boundary(
      symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM, Side::BOTTOM>{}, b,
      grid, time);
  regional_update_boundary(symphas::lib::side_list<Side::BOTTOM, Side::BACK>{},
                           b, grid, time);
  regional_update_boundary(symphas::lib::side_list<Side::BOTTOM, Side::FRONT>{},
                           b, grid, time);
}
#endif
