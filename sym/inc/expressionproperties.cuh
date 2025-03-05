
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
 * PURPOSE: Defines methods of obtaining properties of expressions, in
 * particular about the underlying data of the expressions.
 *
 * ***************************************************************************
 */

#pragma once

#include <tuple>

#include "expressionlib.cuh"

namespace expr {
namespace {

//! Obtains the dimensions from the Grid compatible instance.
template <typename T>
grid::dim_list data_dimensions_cast(GridCUDA<T, 1> const* data) {
  return {data->dims[0]};
}

//! Obtains the dimensions from the Grid compatible instance.
template <typename T>
grid::dim_list data_dimensions_cast(GridCUDA<T, 2> const* data) {
  return {data->dims[0], data->dims[1]};
}

//! Obtains the dimensions from the Grid compatible instance.
template <typename T>
grid::dim_list data_dimensions_cast(GridCUDA<T, 3> const* data) {
  return {data->dims[0], data->dims[1], data->dims[2]};
}
//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto iterable_domain_cast(GridCUDA<T, D> const* data) {
  return grid::get_iterable_domain(*data);
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto iterable_domain_cast(BoundaryGridCUDA<T, D> const* data) {
  return grid::get_iterable_domain(*data);
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto iterable_domain_cast(RegionalGridCUDA<T, D> const* data) {
  region_interval_multiple<D> regions(data->dims, data->region.boundary_size);
  return regions += grid::get_iterable_domain(*data);
}

//! The iterable_domain of a typical data object is 1.
template <typename T>
inline auto iterable_domain_cast(BlockCUDA<T> const* data) {
  return iterable_domain_cast(grid::get_iterable_domain(*data));
}
}  // namespace
}  // namespace expr