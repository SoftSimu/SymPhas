
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
 * PURPOSE: Defines the grid types used in the finite difference and other
 * semi-discrete numerical solvers.
 *
 * ***************************************************************************
 */

#pragma once

#include "grid.h"

#ifdef USING_CUDA

namespace symphas {
namespace cuda {}
}  // namespace symphas

template <typename T>
struct carry_value_cuda;

template <size_t N, typename T>
struct multi_value_cuda;

namespace grid {

template <typename T>
struct selected_entry_cuda;

template <size_t D>
struct select_grid_index_cuda;

template <size_t D>
struct select_region_cuda;
}  // namespace grid

//! Manages an array of values of arbitrary type.
/*!
 * Basic array type object used in constructing the finite difference grid.
 * Values are always initialized to the empty value.
 *
 * \tparam T The value type of the underlying array.
 */
template <typename T>
struct BlockCUDA;

// ***********************************************************************************************

//! Manages an array of values of arbitrary type.
/*!
 * Basic array type object used in constructing the finite difference grid.
 * Values are always initialized to the empty value.
 *
 * \tparam T The value type of the underlying array.
 */
template <size_t N, typename T>
struct MultiBlockCUDA;

// ***********************************************************************************************

//! A grid object of arbitrary dimension and arbitrary value type.
/*!
 * A grid object of arbitrary dimension and arbitrary value type. The grid
 * forms the basis of a finite difference grid. Only manages its own
 * dimensions and list of values. The values are inherited from ::Block, meaning
 * that the data is flattened and is not `D`-dimensional in memory.
 *
 * \tparam T The value type of the underlying array.
 * \tparam D The dimension of the grid.
 */
template <typename T, size_t D>
struct GridCUDA;

//! A grid object of arbitrary dimension and arbitrary value type.
/*!
 * A grid object of arbitrary dimension and arbitrary value type.
 * This grid implementation is an extension of the base ::Grid, but it
 * contains virtual boundary conditions.
 * The grid
 * forms the basis of a finite difference grid. Only manages its own
 * dimensions and list of values. The values are inherited from ::BlockCUDA,
 * meaning that the data is flattened and is not `D`-dimensional in memory.
 *
 * \tparam T The value type of the underlying array.
 * \tparam D The dimension of the grid.
 */
template <typename T, size_t D>
struct BoundaryGridCUDA;

//! A grid object of arbitrary dimension and arbitrary value type.
/*!
 * A grid object of arbitrary dimension and arbitrary value type.
 * This grid implementation is an extension of the base ::Grid, but it
 * is meant to sub-domain a larger domain when there are many fields.
 * The grid
 * forms the basis of a finite difference grid. Only manages its own
 * dimensions and list of values. The values are inherited from ::Block, meaning
 * that the data is flattened and is not `D`-dimensional in memory.
 *
 * \tparam T The value type of the underlying array.
 * \tparam D The dimension of the grid.
 */
template <typename T, size_t D>
struct RegionalGridCUDA;

#endif
