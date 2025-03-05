
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
 * PURPOSE: Defines an iterator that is used for iterating over data so that
 * expressions can be evaluated into the data.
 *
 * ***************************************************************************
 */

#pragma once

#include "dataiterator.h"
#include "gridcuda.cuh"

namespace symphas::internal {

template <typename T, size_t D>
struct data_value_type<RegionalGridCUDA<any_vector_t<T, D>, D>> {
  using type = any_vector_t<T, D>;
  using ref = multi_value_cuda<D, T>;

  ref operator()(RegionalGridCUDA<any_vector_t<T, D>, D>* data, iter_type n) {
    return (*data)[n];
  }
};

template <typename T, size_t D>
struct data_value_type<RegionalGridCUDA<T, D>> {
  using type = T;
  using ref = carry_value_cuda<T>;

  ref operator()(RegionalGridCUDA<T, D>* data, iter_type n) {
    return (*data)[n];
  }
};

template <typename T>
struct data_value_type<BlockCUDA<T>> {
  using type = T;
  using ref = type&;

  ref operator()(BlockCUDA<T>* data, iter_type n) { return (*data)[n]; }
};

template <size_t N, typename T>
struct data_value_type<MultiBlockCUDA<N, T>> {
  using type = any_vector_t<T, N>;
  using ref = multi_value_cuda<N, T>;

  ref operator()(MultiBlockCUDA<N, T>* data, iter_type n) { return (*data)[n]; }
};

template <typename T, size_t D>
struct data_value_type<BoundaryGridCUDA<any_vector_t<T, D>, D>>
    : data_value_type<MultiBlockCUDA<D, T>> {
  using parent_type = data_value_type<MultiBlockCUDA<D, T>>;
  using typename parent_type::ref;
  using typename parent_type::type;
};

template <typename T, size_t D>
struct data_value_type<GridCUDA<any_vector_t<T, D>, D>>
    : data_value_type<MultiBlockCUDA<D, T>> {
  using parent_type = data_value_type<MultiBlockCUDA<D, T>>;
  using typename parent_type::ref;
  using typename parent_type::type;
};

}  // namespace symphas::internal

namespace grid {

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_iterable_domain(GridCUDA<T, D> const& data) {
  return region_size(data.len);
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_iterable_domain(BoundaryGridCUDA<T, D> const& data) {
  region_interval<D> region(data.dims);
  for (iter_type i = 0; i < D; ++i) {
    region[i][0] = BOUNDARY_DEPTH;
    region[i][1] = data.dims[i] - BOUNDARY_DEPTH;
  }
  return region;
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_iterable_domain(RegionalGridCUDA<T, D> const& data) {
  region_interval<D> region(data.dims);
  if (data.region.len > 0) {
    for (iter_type i = 0; i < D; ++i) {
      region[i][0] = data.region.origin[i] + data.region.boundary_size;
      region[i][1] = data.region.origin[i] + data.region.dims[i] -
                     data.region.boundary_size;
    }
  } else {
    for (iter_type i = 0; i < D; ++i) {
      region[i][0] = data.region.origin[i];
      region[i][1] = data.region.origin[i];
    }
  }
  return region;
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_data_domain(GridCUDA<T, D> const& data) {
  return region_size(data.len);
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_data_domain(BoundaryGridCUDA<T, D> const& data) {
  region_interval<D> region(data.dims);
  for (iter_type i = 0; i < D; ++i) {
    region[i][0] = 0;
    region[i][1] = data.dims[i];
  }
  return region;
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_data_domain(RegionalGridCUDA<T, D> const& data) {
  region_interval<D> region(data.dims);
  for (iter_type i = 0; i < D; ++i) {
    region[i][0] = data.region.origin[i];
    region[i][1] = data.region.origin[i] + data.region.dims[i];
  }
  return region;
}

}  // namespace grid
