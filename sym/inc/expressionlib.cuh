
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
 * PURPOSE: Defines the basic functionality of using the symbolic algebra
 * framework.
 *
 * ***************************************************************************
 */

#pragma once

#include "dataiterator.cuh"
#include "expressionlib.h"

template <Axis ax, typename T, size_t D>
struct VectorComponentRegionData<ax, CUDADataType<T> *, D> {
  T *values;
  grid::select_region_cuda<D> region;
  len_type dims[D];
  T empty;

  template <size_t... Is>
  VectorComponentRegionData(T *data, grid::select_region_cuda<D> const &region,
                            const len_type (&dims)[D], T const &empty,
                            std::index_sequence<Is...>)
      : values{data}, region{region}, dims{dims[Is]...}, empty{empty} {}
  VectorComponentRegionData(T *data, grid::select_region_cuda<D> const &region,
                            const len_type (&dims)[D], T const &empty)
      : VectorComponentRegionData(data, region, dims, empty,
                                  std::make_index_sequence<D>{}) {}

  __host__ __device__ auto operator[](iter_type n) const {
    iter_type pos[D]{};
    grid::get_grid_position(pos, dims, n);
    return region(pos, values, dims, empty);
  }

  __host__ __device__ auto operator[](iter_type n) {
    iter_type pos[D]{};
    grid::get_grid_position(pos, dims, n);
    return region(pos, values, dims, empty);
  }
};

template <typename T, size_t D>
struct GridDataCUDA {
  T *values;
  __host__ __device__ GridDataCUDA(T *values) : values{values} {}

  __host__ __device__ const auto &operator[](iter_type i) const {
    return values[i];
  }
  __host__ __device__ auto &operator[](iter_type i) { return values[i]; }
};

template <typename T, size_t D>
struct RegionGridDataCUDA {
  T *values;
  len_type dims[D];
  len_type stride[D];

  template <size_t... Is>
  __host__ __device__ RegionGridDataCUDA(T *values, len_type const (&dims)[D],
                                         len_type const (&stride)[D],
                                         std::index_sequence<Is...>)
      : values{values}, dims{dims[Is]...}, stride{stride[Is]...} {}
  __host__ __device__ RegionGridDataCUDA(T *values, len_type const (&dims)[D],
                                         len_type const (&stride)[D])
      : RegionGridDataCUDA(values, dims, stride,
                           std::make_index_sequence<D>{}) {}

  __host__ __device__ const auto &operator[](iter_type i) const {
    return values[i];
  }
  __host__ __device__ auto &operator[](iter_type i) { return values[i]; }
};

template <Axis ax, size_t N, typename T>
VectorComponentRegionData<ax, CUDADataType<T> *, N>
expr::resolve_axis_component(
    RegionalGridCUDA<any_vector_t<T, N>, N> const &data) {
  return {data.values[symphas::axis_to_index(ax)], data.region, data.dims,
          data.empty[symphas::axis_to_index(ax)]};
}

template <Axis ax, size_t N, typename T>
VectorComponentRegionData<ax, CUDADataType<T> *, N>
expr::resolve_axis_component(RegionalGridCUDA<any_vector_t<T, N>, N> &data) {
  return {data.values[symphas::axis_to_index(ax)], data.region, data.dims,
          data.empty[symphas::axis_to_index(ax)]};
}

template <Axis ax, size_t N, typename T>
VectorComponentData<ax, T *, N> expr::resolve_axis_component(
    MultiBlockCUDA<N, T> const &data) {
  return {data.values[symphas::axis_to_index(ax)]};
}

template <Axis ax, size_t N, typename T>
VectorComponentData<ax, T *, N> expr::resolve_axis_component(
    MultiBlockCUDA<N, T> &data) {
  return {data.values[symphas::axis_to_index(ax)]};
}
