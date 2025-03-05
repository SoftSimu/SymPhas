
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

#include "boundaryupdatedata.h"

namespace symphas::internal {

template <Side... sides, typename T, size_t D>
__host__ __device__ void get_overlap(iter_type (&overlaps)[D][2],
                                     RegionalGridCUDA<T, D> const& grid,
                                     const iter_type (&origin)[D],
                                     const len_type (&dims)[D]) {
  get_overlap<sides...>(overlaps, grid.region, origin, dims);
}

template <typename T, size_t D, Side... sides>
__host__ __device__ void get_overlap(
    iter_type (&overlaps)[D][2], RegionalGridCUDA<T, D> const& grid,
    boundary_interval<D, sides...> const& interval) {
  get_overlap_adjusted<sides...>(overlaps, grid.region, interval.origin,
                                 interval.dims);
}

template <Side... sides, typename T, size_t D>
__host__ __device__ iter_type get_regional_index(iter_type (&pos)[D],
                                                 RegionalGridCUDA<T, D>& grid) {
  return get_regional_index(pos, grid.region, grid.dims);
}
template <Side... sides, typename T, size_t D>
__host__ __device__ bool side_non_wrapping(RegionalGridCUDA<T, D>& grid) {
  return side_non_wrapping<sides...>(grid.region, grid.dims);
}

template <typename T, size_t D, Side... sides>
__global__ void regionalUpdateBoundaryKernel(
    T* values, iter_type m, boundary_interval<D, sides...>* interval,
    grid::select_region<D>* region) {
  int n = blockIdx.x * blockDim.x + threadIdx.x;
  if (n < grid::length<D>(interval->dims)) {
    iter_type pos[D]{};
    grid::get_grid_position_offset(pos, interval->dims, interval->origin, n);
    iter_type index =
        get_regional_index<sides...>(pos, *region, interval->dims);
    values[index - m] = values[index];
  }
}

template <typename T, size_t D, Side... sides>
void regional_update_boundary(symphas::lib::side_list<sides...>,
                              RegionalGridCUDA<T, D>& grid) {
  // to update the LEFT face.
  // first see if the grid region overlaps the left boundary.
  boundary_interval<D, sides...> interval(grid.dims, grid.region.boundary_size);

  iter_type overlaps[D][2];
  get_overlap(overlaps, grid, interval);

  if (is_overlapping(overlaps) && side_non_wrapping<sides...>(grid)) {
    get_dims_and_origin(interval.origin, interval.dims, overlaps);

    // shift origin to region to copy the boundary from.
    periodic_offset<D, BoundaryType::PERIODIC, sides...> offset(
        grid.dims, grid.region.boundary_size);

    for (iter_type i = 0; i < D; ++i) {
      interval.origin[i] += offset[i];
    }

    get_overlap(overlaps, grid, interval);
    if (is_overlapping(overlaps)) {
      get_dims_and_origin(interval.origin, interval.dims, overlaps);

      iter_type m =
          grid::index_from_position(offset.data(), grid.region.stride);

      boundary_interval<D, sides...>* devInterval;
      grid::select_region<D>* devRegion;
      CHECK_CUDA_ERROR(
          cudaMalloc(&devInterval, sizeof(boundary_interval<D, sides...>)));
      CHECK_CUDA_ERROR(cudaMalloc(&devRegion, sizeof(grid::select_region<D>)));
      CHECK_CUDA_ERROR(cudaMemcpy(devInterval, &interval,
                                  sizeof(boundary_interval<D, sides...>),
                                  cudaMemcpyHostToDevice));
      CHECK_CUDA_ERROR(cudaMemcpy(devRegion, &grid.region,
                                  sizeof(grid::select_region<D>),
                                  cudaMemcpyHostToDevice));

      int numBlocks =
          (grid::length<D>(interval.dims) + BLOCK_SIZE - 1) / BLOCK_SIZE;

      regionalUpdateBoundaryKernel CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
          grid.values, m, devInterval, devRegion);
      CHECK_CUDA_ERROR(cudaPeekAtLastError());
      CHECK_CUDA_ERROR(cudaDeviceSynchronize());
      CHECK_CUDA_ERROR(cudaFree(devInterval));
      CHECK_CUDA_ERROR(cudaFree(devRegion));
    }
  }
}

template <typename T, size_t Dm1>
void apply_boundary_function(
    const grid::BoundaryApplied<T, Dm1, BoundaryType::DEFAULT>* bd,
    carry_value_cuda<T> const& value, const double (&pos)[2], double time) {
  bd->update(*value.value, pos[0], pos[1], time);
}

template <typename T, size_t Dm1>
void apply_boundary_function(
    const grid::BoundaryApplied<any_vector_t<T, Dm1 + 1>, Dm1,
                                BoundaryType::DEFAULT>* bd,
    multi_value_cuda<Dm1 + 1, T> const& value, const double (&pos)[2],
    double time) {
  bd->update(value, pos[0], pos[1], time);
}

template <typename T, size_t D, Side... sides>
void regional_update_boundary(symphas::lib::side_list<sides...>,
                              const grid::Boundary<T, D - 1>* b,
                              RegionalGridCUDA<T, D>& grid, double time) {
  // to update the LEFT face.
  // first see if the grid region overlaps the left boundary.
  boundary_interval<D, sides...> interval(grid.dims, grid.region.boundary_size);

  iter_type overlaps[D][2];
  get_overlap(overlaps, grid, interval);

  if (is_overlapping(overlaps)) {
    iter_type origin[D];
    len_type dims[D];
    get_dims_and_origin(origin, dims, overlaps);

    // shift origin to region to copy the boundary from.
    periodic_offset<D, BoundaryType::PERIODIC, sides...> offset(
        grid.dims, grid.region.boundary_size);

    for (iter_type i = 0; i < D; ++i) {
      origin[i] += offset[i];
    }

    get_overlap(overlaps, grid, origin, dims);
    if (is_overlapping(overlaps)) {
      get_dims_and_origin(origin, dims, overlaps);

      iter_type stride[D]{};
      grid::get_stride(stride, grid.dims);

      iter_type m = grid::index_from_position(offset.data(), stride);

      auto* bd = static_cast<
          grid::BoundaryApplied<T, D - 1, BoundaryType::DEFAULT> const*>(b);

      for (iter_type n = 0; n < grid::length<D>(dims); ++n) {
        iter_type pos[D]{};
        grid::get_grid_position_offset(pos, dims, origin, n);
        iter_type index = grid::index_from_position(pos, stride);

        double pos0[2]{};
        set_boundary_position<sides...>(pos0, bd, pos, origin, bd->v);
        apply_boundary_function(bd, grid[index - m], pos0, time);
      }
    }
  }
}
}  // namespace symphas::internal
