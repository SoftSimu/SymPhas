
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

#include <cuda_runtime.h>

#include "boundaryupdatecuda.cuh"

// one dimensional boundaries

template <typename T>
__global__ void applyPeriodicBoundaryLeft1d(T* grid, int gridSize,
                                            int boundaryDepth) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < boundaryDepth) {
    grid[idx] = grid[gridSize - boundaryDepth * 2 + idx];
  }
}

template <typename T>
__global__ void applyPeriodicBoundaryRight1d(T* grid, int gridSize,
                                             int boundaryDepth) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < boundaryDepth) {
    grid[gridSize - boundaryDepth + idx] = grid[idx];
  }
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 0>::
operator()(const grid::Boundary<T, 0>*, GridCUDA<T, 1>& grid) {
  int numBlocks = (BOUNDARY_DEPTH + BLOCK_SIZE - 1) / BLOCK_SIZE;
  applyPeriodicBoundaryLeft1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      grid.values, grid.len, BOUNDARY_DEPTH);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::RIGHT,
                                        0>::operator()(const grid::Boundary<T,
                                                                            0>*,
                                                       GridCUDA<T, 1>& grid) {
  int numBlocks = (BOUNDARY_DEPTH + BLOCK_SIZE - 1) / BLOCK_SIZE;
  applyPeriodicBoundaryRight1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      grid.values, grid.len, BOUNDARY_DEPTH);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

// one dimensional boundaries

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 0>::
operator()(const grid::Boundary<T, 0>*, RegionalGridCUDA<T, 1>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC, Side::RIGHT,
    0>::operator()(const grid::Boundary<T, 0>*, RegionalGridCUDA<T, 1>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, grid);
}

// 1 dimension

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 0>::
operator()(const grid::Boundary<T, 0>* b, GridCUDA<T, 1>& grid, double time) {
  auto* bd =
      static_cast<grid::BoundaryApplied<T, 0, BoundaryType::DEFAULT> const*>(b);
  const double v = bd->v;

  ITER_GRID1_LEFT({ bd->update(grid[INDEX], v, 0, time); });
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 0>::
operator()(const grid::Boundary<T, 0>* b, GridCUDA<T, 1>& grid, double time) {
  auto* bd =
      static_cast<grid::BoundaryApplied<T, 0, BoundaryType::DEFAULT> const*>(b);
  const double v = bd->v;
  iter_type L = grid.dims[0];

  ITER_GRID1_RIGHT({ bd->update(grid[INDEX], v, 0, time); }, L);
}

// 1 dimension

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 0>::
operator()(const grid::Boundary<T, 0>* b, RegionalGridCUDA<T, 1>& grid,
           double time) {
  regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, b, grid,
                           time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 0>::
operator()(const grid::Boundary<T, 0>* b, RegionalGridCUDA<T, 1>& grid,
           double time) {
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, b, grid,
                           time);
}
