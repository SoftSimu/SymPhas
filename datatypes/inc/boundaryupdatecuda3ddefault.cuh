
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

// *********************************************************************
/* DEFAULT BOUNDARY ALGORITHMS
 */

template <Side side0, Side side1, typename T>
void update_default_boundary(symphas::lib::side_list<side0, side1>,
                             const grid::Boundary<T, 2>* b,
                             GridCUDA<T, 3>& grid, double time) {
  auto* bd =
      static_cast<grid::BoundaryApplied<T, 2, BoundaryType::DEFAULT> const*>(b);

  const double* x = bd->v;
  const double* y = bd->v + 2;
  double h[2];

  // go backwards or forwards in iteration depending on the interval
  int fx = (x[0] < x[1]) ? 1 : -1;
  int fy = (y[0] < y[1]) ? 1 : -1;
  h[0] = bd->h[0] * fx;
  h[1] = bd->h[1] * fy;

  //iter_type L = grid.dims[0];
  iter_type M = grid.dims[1];
  //iter_type N = grid.dims[2];

  ///////////////////////////////////////////////////////////////////////////////////////////
  // TODO
  ///////////////////////////////////////////////////////////////////////////////////////////
  T* boundaryHost = new T[BOUNDARY_DEPTH * grid.dims[1]]{};
  T* boundaryDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&boundaryDevice, BOUNDARY_DEPTH * grid.dims[1] * sizeof(T)));

  int numBlocks = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
  // copyLeftBoundaryFromGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
  //     grid.values, boundaryDevice, grid.dims[0], grid.dims[1]);

  CHECK_CUDA_ERROR(cudaMemcpy(boundaryHost, boundaryDevice,
                              BOUNDARY_DEPTH * grid.dims[1] * sizeof(T),
                              cudaMemcpyDeviceToHost));
  ///////////////////////////////////////////////////////////////////////////////////////////

  for (iter_type k = 0; k < BOUNDARY_DEPTH; ++k) {
    // four corners

    for (iter_type i = 0; i < BOUNDARY_DEPTH; ++i) {
      for (iter_type j = 0; j < BOUNDARY_DEPTH; ++j) {
        bd->update(boundaryHost[j * BOUNDARY_DEPTH + i], x[0], y[0], time);
      }
    }
    for (iter_type i = 0; i < BOUNDARY_DEPTH; ++i) {
      for (iter_type j = grid.dims[1] - BOUNDARY_DEPTH; j < grid.dims[1]; ++j) {
        bd->update(boundaryHost[j * BOUNDARY_DEPTH + i], x[0], y[1], time);
      }
    }
    for (iter_type i = grid.dims[1] - BOUNDARY_DEPTH; i < grid.dims[0]; ++i) {
      for (iter_type j = 0; j < BOUNDARY_DEPTH; ++j) {
        bd->update(boundaryHost[j * BOUNDARY_DEPTH + i], x[1], y[0], time);
      }
    }
    for (iter_type i = grid.dims[1] - BOUNDARY_DEPTH; i < grid.dims[0]; ++i) {
      for (iter_type j = grid.dims[1] - BOUNDARY_DEPTH; j < grid.dims[1]; ++j) {
        bd->update(boundaryHost[j * BOUNDARY_DEPTH + i], x[1], y[1], time);
      }
    }

    // edges

    for (iter_type i = 0; i < BOUNDARY_DEPTH; ++i) {
      for (iter_type j = BOUNDARY_DEPTH; j < grid.dims[1] - BOUNDARY_DEPTH;
           ++j) {
        bd->update(boundaryHost[j * BOUNDARY_DEPTH + i], x[0], j * h[1], time);
      }
    }
    for (iter_type i = grid.dims[0] - BOUNDARY_DEPTH; i < grid.dims[0]; ++i) {
      for (iter_type j = BOUNDARY_DEPTH; j < grid.dims[1] - BOUNDARY_DEPTH;
           ++j) {
        bd->update(boundaryHost[j * BOUNDARY_DEPTH + i], x[1], j * h[1], time);
      }
    }
    for (iter_type i = BOUNDARY_DEPTH; i < grid.dims[0] - BOUNDARY_DEPTH; ++i) {
      for (iter_type j = 0; j < BOUNDARY_DEPTH; ++j) {
        bd->update(boundaryHost[j * BOUNDARY_DEPTH + i], i * h[0], y[0], time);
      }
    }
    for (iter_type i = BOUNDARY_DEPTH; i < grid.dims[0] - BOUNDARY_DEPTH; ++i) {
      for (iter_type j = grid.dims[1] - BOUNDARY_DEPTH; j < grid.dims[1]; ++j) {
        bd->update(boundaryHost[j * BOUNDARY_DEPTH + i], i * h[0], y[1], time);
      }
    }

    // face

    for (iter_type i = BOUNDARY_DEPTH; i < grid.dims[0] - BOUNDARY_DEPTH; ++i) {
      for (iter_type j = BOUNDARY_DEPTH; j < grid.dims[1] - BOUNDARY_DEPTH;
           ++j) {
        bd->update(boundaryHost[j * BOUNDARY_DEPTH + i], i * h[0], j * h[1],
                   time);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  // TODO
  ///////////////////////////////////////////////////////////////////////////////////////////
  CHECK_CUDA_ERROR(cudaMemcpy(boundaryDevice, boundaryHost,
                              BOUNDARY_DEPTH * grid.dims[1] * sizeof(T),
                              cudaMemcpyHostToDevice));

  /*copyLeftBoundaryToGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      boundaryDevice, grid.values, grid.dims[0], grid.dims[1]);*/

  CHECK_CUDA_ERROR(cudaFree(boundaryDevice));
  delete[] boundaryHost;
  ///////////////////////////////////////////////////////////////////////////////////////////
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::FRONT, 2>::
operator()(const grid::Boundary<T, 2>* b, GridCUDA<T, 3>& grid, double time) {
  update_default_boundary(symphas::lib::side_list<Side::FRONT, Side::FRONT>{},
                          b, grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::BACK, 2>::
operator()(const grid::Boundary<T, 2>* b, GridCUDA<T, 3>& grid, double time) {
  update_default_boundary(symphas::lib::side_list<Side::BACK, Side::BACK>{}, b,
                          grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 2>::
operator()(const grid::Boundary<T, 2>* b, GridCUDA<T, 3>& grid, double time) {
  update_default_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT>{}, b,
                          grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 2>::
operator()(const grid::Boundary<T, 2>* b, GridCUDA<T, 3>& grid, double time) {
  update_default_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT>{},
                          b, grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::TOP, 2>::
operator()(const grid::Boundary<T, 2>* b, GridCUDA<T, 3>& grid, double time) {
  update_default_boundary(symphas::lib::side_list<Side::TOP, Side::TOP>{}, b,
                          grid, time);
}

template <>
template <typename T>
void symphas::internal::
    update_boundary<BoundaryType::DEFAULT, Side::BOTTOM, 2>::operator()(
        const grid::Boundary<T, 2>* b, GridCUDA<T, 3>& grid, double time) {
  update_default_boundary(symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM>{},
                          b, grid, time);
}

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
