
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
 * PURPOSE: Defines the model used for phase field crystal problems.
 *
 * ***************************************************************************
 */

#pragma once

#ifdef USING_CUDA

#include <cuda_runtime.h>

#include <array>
#include <set>
#include <utility>

#include "solver.h"

namespace symphas::internal {

template <typename T, size_t D, typename grid_type>
void update_max_grid(const SolverSystemFDwSDCUDA<T, D>* _s,
                     grid_type const& s_max, len_type len) {
  RegionalGrid<T, D> work_max(s_max.dims, s_max.empty,
                              s_max.region.boundary_size);
  work_max.adjust(s_max.region.origin, s_max.region.dims);

  CHECK_CUDA_ERROR(cudaMemcpy(work_max.values, s_max.values,
                              s_max.region.len * sizeof(T),
                              cudaMemcpyDeviceToHost));

  auto s_max_it_start =
      symphas::data_iterator_region(work_max, grid::get_iterable_domain(s_max));
  auto s_max_it_end = s_max_it_start + grid::length_interior(s_max);

  for (auto it(s_max_it_start); it < s_max_it_end; ++it) {
    *it = OpVoid{};
  }
  for (iter_type i = 0; i < len; ++i) {
    RegionalGrid<T, D> work(_s[i].dims, _s[i].empty,
                            _s[i].region.boundary_size);
    work.adjust(_s[i].region.origin, _s[i].region.dims);
    CHECK_CUDA_ERROR(cudaMemcpy(work.values, _s[i].values,
                                _s[i].region.len * sizeof(T),
                                cudaMemcpyDeviceToHost));
    auto it0 =
        symphas::data_iterator_region(work, grid::get_iterable_domain(s_max));
    for (auto it(s_max_it_start); it < s_max_it_end; ++it, ++it0) {
      *it = std::max(*it, *it0);
    }
  }
  CHECK_CUDA_ERROR(cudaMemcpy(s_max.values, work_max.values,
                              s_max.region.len * sizeof(T),
                              cudaMemcpyHostToDevice));
}

template <typename T, size_t D, typename grid_type>
void update_max_grid(const SolverSystemFDCUDA<T, D>* _s, grid_type const& s_max,
                     len_type len) {
  BoundaryGrid<T, D> work_max(s_max.dims);

  CHECK_CUDA_ERROR(cudaMemcpy(work_max.values, s_max.values,
                              s_max.len * sizeof(T), cudaMemcpyDeviceToHost));

  auto s_max_it_start =
      symphas::data_iterator_region(work_max, grid::get_iterable_domain(s_max));
  auto s_max_it_end = s_max_it_start + grid::length_interior(s_max);

  for (auto it(s_max_it_start); it < s_max_it_end; ++it) {
    *it = OpVoid{};
  }
  for (iter_type i = 0; i < len; ++i) {
    BoundaryGrid<T, D> work(_s[i].dims);
    CHECK_CUDA_ERROR(cudaMemcpy(work.values, _s[i].values,
                                _s[i].len * sizeof(T), cudaMemcpyDeviceToHost));
    auto it0 =
        symphas::data_iterator_region(work, grid::get_iterable_domain(s_max));
    for (auto it(s_max_it_start); it < s_max_it_end; ++it, ++it0) {
      *it = std::max(*it, *it0);
    }
  }
  CHECK_CUDA_ERROR(cudaMemcpy(s_max.values, work_max.values,
                              s_max.len * sizeof(T), cudaMemcpyHostToDevice));
}

}  // namespace symphas::internal

#endif
