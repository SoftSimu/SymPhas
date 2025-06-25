
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
 * PURPOSE: Defines a system using a grid with boundaries.
 *
 * ***************************************************************************
 */

#pragma once

#include "boundarysystem.h"

#ifdef USING_CUDA

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "boundarygroup.cuh"
#include "expressionlogic.cuh"
#include "gridfunctions.cuh"

//! Representation of a phase field and associated characteristics.
/*!
 * Specialization based on a phase field system in which boundaries must be
 * managed. See ::PhaseFieldSystem.
 *
 * When the phase field is created, it will automatically populate the phase
 * field values from the given initial conditions.
 *
 * The phase field is based on a PersistenSystemData, meaning it has the ability
 * to persist its data to a file if the **io** module is enabled.
 *
 * \tparam T The type of the phase field.
 * \tparam D The dimension of the phase field.
 */
template <typename T, size_t D>
struct PhaseFieldSystem<BoundaryGridCUDA, T, D>
    : PersistentSystemData<BoundaryGridCUDA<T, D>>, BoundaryGroup<T, D> {
  using parent_type = PersistentSystemData<BoundaryGridCUDA<T, D>>;
  using parent_type::info;

  //! Given the system information, generate the phase field.
  /*!
   * Given the system information, generate the phase field. This includes
   * populating the values of the field from the initial conditions.
   * This also includes setting the boundary conditions of the system.
   *
   * \param tdata The information about the initial conditions of the
   * phase field system.
   * \param vdata The interval data about the system.
   * \param bdata The information about the boundary conditions of the phase
   * field system.
   * \param id A special identifier number. This does not have to be unique,
   * but should be unique between different systems in the same model.
   */
  PhaseFieldSystem(symphas::init_data_type const &tdata,
                   symphas::interval_data_type const &vdata,
                   symphas::b_data_type const &bdata, size_t id = 0);

  //! Update the system to prepare for the next solver iteration.
  /*!
   * Update the system to prepare for the next solver iteration. The values
   * on the boundaries are updated based on the solution index and solution
   * time. The way the boundaries are updated is based on the boundary
   * types, and is implemented higher in the class heirarchy.
   *
   * \param index The solution index.
   * \param time The solution time.
   */
  void update(iter_type index = 0, double time = 0);

  //! Implemented for parallelization routines.
  template <typename... Ts>
  static void synchronize(Ts &&...) {}

 protected:
  PhaseFieldSystem();

  symphas::interval_data_type get_extended_intervals(
      symphas::interval_data_type vdata) {
    if (params::extend_boundary) {
      for (auto &[_, interval] : vdata) {
        interval.set_count(interval.get_count() + 2 * BOUNDARY_DEPTH);
        interval.interval_to_domain();
      }
    }
    return vdata;
  }
};

template <typename T, size_t D>
void PhaseFieldSystem<BoundaryGridCUDA, T, D>::update(iter_type index,
                                                      double time) {
  BoundaryGroup<T, D>::update_boundaries(*this, index, time);
}

//! Representation of a phase field and associated characteristics.
/*!
 * Specialization based on a phase field system in which boundaries must be
 * managed. See ::PhaseFieldSystem.
 *
 * When the phase field is created, it will automatically populate the phase
 * field values from the given initial conditions.
 *
 * The phase field is based on a PersistenSystemData, meaning it has the ability
 * to persist its data to a file if the **io** module is enabled.
 *
 * \tparam T The type of the phase field.
 * \tparam D The dimension of the phase field.
 */
template <typename T, size_t D>
struct PhaseFieldSystem<RegionalGridCUDA, T, D>
    : PersistentSystemData<RegionalGridCUDA<T, D>>, BoundaryGroup<T, D> {
  using parent_type = PersistentSystemData<RegionalGridCUDA<T, D>>;
  using parent_type::info;

  //! Given the system information, generate the phase field.
  /*!
   * Given the system information, generate the phase field. This includes
   * populating the values of the field from the initial conditions.
   * This also includes setting the boundary conditions of the system.
   *
   * \param tdata The information about the initial conditions of the
   * phase field system.
   * \param vdata The interval data about the system.
   * \param bdata The information about the boundary conditions of the phase
   * field system.
   * \param id A special identifier number. This does not have to be unique,
   * but should be unique between different systems in the same model.
   */
  PhaseFieldSystem(symphas::init_data_type const &tdata,
                   symphas::interval_data_type const &vdata,
                   symphas::b_data_type const &bdata, size_t id = 0);

  //! Update the system to prepare for the next solver iteration.
  /*!
   * Update the system to prepare for the next solver iteration. The values
   * on the boundaries are updated based on the solution index and solution
   * time. The way the boundaries are updated is based on the boundary
   * types, and is implemented higher in the class heirarchy.
   *
   * \param index The solution index.
   * \param time The solution time.
   */
  void update(iter_type index = 0, double time = 0);

  //! Implemented for parallelization routines.
  template <typename... Ts>
  static void synchronize(Ts &&...) {}

 protected:
  PhaseFieldSystem();

  symphas::interval_data_type get_extended_intervals(
      symphas::interval_data_type vdata) {
    if (params::extend_boundary) {
      for (auto &[_, interval] : vdata) {
        interval.set_count(interval.get_count() + 2 * BOUNDARY_DEPTH);
        interval.interval_to_domain();
      }
    }
    return vdata;
  }

  regional_system_info_type<T> regional_info;
};

namespace grid {
// CUDA kernel to find the minimum value in a device array
template <typename T>
__global__ void findMinKernel(const T *d_array, T *d_min, int size) {
  extern __shared__ __align__(sizeof(T)) unsigned char sdata_raw[];
  T *sdata = reinterpret_cast<T *>(sdata_raw);

  // Each thread loads one element from global to shared memory
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < size) {
    sdata[tid] = d_array[i];
  } else {
    sdata[tid] = grid::getMaxValue<T>();
  }
  __syncthreads();

  // Do reduction in shared memory
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdata[tid] =
          (!compare_cutoff(sdata[tid], sdata[tid + s]))
              ? sdata[tid]
              : sdata[tid + s];
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    d_min[blockIdx.x] = sdata[0];
  }
}

// CUDA kernel to find the minimum value in a device array
template <typename T>
__global__ void findMinKernelVec(const T *d_array0, const T *d_array1,
                                 T *d_min0, T *d_min1, int size) {
  extern __shared__ __align__(sizeof(T)) unsigned char sdata0_raw[];
  T *sdata0 = reinterpret_cast<T *>(sdata0_raw);
  extern __shared__ __align__(sizeof(T)) unsigned char sdata1_raw[];
  T *sdata1 = reinterpret_cast<T *>(sdata1_raw);

  // Each thread loads one element from global to shared memory
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < size) {
    sdata0[tid] = d_array0[i];
    sdata1[tid] = d_array1[i];
  } else {
    sdata0[tid] = grid::getMaxValue<T>();
    sdata1[tid] = grid::getMaxValue<T>();
  }
  __syncthreads();

  // Do reduction in shared memory
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      T value_curr = sdata0[tid] * sdata0[tid] + sdata1[tid] * sdata1[tid];
      T value_s =
          sdata0[tid + s] * sdata0[tid + s] + sdata1[tid + s] * sdata1[tid + s];
      if (value_curr > value_s) {
        sdata0[tid] = sdata0[tid + s];
        sdata1[tid] = sdata1[tid + s];
      }
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    d_min0[blockIdx.x] = sdata0[0];
    d_min1[blockIdx.x] = sdata1[0];
  }
}

// CUDA kernel to find the minimum value in a device array
template <typename T>
__global__ void findMinKernelVec(const T *d_array0, const T *d_array1,
                                 const T *d_array2, T *d_min0, T *d_min1,
                                 T *d_min2, int size) {
  extern __shared__ __align__(sizeof(T)) unsigned char sdata0_raw[];
  T *sdata0 = reinterpret_cast<T *>(sdata0_raw);
  extern __shared__ __align__(sizeof(T)) unsigned char sdata1_raw[];
  T *sdata1 = reinterpret_cast<T *>(sdata1_raw);
  extern __shared__ __align__(sizeof(T)) unsigned char sdata2_raw[];
  T *sdata2 = reinterpret_cast<T *>(sdata2_raw);

  // Each thread loads one element from global to shared memory
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < size) {
    sdata0[tid] = d_array0[i];
    sdata1[tid] = d_array1[i];
    sdata2[tid] = d_array2[i];
  } else {
    sdata0[tid] = grid::getMaxValue<T>();
    sdata1[tid] = grid::getMaxValue<T>();
    sdata2[tid] = grid::getMaxValue<T>();
  }
  __syncthreads();

  // Do reduction in shared memory
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      T value_curr = sdata0[tid] * sdata0[tid] + sdata1[tid] * sdata1[tid] +
                     sdata2[tid] * sdata2[tid];
      T value_s = sdata0[tid + s] * sdata0[tid + s] +
                  sdata1[tid + s] * sdata1[tid + s] +
                  sdata2[tid + s] * sdata2[tid + s];
      if (value_curr > value_s) {
        sdata0[tid] = sdata0[tid + s];
        sdata1[tid] = sdata1[tid + s];
        sdata2[tid] = sdata2[tid + s];
      }
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    d_min0[blockIdx.x] = sdata0[0];
    d_min1[blockIdx.x] = sdata1[0];
    d_min2[blockIdx.x] = sdata2[0];
  }
}

template <size_t D>
struct find_min_val_vec {
  template <typename T>
  void operator()(const T *d_array, T *d_min, int size) {
    int blocksPerGrid = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    findMinKernel CUDA_KERNEL(blocksPerGrid, BLOCK_SIZE,
                              BLOCK_SIZE * sizeof(T))(d_array, d_min, size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
};

template <>
struct find_min_val_vec<1> {
  template <typename T>
  void operator()(const T *(&d_array)[1], T *(&d_min)[1], int size) {
    int blocksPerGrid = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    findMinKernel CUDA_KERNEL(blocksPerGrid, BLOCK_SIZE,
                              BLOCK_SIZE * sizeof(T))(d_array[0], d_min[0],
                                                      size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
};

template <>
struct find_min_val_vec<2> {
  template <typename T>
  void operator()(T *const (&d_array)[2], T *(&d_min)[2], int size) {
    int blocksPerGrid = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    findMinKernelVec CUDA_KERNEL(blocksPerGrid, BLOCK_SIZE,
                                 2 * BLOCK_SIZE * sizeof(T))(
        d_array[0], d_array[1], d_min[0], d_min[1], size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
};

template <>
struct find_min_val_vec<3> {
  template <typename T>
  void operator()(T *const (&d_array)[3], T *(&d_min)[3], int size) {
    int blocksPerGrid = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    findMinKernelVec CUDA_KERNEL(blocksPerGrid, BLOCK_SIZE,
                                 3 * BLOCK_SIZE * sizeof(T))(
        d_array[0], d_array[1], d_array[2], d_min[0], d_min[1], d_min[2], size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
};

template <typename T, size_t D>
auto min_value(RegionalGridCUDA<T, D> const &grid) {
  int blocksPerGrid = (grid.region.len + BLOCK_SIZE - 1) / BLOCK_SIZE;

  // Allocate memory for block results
  T *d_blockMin;
  CHECK_CUDA_ERROR(cudaMalloc(&d_blockMin, blocksPerGrid * sizeof(T)));

  // Launch kernel
  find_min_val_vec<0>{}(grid.values, d_blockMin, grid.region.len);

  // Copy block results to host
  T *h_blockMin = new T[blocksPerGrid];
  CHECK_CUDA_ERROR(cudaMemcpy(h_blockMin, d_blockMin, blocksPerGrid * sizeof(T),
                              cudaMemcpyDeviceToHost));

  // Find the minimum value among block results
  T minVal = std::numeric_limits<T>::max();
  for (int i = 0; i < blocksPerGrid; ++i) {
    minVal = std::min(minVal, h_blockMin[i]);
  }

  // Clean up
  delete[] h_blockMin;
  CHECK_CUDA_ERROR(cudaFree(d_blockMin));

  return minVal;
}

template <typename T, size_t D>
auto min_value(RegionalGridCUDA<any_vector_t<T, D>, D> const &grid) {
  int blocksPerGrid = (grid.region.len + BLOCK_SIZE - 1) / BLOCK_SIZE;

  // Allocate memory for block results
  T *d_blockMin[D];
  for (iter_type i = 0; i < D; ++i) {
    CHECK_CUDA_ERROR(cudaMalloc(&d_blockMin[i], blocksPerGrid * sizeof(T)));
  }

  find_min_val_vec<D>{}(grid.values, d_blockMin, grid.region.len);

  // Copy block results to host
  T *h_blockMin[D];
  for (iter_type i = 0; i < D; ++i) {
    h_blockMin[i] = new T[blocksPerGrid];
    CHECK_CUDA_ERROR(cudaMemcpy(h_blockMin[i], d_blockMin[i],
                                blocksPerGrid * sizeof(T),
                                cudaMemcpyDeviceToHost));
  }

  // Find the minimum value among block results
  T minVal[D]{};

  T sum = std::numeric_limits<T>::max();
  for (int i = 0; i < blocksPerGrid; ++i) {
    T sum_curr = 0;
    for (iter_type d = 0; d < D; ++d) {
      sum_curr += h_blockMin[i][d] * h_blockMin[i][d];
    }
    if (sum_curr < sum) {
      sum = sum_curr;
      for (iter_type d = 0; d < D; ++d) {
        minVal[d] = h_blockMin[i][d];
      }
    }
  }

  // Clean up
  for (iter_type i = 0; i < D; ++i) {
    delete[] h_blockMin[i];
    CHECK_CUDA_ERROR(cudaFree(d_blockMin[i]));
  }

  any_vector_t<T, D> result;
  for (iter_type i = 0; i < D; ++i) {
    result[i] = minVal[i];
  }
  return result;
}
}  // namespace grid

namespace symphas::internal {

template <typename T, size_t D>
void update_regional_system(regional_system_info_type<T> &regional_info,
                            RegionalGridCUDA<T, D> &grid,
                            symphas::grid_info &info, double time) {
  if (regional_info.next_resize == 0) {
    auto min0 = grid::min_value(grid);
    set_value_for_resize(regional_info.cutoff, min0, REGIONAL_GRID_CUTOFF_EPS);
    set_value_for_resize(grid.empty, min0);

    grid::resize_adjust_region(grid, regional_info.cutoff,
                               REGIONAL_GRID_RESIZE_FACTOR);

    for (auto &[axis, interval] : info) {
      iter_type i = symphas::axis_to_index(axis);
      double offset = grid.region.boundary_size * interval.width();
      interval.set_interval(
          grid.region.origin[i] * interval.width() + offset,
          (grid.region.origin[i] + grid.region.dims[i] - 1) * interval.width() -
              offset);
    }
    ++regional_info.next_resize;
  } else if (regional_info.next_resize > 0) {
    iter_type current_resize =
        static_cast<iter_type>(time / regional_info.resize_delta);
    if (current_resize >= regional_info.next_resize) {
      if (regional_info.fixed_resize) {
        grid::adjust_region(grid, regional_info.cutoff);
      } else {
        grid::resize_adjust_region(grid, regional_info.cutoff,
                                   REGIONAL_GRID_RESIZE_FACTOR,
                                   REGIONAL_GRID_RELATIVE_DIMS_EPS);
      }
      regional_info.next_resize = current_resize + 1;
    }
  }
}
}  // namespace symphas::internal

template <typename T, size_t D>
PhaseFieldSystem<BoundaryGridCUDA, T, D>::PhaseFieldSystem()
    : parent_type{}, BoundaryGroup<T, D>{} {}

template <typename T, size_t D>
PhaseFieldSystem<BoundaryGridCUDA, T, D>::PhaseFieldSystem(
    symphas::init_data_type const &tdata,
    symphas::interval_data_type const &vdata, symphas::b_data_type const &bdata,
    size_t id)
    : parent_type{get_extended_intervals(vdata), id},
      BoundaryGroup<T, D>{info.intervals, bdata} {
  grid::region_interval<D> region(parent_type::dims, BOUNDARY_DEPTH);

  Grid<T, D> init_grid(parent_type::dims);
  symphas::internal::populate_tdata(tdata, init_grid, &info, region, id);
  grid::copy(*this, init_grid);
}

template <typename T, size_t D>
void PhaseFieldSystem<RegionalGridCUDA, T, D>::update(iter_type index,
                                                      double time) {
  // BoundaryGroup<T, D>::update_boundaries(*this, index, time);
  symphas::internal::update_regional_system(regional_info, *this, info, time);
}

template <typename T, size_t D>
PhaseFieldSystem<RegionalGridCUDA, T, D>::PhaseFieldSystem()
    : parent_type{}, BoundaryGroup<T, D>{}, regional_info{} {}

template <typename T, size_t D>
PhaseFieldSystem<RegionalGridCUDA, T, D>::PhaseFieldSystem(
    symphas::init_data_type const &tdata,
    symphas::interval_data_type const &vdata, symphas::b_data_type const &bdata,
    size_t id)
    : parent_type{get_extended_intervals(vdata), id},
      BoundaryGroup<T, D>{info.intervals, bdata},
      regional_info{} {
  grid::region_interval<D> region(RegionalGridCUDA<T, D>::region.dims,
                                  RegionalGridCUDA<T, D>::region.boundary_size);

  RegionalGrid<T, D> init_grid(parent_type::dims, parent_type::empty,
                               parent_type::region.boundary_size);
  symphas::internal::populate_tdata(tdata, init_grid, &info, region, id);
  grid::copy(*this, init_grid);

  auto min0 = grid::min_value(*this);
  symphas::internal::set_value_for_resize(regional_info.cutoff, min0,
                                          REGIONAL_GRID_CUTOFF_EPS);
  symphas::internal::set_value_for_resize(parent_type::empty, min0);

  if (grid::has_subdomain(vdata)) {
    grid::resize_adjust_region(*this, vdata);
    regional_info.next_resize = -1;
  }
}

template <typename T, size_t D>
using RegionalSystemCUDA = PhaseFieldSystem<RegionalGridCUDA, T, D>;

template <typename T, size_t D>
using BoundarySystemCUDA = PhaseFieldSystem<BoundaryGridCUDA, T, D>;

DEFINE_SYMBOL_ID((typename T, size_t D), (BoundaryGridCUDA<T, D>),
                 return data.values)
DEFINE_SYMBOL_ID((typename T, size_t D), (RegionalGridCUDA<T, D>),
                 return data.values)

#endif
