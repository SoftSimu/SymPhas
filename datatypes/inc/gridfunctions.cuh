
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
 * PURPOSE: Defines functions for using grids, such as scaling grid values.
 *
 * ***************************************************************************
 */

#pragma once

#include "gridfunctions.h"

#ifdef USING_CUDA

#include <cuda_runtime.h>

#include <limits>

#include "dataiterator.h"
#include "dft.h"

namespace grid {

template <typename T>
__device__ T getMaxValue();

template <>
inline __device__ float getMaxValue<float>() {
  return FLT_MAX;
}

template <>
inline __device__ double getMaxValue<double>() {
  return DBL_MAX;
}

template <>
inline __device__ int getMaxValue<int>() {
  return INT_MAX;
}

template <typename T>
__device__ T getMinValue();

template <>
inline __device__ float getMinValue<float>() {
  return FLT_MIN;
}

template <>
inline __device__ double getMinValue<double>() {
  return DBL_MIN;
}

template <>
inline __device__ int getMinValue<int>() {
  return INT_MIN;
}
//! Copy the interior values of the grid into an array.
/*!
 * The interior values of the given grid are copied into an array.
 * It is assumed that the array has enough space to store all the interior
 * values. For computing the number of interior points, see
 * grid::length_interior(len_type const*). The grid is 1-dimensional.
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(GridCUDA<T, 1> const& input, T* output) {
  copy_interior_cuda_1d(input.values, output, input.dims[0],
                        input.dims[0] - 2 * BOUNDARY_DEPTH);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 2-dimensional Grid, see
 * grid::copy_interior(Grid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(GridCUDA<T, 2> const& input, T* output) {
  copy_interior_cuda_2d(input.values, output, input.dims[0], input.dims[1],
                        input.dims[0] - 2 * BOUNDARY_DEPTH,
                        input.dims[1] - 2 * BOUNDARY_DEPTH);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 3-dimensional Grid, see
 * grid::copy_interior(Grid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(GridCUDA<T, 3> const& input, T* output) {
  copy_interior_cuda_3d(input.values, output, input.dims[0], input.dims[1],
                        input.dims[2], input.dims[0] - 2 * BOUNDARY_DEPTH,
                        input.dims[1] - 2 * BOUNDARY_DEPTH,
                        input.dims[2] - 2 * BOUNDARY_DEPTH);
}

//! Copy the interior values of the grid into an array.
/*!
 * The interior values of the given grid are copied into an array.
 * It is assumed that the array has enough space to store all the interior
 * values. For computing the number of interior points, see
 * grid::length_interior(len_type const*). The grid is 1-dimensional.
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(GridCUDA<any_vector_t<T, 1>, 1> const& input,
                   any_vector_t<T, 1>* output) {
  copy_interior_cuda_1d(input.values, output, input.dims[0],
                        input.dims[0] - 2 * BOUNDARY_DEPTH);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 2-dimensional Grid, see
 * grid::copy_interior(Grid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(GridCUDA<any_vector_t<T, 2>, 2> const& input,
                   any_vector_t<T, 2>* output) {
  copy_interior_cuda_2d(input.values, output, input.dims[0], input.dims[1],
                        input.dims[0] - 2 * BOUNDARY_DEPTH,
                        input.dims[1] - 2 * BOUNDARY_DEPTH);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 3-dimensional Grid, see
 * grid::copy_interior(Grid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(GridCUDA<any_vector_t<T, 3>, 3> const& input,
                   any_vector_t<T, 3>* output) {
  copy_interior_cuda_3d(input.values, output, input.dims[0], input.dims[1],
                        input.dims[2], input.dims[0] - 2 * BOUNDARY_DEPTH,
                        input.dims[1] - 2 * BOUNDARY_DEPTH,
                        input.dims[2] - 2 * BOUNDARY_DEPTH);
}

//! Copy the interior values of the regional grid into an array.
/*!
 * The interior values of the given grid are copied into an array.
 * It is assumed that the array has enough space to store all the interior
 * values. For computing the number of interior points, see
 * grid::length_interior(len_type const*). The grid is 1-dimensional.
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGridCUDA<T, 1> const& input, T* output) {
  copy_interior_cuda_1d(input.values, output, input.dims[0],
                        input.dims[0] - 2 * BOUNDARY_DEPTH);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 2-dimensional Grid, see
 * grid::copy_interior(RegionalGrid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGridCUDA<T, 2> const& input, T* output) {
  copy_interior_cuda_2d(input.values, output, input.dims[0], input.dims[1],
                        input.dims[0] - 2 * BOUNDARY_DEPTH,
                        input.dims[1] - 2 * BOUNDARY_DEPTH);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 3-dimensional Grid, see
 * grid::copy_interior(RegionalGrid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGridCUDA<T, 3> const& input, T* output) {
  copy_interior_cuda_3d(input.values, output, input.dims[0], input.dims[1],
                        input.dims[2], input.dims[0] - 2 * BOUNDARY_DEPTH,
                        input.dims[1] - 2 * BOUNDARY_DEPTH,
                        input.dims[2] - 2 * BOUNDARY_DEPTH);
}

//! Copy the interior values of the grid into an array.
/*!
 * The interior values of the given grid are copied into an array.
 * It is assumed that the array has enough space to store all the interior
 * values. For computing the number of interior points, see
 * grid::length_interior(len_type const*). The grid is 1-dimensional.
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGridCUDA<any_vector_t<T, 1>, 1> const& input,
                   any_vector_t<T, 1>* output) {
  copy_interior_cuda_1d(input.values, output, input.dims[0],
                        input.dims[0] - 2 * BOUNDARY_DEPTH);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 2-dimensional Grid, see
 * grid::copy_interior(RegionalGrid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGridCUDA<any_vector_t<T, 2>, 2> const& input,
                   any_vector_t<T, 2>* output) {
  copy_interior_cuda_2d(input.values, output, input.dims[0], input.dims[1],
                        input.dims[0] - 2 * BOUNDARY_DEPTH,
                        input.dims[1] - 2 * BOUNDARY_DEPTH);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 3-dimensional Grid, see
 * grid::copy_interior(RegionalGrid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGridCUDA<any_vector_t<T, 3>, 3> const& input,
                   any_vector_t<T, 3>* output) {
  copy_interior_cuda_3d(input.values, output, input.dims[0], input.dims[1],
                        input.dims[2], input.dims[0] - 2 * BOUNDARY_DEPTH,
                        input.dims[1] - 2 * BOUNDARY_DEPTH,
                        input.dims[2] - 2 * BOUNDARY_DEPTH);
}

}  // namespace grid

// 1D Kernel
template <typename T>
__global__ void kernelCopyRegion1D(const T* smallGrid, T* largeGrid,
                                   iter_type start, iter_type end,
                                   iter_type largeLength) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  len_type length = end - start;
  if (idx < length) {
    smallGrid[idx] = largeGrid[start + idx];
  }
}

// 2D Kernel
template <typename T>
__global__ void kernelCopyRegion2D(T* smallGrid, const T* largeGrid,
                                   iter_type startX, iter_type endX,
                                   iter_type startY, iter_type endY,
                                   len_type largeWidth, len_type largeHeight) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  len_type width = endX - startX;
  len_type height = endY - startY;
  if (x < width && y < height) {
    smallGrid[y * width + x] =
        largeGrid[(startY + y) * largeWidth + (startX + x)];
  }
}

// 3D Kernel
template <typename T>
__global__ void kernelCopyRegion3D(T* smallGrid, const T* largeGrid,
                                   iter_type startX, iter_type endX,
                                   iter_type startY, iter_type endY,
                                   iter_type startZ, iter_type endZ,
                                   len_type largeWidth, len_type largeHeight,
                                   len_type largeDepth) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;
  iter_type width = endX - startX;
  iter_type height = endY - startY;
  iter_type depth = endZ - startZ;
  if (x < width && y < height && z < depth) {
    smallGrid[(z * height + y) * width + x] =
        largeGrid[((startZ + z) * largeHeight + (startY + y)) * largeWidth +
                  (startX + x)];
  }
}

// 1D Kernel
template <typename T>
__global__ void kernelCopyRegionVec1D(T* smallGrid0, const T* largeGrid0,
                                      iter_type start, iter_type end,
                                      iter_type largeLength) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  len_type length = end - start;
  if (idx < length) {
    smallGrid0[idx] = largeGrid0[start + idx];
  }
}

// 2D Kernel
template <typename T>
__global__ void kernelCopyRegionVec2D(T* smallGrid0, T* smallGrid1,
                                      const T* largeGrid0, const T* largeGrid1,
                                      iter_type startX, iter_type endX,
                                      iter_type startY, iter_type endY,
                                      iter_type largeWidth) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  len_type width = endX - startX;
  len_type height = endY - startY;
  if (x < width && y < height) {
    iter_type largeGridIndex = (startY + y) * largeWidth + (startX + x);
    iter_type smallGridIndex = y * width + x;
    smallGrid0[smallGridIndex] = largeGrid0[largeGridIndex];
    smallGrid1[smallGridIndex] = largeGrid1[largeGridIndex];
  }
}

// 3D Kernel
template <typename T>
__global__ void kernelCopyRegionVec3D(
    T* smallGrid0, T* smallGrid1, T* smallGrid2, const T* largeGrid0,
    const T* largeGrid1, const T* largeGrid2, iter_type startX, iter_type endX,
    iter_type startY, iter_type endY, iter_type startZ, iter_type endZ,
    iter_type largeWidth, iter_type largeHeight) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;
  iter_type width = endX - startX;
  iter_type height = endY - startY;
  iter_type depth = endZ - startZ;
  if (x < width && y < height && z < depth) {
    iter_type largeGridIndex =
        ((startZ + z) * largeHeight + (startY + y)) * largeWidth + (startX + x);
    iter_type smallGridIndex = (z * height + y) * width + x;
    smallGrid0[smallGridIndex] = largeGrid0[largeGridIndex];
    smallGrid1[smallGridIndex] = largeGrid1[largeGridIndex];
    smallGrid2[smallGridIndex] = largeGrid2[largeGridIndex];
  }
}

// 1D Kernel
template <typename T>
__global__ void kernelCopyRegionVec1D(any_vector_t<T, 1>* smallGrid,
                                      const T* largeGrid0, iter_type start,
                                      iter_type end, iter_type largeLength) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  len_type length = end - start;
  if (idx < length) {
    smallGrid[idx][0] = largeGrid0[start + idx];
  }
}

// 2D Kernel
template <typename T>
__global__ void kernelCopyRegionVec2D(any_vector_t<T, 2>* smallGrid,
                                      const T* largeGrid0, const T* largeGrid1,
                                      iter_type startX, iter_type endX,
                                      iter_type startY, iter_type endY,
                                      iter_type largeWidth) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  len_type width = endX - startX;
  len_type height = endY - startY;
  if (x < width && y < height) {
    iter_type largeGridIndex = (startY + y) * largeWidth + (startX + x);
    iter_type smallGridIndex = y * width + x;
    smallGrid[smallGridIndex][0] = largeGrid0[largeGridIndex];
    smallGrid[smallGridIndex][1] = largeGrid1[largeGridIndex];
  }
}

// 3D Kernel
template <typename T>
__global__ void kernelCopyRegionVec3D(any_vector_t<T, 3>* smallGrid,
                                      const T* largeGrid0, const T* largeGrid1,
                                      const T* largeGrid2, iter_type startX,
                                      iter_type endX, iter_type startY,
                                      iter_type endY, iter_type startZ,
                                      iter_type endZ, iter_type largeWidth,
                                      iter_type largeHeight) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;
  iter_type width = endX - startX;
  iter_type height = endY - startY;
  iter_type depth = endZ - startZ;
  if (x < width && y < height && z < depth) {
    iter_type largeGridIndex =
        ((startZ + z) * largeHeight + (startY + y)) * largeWidth + (startX + x);
    iter_type smallGridIndex = (z * height + y) * width + x;
    smallGrid[smallGridIndex][0] = largeGrid0[largeGridIndex];
    smallGrid[smallGridIndex][1] = largeGrid1[largeGridIndex];
    smallGrid[smallGridIndex][2] = largeGrid2[largeGridIndex];
  }
}

template <typename T, size_t D>
struct initiate_region_copy;

template <typename T>
struct initiate_region_copy<T, 1> {
  void operator()(grid::region_interval<1> const& interval,
                  const T* inputDevice, T* outputDevice) {
    len_type len = grid::length<1>(interval.dims);
    int numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelCopyRegion1D CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        inputDevice, outputDevice, interval.intervals[0][0],
        interval.intervals[0][1], interval.dims[0]);
  }
};

template <typename T>
struct initiate_region_copy<T, 2> {
  void operator()(grid::region_interval<2> const& interval,
                  const T* inputDevice, T* outputDevice) {
    dim3 blockDim(32, 32);
    dim3 gridDim((interval.dims[0] + blockDim.x - 1) / blockDim.x,
                 (interval.dims[1] + blockDim.y - 1) / blockDim.y);

    kernelCopyRegion2D CUDA_KERNEL(gridDim, blockDim)(
        outputDevice, inputDevice, interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.dims[0], interval.dims[1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
};

template <typename T>
struct initiate_region_copy<T, 3> {
  void operator()(grid::region_interval<3> const& interval,
                  const T* inputDevice, T* outputDevice) {
    dim3 blockDim(8, 8, 8);
    dim3 gridDim((interval.dims[0] + blockDim.x - 1) / blockDim.x,
                 (interval.dims[1] + blockDim.y - 1) / blockDim.y,
                 (interval.dims[2] + blockDim.z - 1) / blockDim.z);

    kernelCopyRegion3D CUDA_KERNEL(gridDim, blockDim)(
        outputDevice, inputDevice, interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1], interval.dims[0], interval.dims[1],
        interval.dims[2]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
};

template <typename T>
struct initiate_region_copy<any_vector_t<T, 1>, 1> {
  void operator()(grid::region_interval<1> const& interval,
                  T* const (&inputDevice)[1], T* (&outputDevice)[1]) {
    len_type len = grid::length<1>(interval.dims);
    int numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelCopyRegionVec1D CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        outputDevice[0], inputDevice[0], interval.intervals[0][0],
        interval.intervals[0][1], interval.dims[0]);
  }
  void operator()(grid::region_interval<1> const& interval,
                  T* const (&inputDevice)[1],
                  any_vector_t<T, 1>* outputDevice) {
    len_type len = grid::length<1>(interval.dims);
    int numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    kernelCopyRegionVec1D CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        outputDevice, inputDevice[0], interval.intervals[0][0],
        interval.intervals[0][1], interval.dims[0]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
};

template <typename T>
struct initiate_region_copy<any_vector_t<T, 2>, 2> {
  void operator()(grid::region_interval<2> const& interval,
                  T* const (&inputDevice)[2], T* (&outputDevice)[2]) {
    dim3 blockDim(32, 32);
    dim3 gridDim((interval.dims[0] + blockDim.x - 1) / blockDim.x,
                 (interval.dims[1] + blockDim.y - 1) / blockDim.y);
    kernelCopyRegionVec2D CUDA_KERNEL(gridDim, blockDim)(
        outputDevice[0], outputDevice[1], inputDevice[0], inputDevice[1],
        interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1], interval.dims[0]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
  void operator()(grid::region_interval<2> const& interval,
                  T* const (&inputDevice)[2],
                  any_vector_t<T, 2>* outputDevice) {
    dim3 blockDim(32, 32);
    dim3 gridDim((interval.dims[0] + blockDim.x - 1) / blockDim.x,
                 (interval.dims[1] + blockDim.y - 1) / blockDim.y);
    kernelCopyRegionVec2D CUDA_KERNEL(gridDim, blockDim)(
        outputDevice, inputDevice[0], inputDevice[1], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.dims[0]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
};

template <typename T>
struct initiate_region_copy<any_vector_t<T, 3>, 3> {
  void operator()(grid::region_interval<3> const& interval,
                  T* const (&inputDevice)[3], T* (&outputDevice)[3]) {
    dim3 blockDim(8, 8, 8);
    dim3 gridDim((interval.dims[0] + blockDim.x - 1) / blockDim.x,
                 (interval.dims[1] + blockDim.y - 1) / blockDim.y,
                 (interval.dims[2] + blockDim.z - 1) / blockDim.z);
    kernelCopyRegionVec3D CUDA_KERNEL(gridDim, blockDim)(
        outputDevice[0], outputDevice[1], outputDevice[2], inputDevice[0],
        inputDevice[1], inputDevice[2], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1], interval.dims[0], interval.dims[1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
  void operator()(grid::region_interval<3> const& interval,
                  T* const (&inputDevice)[3],
                  any_vector_t<T, 3>* outputDevice) {
    dim3 blockDim(8, 8, 8);
    dim3 gridDim((interval.dims[0] + blockDim.x - 1) / blockDim.x,
                 (interval.dims[1] + blockDim.y - 1) / blockDim.y,
                 (interval.dims[2] + blockDim.z - 1) / blockDim.z);
    kernelCopyRegionVec3D CUDA_KERNEL(gridDim, blockDim)(
        outputDevice, inputDevice[0], inputDevice[1], inputDevice[2],
        interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1],
        interval.intervals[2][0], interval.intervals[2][1], interval.dims[0],
        interval.dims[1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
};

namespace grid {
template <typename T, size_t D>
void copy_region(RegionalGridCUDA<T, D> const& input, T* output) {
  auto interval = get_iterable_domain(input);
  T* outputDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&outputDevice, grid::length<D>(interval) * sizeof(T)));
  initiate_region_copy<T, D>{}(interval, input.values, outputDevice);
  CHECK_CUDA_ERROR(cudaMemcpy(output, outputDevice,
                              grid::length<D>(interval) * sizeof(T),
                              cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(outputDevice));
}

template <typename T, size_t D>
void copy_region(RegionalGridCUDA<any_vector_t<T, D>, D> const& input,
                 T* (&output)[D]) {
  auto interval = get_iterable_domain(input);
  T* outputDevice[D];
  for (iter_type i = 0; i < D; ++i) {
    CHECK_CUDA_ERROR(
        cudaMalloc(&outputDevice[i], grid::length<D>(input.dims) * sizeof(T)));
  }
  initiate_region_copy<any_vector_t<T, D>, D>{}(interval, input.values,
                                                outputDevice);
  for (iter_type i = 0; i < D; ++i) {
    CHECK_CUDA_ERROR(cudaMemcpy(output[i], outputDevice[i],
                                grid::length<D>(interval) * sizeof(T),
                                cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaFree(outputDevice[i]));
  }
}

template <typename T, size_t D>
void copy_region(RegionalGridCUDA<any_vector_t<T, D>, D> const& input,
                 any_vector_t<T, D>* output) {
  auto interval = get_iterable_domain(input);
  any_vector_t<T, D>* outputDevice;
  CHECK_CUDA_ERROR(cudaMalloc(
      &outputDevice, grid::length<D>(input.dims) * sizeof(any_vector_t<T, D>)));
  initiate_region_copy<any_vector_t<T, D>, D>{}(interval, input.values,
                                                outputDevice);
  CHECK_CUDA_ERROR(cudaMemcpy(output, outputDevice,
                              grid::length<D>(interval) * sizeof(T) * D,
                              cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaFree(outputDevice));
}

}  // namespace grid

template <size_t D>
__global__ void kernelAdjustRegionToFrom(scalar_t* newValuesDevice,
                                         const scalar_t* oldValuesDevice,
                                         const grid::RegionAdjustParams<D>* rp,
                                         len_type boundarySize) {
  for (iter_type n = 0; n < grid::length<D>(rp->new_interior_dims); ++n) {
    iter_type pos[D];

    grid::get_grid_position(pos, rp->new_interior_dims, n);
    iter_type new_index =
        grid::index_from_position(pos, rp->new_stride, boundarySize);
    rp->compute_old_position(pos, n);

    bool is_in_old = grid::is_in_region(pos, rp->old_interior_dims);
    if (is_in_old) {
      iter_type old_index =
          grid::index_from_position(pos, rp->old_stride, boundarySize);
      newValuesDevice[new_index] = oldValuesDevice[old_index];
    }
  }
}

// Kernel to find the minimal region with values greater than the cutoff
template <typename T, size_t D>
__global__ void findMinimalRegion(const T* values, T cutoff_value,
                                  grid::MinimalRegionParams<D>* params,
                                  iter_type* min_indices,
                                  iter_type* max_indices, len_type total_size) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= total_size) return;
  /*compare_cutoff(grid.values, index, cutoff_value) &&
      !grid::is_in_region(pos, intervals) */
  if (grid::compare_cutoff(values, idx, cutoff_value)) {
    for (size_t d = 0; d < D; ++d) {
      iter_type coord = (idx / params->stride[d]) % params->dims[d];
      atomicMin(&min_indices[d], coord);
      atomicMax(&max_indices[d], coord);
    }
  }
}

// Kernel to find the minimal region with values greater than the cutoff
template <typename T>
__global__ void findMinimalRegionVec1D(const T* values0,
                                       any_vector_t<T, 1> cutoff_value,
                                       grid::MinimalRegionParams<1>* params,
                                       iter_type* min_indices,
                                       iter_type* max_indices,
                                       len_type total_size) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= total_size) return;
  /*compare_cutoff(grid.values, index, cutoff_value) &&
      !grid::is_in_region(pos, intervals) */
  if (grid::compare_cutoff(values0, idx, cutoff_value)) {
    iter_type coord = (idx / params->stride[0]) % params->dims[0];
    atomicMin(&min_indices[0], coord);
    atomicMax(&max_indices[0], coord);
  }
}

// Kernel to find the minimal region with values greater than the cutoff
template <typename T>
__global__ void findMinimalRegionVec2D(const T* values0, const T* values1,
                                       any_vector_t<T, 2> cutoff_value,
                                       grid::MinimalRegionParams<2>* params,
                                       iter_type* min_indices,
                                       iter_type* max_indices,
                                       len_type total_size) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= total_size) return;
  /*compare_cutoff(grid.values, index, cutoff_value) &&
      !grid::is_in_region(pos, intervals) */
  T value[]{values0[idx], values1[idx]};
  if (grid::compare_cutoff(value, cutoff_value)) {
    for (size_t d = 0; d < 2; ++d) {
      iter_type coord = (idx / params->stride[d]) % params->dims[d];
      atomicMin(&min_indices[d], coord);
      atomicMax(&max_indices[d], coord);
    }
  }
}

// Kernel to find the minimal region with values greater than the cutoff
template <typename T>
__global__ void findMinimalRegionVec3D(
    const T* values0, const T* values1, const T* values2,
    any_vector_t<T, 3> cutoff_value, grid::MinimalRegionParams<3>* params,
    iter_type* min_indices, iter_type* max_indices, len_type total_size) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= total_size) return;
  /*compare_cutoff(grid.values, index, cutoff_value) &&
      !grid::is_in_region(pos, intervals) */
  T value[]{values0[idx], values1[idx], values2[idx]};
  if (grid::compare_cutoff(value, cutoff_value)) {
    for (size_t d = 0; d < 3; ++d) {
      iter_type coord = (idx / params->stride[d]) % params->dims[d];
      atomicMin(&min_indices[d], coord);
      atomicMax(&max_indices[d], coord);
    }
  }
}

template <size_t D>
struct run_find_minimal_region;

template <>
struct run_find_minimal_region<1> {
  template <typename T>
  void operator()(T* const (&values)[1], any_vector_t<T, 1> const& cutoff_value,
                  grid::MinimalRegionParams<1>* params, iter_type* min_indices,
                  iter_type* max_indices, len_type total_size) {
    int numBlocks = (total_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    findMinimalRegionVec1D CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        values[0], cutoff_value, params, min_indices, max_indices, total_size);
  }

  template <typename T>
  void operator()(const T* values, T cutoff_value,
                  grid::MinimalRegionParams<1>* params, iter_type* min_indices,
                  iter_type* max_indices, len_type total_size) {
    int numBlocks = (total_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    findMinimalRegion CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        values, cutoff_value, params, min_indices, max_indices, total_size);
  }
};

template <>
struct run_find_minimal_region<2> {
  template <typename T>
  void operator()(T* const (&values)[2], any_vector_t<T, 2> const& cutoff_value,
                  grid::MinimalRegionParams<2>* params, iter_type* min_indices,
                  iter_type* max_indices, len_type total_size) {
    int numBlocks = (total_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    findMinimalRegionVec2D CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        values[0], values[1], cutoff_value, params, min_indices, max_indices,
        total_size);
  }

  template <typename T>
  void operator()(const T* values, T cutoff_value,
                  grid::MinimalRegionParams<2>* params, iter_type* min_indices,
                  iter_type* max_indices, len_type total_size) {
    int numBlocks = (total_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    findMinimalRegion CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        values, cutoff_value, params, min_indices, max_indices, total_size);
  }
};

template <>
struct run_find_minimal_region<3> {
  template <typename T>
  void operator()(T* const (&values)[3], any_vector_t<T, 3> const& cutoff_value,
                  grid::MinimalRegionParams<3>* params, iter_type* min_indices,
                  iter_type* max_indices, len_type total_size) {
    int numBlocks = (total_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    findMinimalRegionVec3D CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        values[0], values[1], values[2], cutoff_value, params, min_indices,
        max_indices, total_size);
  }

  template <typename T>
  void operator()(const T* values, T cutoff_value,
                  grid::MinimalRegionParams<3>* params, iter_type* min_indices,
                  iter_type* max_indices, len_type total_size) {
    int numBlocks = (total_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    findMinimalRegion CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        values, cutoff_value, params, min_indices, max_indices, total_size);
  }
};
template <typename T>
__global__ void computeEnclosingInterval1d(
    const T* values, iter_type* intervals_0, iter_type* intervals_1,
    T cutoff_value, iter_type total_length, iter_type start, len_type dimx) {
  extern __shared__ iter_type shared_intervals[];

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int tidx = threadIdx.x;

  // Initialize shared memory
  shared_intervals[0] = grid::getMaxValue<iter_type>();
  shared_intervals[1] = grid::getMinValue<iter_type>();
  __syncthreads();

  if (idx < dimx) {
    iter_type idxn = start + idx;
    iter_type n = idxn % dimx;
    if (values[n] > cutoff_value) {
      atomicMin(&shared_intervals[0], idxn);
      atomicMax(&shared_intervals[1], idxn);
    }
  }
  __syncthreads();

  if (tidx == 0) {
    atomicMin(&intervals_0[0], shared_intervals[0]);
    atomicMax(&intervals_1[0], shared_intervals[1]);
  }
}

template <typename T>
__global__ void computeEnclosingInterval2d(
    const T* values, iter_type* intervals_0, iter_type* intervals_1,
    T cutoff_value, iter_type total_length, iter_type posx, iter_type posy,
    len_type dimx, len_type dimy, iter_type startx, iter_type starty,
    len_type stride) {
  extern __shared__ iter_type shared_intervals[];
  iter_type* shared_intervals_x = &shared_intervals[0];
  iter_type* shared_intervals_y = &shared_intervals_x[blockDim.x * 2];

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y + threadIdx.y;
  int tidx = threadIdx.x;
  int tidy = threadIdx.y;
  // Initialize shared memory
  if (tidy == 0) {
    shared_intervals_x[tidx * 2] = grid::getMaxValue<iter_type>();
    shared_intervals_x[tidx * 2 + 1] = grid::getMinValue<iter_type>();
  }
  if (tidx == 0) {
    shared_intervals_y[tidy * 2] = grid::getMaxValue<iter_type>();
    shared_intervals_y[tidy * 2 + 1] = grid::getMinValue<iter_type>();
  }
  __syncthreads();

  if (idx < dimx && idy < dimy) {
    iter_type idxn = posx + idx;
    iter_type idyn = posy + idy;

    iter_type xn = ((idxn >= dimx) ? idxn - dimx : idxn) + startx;
    iter_type yn = ((idyn >= dimy) ? idyn - dimy : idyn) + starty;
    iter_type n = xn + yn * stride;
    if (values[n] > cutoff_value) {
      atomicMin(&shared_intervals_x[tidx * 2], idxn);
      atomicMax(&shared_intervals_x[tidx * 2 + 1], idxn);
      atomicMin(&shared_intervals_y[tidy * 2], idyn);
      atomicMax(&shared_intervals_y[tidy * 2 + 1], idyn);
    }
  }
  __syncthreads();

  iter_type* start_interval_x = intervals_0;
  iter_type* start_interval_y = intervals_0 + dimx;
  iter_type* end_interval_x = intervals_1;
  iter_type* end_interval_y = intervals_1 + dimx;

  // Write local minimums to global memory
  if (idx < dimx && idy < dimy) {
    if (tidy == 0) {
      atomicMin(&start_interval_x[idx], shared_intervals_x[tidx * 2]);
      atomicMax(&end_interval_x[idx], shared_intervals_x[tidx * 2 + 1]);
    }
    if (tidx == 0) {
      atomicMin(&start_interval_y[idy], shared_intervals_y[tidy * 2]);
      atomicMax(&end_interval_y[idy], shared_intervals_y[tidy * 2 + 1]);
    }
  }
}

template <typename T>
__global__ void computeEnclosingInterval2d(
    const T* values0, const T* values1, iter_type* intervals_0,
    iter_type* intervals_1, any_vector_t<T, 2> cutoff_value,
    iter_type total_length, iter_type posx, iter_type posy, len_type dimx,
    len_type dimy, iter_type startx, iter_type starty, len_type stride) {
  extern __shared__ iter_type shared_intervals[];
  iter_type* shared_intervals_x = &shared_intervals[0];
  iter_type* shared_intervals_y = &shared_intervals_x[blockDim.x * 2];

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y + threadIdx.y;
  int tidx = threadIdx.x;
  int tidy = threadIdx.y;
  // Initialize shared memory
  if (tidy == 0) {
    shared_intervals_x[tidx * 2] = grid::getMaxValue<iter_type>();
    shared_intervals_x[tidx * 2 + 1] = grid::getMinValue<iter_type>();
  }
  if (tidx == 0) {
    shared_intervals_y[tidy * 2] = grid::getMaxValue<iter_type>();
    shared_intervals_y[tidy * 2 + 1] = grid::getMinValue<iter_type>();
  }
  __syncthreads();

  if (idx < dimx && idy < dimy) {
    iter_type idxn = posx + idx;
    iter_type idyn = posy + idy;

    iter_type xn = ((idxn >= dimx) ? idxn - dimx : idxn) + startx;
    iter_type yn = ((idyn >= dimy) ? idyn - dimy : idyn) + starty;
    iter_type n = xn + yn * stride;
    T value[]{values0[n], values1[n]};
    if (grid::compare_cutoff(value, cutoff_value)) {
      atomicMin(&shared_intervals_x[tidx * 2], idxn);
      atomicMax(&shared_intervals_x[tidx * 2 + 1], idxn);
      atomicMin(&shared_intervals_y[tidy * 2], idyn);
      atomicMax(&shared_intervals_y[tidy * 2 + 1], idyn);
    }
  }
  __syncthreads();

  iter_type* start_interval_x = intervals_0;
  iter_type* start_interval_y = intervals_0 + dimx;
  iter_type* end_interval_x = intervals_1;
  iter_type* end_interval_y = intervals_1 + dimx;

  // Write local minimums to global memory
  if (idx < dimx && idy < dimy) {
    if (tidy == 0) {
      atomicMin(&start_interval_x[idx], shared_intervals_x[tidx * 2]);
      atomicMax(&end_interval_x[idx], shared_intervals_x[tidx * 2 + 1]);
    }
    if (tidx == 0) {
      atomicMin(&start_interval_y[idy], shared_intervals_y[tidy * 2]);
      atomicMax(&end_interval_y[idy], shared_intervals_y[tidy * 2 + 1]);
    }
  }
}

template <typename T>
__global__ void computeEnclosingInterval3d(
    const T* values, iter_type* intervals_0, iter_type* intervals_1,
    T cutoff_value, iter_type total_length, iter_type posx, iter_type posy,
    iter_type posz, len_type dimx, len_type dimy, len_type dimz,
    iter_type startx, iter_type starty, iter_type startz, len_type stridey,
    len_type stridez) {
  extern __shared__ iter_type shared_intervals[];
  iter_type* shared_intervals_x = &shared_intervals[0];
  iter_type* shared_intervals_y = &shared_intervals_x[0] + blockDim.x * 2;
  iter_type* shared_intervals_z = &shared_intervals_y[0] + blockDim.y * 2;

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y + threadIdx.y;
  int idz = blockIdx.z * blockDim.z + threadIdx.z;
  int tidx = threadIdx.x;
  int tidy = threadIdx.y;
  int tidz = threadIdx.z;

  // Initialize shared memory
  shared_intervals_x[tidx * 2] = grid::getMaxValue<iter_type>();
  shared_intervals_x[tidx * 2 + 1] = grid::getMinValue<iter_type>();
  shared_intervals_y[tidy * 2] = grid::getMaxValue<iter_type>();
  shared_intervals_y[tidy * 2 + 1] = grid::getMinValue<iter_type>();
  shared_intervals_z[tidy * 2] = grid::getMaxValue<iter_type>();
  shared_intervals_z[tidy * 2 + 1] = grid::getMinValue<iter_type>();
  __syncthreads();

  if (idx < dimx && idy < dimy) {
    iter_type idxn = posx + idx;
    iter_type idyn = posy + idy;
    iter_type idzn = posz + idz;

    iter_type xn = ((idxn > dimx) ? idxn - dimx : idxn) + startx;
    iter_type yn = ((idyn > dimy) ? idyn - dimy : idyn) + starty;
    iter_type zn = ((idzn > dimz) ? idzn - dimz : idzn) + startz;
    iter_type n = xn + yn * stridey + zn * stridez;
    if (values[n] > cutoff_value) {
      atomicMin(&shared_intervals_x[tidx * 2], idxn);
      atomicMax(&shared_intervals_x[tidx * 2 + 1], idxn);
      atomicMin(&shared_intervals_y[tidy * 2], idyn);
      atomicMax(&shared_intervals_y[tidy * 2 + 1], idyn);
      atomicMin(&shared_intervals_z[tidz * 2], idzn);
      atomicMax(&shared_intervals_z[tidz * 2 + 1], idzn);
    }
  }
  __syncthreads();

  iter_type* start_interval_x = &intervals_0[0];
  iter_type* start_interval_y = &start_interval_x[0] + dimx;
  iter_type* start_interval_z = &start_interval_y[0] + dimy;
  iter_type* end_interval_x = &intervals_1[0];
  iter_type* end_interval_y = &end_interval_x[0] + dimx;
  iter_type* end_interval_z = &end_interval_y[0] + dimy;

  // Write local minimums to global memory
  if (tidy == 0 && tidz == 0) {
    atomicMin(&start_interval_x[idx], shared_intervals_x[tidx * 2]);
    atomicMax(&end_interval_x[idx], shared_intervals_x[tidx * 2 + 1]);
  }
  if (tidx == 0 && tidz == 0) {
    atomicMin(&start_interval_y[idy], shared_intervals_y[tidy * 2]);
    atomicMax(&end_interval_y[idy], shared_intervals_y[tidy * 2 + 1]);
  }
  if (tidx == 0 && tidy == 0) {
    atomicMin(&start_interval_z[idz], shared_intervals_z[tidz * 2]);
    atomicMax(&end_interval_z[idz], shared_intervals_z[tidz * 2 + 1]);
  }
}

template <typename T>
__global__ void computeEnclosingInterval3d(
    const T* values0, const T* values1, const T* values2,
    iter_type* intervals_0, iter_type* intervals_1,
    any_vector_t<T, 3> cutoff_value, iter_type total_length, iter_type posx,
    iter_type posy, iter_type posz, len_type dimx, len_type dimy, len_type dimz,
    iter_type startx, iter_type starty, iter_type startz, len_type stridey,
    len_type stridez) {
  extern __shared__ iter_type shared_intervals[];
  iter_type* shared_intervals_x = &shared_intervals[0];
  iter_type* shared_intervals_y = &shared_intervals_x[0] + blockDim.x * 2;
  iter_type* shared_intervals_z = &shared_intervals_y[0] + blockDim.y * 2;

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y + threadIdx.y;
  int idz = blockIdx.z * blockDim.z + threadIdx.z;
  int tidx = threadIdx.x;
  int tidy = threadIdx.y;
  int tidz = threadIdx.z;

  // Initialize shared memory
  shared_intervals_x[tidx * 2] = grid::getMaxValue<iter_type>();
  shared_intervals_x[tidx * 2 + 1] = grid::getMinValue<iter_type>();
  shared_intervals_y[tidy * 2] = grid::getMaxValue<iter_type>();
  shared_intervals_y[tidy * 2 + 1] = grid::getMinValue<iter_type>();
  shared_intervals_z[tidy * 2] = grid::getMaxValue<iter_type>();
  shared_intervals_z[tidy * 2 + 1] = grid::getMinValue<iter_type>();
  __syncthreads();

  if (idx < dimx && idy < dimy) {
    iter_type idxn = posx + idx;
    iter_type idyn = posy + idy;
    iter_type idzn = posz + idz;

    iter_type xn = ((idxn > dimx) ? idxn - dimx : idxn) + startx;
    iter_type yn = ((idyn > dimy) ? idyn - dimy : idyn) + starty;
    iter_type zn = ((idzn > dimz) ? idzn - dimz : idzn) + startz;
    iter_type n = xn + yn * stridey + zn * stridez;
    T value[]{values0[n], values1[n], values2[n]};
    if (grid::compare_cutoff(value, cutoff_value)) {
      atomicMin(&shared_intervals_x[tidx * 2], idxn);
      atomicMax(&shared_intervals_x[tidx * 2 + 1], idxn);
      atomicMin(&shared_intervals_y[tidy * 2], idyn);
      atomicMax(&shared_intervals_y[tidy * 2 + 1], idyn);
      atomicMin(&shared_intervals_z[tidz * 2], idzn);
      atomicMax(&shared_intervals_z[tidz * 2 + 1], idzn);
    }
  }
  __syncthreads();

  iter_type* start_interval_x = &intervals_0[0];
  iter_type* start_interval_y = &start_interval_x[0] + dimx;
  iter_type* start_interval_z = &start_interval_y[0] + dimy;
  iter_type* end_interval_x = &intervals_1[0];
  iter_type* end_interval_y = &end_interval_x[0] + dimx;
  iter_type* end_interval_z = &end_interval_y[0] + dimy;

  // Write local minimums to global memory
  if (tidy == 0 && tidz == 0) {
    atomicMin(&start_interval_x[idx], shared_intervals_x[tidx * 2]);
    atomicMax(&end_interval_x[idx], shared_intervals_x[tidx * 2 + 1]);
  }
  if (tidx == 0 && tidz == 0) {
    atomicMin(&start_interval_y[idy], shared_intervals_y[tidy * 2]);
    atomicMax(&end_interval_y[idy], shared_intervals_y[tidy * 2 + 1]);
  }
  if (tidx == 0 && tidy == 0) {
    atomicMin(&start_interval_z[idz], shared_intervals_z[tidz * 2]);
    atomicMax(&end_interval_z[idz], shared_intervals_z[tidz * 2 + 1]);
  }
}

template <typename T>
__global__ void reduceIntervals2d(iter_type* intervals_0,
                                  iter_type* intervals_1, len_type dimx,
                                  len_type dimy) {
  extern __shared__ iter_type sdata[];

  iter_type* shared_intervals_x_start = &sdata[0];
  iter_type* shared_intervals_y_start = &shared_intervals_x_start[blockDim.x];
  iter_type* shared_intervals_x_end = &shared_intervals_y_start[blockDim.x];
  iter_type* shared_intervals_y_end = &shared_intervals_x_end[blockDim.x];

  iter_type* intervals_x_start = &intervals_0[0];
  iter_type* intervals_y_start = &intervals_0[dimx];
  iter_type* intervals_x_end = &intervals_1[0];
  iter_type* intervals_y_end = &intervals_1[dimx];

  int tid = threadIdx.x;
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  // Load data into shared memory
  shared_intervals_x_start[tid] =
      (idx < dimx) ? intervals_x_start[idx] : grid::getMaxValue<iter_type>();
  shared_intervals_y_start[tid] =
      (idx < dimy) ? intervals_y_start[idx] : grid::getMaxValue<iter_type>();
  shared_intervals_x_end[tid] =
      (idx < dimx) ? intervals_x_end[idx] : grid::getMinValue<iter_type>();
  shared_intervals_y_end[tid] =
      (idx < dimy) ? intervals_y_end[idx] : grid::getMinValue<iter_type>();
  __syncthreads();

  // Perform reduction to find max
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      shared_intervals_x_start[tid] =
          min(shared_intervals_x_start[tid], shared_intervals_x_start[tid + s]);
      shared_intervals_y_start[tid] =
          min(shared_intervals_y_start[tid], shared_intervals_y_start[tid + s]);
      shared_intervals_x_end[tid] =
          max(shared_intervals_x_end[tid], shared_intervals_x_end[tid + s]);
      shared_intervals_y_end[tid] =
          max(shared_intervals_y_end[tid], shared_intervals_y_end[tid + s]);
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    atomicMin(&intervals_x_start[0], shared_intervals_x_start[0]);
    atomicMin(&intervals_y_start[0], shared_intervals_y_start[0]);
    atomicMax(&intervals_x_end[0], shared_intervals_x_end[0]);
    atomicMax(&intervals_y_end[0], shared_intervals_y_end[0]);
  }
}

template <typename T>
__global__ void reduceIntervals3d(iter_type* intervals_0,
                                  iter_type* intervals_1, len_type dimx,
                                  len_type dimy, len_type dimz) {
  extern __shared__ iter_type sdata[];

  iter_type* shared_intervals_x_start = &sdata[0];
  iter_type* shared_intervals_y_start = &shared_intervals_x_start[blockDim.x];
  iter_type* shared_intervals_z_start = &shared_intervals_y_start[blockDim.x];
  iter_type* shared_intervals_x_end = &shared_intervals_z_start[blockDim.x];
  iter_type* shared_intervals_y_end = &shared_intervals_x_end[blockDim.x];
  iter_type* shared_intervals_z_end = &shared_intervals_y_end[blockDim.x];

  iter_type* intervals_x_start = &intervals_0[0];
  iter_type* intervals_y_start = &intervals_x_start[dimx];
  iter_type* intervals_z_start = &intervals_y_start[dimy];
  iter_type* intervals_x_end = &intervals_1[0];
  iter_type* intervals_y_end = &intervals_x_end[dimx];
  iter_type* intervals_z_end = &intervals_y_end[dimy];

  int tid = threadIdx.x;
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  // Load data into shared memory
  shared_intervals_x_start[tid] =
      (idx < dimx) ? intervals_x_start[idx] : grid::getMaxValue<iter_type>();
  shared_intervals_y_start[tid] =
      (idx < dimy) ? intervals_y_start[idx] : grid::getMaxValue<iter_type>();
  shared_intervals_z_start[tid] =
      (idx < dimz) ? intervals_z_start[idx] : grid::getMaxValue<iter_type>();
  shared_intervals_x_end[tid] =
      (idx < dimx) ? intervals_x_end[idx] : grid::getMinValue<iter_type>();
  shared_intervals_y_end[tid] =
      (idx < dimy) ? intervals_y_end[idx] : grid::getMinValue<iter_type>();
  shared_intervals_z_end[tid] =
      (idx < dimz) ? intervals_z_end[idx] : grid::getMinValue<iter_type>();
  __syncthreads();

  // Perform reduction to find max
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      shared_intervals_x_start[tid] =
          min(shared_intervals_x_start[tid], shared_intervals_x_start[tid + s]);
      shared_intervals_z_start[tid] =
          min(shared_intervals_z_start[tid], shared_intervals_z_start[tid + s]);
      shared_intervals_y_start[tid] =
          min(shared_intervals_y_start[tid], shared_intervals_y_start[tid + s]);
      shared_intervals_x_end[tid] =
          max(shared_intervals_x_end[tid], shared_intervals_x_end[tid + s]);
      shared_intervals_y_end[tid] =
          max(shared_intervals_y_end[tid], shared_intervals_y_end[tid + s]);
      shared_intervals_z_end[tid] =
          max(shared_intervals_z_end[tid], shared_intervals_z_end[tid + s]);
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    atomicMin(&intervals_x_start[0], shared_intervals_x_start[0]);
    atomicMin(&intervals_y_start[0], shared_intervals_y_start[0]);
    atomicMin(&intervals_z_start[0], shared_intervals_z_start[0]);
    atomicMax(&intervals_x_end[0], shared_intervals_x_end[0]);
    atomicMax(&intervals_y_end[0], shared_intervals_y_end[0]);
    atomicMax(&intervals_z_end[0], shared_intervals_z_end[0]);
  }
}

template <size_t D>
struct run_find_enclosing_intervals;

template <>
struct run_find_enclosing_intervals<1> {
  template <typename T>
  void call_kernel(iter_type* d_intervals_0, iter_type* d_intervals_1,
                   const T* values, T cutoff_value, iter_type total_length,
                   const iter_type (&dims)[1], const iter_type (&stride)[1],
                   const iter_type (&start_pos)[1]) {
    int numBlocks = (total_length + BLOCK_SIZE - 1) / BLOCK_SIZE;
    computeEnclosingInterval1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        values, d_intervals_0, d_intervals_1, cutoff_value, start_pos[0],
        dims[0], total_length);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
  template <typename T>
  void call_kernel(iter_type* d_intervals_0, iter_type* d_intervals_1,
                   T* const (&values)[1], any_vector_t<T, 1> cutoff_value,
                   iter_type total_length, const iter_type (&dims)[1],
                   const iter_type (&stride)[1],
                   const iter_type (&start_pos)[1]) {
    int numBlocks = (total_length + BLOCK_SIZE - 1) / BLOCK_SIZE;
    computeEnclosingInterval1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        values[0], d_intervals_0, d_intervals_1, cutoff_value, start_pos[0],
        dims[0], total_length);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename S, typename T>
  void operator()(S&& values, T const& cutoff_value, iter_type total_length,
                  const iter_type (&dims)[1], const iter_type (&stride)[1],
                  const iter_type (&offset)[1], const iter_type (&start_pos)[1],
                  iter_type (&intervals)[1][2]) {
    iter_type* d_intervals_0;
    iter_type* d_intervals_1;
    CHECK_CUDA_ERROR(cudaMalloc(&d_intervals_0, sizeof(iter_type)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_intervals_1, sizeof(iter_type)));
    symphas::cuda::initializeArray CUDA_KERNEL(1, BLOCK_SIZE)(
        d_intervals_0, std::numeric_limits<iter_type>::max(), 1);
    symphas::cuda::initializeArray CUDA_KERNEL(1, BLOCK_SIZE)(d_intervals_1, 0,
                                                              1);
    call_kernel(d_intervals_0, d_intervals_1, std::forward<S>(values),
                cutoff_value, total_length, dims, stride, offset, start_pos);
    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[0][0], &d_intervals_0[0],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[0][1], &d_intervals_1[0],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaFree(d_intervals_0));
    CHECK_CUDA_ERROR(cudaFree(d_intervals_1));
  }
};

template <>
struct run_find_enclosing_intervals<2> {
  template <typename T>
  void call_kernel(iter_type* d_intervals_0, iter_type* d_intervals_1,
                   const T* values, T cutoff_value, iter_type total_length,
                   const iter_type (&dims)[2], const iter_type (&stride)[2],
                   const iter_type (&offset)[2],
                   const iter_type (&start_pos)[2]) {
    dim3 blockDim(32, 32);
    dim3 gridDim((dims[0] + blockDim.x - 1) / blockDim.x,
                 (dims[1] + blockDim.y - 1) / blockDim.y);

    size_t sharedMemSize =
        (2 * blockDim.x + 2 * blockDim.y) * sizeof(iter_type);
    computeEnclosingInterval2d CUDA_KERNEL(gridDim, blockDim, sharedMemSize)(
        values, d_intervals_0, d_intervals_1, cutoff_value, total_length,
        start_pos[0], start_pos[1], dims[0], dims[1], offset[0], offset[1],
        stride[1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
  template <typename T>
  void call_kernel(iter_type* d_intervals_0, iter_type* d_intervals_1,
                   T* const (&values)[2], any_vector_t<T, 2> cutoff_value,
                   iter_type total_length, const iter_type (&dims)[2],
                   const iter_type (&stride)[2], const iter_type (&offset)[2],
                   const iter_type (&start_pos)[2]) {
    dim3 blockDim(32, 32);
    dim3 gridDim((dims[0] + blockDim.x - 1) / blockDim.x,
                 (dims[1] + blockDim.y - 1) / blockDim.y);

    size_t sharedMemSize =
        (2 * blockDim.x + 2 * blockDim.y) * sizeof(iter_type);
    computeEnclosingInterval2d CUDA_KERNEL(gridDim, blockDim, sharedMemSize)(
        values[0], values[1], d_intervals_0, d_intervals_1, cutoff_value,
        total_length, start_pos[0], start_pos[1], dims[0], dims[1], offset[0],
        offset[1], stride[1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename S, typename T>
  void operator()(S&& values, T const& cutoff_value, iter_type total_length,
                  const iter_type (&dims)[2], const iter_type (&stride)[2],
                  const iter_type (&offset)[2], const iter_type (&start_pos)[2],
                  iter_type (&intervals)[2][2]) {
    iter_type* d_intervals_0;
    iter_type* d_intervals_1;
    CHECK_CUDA_ERROR(
        cudaMalloc(&d_intervals_0, (dims[0] + dims[1]) * sizeof(iter_type)));
    CHECK_CUDA_ERROR(
        cudaMalloc(&d_intervals_1, (dims[0] + dims[1]) * sizeof(iter_type)));
    int numBlocks = (dims[0] + dims[1] + BLOCK_SIZE - 1) / BLOCK_SIZE;
    symphas::cuda::initializeArray CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        d_intervals_0, std::numeric_limits<iter_type>::max(),
        dims[0] + dims[1]);
    symphas::cuda::initializeArray CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        d_intervals_1, 0, dims[0] + dims[1]);

    call_kernel(d_intervals_0, d_intervals_1, std::forward<S>(values),
                cutoff_value, total_length, dims, stride, offset, start_pos);

    numBlocks = (total_length + BLOCK_SIZE - 1) / BLOCK_SIZE;
    reduceIntervals2d<T> CUDA_KERNEL(numBlocks, BLOCK_SIZE,
                                     4 * BLOCK_SIZE * sizeof(iter_type))(
        d_intervals_0, d_intervals_1, dims[0], dims[1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[0][0], &d_intervals_0[0],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[0][1], &d_intervals_1[0],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[1][0], &d_intervals_0[dims[0]],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[1][1], &d_intervals_1[dims[0]],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaFree(d_intervals_0));
    CHECK_CUDA_ERROR(cudaFree(d_intervals_1));
  }
};

template <>
struct run_find_enclosing_intervals<3> {
  template <typename T>
  void call_kernel(iter_type* d_intervals_0, iter_type* d_intervals_1,
                   const T* values, T cutoff_value, iter_type total_length,
                   const iter_type (&dims)[3], const iter_type (&stride)[3],
                   const iter_type (&offset)[3],
                   const iter_type (&start_pos)[3]) {
    dim3 blockDim(8, 8, 8);  // 16x16 threads per block
    dim3 gridDim((dims[0] + blockDim.x - 1) / blockDim.x,
                 (dims[1] + blockDim.y - 1) / blockDim.y,
                 (dims[2] + blockDim.z - 1) / blockDim.z);

    size_t sharedMemSize =
        (2 * blockDim.x + 2 * blockDim.y + 2 * blockDim.z) * sizeof(iter_type);

    computeEnclosingInterval3d CUDA_KERNEL(gridDim, blockDim, sharedMemSize)(
        values, d_intervals_0, d_intervals_1, cutoff_value, total_length,
        start_pos[0], start_pos[1], start_pos[2], dims[0], dims[1], dims[2],
        offset[0], offset[1], offset[2], stride[1], stride[2]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }
  template <typename T>
  void call_kernel(iter_type* d_intervals_0, iter_type* d_intervals_1,
                   T* const (&values)[3], any_vector_t<T, 3> cutoff_value,
                   iter_type total_length, const iter_type (&dims)[3],
                   const iter_type (&stride)[3], const iter_type (&offset)[3],
                   const iter_type (&start_pos)[3]) {
    dim3 blockDim(8, 8, 8);  // 16x16 threads per block
    dim3 gridDim((dims[0] + blockDim.x - 1) / blockDim.x,
                 (dims[1] + blockDim.y - 1) / blockDim.y,
                 (dims[2] + blockDim.z - 1) / blockDim.z);

    size_t sharedMemSize =
        (2 * blockDim.x + 2 * blockDim.y + 2 * blockDim.z) * sizeof(iter_type);

    computeEnclosingInterval3d CUDA_KERNEL(gridDim, blockDim, sharedMemSize)(
        values[0], values[1], values[2], d_intervals_0, d_intervals_1,
        cutoff_value, total_length, start_pos[0], start_pos[1], start_pos[2],
        dims[0], dims[1], dims[2], offset[0], offset[1], offset[2], stride[1],
        stride[2]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename S, typename T>
  void operator()(S&& values, T const& cutoff_value, iter_type total_length,
                  const iter_type (&dims)[3], const iter_type (&stride)[3],
                  const iter_type (&offset)[3], const iter_type (&start_pos)[3],
                  iter_type (&intervals)[3][2]) {
    iter_type* d_intervals_0;
    iter_type* d_intervals_1;
    CHECK_CUDA_ERROR(cudaMalloc(
        &d_intervals_0, (dims[0] + dims[1] + dims[2]) * sizeof(iter_type)));
    CHECK_CUDA_ERROR(cudaMalloc(
        &d_intervals_1, (dims[0] + dims[1] + dims[2]) * sizeof(iter_type)));
    CHECK_CUDA_ERROR(
        cudaMemset(d_intervals_0, std::numeric_limits<iter_type>::max(),
                   (dims[0] + dims[1] + dims[2]) * sizeof(iter_type)));
    CHECK_CUDA_ERROR(cudaMemset(
        d_intervals_1, 0, (dims[0] + dims[1] + dims[2]) * sizeof(iter_type)));

    call_kernel(d_intervals_0, d_intervals_1, std::forward<S>(values),
                cutoff_value, total_length, dims, stride, offset, start_pos);

    int numBlocks = (total_length + BLOCK_SIZE - 1) / BLOCK_SIZE;
    reduceIntervals3d<T> CUDA_KERNEL(numBlocks, BLOCK_SIZE,
                                     6 * BLOCK_SIZE * sizeof(iter_type))(
        d_intervals_0, d_intervals_1, dims[0], dims[1], dims[2]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[0][0], &d_intervals_0[0],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[0][1], &d_intervals_1[0],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[1][0], &d_intervals_0[dims[0]],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[1][1], &d_intervals_1[dims[0]],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[2][0],
                                &d_intervals_0[dims[0] + dims[1]],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&intervals[2][1],
                                &d_intervals_1[dims[0] + dims[1]],
                                sizeof(iter_type), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaFree(d_intervals_0));
    CHECK_CUDA_ERROR(cudaFree(d_intervals_1));
  }
};

template <typename T>
__global__ void computeStartingPosKernel1d(const T* values,
                                           iter_type* starting_pos,
                                           T cutoff_value,
                                           iter_type total_length,
                                           iter_type start) {
  extern __shared__ iter_type shared_min_pos[];
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  shared_min_pos[threadIdx.x] = grid::getMaxValue<iter_type>();
  __syncthreads();

  iter_type n = idx + start;
  if (idx < total_length && values[n] <= cutoff_value) {
    atomicMin(&shared_min_pos[0], idx);
  }
  __syncthreads();
  // Write local minimums to global memory
  if (threadIdx.x == 0) {
    atomicMin(&starting_pos[0], shared_min_pos[0]);
  }
}

template <typename T>
__global__ void computeStartingPosKernel1d(const T* values0,
                                           iter_type* starting_pos,
                                           any_vector_t<T, 1> cutoff_value,
                                           iter_type total_length,
                                           iter_type start) {
  extern __shared__ iter_type shared_min_pos[];
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  shared_min_pos[threadIdx.x] = grid::getMaxValue<iter_type>();
  __syncthreads();

  iter_type n = idx + start;
  T value[]{values0[n]};
  if (idx < total_length && grid::compare_cutoff(value, n, cutoff_value)) {
    atomicMin(&shared_min_pos[0], idx);
  }
  __syncthreads();
  // Write local minimums to global memory
  if (threadIdx.x == 0) {
    atomicMin(&starting_pos[0], shared_min_pos[0]);
  }
}

template <typename T>
__global__ void computeStartingPosKernel2d(
    const T* values, iter_type* local_min_pos_x, iter_type* local_min_pos_y,
    T cutoff_value, iter_type total_length, len_type dimx, len_type dimy,
    iter_type startx, iter_type starty, len_type stride) {
  extern __shared__ iter_type shared_min_pos[];

  iter_type* shared_min_pos_x = &shared_min_pos[0];
  iter_type* shared_min_pos_y = &shared_min_pos_x[blockDim.x];

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y + threadIdx.y;
  int tidx = threadIdx.x;
  int tidy = threadIdx.y;

  // Initialize shared memory
  shared_min_pos_x[tidx] = grid::getMaxValue<iter_type>();
  shared_min_pos_y[tidy] = grid::getMaxValue<iter_type>();
  __syncthreads();

  iter_type n = (startx + idx) + (starty + idy) * stride;
  if (idx < dimx && idy < dimy && values[n] <= cutoff_value) {
    atomicMin(&local_min_pos_x[idx], idy);
    atomicMin(&local_min_pos_y[idy], idx);
  }
  //__syncthreads();

  //// Write local minimums to global memory
  // if (tidy == 0) {
  //   atomicMin(&local_min_pos_x[idx], shared_min_pos_x[tidx]);
  // }
  // if (tidx == 0) {
  //   atomicMin(&local_min_pos_y[idy], shared_min_pos_y[tidy]);
  // }
}

template <typename T>
__global__ void computeStartingPosKernel2d(
    const T* values0, const T* values1, iter_type* local_min_pos_x,
    iter_type* local_min_pos_y, any_vector_t<T, 2> cutoff_value,
    iter_type total_length, len_type dimx, len_type dimy, iter_type startx,
    iter_type starty, len_type stride) {
  extern __shared__ iter_type shared_min_pos[];

  iter_type* shared_min_pos_x = &shared_min_pos[0];
  iter_type* shared_min_pos_y = &shared_min_pos_x[blockDim.x];

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y + threadIdx.y;
  int tidx = threadIdx.x;
  int tidy = threadIdx.y;

  // Initialize shared memory
  shared_min_pos_x[tidx] = grid::getMaxValue<iter_type>();
  shared_min_pos_y[tidy] = grid::getMaxValue<iter_type>();
  __syncthreads();

  iter_type n = (startx + idx) + (starty + idy) * stride;
  T value[]{values0[n], values1[n]};
  if (idx < dimx && idy < dimy && grid::compare_cutoff(value, cutoff_value)) {
    atomicMin(&local_min_pos_x[idx], idy);
    atomicMin(&local_min_pos_y[idy], idx);
  }
  //__syncthreads();

  //// Write local minimums to global memory
  // if (tidy == 0) {
  //   atomicMin(&local_min_pos_x[idx], shared_min_pos_x[tidx]);
  // }
  // if (tidx == 0) {
  //   atomicMin(&local_min_pos_y[idy], shared_min_pos_y[tidy]);
  // }
}

template <typename T>
__global__ void computeStartingPosKernel3d(
    const T* values, iter_type* local_min_pos_xy, iter_type* local_min_pos_xz,
    iter_type* local_min_pos_yz, T cutoff_value, iter_type total_length,
    len_type dimx, len_type dimy, len_type dimz, iter_type startx,
    iter_type starty, iter_type startz, len_type stridey, len_type stridez) {
  extern __shared__ iter_type shared_min_pos[];

  iter_type* shared_min_pos_xy = &shared_min_pos[0];
  iter_type* shared_min_pos_xz = &shared_min_pos_xy[blockDim.x * blockDim.y];
  iter_type* shared_min_pos_yz = &shared_min_pos_xz[blockDim.y * blockDim.z];

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y + threadIdx.y;
  int idz = blockIdx.z * blockDim.z + threadIdx.z;
  int tidx = threadIdx.x;
  int tidy = threadIdx.y;
  int tidz = threadIdx.z;

  // Initialize shared memory
  shared_min_pos_xy[tidx + tidy * blockDim.x] = grid::getMaxValue<iter_type>();
  shared_min_pos_xz[tidx + tidz * blockDim.x] = grid::getMaxValue<iter_type>();
  shared_min_pos_yz[tidy + tidz * blockDim.x] = grid::getMaxValue<iter_type>();
  __syncthreads();

  iter_type n =
      (startx + idx) + (starty + idy) * stridey + (startz + idz) * stridez;
  if (idx < dimx && idy < dimy && idz < dimz && values[n] <= cutoff_value) {
    atomicMin(&local_min_pos_xy[idx + idy * dimx], idz);
    atomicMin(&local_min_pos_xz[idx + idz * dimx], idy);
    atomicMin(&local_min_pos_yz[idy + idz * dimy], idx);
  }
  //__syncthreads();

  //// Write local minimums to global memory
  // if (tidz == 0) {
  //   atomicMin(&local_min_pos_xy[idx + idy * dimx],
  //             shared_min_pos_xy[tidx + tidy * blockDim.x]);
  // }
  // if (tidy == 0) {
  //   atomicMin(&local_min_pos_xz[idx + idz * dimx],
  //             shared_min_pos_xz[tidx + tidz * blockDim.x]);
  // }
  // if (tidx == 0) {
  //   atomicMin(&local_min_pos_yz[idy + idz * dimy],
  //             shared_min_pos_yz[tidy + tidz * blockDim.y]);
  // }
}

template <typename T>
__global__ void computeStartingPosKernel3d(
    const T* values0, const T* values1, const T* values2,
    iter_type* local_min_pos_xy, iter_type* local_min_pos_xz,
    iter_type* local_min_pos_yz, any_vector_t<T, 3> cutoff_value,
    iter_type total_length, len_type dimx, len_type dimy, len_type dimz,
    iter_type startx, iter_type starty, iter_type startz, len_type stridey,
    len_type stridez) {
  extern __shared__ iter_type shared_min_pos[];

  iter_type* shared_min_pos_xy = &shared_min_pos[0];
  iter_type* shared_min_pos_xz = &shared_min_pos_xy[blockDim.x * blockDim.y];
  iter_type* shared_min_pos_yz = &shared_min_pos_xz[blockDim.y * blockDim.z];

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y + threadIdx.y;
  int idz = blockIdx.z * blockDim.z + threadIdx.z;
  int tidx = threadIdx.x;
  int tidy = threadIdx.y;
  int tidz = threadIdx.z;

  // Initialize shared memory
  shared_min_pos_xy[tidx + tidy * blockDim.x] = grid::getMaxValue<iter_type>();
  shared_min_pos_xz[tidx + tidz * blockDim.x] = grid::getMaxValue<iter_type>();
  shared_min_pos_yz[tidy + tidz * blockDim.x] = grid::getMaxValue<iter_type>();
  __syncthreads();

  iter_type n =
      (startx + idx) + (starty + idy) * stridey + (startz + idz) * stridez;
  T value[]{values0[n], values1[n], values2[n]};
  if (idx < dimx && idy < dimy && idz < dimz &&
      grid::compare_cutoff(value, n, cutoff_value)) {
    atomicMin(&local_min_pos_xy[idx + idy * dimx], idz);
    atomicMin(&local_min_pos_xz[idx + idz * dimx], idy);
    atomicMin(&local_min_pos_yz[idy + idz * dimy], idx);
  }

  // copy from top
}

template <typename T>
__global__ void findMaxStartingPos2d(iter_type* starting_pos_x,
                                     iter_type* starting_pos_y, len_type dimx,
                                     len_type dimy) {
  extern __shared__ iter_type sdata[];

  iter_type* sdata_x = &sdata[0];
  iter_type* sdata_y = &sdata[blockDim.x];

  int tid = threadIdx.x;
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  // Load data into shared memory
  sdata_x[tid] = (idx < dimx) ? starting_pos_x[idx] : -FLT_MAX;
  sdata_y[tid] = (idx < dimy) ? starting_pos_y[idx] : -FLT_MAX;
  __syncthreads();

  // Perform reduction to find max
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdata_x[tid] = max(sdata_x[tid], sdata_x[tid + s]);
      sdata_y[tid] = max(sdata_y[tid], sdata_y[tid + s]);
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    atomicMax(&starting_pos_x[0], sdata_x[0]);
    atomicMax(&starting_pos_y[0], sdata_y[0]);
  }
}

template <typename T>
__global__ void findMaxStartingPos3d(iter_type* starting_pos_xy,
                                     iter_type* starting_pos_xz,
                                     iter_type* starting_pos_yz, len_type dimx,
                                     len_type dimy, len_type dimz) {
  extern __shared__ iter_type sdata[];

  iter_type* sdata_xy = &sdata[0];
  iter_type* sdata_xz = &sdata_xy[blockDim.x];
  iter_type* sdata_yz = &sdata_xz[blockDim.x];

  int tid = threadIdx.x;
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  // Load data into shared memory
  sdata_xy[tid] = (idx < dimx) ? starting_pos_xy[idx] : -FLT_MAX;
  sdata_xz[tid] = (idx < dimy) ? starting_pos_xz[idx] : -FLT_MAX;
  sdata_yz[tid] = (idx < dimz) ? starting_pos_yz[idx] : -FLT_MAX;
  __syncthreads();

  // Perform reduction to find max
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdata_xy[tid] = max(sdata_xy[tid], sdata_xy[tid + s]);
      sdata_xz[tid] = max(sdata_xz[tid], sdata_xz[tid + s]);
      sdata_yz[tid] = max(sdata_yz[tid], sdata_yz[tid + s]);
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    atomicMax(&starting_pos_xy[0], sdata_xy[0]);
    atomicMax(&starting_pos_xz[0], sdata_xz[0]);
    atomicMax(&starting_pos_yz[0], sdata_yz[0]);
  }
}

template <size_t D>
struct run_compute_starting_pos;

template <>
struct run_compute_starting_pos<1> {
  template <typename T>
  void call_kernel(iter_type* starting_pos, const T* values, T cutoff_value,
                   iter_type total_length, const iter_type (&dims)[1],
                   const iter_type (&stride)[1], const iter_type (&offset)[1]) {
    int numBlocks = (total_length + BLOCK_SIZE - 1) / BLOCK_SIZE;
    computeStartingPosKernel1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        values, starting_pos, cutoff_value, total_length, offset[0]);

    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename T>
  void call_kernel(iter_type* starting_pos, T* const (&values)[1],
                   any_vector_t<T, 1> const& cutoff_value,
                   iter_type total_length, const iter_type (&dims)[1],
                   const iter_type (&stride)[1], const iter_type (&offset)[1]) {
    int numBlocks = (total_length + BLOCK_SIZE - 1) / BLOCK_SIZE;
    computeStartingPosKernel1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        values, starting_pos, cutoff_value, total_length, offset[0]);

    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename S, typename T>
  void operator()(S&& values, T const& cutoff_value, iter_type total_length,
                  const iter_type (&dims)[1], const iter_type (&stride)[1],
                  const iter_type (&offset)[1], iter_type (&pos)[1]) {
    iter_type* starting_pos;
    CHECK_CUDA_ERROR(cudaMalloc(&starting_pos, sizeof(iter_type)));

    call_kernel(starting_pos, std::forward<S>(values), cutoff_value,
                total_length, dims, stride, offset);

    CHECK_CUDA_ERROR(cudaMemcpy(&pos[0], starting_pos, sizeof(iter_type),
                                cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaFree(starting_pos));
  }
};

template <>
struct run_compute_starting_pos<2> {
  template <typename T>
  void call_kernel(iter_type* starting_pos_x, iter_type* starting_pos_y,
                   const T* values, T cutoff_value, iter_type total_length,
                   const iter_type (&dims)[2], const iter_type (&stride)[2],
                   const iter_type (&offset)[2]) {
    dim3 blockDim(16, 16);  // 16x16 threads per block
    dim3 gridDim((dims[0] + blockDim.x - 1) / blockDim.x,
                 (dims[1] + blockDim.y - 1) / blockDim.y);
    computeStartingPosKernel2d<T> CUDA_KERNEL(
        gridDim, blockDim, (blockDim.x + blockDim.y) * sizeof(iter_type))(
        values, starting_pos_x, starting_pos_y, cutoff_value, total_length,
        dims[0], dims[1], offset[0], offset[1], stride[1]);

    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename T>
  void call_kernel(iter_type* starting_pos_x, iter_type* starting_pos_y,
                   T* const (&values)[2],
                   any_vector_t<T, 2> const& cutoff_value,
                   iter_type total_length, const iter_type (&dims)[2],
                   const iter_type (&stride)[2], const iter_type (&offset)[2]) {
    dim3 blockDim(16, 16);  // 16x16 threads per block
    dim3 gridDim((dims[0] + blockDim.x - 1) / blockDim.x,
                 (dims[1] + blockDim.y - 1) / blockDim.y);
    computeStartingPosKernel2d<T> CUDA_KERNEL(
        gridDim, blockDim, (blockDim.x + blockDim.y) * sizeof(iter_type))(
        values[0], values[1], starting_pos_x, starting_pos_y, cutoff_value,
        total_length, dims[0], dims[1], offset[0], offset[1], stride[1]);

    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename S, typename T>
  void operator()(S&& values, T const& cutoff_value, iter_type total_length,
                  const iter_type (&dims)[2], const iter_type (&stride)[2],
                  const iter_type (&offset)[2], iter_type (&pos)[2]) {
    iter_type* starting_pos;

    CHECK_CUDA_ERROR(
        cudaMalloc(&starting_pos, (dims[0] + dims[1]) * sizeof(iter_type)));
    int numBlocks = (dims[0] + dims[1] + BLOCK_SIZE - 1) / BLOCK_SIZE;
    symphas::cuda::initializeArray CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        starting_pos, std::numeric_limits<iter_type>::max(), dims[0] + dims[1]);

    iter_type* starting_pos_host = new iter_type[dims[0] + dims[1]];

    iter_type* starting_pos_x = &starting_pos[0];
    iter_type* starting_pos_y = &starting_pos[dims[0]];

    CHECK_CUDA_ERROR(cudaMemcpy(starting_pos_host, starting_pos,
                                (dims[0] + dims[1]) * sizeof(iter_type),
                                cudaMemcpyDeviceToHost));

    call_kernel(starting_pos_x, starting_pos_y, std::forward<S>(values),
                cutoff_value, total_length, dims, stride, offset);

    len_type max_dim = std::max(dims[0], dims[1]);
    int num_blocks = (max_dim + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_CUDA_ERROR(cudaMemcpy(starting_pos_host, starting_pos,
                                (dims[0] + dims[1]) * sizeof(iter_type),
                                cudaMemcpyDeviceToHost));

    findMaxStartingPos2d<T> CUDA_KERNEL(num_blocks, BLOCK_SIZE,
                                        2 * BLOCK_SIZE * sizeof(iter_type))(
        starting_pos_x, starting_pos_y, dims[0], dims[1]);

    // Copy results to host
    CHECK_CUDA_ERROR(cudaMemcpy(&pos[1], &starting_pos_x[0], sizeof(iter_type),
                                cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&pos[0], &starting_pos_y[0], sizeof(iter_type),
                                cudaMemcpyDeviceToHost));

    // Free allocated memory
    delete[] starting_pos_host;
    CHECK_CUDA_ERROR(cudaFree(starting_pos));
  }
};

template <>
struct run_compute_starting_pos<3> {
  template <typename T>
  void call_kernel(iter_type* starting_pos_xy, iter_type* starting_pos_xz,
                   iter_type* starting_pos_yz, const T* values, T cutoff_value,
                   iter_type total_length, const iter_type (&dims)[3],
                   const iter_type (&stride)[3], const iter_type (&offset)[3]) {
    dim3 blockDim(8, 8, 8);  // 8x8x8 threads per block
    dim3 gridDim((dims[0] + blockDim.x - 1) / blockDim.x,
                 (dims[1] + blockDim.y - 1) / blockDim.y,
                 (dims[2] + blockDim.z - 1) / blockDim.z);

    computeStartingPosKernel3d<T> CUDA_KERNEL(
        num_blocks, BLOCK_SIZE,
        (blockDim.x * blockDim.y + blockDim.x * blockDim.z +
         blockDim.y * blockDim.z) *
            sizeof(iter_type))(values, starting_pos_xy, starting_pos_xz,
                               starting_pos_yz, cutoff_value, total_length,
                               dims[0], dims[1], dims[2], offset[0], offset[1],
                               offset[2], stride[1], stride[2]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename T>
  void call_kernel(iter_type* starting_pos_xy, iter_type* starting_pos_xz,
                   iter_type* starting_pos_yz, T* const (&values)[3],
                   any_vector_t<T, 3> const& cutoff_value,
                   iter_type total_length, const iter_type (&dims)[3],
                   const iter_type (&stride)[3], const iter_type (&offset)[3]) {
    dim3 blockDim(8, 8, 8);  // 8x8x8 threads per block
    dim3 gridDim((dims[0] + blockDim.x - 1) / blockDim.x,
                 (dims[1] + blockDim.y - 1) / blockDim.y,
                 (dims[2] + blockDim.z - 1) / blockDim.z);

    computeStartingPosKernel3d<T> CUDA_KERNEL(
        num_blocks, BLOCK_SIZE,
        (blockDim.x * blockDim.y + blockDim.x * blockDim.z +
         blockDim.y * blockDim.z) *
            sizeof(iter_type))(
        values[0], values[1], values[2], starting_pos_xy, starting_pos_xz,
        starting_pos_yz, cutoff_value, total_length, dims[0], dims[1], dims[2],
        offset[0], offset[1], offset[2], stride[1], stride[2]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename S, typename T>
  void operator()(S&& values, T const& cutoff_value, iter_type total_length,
                  const iter_type (&dims)[3], const iter_type (&stride)[3],
                  const iter_type (&offset)[3], iter_type (&pos)[3]) {
    iter_type* starting_pos;
    iter_type num_blocks = (total_length + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CHECK_CUDA_ERROR(
        cudaMalloc(&starting_pos,
                   (dims[0] * dims[1] + dims[0] * dims[2] + dims[1] * dims[2]) *
                       sizeof(iter_type)));
    int numBlocks =
        ((dims[0] * dims[1] + dims[0] * dims[2] + dims[1] * dims[2]) +
         BLOCK_SIZE - 1) /
        BLOCK_SIZE;
    symphas::cuda::initializeArray CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        starting_pos, std::numeric_limits<iter_type>::max(),
        (dims[0] * dims[1] + dims[0] * dims[2] + dims[1] * dims[2]));

    iter_type* starting_pos_xy = &starting_pos[0];
    iter_type* starting_pos_xz = &starting_pos[dims[0] * dims[1]];
    iter_type* starting_pos_yz =
        &starting_pos[dims[0] * dims[1] + dims[0] * dims[2]];

    call_kernel(starting_pos_xy, starting_pos_xz, starting_pos_yz,
                std::forward<S>(values), cutoff_value, total_length, dims,
                stride, offset);

    len_type max_dim = std::max(dims[0] * dims[1],
                                std::max(dims[0] * dims[2], dims[1] * dims[2]));
    num_blocks = (max_dim + BLOCK_SIZE - 1) / BLOCK_SIZE;
    findMaxStartingPos3d<T> CUDA_KERNEL(num_blocks, BLOCK_SIZE,
                                        3 * BLOCK_SIZE * sizeof(iter_type))(
        starting_pos_xy, starting_pos_xz, starting_pos_yz, dims[0] * dims[1],
        dims[0] * dims[2], dims[1] * dims[2]);

    // Copy results to host
    CHECK_CUDA_ERROR(cudaMemcpy(&pos[2], starting_pos_xy, sizeof(iter_type),
                                cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&pos[1], starting_pos_xz, sizeof(iter_type),
                                cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(&pos[0], starting_pos_yz, sizeof(iter_type),
                                cudaMemcpyDeviceToHost));

    // Free allocated memory
    CHECK_CUDA_ERROR(cudaFree(starting_pos));
  }
};
//
// template <size_t D>
// struct run_update_intervals;
//
// template <>
// struct run_update_intervals<1> {
//  template <typename T>
//  void operator()(const T* values, iter_type* intervals, T cutoff_value,
//                  iter_type total_length, iter_type (&stride)[1]) {
//    int numBlocks = (total_length + BLOCK_SIZE - 1) / BLOCK_SIZE;
//    updateIntervalsKernel1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
//        values, intervals, cutoff_value, total_length);
//    CHECK_CUDA_ERROR(cudaPeekAtLastError());
//    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
//  }
//};
//
// template <>
// struct run_update_intervals<2> {
//  template <typename T>
//  void operator()(const T* values, iter_type* intervals, T cutoff_value,
//                  iter_type total_length, iter_type (&stride)[2]) {
//    int numBlocks = (total_length + BLOCK_SIZE - 1) / BLOCK_SIZE;
//    updateIntervalsKernel2d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
//        values, intervals, cutoff_value, total_length, stride[0]);
//    CHECK_CUDA_ERROR(cudaPeekAtLastError());
//    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
//  }
//};
//
// template <>
// struct run_update_intervals<3> {
//  template <typename T>
//  void operator()(const T* values, iter_type* intervals, T cutoff_value,
//                  iter_type total_length, iter_type (&stride)[3]) {
//    int numBlocks = (total_length + BLOCK_SIZE - 1) / BLOCK_SIZE;
//    updateIntervalsKernel3d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
//        values, intervals, cutoff_value, total_length, stride[0],
//        stride[1]);
//    CHECK_CUDA_ERROR(cudaPeekAtLastError());
//    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
//  }
//};

namespace grid {

//! Adjust the floating region inside the grid to a new position and
//! dimensions
/*!
 * Adjust the floating region inside the grid to a new position and
 * dimensions.
 *
 * \param new_values The new region to fit the old region inside. This needs
 * to have as much space as the given dimensions with an added boundary
 * layer.
 *
 * \param new_origin The new starting point of the new region inside the
 * global domain. This origin point includes the boundary layer as well,
 * meaning that the interior data starts after this origin point.
 *
 * \param new_dimensions The dimensions of the new region, also counting the
 * boundary layer.
 *
 * \param old_values The old region which contains data, from which the data
 * that is overlapping with the new region will be copied.
 *
 * \param old_origin The starting point of the old region inside the global
 * domain. The origin point includes the boundary layer as well, meaning
 * that the interior data starts after this origin point.
 *
 * \param old_dims The dimensions of the old region, also counting the
 * boundary.
 *
 * \param global_dims The dimensions of the global domain.
 *
 * \param boundary_size The width of the boundary region inside both the old
 * and new regions.
 */
template <typename T, size_t D>
void adjust_region_to_from_cuda(
    T* new_values_device, const iter_type (&new_origin)[D],
    const len_type (&new_dims)[D], const T* old_values_device,
    const iter_type (&old_origin)[D], const len_type (&old_dims)[D],
    const len_type (&global_dims)[D], const T empty, len_type boundary_size) {
  RegionAdjustParams<D> rp{new_origin, new_dims,    old_origin,
                           old_dims,   global_dims, boundary_size};
  RegionAdjustParams<D>* rpDev;
  CHECK_CUDA_ERROR(cudaMalloc(&rpDev, sizeof(RegionAdjustParams<D>)));
  CHECK_CUDA_ERROR(cudaMemcpy(rpDev, &rp, sizeof(RegionAdjustParams<D>),
                              cudaMemcpyHostToDevice));
  int numBlocks =
      (grid::length<D>(rp.new_interior_dims) + BLOCK_SIZE - 1) / BLOCK_SIZE;

  kernelAdjustRegionToFrom CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      new_values_device, old_values_device, rpDev, boundary_size);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(rpDev));
}

//! Adjust the floating region inside the grid to a new position and
//! dimensions
/*!
 * Adjust the floating region inside the grid to a new position and
 * dimensions.
 *
 * \param values The existing region which will have its values moved in
 * place to accommodate a new origin.
 *
 * \param new_origin The new starting point of the new region inside the
 * global domain. This origin point includes the boundary layer as well,
 * meaning that the interior data starts after this origin point.
 *
 * \param old_origin The starting point of the old region inside the global
 * domain. The origin point includes the boundary layer as well, meaning
 * that the interior data starts after this origin point.
 *
 * \param dims The dimensions of the existing region, counting the boundary
 * layer.
 *
 * \param global_dims The dimensions of the global domain.
 *
 * \param boundary_size The width of the boundary region inside both the old
 * and new regions.
 */
template <typename T, size_t D>
void adjust_origin_to_from_replace_cuda(T*(&values_device),
                                        const iter_type (&new_origin)[D],
                                        const iter_type (&old_origin)[D],
                                        const len_type (&dims)[D],
                                        const len_type (&global_dims)[D],
                                        T empty, len_type boundary_size) {
  T* newValuesDevice;
  CHECK_CUDA_ERROR(
      cudaMalloc(&newValuesDevice, grid::length<D>(dims) * sizeof(T)));

  RegionAdjustParams<D> rp{new_origin, dims,        old_origin,
                           dims,       global_dims, boundary_size};
  RegionAdjustParams<D>* rpDev;
  CHECK_CUDA_ERROR(cudaMalloc(&rpDev, sizeof(RegionAdjustParams<D>)));
  CHECK_CUDA_ERROR(cudaMemcpy(rpDev, &rp, sizeof(RegionAdjustParams<D>),
                              cudaMemcpyHostToDevice));
  int numBlocks =
      (grid::length<D>(rp.new_interior_dims) + BLOCK_SIZE - 1) / BLOCK_SIZE;

  symphas::cuda::initializeArray CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      newValuesDevice, empty, grid::length<D>(dims));
  kernelAdjustRegionToFrom CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      newValuesDevice, values_device, rpDev, boundary_size);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaFree(rpDev));

  std::swap(newValuesDevice, values_device);
  CHECK_CUDA_ERROR(cudaFree(newValuesDevice));
}

//! Adjust the floating region inside the grid to a new position and
//! dimensions
/*!
 * Adjust the floating region inside the grid to a new position and
 * dimensions.
 *
 * \param values The existing region which will have its values moved in
 * place to accommodate a new origin.
 *
 * \param new_origin The new starting point of the new region inside the
 * global domain. This origin point includes the boundary layer as well,
 * meaning that the interior data starts after this origin point.
 *
 * \param old_origin The starting point of the old region inside the global
 * domain. The origin point includes the boundary layer as well, meaning
 * that the interior data starts after this origin point.
 *
 * \param dims The dimensions of the existing region, counting the boundary
 * layer.
 *
 * \param global_dims The dimensions of the global domain.
 *
 * \param boundary_size The width of the boundary region inside both the old
 * and new regions.
 */
template <typename T, size_t D>
void adjust_origin_to_from_cuda(T*(&values), const iter_type (&new_origin)[D],
                                const iter_type (&old_origin)[D],
                                const len_type (&dims)[D],
                                const len_type (&global_dims)[D], T empty,
                                len_type boundary_size) {
  adjust_origin_to_from_replace_cuda(values, new_origin, old_origin, dims,
                                     global_dims, empty, boundary_size);
}

//! Construct a new view around a region with values greater than a cutoff.
/*!
 * Given the regional grid, populate the origin and dimensions for the
 * minimal new region where all values outside are less than cutoff. The
 * regional grid dimensions give the total region size including boundaries,
 * and this algorithm will compute the dimensions without extending them to
 * boundaries.
 *
 * \param grid The regional grid from which to find a new view.
 * \param origin The values of the new region origin are populated here.
 * \param dims The values of the dimensions of the new region are populated
 * here. \param boundary_size The width of the boundary used in the
 * regional_grid.
 */
template <typename T, size_t D>
void get_view_resized_cuda(RegionalGridCUDA<T, D>& grid, T cutoff_value,
                           iter_type (&origin)[D], len_type (&dims)[D]) {
  MinimalRegionParams<D> params(grid.region.dims);

  // Allocate device memory for min and max indices
  iter_type* d_min_indices;
  iter_type* d_max_indices;
  MinimalRegionParams<D>* d_params;
  CHECK_CUDA_ERROR(cudaMalloc(&d_min_indices, D * sizeof(iter_type)));
  CHECK_CUDA_ERROR(cudaMalloc(&d_max_indices, D * sizeof(iter_type)));
  CHECK_CUDA_ERROR(cudaMalloc(&d_params, sizeof(MinimalRegionParams<D>)));

  CHECK_CUDA_ERROR(cudaMemcpy(d_params, &params, sizeof(MinimalRegionParams<D>),
                              cudaMemcpyHostToDevice));

  // Initialize min and max indices
  iter_type h_min_indices[D];
  iter_type h_max_indices[D];
  for (size_t i = 0; i < D; ++i) {
    h_min_indices[i] = grid.region.dims[i] - grid.region.boundary_size;
    h_max_indices[i] = grid.region.boundary_size;
  }
  CHECK_CUDA_ERROR(cudaMemcpy(d_min_indices, h_min_indices,
                              D * sizeof(iter_type), cudaMemcpyHostToDevice));
  CHECK_CUDA_ERROR(cudaMemcpy(d_max_indices, h_max_indices,
                              D * sizeof(iter_type), cudaMemcpyHostToDevice));

  len_type len = grid::length<D>(grid.region.dims);

  // Launch kernel to find the minimal region
  run_find_minimal_region<D>{}(grid.values, cutoff_value, d_params,
                               d_min_indices, d_max_indices, len);

  // Copy results back to host
  CHECK_CUDA_ERROR(cudaMemcpy(h_min_indices, d_min_indices,
                              D * sizeof(iter_type), cudaMemcpyDeviceToHost));
  CHECK_CUDA_ERROR(cudaMemcpy(h_max_indices, d_max_indices,
                              D * sizeof(iter_type), cudaMemcpyDeviceToHost));

  // Set the origin and dimensions of the new region
  for (size_t i = 0; i < D; ++i) {
    origin[i] =
        h_min_indices[i] - grid.region.boundary_size + grid.region.origin[i];
    dims[i] = std::max(0, h_max_indices[i] - h_min_indices[i] + 1 +
                              grid.region.boundary_size * 2);
  }

  // Free device memory
  CHECK_CUDA_ERROR(cudaFree(d_min_indices));
  CHECK_CUDA_ERROR(cudaFree(d_max_indices));
  CHECK_CUDA_ERROR(cudaFree(d_params));
}

//! Construct a new view around a region with values greater than a cutoff.
/*!
 * Given the regional grid, populate the origin for the new region
 * where all values outside are less than cutoff. The dimensions of the
 * regional grid will be the same, since only the position of the view
 * itself is adjusted.
 *
 * \param grid The regional grid from which to find a new view.
 * \param origin The values of the new region origin are populated here.
 * \param boundary_size The width of the boundary used in the regional_grid.
 */
template <typename T, size_t D>
void get_view(RegionalGridCUDA<T, D>& grid, T cutoff_value,
              iter_type (&origin)[D]) {
  iter_type dims[D];

  for (iter_type i = 0; i < D; ++i) {
    dims[i] = grid.region.dims[i];
  }

  get_view_resized_cuda(grid, cutoff_value, origin, dims);

  for (iter_type i = 0; i < D; ++i) {
    iter_type delta = grid.region.dims[i] - dims[i];
    origin[i] = grid.region.origin[i] + (origin[i] - (delta / 2));
  }
}

/**
 * @brief Identifies a sub-region within a multi-dimensional grid where the
 * values exceed a specified cutoff value. Adjusts the view's origin and
 * dimensions to encompass this sub-region, considering periodic boundary
 * conditions.
 *
 * @tparam T The data type of the grid values.
 * @tparam D The dimensionality of the grid.
 * @param grid A reference to a RegionalGrid object containing the grid data
 * and metadata.
 * @param cutoff_value The threshold value used to determine significant
 * grid values.
 * @param origin An array to store the origin of the resized view.
 * @param dims An array to store the dimensions of the resized view.
 *
 * This function performs the following steps:
 * 1. Initializes `offset` and `stride` arrays to store the offset and
 * stride for each dimension.
 * 2. Adjusts the dimensions to exclude the boundary size and sets the
 * offset to the boundary size.
 * 3. Initializes the stride based on the grid dimensions.
 * 4. Iterates through the grid to find the first position where the value
 * exceeds the cutoff.
 * 5. If no values exceed the cutoff, sets all dimensions to zero.
 * 6. Initializes intervals to track the region of significant values.
 * 7. Iteratively expands the intervals to include all significant values.
 * 8. Updates the origin and dimensions based on the final intervals.
 *
 * Example usage:
 * @code
 * RegionalGrid<double, 3> grid = ...; // Assume this is initialized
 * double cutoff_value = 0.5;
 * iter_type origin[3];
 * len_type dims[3];
 *
 * get_view_resized_periodic(grid, cutoff_value, origin, dims);
 * // Now, origin and dims contain the starting indices and sizes of the
 * resized view
 * @endcode
 *
 * @note The function uses periodic boundary conditions, meaning it wraps
 * around the grid boundaries.
 * @note The function is templated to work with any data type `T` and any
 * dimensionality `D`.
 *
 * @see init_region_stride
 * @see get_grid_position_offset
 * @see index_from_position
 * @see compare_cutoff
 * @see is_fully_overlapping
 */
template <typename T, size_t D>
void get_view_resized_periodic_cuda(RegionalGridCUDA<T, D>& grid,
                                    T cutoff_value, iter_type (&origin)[D],
                                    len_type (&dims)[D]) {
  iter_type stride[D]{};
  iter_type total_length = 1;
  iter_type offset[D]{};
  iter_type pos[D]{};

  for (iter_type i = 0; i < D; ++i) {
    dims[i] = grid.region.dims[i] - grid.region.boundary_size * 2;
    offset[i] = grid.region.boundary_size;
    total_length *= dims[i];
  }
  init_region_stride(stride, grid.region.dims);

  run_compute_starting_pos<D>{}(grid.values, cutoff_value, total_length, dims,
                                stride, offset, pos);

  if (grid::index_from_position(pos, stride) >= total_length) {
    for (iter_type i = 0; i < D; ++i) {
      dims[i] = 0;
    }
    return;
  }

  iter_type intervals[D][2]{};
  run_find_enclosing_intervals<D>{}(grid.values, cutoff_value, total_length,
                                    dims, stride, offset, pos, intervals);

  for (iter_type i = 0; i < D; ++i) {
    origin[i] = intervals[i][0] + grid.region.origin[i];
    dims[i] = std::max(0, intervals[i][1] - intervals[i][0] + 1 +
                              grid.region.boundary_size * 2);
  }
}

template <typename T, size_t D>
void adjust_region(RegionalGridCUDA<T, D>& grid, T cutoff) {
  iter_type origin[D];
  get_view(grid, cutoff, origin);

  for (iter_type i = 0; i < D; ++i) {
    len_type delta = (grid.dims[i] - grid.region.boundary_size * 2);
    origin[i] += (origin[i] < 0) ? delta : (origin[i] >= delta) ? -delta : 0;
  }

  grid.adjust(origin);
}

//! Adjust the region by finding a new origin and dimensions.
/*!
 * Adjust the region by finding a new origin and dimensions.
 *
 * \param grid The grid to resize.
 * \param cutoff The cutoff value to use when determining the region.
 * \param padding_factor Increases the smallest determined dimensions by the
 * given ratio. \param dims_relative_eps Computes the relative distance
 * between new and current dimensions, and uses the new dimensions if the
 * relative distance exceeds this value. The relative distance is defined
 * as:
 * $\f(\\text{dims}_0 - \\text{dims}_1) / \\text{dims}_0\f$, where subscript
 * 0 represents current dimensions and subscript 1 represents the newly
 * determined dimensions.
 */
template <typename T, size_t D>
void resize_adjust_region(RegionalGridCUDA<T, D>& grid, T cutoff,
                          double padding_factor = 1.0,
                          double dims_relative_eps = 0.0) {
  iter_type origin[D]{};
  len_type dims[D]{};

  len_type dims_set[D]{};
  len_type origin_set[D]{};

  if (grid::is_same_point(grid.dims, grid.region.dims)) {
    get_view_resized_periodic_cuda(grid, cutoff, origin, dims);
  } else {
    get_view_resized_cuda(grid, cutoff, origin, dims);
  }

  bool same_dims_flag = true;
  for (iter_type i = 0; i < D; ++i) {
    auto dim0 = std::min(grid.dims[i], iter_type(dims[i] * padding_factor));
    origin_set[i] = origin[i] - (dim0 - dims[i]) / 2;
    origin_set[i] += (origin_set[i] < 0)
                         ? (grid.dims[i] - grid.region.boundary_size * 2)
                         : 0;
    dims_set[i] = dim0;

    // if the new dimensions are close to the current dimensions, don't
    // change it
    if (std::abs(dims_set[i] - grid.region.dims[i]) /
            double(grid.region.dims[i]) >
        dims_relative_eps) {
      same_dims_flag = false;
    }
  }

  if (same_dims_flag) {
    for (iter_type i = 0; i < D; ++i) {
      origin_set[i] = origin[i] - (grid.region.dims[i] - dims[i]) / 2;
      origin_set[i] += (origin_set[i] < 0)
                           ? (grid.dims[i] - grid.region.boundary_size * 2)
                           : 0;
    }
    grid.adjust(grid::select_region_cuda<D>(origin_set, grid.region.dims));
  } else {
    grid.adjust(grid::select_region_cuda<D>(origin_set, dims_set));
  }
}

template <typename T, size_t D>
void resize_adjust_region(RegionalGridCUDA<T, D>& grid,
                          symphas::grid_info const& info) {
  len_type origin[D]{};
  len_type dims[D]{};
  for (auto const& [axis, interval] : info) {
    dims[symphas::axis_to_index(axis)] = interval.get_interval_count();
    origin[symphas::axis_to_index(axis)] =
        len_type((interval.left() - interval.domain_left()) / interval.width());
  }
  grid.adjust(grid::select_region_cuda<D>(origin, dims));
}

template <typename T, size_t D>
void resize_adjust_region(RegionalGridCUDA<T, D>& grid,
                          const len_type (&intervals)[D][2]) {
  len_type origin[D]{};
  len_type dims[D]{};
  for (iter_type i = 0; i < D; ++i) {
    dims[i] = intervals[i][1] - intervals[i][0];
    origin[i] = intervals[i][0];
  }
  grid.adjust(grid::select_region_cuda<D>(origin, dims));
}

template <typename T, size_t D>
void resize_adjust_region(RegionalGridCUDA<T, D>& grid,
                          grid::region_interval<D> const& interval) {
  resize_adjust_region(grid, interval.intervals);
}

template <typename T, size_t D>
void resize_adjust_region(RegionalGridCUDA<T, D>& grid,
                          grid::region_interval_multiple<D> const& regions) {
  len_type intervals[D][2]{};

  for (iter_type i = 0; i < D; ++i) {
    intervals[i][0] = grid.dims[i];
    intervals[i][1] = 0;
  }

  for (grid::region_interval<D> region : regions) {
    for (iter_type i = 0; i < D; ++i) {
      intervals[i][0] = std::min(intervals[i][0], region[i][0]);
      intervals[i][1] = std::max(intervals[i][1], region[i][1]);
    }
  }
  resize_adjust_region(grid, intervals);
}

template <typename T, size_t D>
void resize_adjust_region(RegionalGridCUDA<T, D>& grid, T cutoff,
                          const len_type (&minimum_dims)[D]) {
  iter_type origin[D];
  len_type dims[D];
  if (grid::is_same_point(grid.dims, grid.region.dims)) {
    get_view_resized_periodic_cuda(grid, cutoff, origin, dims);
  } else {
    get_view_resized_cuda(grid, cutoff, origin, dims);
  }

  for (iter_type i = 0; i < D; ++i) {
    iter_type delta = minimum_dims[i] - dims[i];
    if (delta > 0) {
      origin[i] -= delta / 2;
      dims[i] += delta - (delta / 2);
    }
  }
  grid.adjust(grid::select_region_cuda<D>(origin, dims));
}

}  // namespace grid
#endif
