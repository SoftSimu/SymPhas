
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
 * MODULE:  sym
 * PURPOSE: Defines all the expression types used in the symbolic
 * algebra for cuda functionality.
 *
 * ***************************************************************************
 */

#include <boundarysystem.h>

#include "expressionderivatives.cuh"
#include "expressionsymbolscuda.cuh"
#include "expressiontypeincludes.h"
#include "expressiontypenoise.cuh"

#ifdef USING_CUDA

#include <cuda_runtime.h>

namespace symphas::internal {

template <typename T>
struct coeff_data_cuda {
  coeff_data_cuda(coeff_data<T> const &data)
      : coeff_data_cuda(data.data, data.len) {}

  coeff_data_cuda(const T **data, len_type len) : data{}, len{len} {
    if (len > 0) {
      CHECK_CUDA_ERROR(cudaMalloc(&this->data, sizeof(T) * len));
      CHECK_CUDA_ERROR(cudaMemcpy(this->data, *data, sizeof(T) * len,
                                  cudaMemcpyDeviceToDevice));
    }
  }

  coeff_data_cuda(const T *data, len_type len) : data{}, len{len} {
    if (len > 0) {
      CHECK_CUDA_ERROR(cudaMalloc(&this->data, sizeof(T) * len));
      CHECK_CUDA_ERROR(cudaMemcpy(this->data, data, sizeof(T) * len,
                                  cudaMemcpyHostToDevice));
    }
  }

  coeff_data_cuda(coeff_data_cuda<T> const &other)
      : coeff_data_cuda(&other.data, other.len) {}
  coeff_data_cuda(coeff_data_cuda<T> &&other) : coeff_data_cuda() {
    swap(*this, other);
  }

  auto operator=(coeff_data_cuda<T> other) { swap(*this, other); }

  friend void swap(coeff_data_cuda<T> &first, coeff_data_cuda<T> &second) {
    using std::swap;
    swap(first.data, second.data);
    swap(first.len, second.len);
  }

  const auto &operator[](iter_type i) const { return data[i]; }

  auto &operator[](iter_type i) { return data[i]; }

  T *data;
  len_type len;

  ~coeff_data_cuda() { CHECK_CUDA_ERROR(cudaFree(data)); }
};
}  // namespace symphas::internal

template <size_t D>
__device__ bool is_in_region(const iter_type (&intervals)[D][2],
                             const len_type (&dimensions)[D], iter_type n) {
  int pos[D]{};
  int stride[D]{};

  len_type len = 1;

  for (iter_type i = 0; i < D; ++i) {
    stride[i] = 1;
    len *= dimensions[i];
    for (iter_type j = 0; j < i; ++j) {
      stride[i] *= dimensions[j];
    }
  }

  if (n >= len) return false;

  for (iter_type i = 0; i < D; ++i) {
    pos[i] = (n / stride[i]) % dimensions[i];
  }

  bool in_region = true;
  for (iter_type i = 0; i < D; ++i) {
    in_region =
        in_region && (pos[i] >= intervals[i][0]) && (pos[i] < intervals[i][1]);
  }
  return in_region;
}

template <typename T, typename E>
__global__ void evaluateCudaExpr(E *e, T *values) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    values[0] = e->eval(0);
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExpr(E *e, T *values, len_type len) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n < len) {
    values[n] = e->eval(n);
  }
}

template <typename T, typename E, typename... Us>
__global__ void evaluateCudaExpr1d(E *e, T *values, len_type dim0,
                                   iter_type interval00, iter_type interval01) {
  iter_type intervals[1][2]{{interval00, interval01}};
  iter_type dimensions[1]{dim0};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    values[n] = e->eval(n);
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExpr2d(E *e, T *values, len_type dim0,
                                   len_type dim1, iter_type interval00,
                                   iter_type interval01, iter_type interval10,
                                   iter_type interval11) {
  iter_type intervals[2][2]{{interval00, interval01}, {interval10, interval11}};
  iter_type dimensions[2]{dim0, dim1};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    values[n] = e->eval(n);
  }
}

template <typename T, typename E, typename... Us>
__global__ void evaluateCudaExpr3d(E *e, T *values, len_type dim0,
                                   len_type dim1, len_type dim2,
                                   iter_type interval00, iter_type interval01,
                                   iter_type interval10, iter_type interval11,
                                   iter_type interval20, iter_type interval21) {
  iter_type intervals[3][2]{{interval00, interval01},
                            {interval10, interval11},
                            {interval20, interval21}};
  iter_type dimensions[3]{dim0, dim1, dim2};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    values[n] = e->eval(n);
  }
}

template <typename T, typename E, typename... Us>
__global__ void evaluateCudaExpr1d(E *e, T *values, len_type dim0,
                                   iter_type interval00, iter_type interval01,
                                   len_type region_dim0,
                                   iter_type region_interval00,
                                   iter_type region_interval01,
                                   len_type boundary_size) {
  iter_type intervals[1][2]{{interval00, interval01}};
  iter_type region_intervals[1][2]{{region_interval00, region_interval01}};
  iter_type dimensions[]{dim0};
  iter_type region_dimensions[]{region_dim0};
  auto xi = threadIdx.x + blockIdx.x * blockDim.x;

  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }

  if (xi < intervals[0][1] - intervals[0][0]) {
    unsigned pos[]{xi + intervals[0][0]};
    iter_type n = intervals[0][0] - region_intervals[0][0] + xi;
    values[n] = e->eval(position_type<1>(dimensions, pos));
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExpr2d(
    E *e, T *values, len_type dim0, len_type dim1, iter_type interval00,
    iter_type interval01, iter_type interval10, iter_type interval11,
    len_type region_dim0, len_type region_dim1, iter_type region_interval00,
    iter_type region_interval01, iter_type region_interval10,
    iter_type region_interval11, len_type boundary_size) {
  iter_type intervals[2][2]{{interval00, interval01}, {interval10, interval11}};
  iter_type region_intervals[2][2]{{region_interval00, region_interval01},
                                   {region_interval10, region_interval11}};
  iter_type dimensions[]{dim0, dim1};
  iter_type region_dimensions[]{region_dim0, region_dim1};

  // normalize
  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }
  if (interval10 < region_interval10) {
    region_intervals[1][0] -= dim1 - 2 * boundary_size;
    region_intervals[1][1] -= dim1 - 2 * boundary_size;
  }

  auto xi = threadIdx.x + blockIdx.x * blockDim.x;
  auto yi = threadIdx.y + blockIdx.y * blockDim.y;

  unsigned pos0[]{xi, yi};
  if (xi < intervals[0][1] - intervals[0][0] &&
      yi < intervals[1][1] - intervals[1][0]) {
    unsigned pos[]{xi + intervals[0][0], yi + intervals[1][0]};
    iter_type n = 0;
    len_type stride = 1;
    for (iter_type i = 0; i < 2; ++i) {
      iter_type delta = intervals[i][0] - region_intervals[i][0];
      n += (delta + pos0[i]) * stride;
      stride *= region_dimensions[i];
    }
    auto r = e->eval(position_type<2>(dimensions, pos));
    if (r != 0) printf("got %lf for %d %d\n", r, pos[0], pos[1]);
    values[n] = r;
  }
}

template <typename T, typename E, typename... Us>
__global__ void evaluateCudaExpr3d(
    E *e, T *values, len_type dim0, len_type dim1, len_type dim2,
    iter_type interval00, iter_type interval01, iter_type interval10,
    iter_type interval11, iter_type interval20, iter_type interval21,
    len_type region_dim0, len_type region_dim1, len_type region_dim2,
    iter_type region_interval00, iter_type region_interval01,
    iter_type region_interval10, iter_type region_interval11,
    iter_type region_interval20, iter_type region_interval21,
    len_type boundary_size) {
  iter_type intervals[3][2]{{interval00, interval01},
                            {interval10, interval11},
                            {interval20, interval21}};
  iter_type region_intervals[3][2]{{region_interval00, region_interval01},
                                   {region_interval10, region_interval11},
                                   {region_interval20, region_interval21}};
  iter_type dimensions[]{dim0, dim1, dim2};
  iter_type region_dimensions[]{region_dim0, region_dim1, region_dim2};

  // normalize
  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }
  if (interval10 < region_interval10) {
    region_intervals[1][0] -= dim1 - 2 * boundary_size;
    region_intervals[1][1] -= dim1 - 2 * boundary_size;
  }
  if (interval20 < region_interval20) {
    region_intervals[2][0] -= dim2 - 2 * boundary_size;
    region_intervals[2][1] -= dim2 - 2 * boundary_size;
  }

  auto xi = threadIdx.x + blockIdx.x * blockDim.x;
  auto yi = threadIdx.y + blockIdx.y * blockDim.y;
  auto zi = threadIdx.z + blockIdx.z * blockDim.z;

  unsigned pos0[]{xi, yi, zi};
  if (xi < intervals[0][1] - intervals[0][0] &&
      yi < intervals[1][1] - intervals[1][0] &&
      zi < intervals[2][1] - intervals[2][0]) {
    unsigned pos[]{xi + intervals[0][0], yi + intervals[1][0],
                   zi + intervals[2][0]};

    iter_type n = 0;
    len_type stride = 1;
    for (iter_type i = 0; i < 3; ++i) {
      iter_type delta = intervals[i][0] - region_intervals[i][0];
      n += (delta + pos0[i]) * stride;
      stride *= region_dimensions[i];
    }
    values[n] = e->eval(position_type<3>(dimensions, pos));
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExpr(E *e, T *values) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    values[0] += e->eval(0);
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExpr(E *e, T *values, len_type len) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n < len) {
    values[n] += e->eval(n);
  }
}

template <typename T, typename E, typename... Us>
__global__ void accumulateCudaExpr1d(E *e, T *values, len_type dim0,
                                     iter_type interval00,
                                     iter_type interval01) {
  iter_type intervals[1][2]{{interval00, interval01}};
  iter_type dimensions[1]{dim0};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    values[n] += e->eval(n);
  }
}

template <typename T, typename E, typename... Us>
__global__ void accumulateCudaExpr2d(E *e, T *values, len_type dim0,
                                     len_type dim1, iter_type interval00,
                                     iter_type interval01, iter_type interval10,
                                     iter_type interval11) {
  iter_type intervals[2][2]{{interval00, interval01}, {interval10, interval11}};
  iter_type dimensions[2]{dim0, dim1};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    values[n] += e->eval(n);
  }
}

template <typename T, typename E, typename... Us>
__global__ void accumulateCudaExpr3d(E *e, T *values, len_type dim0,
                                     len_type dim1, len_type dim2,
                                     iter_type interval00, iter_type interval01,
                                     iter_type interval10, iter_type interval11,
                                     iter_type interval20,
                                     iter_type interval21) {
  iter_type intervals[3][2]{{interval00, interval01},
                            {interval10, interval11},
                            {interval20, interval21}};
  iter_type dimensions[3]{dim0, dim1, dim2};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    values[n] += e->eval(n);
  }
}

template <typename T, typename E, typename... Us>
__global__ void accumulateCudaExpr1d(E *e, T *values, len_type dim0,
                                     iter_type interval00, iter_type interval01,
                                     len_type region_dim0,
                                     iter_type region_interval00,
                                     iter_type region_interval01,
                                     len_type boundary_size) {
  iter_type intervals[1][2]{{interval00, interval01}};
  iter_type region_intervals[1][2]{{region_interval00, region_interval01}};
  iter_type dimensions[]{dim0};
  iter_type region_dimensions[]{region_dim0};

  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }

  auto xi = threadIdx.x + blockIdx.x * blockDim.x;

  if (xi < intervals[0][1] - intervals[0][0]) {
    unsigned pos[]{xi + intervals[0][0]};
    iter_type n = intervals[0][0] - region_intervals[0][0] + xi;
    values[n] += e->eval(position_type<1>(dimensions, pos));
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExpr2d(
    E *e, T *values, len_type dim0, len_type dim1, iter_type interval00,
    iter_type interval01, iter_type interval10, iter_type interval11,
    len_type region_dim0, len_type region_dim1, iter_type region_interval00,
    iter_type region_interval01, iter_type region_interval10,
    iter_type region_interval11, len_type boundary_size) {
  iter_type intervals[2][2]{{interval00, interval01}, {interval10, interval11}};
  iter_type region_intervals[2][2]{{region_interval00, region_interval01},
                                   {region_interval10, region_interval11}};
  iter_type dimensions[]{dim0, dim1};
  iter_type region_dimensions[]{region_dim0, region_dim1};

  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }
  if (interval10 < region_interval10) {
    region_intervals[1][0] -= dim1 - 2 * boundary_size;
    region_intervals[1][1] -= dim1 - 2 * boundary_size;
  }

  auto xi = threadIdx.x + blockIdx.x * blockDim.x;
  auto yi = threadIdx.y + blockIdx.y * blockDim.y;

  unsigned pos0[]{xi, yi};
  if (xi < intervals[0][1] - intervals[0][0] &&
      yi < intervals[1][1] - intervals[1][0]) {
    unsigned pos[]{xi + intervals[0][0], yi + intervals[1][0]};
    iter_type n = 0;
    len_type stride = 1;
    for (iter_type i = 0; i < 2; ++i) {
      iter_type delta = intervals[i][0] - region_intervals[i][0];
      n += (delta + pos0[i]) * stride;
      stride *= region_dimensions[i];
    }
    auto r = e->eval(position_type<2>(dimensions, pos));
    // printf("got %lf for %d %d\n", r, pos[0], pos[1]);
    values[n] += r;
  }
}

template <typename T, typename E, typename... Us>
__global__ void accumulateCudaExpr3d(
    E *e, T *values, len_type dim0, len_type dim1, len_type dim2,
    iter_type interval00, iter_type interval01, iter_type interval10,
    iter_type interval11, iter_type interval20, iter_type interval21,
    len_type region_dim0, len_type region_dim1, len_type region_dim2,
    iter_type region_interval00, iter_type region_interval01,
    iter_type region_interval10, iter_type region_interval11,
    iter_type region_interval20, iter_type region_interval21,
    len_type boundary_size) {
  iter_type intervals[3][2]{{interval00, interval01},
                            {interval10, interval11},
                            {interval20, interval21}};
  iter_type region_intervals[3][2]{{region_interval00, region_interval01},
                                   {region_interval10, region_interval11},
                                   {region_interval20, region_interval21}};
  iter_type dimensions[]{dim0, dim1, dim2};
  iter_type region_dimensions[]{region_dim0, region_dim1, region_dim2};

  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }
  if (interval10 < region_interval10) {
    region_intervals[1][0] -= dim1 - 2 * boundary_size;
    region_intervals[1][1] -= dim1 - 2 * boundary_size;
  }
  if (interval20 < region_interval20) {
    region_intervals[2][0] -= dim2 - 2 * boundary_size;
    region_intervals[2][1] -= dim2 - 2 * boundary_size;
  }

  auto xi = threadIdx.x + blockIdx.x * blockDim.x;
  auto yi = threadIdx.y + blockIdx.y * blockDim.y;
  auto zi = threadIdx.z + blockIdx.z * blockDim.z;

  unsigned pos0[]{xi, yi, zi};
  if (xi < intervals[0][1] - intervals[0][0] &&
      yi < intervals[1][1] - intervals[1][0] &&
      zi < intervals[2][1] - intervals[2][0]) {
    unsigned pos[]{xi + intervals[0][0], yi + intervals[1][0],
                   zi + intervals[2][0]};

    iter_type n = 0;
    len_type stride = 1;
    for (iter_type i = 0; i < 3; ++i) {
      iter_type delta = intervals[i][0] - region_intervals[i][0];
      n += (delta + pos0[i]) * stride;
      stride *= region_dimensions[i];
    }
    values[n] += e->eval(position_type<3>(dimensions, pos));
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector(E *e, T *values0) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    auto result = e->eval(0);
    values0 = result[0];
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector(E *e, T *values0, T *values1) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    auto result = e->eval(0);
    values0 = result[0];
    values1 = result[1];
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector(E *e, T *values0, T *values1,
                                       T *values2) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    auto result = e->eval(0);
    values0 = result[0];
    values1 = result[1];
    values2 = result[2];
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector(E *e, T *values0, len_type len) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    auto result = e->eval(n);
    values0 = result[0];
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector(E *e, T *values0, T *values1,
                                       len_type len) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n < len) {
    auto result = e->eval(n);
    values0 = e->eval(0);
    values1 = e->eval(0);
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector(E *e, T *values0, T *values1, T *values2,
                                       len_type len) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n < len) {
    auto result = e->eval(n);
    values0 = result[0];
    values1 = result[1];
    values2 = result[2];
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector1d(E *e, T *values0, len_type dim0,
                                         iter_type interval00,
                                         iter_type interval01) {
  iter_type intervals[1][2]{{interval00, interval01}};
  iter_type dimensions[1]{dim0};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    auto result = e->eval(n);
    values0[n] = result[0];
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector2d(E *e, T *values0, T *values1,
                                         len_type dim0, len_type dim1,
                                         iter_type interval00,
                                         iter_type interval01,
                                         iter_type interval10,
                                         iter_type interval11) {
  iter_type intervals[2][2]{{interval00, interval01}, {interval10, interval11}};
  iter_type dimensions[2]{dim0, dim1};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    auto result = e->eval(n);
    values0[n] = result[0];
    values1[n] = result[1];
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector3d(
    E *e, T *values0, T *values1, T *values2, len_type dim0, len_type dim1,
    len_type dim2, iter_type interval00, iter_type interval01,
    iter_type interval10, iter_type interval11, iter_type interval20,
    iter_type interval21) {
  iter_type intervals[3][2]{{interval00, interval01},
                            {interval10, interval11},
                            {interval20, interval21}};
  iter_type dimensions[3]{dim0, dim1, dim2};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    auto result = e->eval(n);
    values0[n] = result[0];
    values1[n] = result[1];
    values2[n] = result[2];
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector1d(
    E *e, T *values0, len_type dim0, iter_type interval00, iter_type interval01,
    len_type region_dim0, iter_type region_interval00,
    iter_type region_interval01, len_type boundary_size) {
  iter_type intervals[1][2]{{interval00, interval01}};
  iter_type region_intervals[1][2]{{region_interval00, region_interval01}};
  iter_type dimensions[]{dim0};
  iter_type region_dimensions[]{region_dim0};

  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }

  auto xi = threadIdx.x + blockIdx.x * blockDim.x;

  if (xi < intervals[0][1] - intervals[0][0]) {
    unsigned pos[]{xi + intervals[0][0]};

    iter_type n = 0;
    len_type stride = 1;
    iter_type delta = intervals[0][0] - region_intervals[0][0] + xi;
    n += delta;

    auto result = e->eval(position_type<1>(dimensions, pos));
    values0[n] = result[0];
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector2d(
    E *e, T *values0, T *values1, len_type dim0, len_type dim1,
    iter_type interval00, iter_type interval01, iter_type interval10,
    iter_type interval11, len_type region_dim0, len_type region_dim1,
    iter_type region_interval00, iter_type region_interval01,
    iter_type region_interval10, iter_type region_interval11,
    len_type boundary_size) {
  iter_type intervals[2][2]{{interval00, interval01}, {interval10, interval11}};
  iter_type region_intervals[2][2]{{region_interval00, region_interval01},
                                   {region_interval10, region_interval11}};
  iter_type dimensions[]{dim0, dim1};
  iter_type region_dimensions[]{region_dim0, region_dim1};

  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }
  if (interval10 < region_interval10) {
    region_intervals[1][0] -= dim1 - 2 * boundary_size;
    region_intervals[1][1] -= dim1 - 2 * boundary_size;
  }

  auto xi = threadIdx.x + blockIdx.x * blockDim.x;
  auto yi = threadIdx.y + blockIdx.y * blockDim.y;

  unsigned pos0[]{xi, yi};
  if (xi < intervals[0][1] - intervals[0][0] &&
      yi < intervals[1][1] - intervals[1][0]) {
    unsigned pos[]{xi + intervals[0][0], yi + intervals[1][0]};

    iter_type n = 0;
    len_type stride = 1;
    for (iter_type i = 0; i < 2; ++i) {
      iter_type delta = intervals[i][0] - region_intervals[i][0];
      n += (pos0[i] + delta) * stride;
      stride *= region_dimensions[i];
    }
    auto result = e->eval(position_type<2>(dimensions, pos));
    values0[n] = result[0];
    values1[n] = result[1];
  }
}

template <typename T, typename E>
__global__ void evaluateCudaExprVector3d(
    E *e, T *values0, T *values1, T *values2, len_type dim0, len_type dim1,
    len_type dim2, iter_type interval00, iter_type interval01,
    iter_type interval10, iter_type interval11, iter_type interval20,
    iter_type interval21, len_type region_dim0, len_type region_dim1,
    len_type region_dim2, iter_type region_interval00,
    iter_type region_interval01, iter_type region_interval10,
    iter_type region_interval11, iter_type region_interval20,
    iter_type region_interval21, len_type boundary_size) {
  iter_type intervals[3][2]{{interval00, interval01},
                            {interval10, interval11},
                            {interval20, interval21}};

  iter_type region_intervals[3][2]{{region_interval00, region_interval01},
                                   {region_interval10, region_interval11},
                                   {region_interval20, region_interval21}};
  iter_type dimensions[]{dim0, dim1, dim2};
  iter_type region_dimensions[]{region_dim0, region_dim1, region_dim2};

  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }
  if (interval10 < region_interval10) {
    region_intervals[1][0] -= dim1 - 2 * boundary_size;
    region_intervals[1][1] -= dim1 - 2 * boundary_size;
  }
  if (interval20 < region_interval20) {
    region_intervals[2][0] -= dim2 - 2 * boundary_size;
    region_intervals[2][1] -= dim2 - 2 * boundary_size;
  }

  auto xi = threadIdx.x + blockIdx.x * blockDim.x;
  auto yi = threadIdx.y + blockIdx.y * blockDim.y;
  auto zi = threadIdx.z + blockIdx.z * blockDim.z;

  unsigned pos0[]{xi, yi, zi};
  if (xi < intervals[0][1] - intervals[0][0] &&
      yi < intervals[1][1] - intervals[1][0] &&
      zi < intervals[2][1] - intervals[2][0]) {
    unsigned pos[]{xi + intervals[0][0], yi + intervals[1][0],
                   zi + intervals[2][0]};

    iter_type n = 0;
    len_type stride = 1;
    for (iter_type i = 0; i < 3; ++i) {
      iter_type delta = intervals[i][0] - region_intervals[i][0];
      n += (delta + pos0[i]) * stride;
      stride *= region_dimensions[i];
    }
    auto result = e->eval(position_type<3>(dimensions, pos));
    values0[n] = result[0];
    values1[n] = result[1];
    values2[n] = result[2];
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector(E *e, T *values0) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    values0 += e->eval(0)[0];
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector(E *e, T *values0, T *values1) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    auto result = e->eval(0);
    values0 += result[0];
    values1 += result[1];
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector(E *e, T *values0, T *values1,
                                         T *values2) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    auto result = e->eval(0);
    values0 += result[0];
    values1 += result[1];
    values2 += result[2];
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector(E *e, T *values0, len_type len) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    auto result = e->eval(n);
    values0 += result[0];
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector(E *e, T *values0, T *values1,
                                         len_type len) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n < len) {
    auto result = e->eval(n);
    values0 += e->eval(0);
    values1 += e->eval(0);
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector(E *e, T *values0, T *values1,
                                         T *values2, len_type len) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n < len) {
    auto result = e->eval(n);
    values0 += result[0];
    values1 += result[1];
    values2 += result[2];
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector1d(E *e, T *values0, len_type dim0,
                                           iter_type interval00,
                                           iter_type interval01) {
  iter_type intervals[1][2]{{interval00, interval01}};
  iter_type dimensions[1]{dim0};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    auto result = e->eval(n);
    values0[n] += result[0];
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector2d(E *e, T *values0, T *values1,
                                           len_type dim0, len_type dim1,
                                           iter_type interval00,
                                           iter_type interval01,
                                           iter_type interval10,
                                           iter_type interval11) {
  iter_type intervals[2][2]{{interval00, interval01}, {interval10, interval11}};
  iter_type dimensions[2]{dim0, dim1};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    auto result = e->eval(n);
    values0[n] += result[0];
    values1[n] += result[1];
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector3d(
    E *e, T *values0, T *values1, T *values2, len_type dim0, len_type dim1,
    len_type dim2, iter_type interval00, iter_type interval01,
    iter_type interval10, iter_type interval11, iter_type interval20,
    iter_type interval21) {
  iter_type intervals[3][2]{{interval00, interval01},
                            {interval10, interval11},
                            {interval20, interval21}};
  iter_type dimensions[3]{dim0, dim1, dim2};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (is_in_region(intervals, dimensions, n)) {
    auto result = e->eval(n);
    values0[n] += result[0];
    values1[n] += result[1];
    values2[n] += result[2];
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector1d(
    E *e, T *values0, len_type dim0, iter_type interval00, iter_type interval01,
    len_type region_dim0, iter_type region_interval00,
    iter_type region_interval01, len_type boundary_size) {
  iter_type intervals[1][2]{{interval00, interval01}};
  iter_type region_intervals[1][2]{{region_interval00, region_interval01}};
  iter_type dimensions[]{dim0};
  iter_type region_dimensions[]{region_dim0};

  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }
  auto xi = threadIdx.x + blockIdx.x * blockDim.x;

  if (xi < intervals[0][1] - intervals[0][0]) {
    unsigned pos[]{xi + intervals[0][0]};

    iter_type n = 0;
    len_type stride = 1;
    iter_type delta = intervals[0][0] - region_intervals[0][0];
    n += delta;

    auto result = e->eval(position_type<1>(dimensions, pos));
    values0[n] += result[0];
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector2d(
    E *e, T *values0, T *values1, len_type dim0, len_type dim1,
    iter_type interval00, iter_type interval01, iter_type interval10,
    iter_type interval11, len_type region_dim0, len_type region_dim1,
    iter_type region_interval00, iter_type region_interval01,
    iter_type region_interval10, iter_type region_interval11,
    len_type boundary_size) {
  iter_type intervals[2][2]{{interval00, interval01}, {interval10, interval11}};
  iter_type region_intervals[2][2]{{region_interval00, region_interval01},
                                   {region_interval10, region_interval11}};
  iter_type dimensions[]{dim0, dim1};
  iter_type region_dimensions[]{region_dim0, region_dim1};

  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }
  if (interval10 < region_interval10) {
    region_intervals[1][0] -= dim1 - 2 * boundary_size;
    region_intervals[1][1] -= dim1 - 2 * boundary_size;
  }

  auto xi = threadIdx.x + blockIdx.x * blockDim.x;
  auto yi = threadIdx.y + blockIdx.y * blockDim.y;

  if (xi < intervals[0][1] - intervals[0][0] &&
      yi < intervals[1][1] - intervals[1][0]) {
    unsigned pos[]{xi + intervals[0][0], yi + intervals[1][0]};

    iter_type n = 0;
    len_type stride = 1;
    for (iter_type i = 0; i < 2; ++i) {
      iter_type delta = intervals[i][0] - region_intervals[i][0];
      n += delta * stride;
      stride *= region_dimensions[i];
    }
    auto result = e->eval(position_type<2>(dimensions, pos));
    values0[n] += result[0];
    values1[n] += result[1];
  }
}

template <typename T, typename E>
__global__ void accumulateCudaExprVector3d(
    E *e, T *values0, T *values1, T *values2, len_type dim0, len_type dim1,
    len_type dim2, iter_type interval00, iter_type interval01,
    iter_type interval10, iter_type interval11, iter_type interval20,
    iter_type interval21, len_type region_dim0, len_type region_dim1,
    len_type region_dim2, iter_type region_interval00,
    iter_type region_interval01, iter_type region_interval10,
    iter_type region_interval11, iter_type region_interval20,
    iter_type region_interval21, len_type boundary_size) {
  iter_type intervals[3][2]{{interval00, interval01},
                            {interval10, interval11},
                            {interval20, interval21}};

  iter_type region_intervals[3][2]{{region_interval00, region_interval01},
                                   {region_interval10, region_interval11},
                                   {region_interval20, region_interval21}};
  iter_type dimensions[]{dim0, dim1, dim2};
  iter_type region_dimensions[]{region_dim0, region_dim1, region_dim2};

  if (interval00 < region_interval00) {
    region_intervals[0][0] -= dim0 - 2 * boundary_size;
    region_intervals[0][1] -= dim0 - 2 * boundary_size;
  }
  if (interval10 < region_interval10) {
    region_intervals[1][0] -= dim1 - 2 * boundary_size;
    region_intervals[1][1] -= dim1 - 2 * boundary_size;
  }
  if (interval20 < region_interval20) {
    region_intervals[2][0] -= dim2 - 2 * boundary_size;
    region_intervals[2][1] -= dim2 - 2 * boundary_size;
  }

  auto xi = threadIdx.x + blockIdx.x * blockDim.x;
  auto yi = threadIdx.y + blockIdx.y * blockDim.y;
  auto zi = threadIdx.z + blockIdx.z * blockDim.z;

  if (xi < intervals[0][1] - intervals[0][0] &&
      yi < intervals[1][1] - intervals[1][0] &&
      zi < intervals[2][1] - intervals[2][0]) {
    unsigned pos[]{xi + intervals[0][0], yi + intervals[1][0],
                   zi + intervals[2][0]};

    iter_type n = 0;
    len_type stride = 1;
    for (iter_type i = 0; i < 3; ++i) {
      iter_type delta = intervals[i][0] - region_intervals[i][0];
      n += delta * stride;
      stride *= region_dimensions[i];
    }
    auto result = e->eval(position_type<3>(dimensions, pos));
    values0[n] += result[0];
    values1[n] += result[1];
    values2[n] += result[2];
  }
}

template <typename grid_type>
struct eval_submit {
  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<0> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExpr CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values));
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<1> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExpr1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.intervals[0][0], interval.intervals[0][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<2> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExpr2d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.dims[1], interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<3> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExpr3d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.dims[1], interval.dims[2], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<1> const &interval,
            grid::region_interval<1> const &region_interval,
            len_type boundary_size) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExpr1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.intervals[0][0], interval.intervals[0][1],
        region_interval.dims[0], region_interval.intervals[0][0],
        region_interval.intervals[0][1], boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<2> const &interval,
            grid::region_interval<2> const &region_interval,
            len_type boundary_size) {
    dim3 blockDim(16, 16);
    dim3 gridDim(
        (interval.intervals[0][1] - interval.intervals[0][0] + blockDim.x - 1) /
            blockDim.x,
        (interval.intervals[1][1] - interval.intervals[1][0] + blockDim.y - 1) /
            blockDim.y);

    evaluateCudaExpr2d CUDA_KERNEL(gridDim, blockDim)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.dims[1], interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1],
        region_interval.dims[0], region_interval.dims[1],
        region_interval.intervals[0][0], region_interval.intervals[0][1],
        region_interval.intervals[1][0], region_interval.intervals[1][1],
        boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<3> const &interval,
            grid::region_interval<3> const &region_interval,
            len_type boundary_size) {
    dim3 blockDim(8, 8, 8);
    dim3 gridDim(
        (interval.intervals[0][1] - interval.intervals[0][0] + blockDim.x - 1) /
            blockDim.x,
        (interval.intervals[1][1] - interval.intervals[1][0] + blockDim.y - 1) /
            blockDim.y,
        (interval.intervals[2][1] - interval.intervals[2][0] + blockDim.z - 1) /
            blockDim.z);

    evaluateCudaExpr3d CUDA_KERNEL(gridDim, blockDim)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.dims[1], interval.dims[2], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1], region_interval.dims[0],
        region_interval.dims[1], region_interval.dims[2],
        region_interval.intervals[0][0], region_interval.intervals[0][1],
        region_interval.intervals[1][0], region_interval.intervals[1][1],
        region_interval.intervals[2][0], region_interval.intervals[2][1],
        boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values, len_type len0) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExpr CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), len0);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values, len_type len0) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExpr CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), len0);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<0> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExpr CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values));
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<1> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExpr1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.intervals[0][0], interval.intervals[0][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<2> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExpr2d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.dims[1], interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<3> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExpr3d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.dims[1], interval.dims[2], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<1> const &interval,
           grid::region_interval<1> const &region_interval,
           len_type boundary_size) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExpr1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.intervals[0][0], interval.intervals[0][1],
        region_interval.dims[0], region_interval.intervals[0][0],
        region_interval.intervals[0][1], boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<2> const &interval,
           grid::region_interval<2> const &region_interval,
           len_type boundary_size) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;

    dim3 blockDim(16, 16);
    dim3 gridDim(
        (interval.intervals[0][1] - interval.intervals[0][0] + blockDim.x - 1) /
            blockDim.x,
        (interval.intervals[1][1] - interval.intervals[1][0] + blockDim.y - 1) /
            blockDim.y);

    accumulateCudaExpr2d CUDA_KERNEL(gridDim, blockDim)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.dims[1], interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1],
        region_interval.dims[0], region_interval.dims[1],
        region_interval.intervals[0][0], region_interval.intervals[0][1],
        region_interval.intervals[1][0], region_interval.intervals[1][1],
        boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<3> const &interval,
           grid::region_interval<3> const &region_interval,
           len_type boundary_size) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;

    dim3 blockDim(8, 8, 8);
    dim3 gridDim(
        (interval.intervals[0][1] - interval.intervals[0][0] + blockDim.x - 1) /
            blockDim.x,
        (interval.intervals[1][1] - interval.intervals[1][0] + blockDim.y - 1) /
            blockDim.y,
        (interval.intervals[2][1] - interval.intervals[2][0] + blockDim.z - 1) /
            blockDim.z);

    accumulateCudaExpr3d CUDA_KERNEL(gridDim, blockDim)(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.dims[1], interval.dims[2], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1], region_interval.dims[0],
        region_interval.dims[1], region_interval.dims[2],
        region_interval.intervals[0][0], region_interval.intervals[0][1],
        region_interval.intervals[1][0], region_interval.intervals[1][1],
        region_interval.intervals[2][0], region_interval.intervals[2][1],
        boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }
};

template <typename grid_type, size_t dimension>
struct eval_submit<any_vector_t<grid_type, dimension>> {
  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<0> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExpr1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values)[0]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<1> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExprVector1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values)[0], interval.dims[0],
        interval.intervals[0][0], interval.intervals[0][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<2> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExprVector2d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values)[0],
        std::forward<assign_type>(values)[1], interval.dims[0],
        interval.dims[1], interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<3> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExprVector3d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values)[0],
        std::forward<assign_type>(values)[1],
        std::forward<assign_type>(values)[2], interval.dims[0],
        interval.dims[1], interval.dims[2], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<1> const &interval,
            grid::region_interval<1> const &region_interval,
            len_type boundary_size) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExprVector1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values)[0], interval.dims[0],
        interval.intervals[0][0], interval.intervals[0][1],
        region_interval.dims[0], region_interval.intervals[0][0],
        region_interval.intervals[0][1], boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<2> const &interval,
            grid::region_interval<2> const &region_interval,
            len_type boundary_size) {
    dim3 blockDim(16, 16);
    dim3 gridDim(
        (interval.intervals[0][1] - interval.intervals[0][0] + blockDim.x - 1) /
            blockDim.x,
        (interval.intervals[1][1] - interval.intervals[1][0] + blockDim.y - 1) /
            blockDim.y);

    evaluateCudaExprVector2d CUDA_KERNEL(gridDim, blockDim)(
        e, std::forward<assign_type>(values)[0],
        std::forward<assign_type>(values)[1], interval.dims[0],
        interval.dims[1], interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1],
        region_interval.dims[0], region_interval.dims[1],
        region_interval.intervals[0][0], region_interval.intervals[0][1],
        region_interval.intervals[1][0], region_interval.intervals[1][1],
        boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<3> const &interval,
            grid::region_interval<3> const &region_interval,
            len_type boundary_size) {
    dim3 blockDim(8, 8, 8);
    dim3 gridDim(
        (interval.intervals[0][1] - interval.intervals[0][0] + blockDim.x - 1) /
            blockDim.x,
        (interval.intervals[1][1] - interval.intervals[1][0] + blockDim.y - 1) /
            blockDim.y,
        (interval.intervals[2][1] - interval.intervals[2][0] + blockDim.z - 1) /
            blockDim.z);

    evaluateCudaExprVector3d CUDA_KERNEL(gridDim, blockDim)(
        e, std::forward<assign_type>(values)[0],
        std::forward<assign_type>(values)[1],
        std::forward<assign_type>(values)[2], interval.dims[0],
        interval.dims[1], interval.dims[2], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1], region_interval.dims[0],
        region_interval.dims[1], region_interval.dims[2],
        region_interval.intervals[0][0], region_interval.intervals[0][1],
        region_interval.intervals[1][0], region_interval.intervals[1][1],
        region_interval.intervals[2][0], region_interval.intervals[2][1],
        boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values, len_type len0) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    evaluateCudaExprVector CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), len0);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values, len_type len0) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExprVector CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), len0);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<0> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExpr1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values)[0]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<1> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExprVector1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values)[0], interval.dims[0],
        interval.intervals[0][0], interval.intervals[0][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<2> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExprVector2d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values)[0],
        std::forward<assign_type>(values)[1], interval.dims[0],
        interval.dims[1], interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<3> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExprVector3d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values)[0],
        std::forward<assign_type>(values)[1],
        std::forward<assign_type>(values)[2], interval.dims[0],
        interval.dims[1], interval.dims[2], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<1> const &interval,
           grid::region_interval<1> const &region_interval,
           len_type boundary_size) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    accumulateCudaExprVector1d CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values)[0], interval.dims[0],
        interval.intervals[0][0], interval.intervals[0][1],
        region_interval.dims[0], region_interval.intervals[0][0],
        region_interval.intervals[0][1], boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<2> const &interval,
           grid::region_interval<2> const &region_interval,
           len_type boundary_size) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;

    dim3 blockDim(16, 16);
    dim3 gridDim(
        (interval.intervals[0][1] - interval.intervals[0][0] + blockDim.x - 1) /
            blockDim.x,
        (interval.intervals[1][1] - interval.intervals[1][0] + blockDim.y - 1) /
            blockDim.y);

    accumulateCudaExprVector2d CUDA_KERNEL(gridDim, blockDim)(
        e, std::forward<assign_type>(values)[0],
        std::forward<assign_type>(values)[1], interval.dims[0],
        interval.dims[1], interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1],
        region_interval.dims[0], region_interval.dims[1],
        region_interval.intervals[0][0], region_interval.intervals[0][1],
        region_interval.intervals[1][0], region_interval.intervals[1][1],
        boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void acc(len_type len, E *e, assign_type &&values,
           grid::region_interval<3> const &interval,
           grid::region_interval<3> const &region_interval,
           len_type boundary_size) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;

    dim3 blockDim(8, 8, 8);
    dim3 gridDim(
        (interval.intervals[0][1] - interval.intervals[0][0] + blockDim.x - 1) /
            blockDim.x,
        (interval.intervals[1][1] - interval.intervals[1][0] + blockDim.y - 1) /
            blockDim.y,
        (interval.intervals[2][1] - interval.intervals[2][0] + blockDim.z - 1) /
            blockDim.z);

    accumulateCudaExprVector3d CUDA_KERNEL(gridDim, blockDim)(
        e, std::forward<assign_type>(values)[0],
        std::forward<assign_type>(values)[1],
        std::forward<assign_type>(values)[2], interval.dims[0],
        interval.dims[1], interval.dims[2], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1], region_interval.dims[0],
        region_interval.dims[1], region_interval.dims[2],
        region_interval.intervals[0][0], region_interval.intervals[0][1],
        region_interval.intervals[1][0], region_interval.intervals[1][1],
        region_interval.intervals[2][0], region_interval.intervals[2][1],
        boundary_size);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }
};

template <typename T, typename E>
__global__ void sumCudaExpr(E *e, T *output) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    output[0] = e->eval(0);
  }
}

template <typename T, typename E>
__global__ void sumCudaExpr(E *e, T *output, len_type len) {
  extern __shared__ T sdata[];
  auto n = threadIdx.x + blockIdx.x * blockDim.x;

  // Each thread loads one element from global to shared memory
  unsigned int tid = threadIdx.x;
  sdata[tid] = (n < len) ? e->eval(n) : T{};
  __syncthreads();

  // Do reduction in shared memory
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdata[tid] += sdata[tid + s];
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    output[blockIdx.x] = sdata[0];
  }
}

template <typename T, typename E>
__global__ void sumCudaExpr1d(E *e, T *output, len_type dim0,
                              iter_type interval00, iter_type interval01) {
  extern __shared__ T sdata[];
  iter_type intervals[1][2]{{interval00, interval01}};
  iter_type dimensions[1]{dim0};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;

  // Each thread loads one element from global to shared memory
  unsigned int tid = threadIdx.x;
  sdata[tid] = (is_in_region(intervals, dimensions, n)) ? e->eval(n) : T{};
  __syncthreads();

  // Do reduction in shared memory
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdata[tid] += sdata[tid + s];
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    output[blockIdx.x] = sdata[0];
  }
}

template <typename T, typename E>
__global__ void sumCudaExpr2d(E *e, T *output, len_type dim0, len_type dim1,
                              iter_type interval00, iter_type interval01,
                              iter_type interval10, iter_type interval11) {
  extern __shared__ T sdata0[];

  iter_type intervals[2][2]{{interval00, interval01}, {interval10, interval11}};
  iter_type dimensions[2]{dim0, dim1};

  auto n = threadIdx.x + blockIdx.x * blockDim.x;

  // Each thread loads one element from global to shared memory
  unsigned int tid = threadIdx.x;
  sdata0[tid] = (is_in_region(intervals, dimensions, n)) ? e->eval(n) : T{};
  __syncthreads();

  // Do reduction in shared memory
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdata0[tid] += sdata0[tid + s];
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    output[blockIdx.x] = sdata0[0];
  }
}

template <typename T, typename E>
__global__ void sumCudaExpr3d(E *e, T *output, len_type dim0, len_type dim1,
                              len_type dim2, iter_type interval00,
                              iter_type interval01, iter_type interval10,
                              iter_type interval11, iter_type interval20,
                              iter_type interval21) {
  extern __shared__ T sdata0[];
  iter_type intervals[3][2]{{interval00, interval01},
                            {interval10, interval11},
                            {interval20, interval21}};
  iter_type dimensions[3]{dim0, dim1, dim2};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;

  // Each thread loads one element from global to shared memory
  unsigned int tid = threadIdx.x;
  sdata0[tid] = (is_in_region(intervals, dimensions, n)) ? e->eval(n) : T{};
  __syncthreads();

  // Do reduction in shared memory
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdata0[tid] += sdata0[tid + s];
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    output[blockIdx.x] = sdata0[0];
  }
}

template <typename T, typename E>
__global__ void sumCudaExprVector(E *e, T *output0) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    auto result = e->eval(0);
    output0[0] = result[0];
  }
}

template <typename T, typename E>
__global__ void sumCudaExprVector(E *e, T *output0, T *output1) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    auto result = e->eval(0);
    output0[0] = result[0];
    output1[0] = result[1];
  }
}

template <typename T, typename E>
__global__ void sumCudaExprVector(E *e, T *output0, T *output1, T *output2) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    auto result = e->eval(0);
    output0[0] = result[0];
    output1[0] = result[1];
    output2[0] = result[2];
  }
}

template <typename T, typename E>
__global__ void sumCudaExprVector(E *e, T *output0, len_type len) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n == 0) {
    auto result = e->eval(n);
    output0[n] = result[0];
  }
}

template <typename T, typename E>
__global__ void sumCudaExprVector(E *e, T *output0, T *output1, len_type len) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n < len) {
    auto result = e->eval(n);
    output0[n] = e->eval(0);
    output1[n] = e->eval(0);
  }
}

template <typename T, typename E>
__global__ void sumCudaExprVector(E *e, T *output0, T *output1, T *output2,
                                  len_type len) {
  auto n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n < len) {
    auto result = e->eval(n);
    output0[n] = result[0];
    output1[n] = result[1];
    output2[n] = result[2];
  }
}

template <typename T, typename E>
__global__ void sumCudaExprVector1d(E *e, T *output0, len_type dim0,
                                    iter_type interval00,
                                    iter_type interval01) {
  iter_type intervals[1][2]{{interval00, interval01}};
  iter_type dimensions[1]{dim0};

  auto n = threadIdx.x + blockIdx.x * blockDim.x;

  extern __shared__ T sdata[];
  T *sdatax = &sdata[0];

  // Each thread loads one element from global to shared memory
  unsigned int tid = threadIdx.x;
  if (is_in_region(intervals, dimensions, n)) {
    auto result = e->eval(n);
    sdatax[tid] = result[0];
  } else {
    sdatax[tid] = T{};
  }
  __syncthreads();

  // Do reduction in shared memory
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdatax[tid] += sdatax[tid + s];
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    output0[blockIdx.x] = sdatax[0];
  }
}

template <typename T, typename E>
__global__ void sumCudaExprVector2d(E *e, T *output0, T *output1, len_type dim0,
                                    len_type dim1, iter_type interval00,
                                    iter_type interval01, iter_type interval10,
                                    iter_type interval11) {
  iter_type intervals[2][2]{{interval00, interval01}, {interval10, interval11}};
  iter_type dimensions[2]{dim0, dim1};

  auto n = threadIdx.x + blockIdx.x * blockDim.x;

  extern __shared__ T sdata0[];
  T *sdatax = &sdata0[0];
  T *sdatay = &sdatax[blockDim.x];

  // Each thread loads one element from global to shared memory
  unsigned int tid = threadIdx.x;
  if (is_in_region(intervals, dimensions, n)) {
    auto result = e->eval(n);
    sdatax[tid] = result[0];
    sdatay[tid] = result[1];
  } else {
    sdatax[tid] = T{};
    sdatay[tid] = T{};
  }
  __syncthreads();

  // Do reduction in shared memory
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdatax[tid] += sdatax[tid + s];
      sdatay[tid] += sdatay[tid + s];
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    output0[blockIdx.x] = sdatax[0];
    output1[blockIdx.x] = sdatay[0];
  }
}

template <typename T, typename E>
__global__ void sumCudaExprVector3d(E *e, T *output0, T *output1, T *output2,
                                    len_type dim0, len_type dim1, len_type dim2,
                                    iter_type interval00, iter_type interval01,
                                    iter_type interval10, iter_type interval11,
                                    iter_type interval20,
                                    iter_type interval21) {
  iter_type intervals[3][2]{{interval00, interval01},
                            {interval10, interval11},
                            {interval20, interval21}};
  iter_type dimensions[3]{dim0, dim1, dim2};
  auto n = threadIdx.x + blockIdx.x * blockDim.x;

  extern __shared__ T sdata[];
  T *sdatax = &sdata[0];
  T *sdatay = &sdatax[blockDim.x];
  T *sdataz = &sdatay[blockDim.x];

  // Each thread loads one element from global to shared memory
  unsigned int tid = threadIdx.x;
  if (is_in_region(intervals, dimensions, n)) {
    auto result = e->eval(n);
    sdatax[tid] = result[0];
    sdatay[tid] = result[1];
    sdataz[tid] = result[2];
  } else {
    sdatax[tid] = T{};
    sdatay[tid] = T{};
    sdataz[tid] = T{};
  }
  __syncthreads();

  // Do reduction in shared memory
  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      sdatax[tid] += sdatax[tid + s];
      sdatay[tid] += sdatay[tid + s];
      sdataz[tid] += sdataz[tid + s];
    }
    __syncthreads();
  }

  // Write the result for this block to global memory
  if (tid == 0) {
    output0[blockIdx.x] = sdatax[0];
    output1[blockIdx.x] = sdatay[0];
    output2[blockIdx.x] = sdataz[0];
  }
}

template <typename grid_type>
struct sum_submit {
  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<0> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    sumCudaExpr CUDA_KERNEL(numBlocks, BLOCK_SIZE,
                            BLOCK_SIZE * sizeof(assign_type))(
        e, std::forward<assign_type>(values));
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<1> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    sumCudaExpr1d CUDA_KERNEL(numBlocks, BLOCK_SIZE,
                              BLOCK_SIZE * sizeof(assign_type))(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.intervals[0][0], interval.intervals[0][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<2> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    sumCudaExpr2d CUDA_KERNEL(numBlocks, BLOCK_SIZE,
                              BLOCK_SIZE * sizeof(assign_type))(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.dims[1], interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<3> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    sumCudaExpr3d CUDA_KERNEL(numBlocks, BLOCK_SIZE,
                              BLOCK_SIZE * sizeof(assign_type))(
        e, std::forward<assign_type>(values), interval.dims[0],
        interval.dims[1], interval.dims[2], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type, typename interval_type>
  void eval(len_type len, E *e, assign_type &&values,
            interval_type const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    sumCudaExpr CUDA_KERNEL(numBlocks, BLOCK_SIZE,
                            BLOCK_SIZE * sizeof(assign_type))(
        e, std::forward<assign_type>(values), len);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }
};

template <typename grid_type, size_t dimension>
struct sum_submit<any_vector_t<grid_type, dimension>> {
  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<0> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    sumCudaExpr1d CUDA_KERNEL(numBlocks, BLOCK_SIZE,
                              BLOCK_SIZE * sizeof(assign_type))(
        e, std::forward<assign_type>(values)[0]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<1> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    sumCudaExprVector1d CUDA_KERNEL(numBlocks, BLOCK_SIZE,
                                    BLOCK_SIZE * sizeof(assign_type))(
        e, std::forward<assign_type>(values)[0], interval.dims[0],
        interval.intervals[0][0], interval.intervals[0][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<2> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    sumCudaExprVector2d CUDA_KERNEL(numBlocks, BLOCK_SIZE,
                                    2 * BLOCK_SIZE * sizeof(assign_type))(
        e, std::forward<assign_type>(values)[0],
        std::forward<assign_type>(values)[1], interval.dims[0],
        interval.dims[1], interval.intervals[0][0], interval.intervals[0][1],
        interval.intervals[1][0], interval.intervals[1][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type>
  void eval(len_type len, E *e, assign_type &&values,
            grid::region_interval<3> const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    sumCudaExprVector3d CUDA_KERNEL(numBlocks, BLOCK_SIZE,
                                    3 * BLOCK_SIZE * sizeof(assign_type))(
        e, std::forward<assign_type>(values)[0],
        std::forward<assign_type>(values)[1],
        std::forward<assign_type>(values)[2], interval.dims[0],
        interval.dims[1], interval.dims[2], interval.intervals[0][0],
        interval.intervals[0][1], interval.intervals[1][0],
        interval.intervals[1][1], interval.intervals[2][0],
        interval.intervals[2][1]);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }

  template <typename E, typename assign_type, typename interval_type>
  void eval(len_type len, E *e, assign_type &&values,
            interval_type const &interval) {
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    sumCudaExprVector CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        e, std::forward<assign_type>(values), interval);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }
};

template <typename T>
struct inspect_result_for {
  template <typename E, size_t D>
  auto operator()(OpEvaluable<E> const &e,
                  grid::region_interval<D> const &interval) {
    len_type len = grid::length<D>(interval.dims);
    T *values;
    CHECK_CUDA_ERROR(cudaMalloc(&values, len * sizeof(T)));

    auto ec = expr::to_cuda_expr(*static_cast<E const *>(&e));
    using cuda_expr_type = decltype(ec);
    cuda_expr_type *ecDevPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&ecDevPtr, sizeof(cuda_expr_type)));
    CHECK_CUDA_ERROR(cudaMemcpy(ecDevPtr, &ec, sizeof(cuda_expr_type),
                                cudaMemcpyHostToDevice));

    eval_submit<T>{}.eval(len, ecDevPtr, values, interval);

    CHECK_CUDA_ERROR(cudaFree(ecDevPtr));
    T *hostValues = new T[len];
    CHECK_CUDA_ERROR(cudaMemcpy(hostValues, values, len * sizeof(T),
                                cudaMemcpyDeviceToHost));

    delete[] hostValues;
    CHECK_CUDA_ERROR(cudaFree(values));
  }
};

template <typename T, size_t D>
struct inspect_result_for<any_vector_t<T, D>> {
  template <typename E>
  auto operator()(OpEvaluable<E> const &e,
                  grid::region_interval<D> const &interval) {
    len_type len = grid::length<D>(interval.dims);
    T *values[D];
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMalloc(&values[i], len * sizeof(T)));
    }

    auto ec = expr::to_cuda_expr(*static_cast<E const *>(&e));
    using cuda_expr_type = decltype(ec);
    cuda_expr_type *ecDevPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&ecDevPtr, sizeof(cuda_expr_type)));
    CHECK_CUDA_ERROR(cudaMemcpy(ecDevPtr, &ec, sizeof(cuda_expr_type),
                                cudaMemcpyHostToDevice));

    eval_submit<any_vector_t<T, D>>{}.eval(len, ecDevPtr, values, interval);

    CHECK_CUDA_ERROR(cudaFree(ecDevPtr));
    T *hostValues[D];
    for (iter_type i = 0; i < D; ++i) {
      hostValues[i] = new T[len];
      CHECK_CUDA_ERROR(cudaMemcpy(hostValues[i], values[i], len * sizeof(T),
                                  cudaMemcpyDeviceToHost));
    }

    for (iter_type i = 0; i < D; ++i) {
      delete[] hostValues[i];
      CHECK_CUDA_ERROR(cudaFree(values[i]));
    }
  }
};

template <typename E, size_t D>
void inspect_result(OpEvaluable<E> const &e,
                    grid::region_interval<D> const &interval) {
  using eval_type = expr::eval_type_t<E>;
  inspect_result_for<eval_type>{}(*static_cast<E const *>(&e), interval);
}

template <typename A, typename B, size_t D>
void inspect_result(OpBinaryMul<A, B> const &mul,
                    grid::region_interval<D> const &interval) {
  inspect_result(mul.a, interval);
  inspect_result(mul.b, interval);
  using eval_type = expr::eval_type_t<OpBinaryMul<A, B>>;
  inspect_result_for<eval_type>{}(mul, interval);
}

template <typename... Es, size_t D, size_t... Is>
void inspect_result(OpAdd<Es...> const &add,
                    grid::region_interval<D> const &interval,
                    std::index_sequence<Is...>) {
  ((inspect_result(expr::get<Is>(add), interval), ...));
  using eval_type = expr::eval_type_t<OpAdd<Es...>>;
  inspect_result_for<eval_type>{}(add, interval);
}

template <typename... Es, size_t D>
void inspect_result(OpAdd<Es...> const &add,
                    grid::region_interval<D> const &interval) {
  inspect_result(add, interval, std::make_index_sequence<sizeof...(Es)>{});
}

template <typename T, size_t D>
struct expr::evaluate_expression_trait<GridCUDA<T, D>> {
  template <typename E, typename assign_type, typename interval_type>
  auto submit_expr(len_type len, E const &ec, assign_type &&values,
                   interval_type &&interval) {
    E *ecDevPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&ecDevPtr, sizeof(E)));
    CHECK_CUDA_ERROR(
        cudaMemcpy(ecDevPtr, &ec, sizeof(E), cudaMemcpyHostToDevice));

    eval_submit<T>{}.eval(len, ecDevPtr, std::forward<assign_type>(values),
                          std::forward<interval_type>(interval));
    CHECK_CUDA_ERROR(cudaFree(ecDevPtr));
  }

  template <typename E, typename assign_type>
  auto submit_expr(len_type len, E const &ec, assign_type &&values,
                   grid::region_interval<D> const &interval,
                   grid::region_interval<D> const &region_interval,
                   len_type boundary_size) {
    E *ecDevPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&ecDevPtr, sizeof(E)));
    CHECK_CUDA_ERROR(
        cudaMemcpy(ecDevPtr, &ec, sizeof(E), cudaMemcpyHostToDevice));

    eval_submit<T>{}.eval(len, ecDevPtr, std::forward<assign_type>(values),
                          interval, region_interval, boundary_size);
    CHECK_CUDA_ERROR(cudaFree(ecDevPtr));
  }

  template <typename E>
  evaluate_expression_trait(OpEvaluable<E> const &e, GridCUDA<T, D> &data,
                            grid::region_interval<D> const &interval,
                            bool synchronize = true) {
    auto ec = to_cuda_expr(*static_cast<E const *>(&e));
    iter_type len = 1;
    for (iter_type i = 0; i < D; ++i) {
      len *= interval.dims[i];
    }

    submit_expr(len, ec, data.values, interval);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename E>
  evaluate_expression_trait(OpEvaluable<E> const &e, GridCUDA<T, D> &data,
                            grid::region_interval<0> const &interval,
                            bool synchronize = true) {
    data[0] = static_cast<E const *>(&e)->eval(0);
  }

  template <typename E>
  evaluate_expression_trait(OpEvaluable<E> const &e, GridCUDA<T, D> &data,
                            len_type len, bool synchronize = true) {
    auto ec = to_cuda_expr(*static_cast<E const *>(&e));
    submit_expr(len, ec, data.values, len);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename E>
  evaluate_expression_trait(OpEvaluable<E> const &e,
                            RegionalGridCUDA<T, D> &data,
                            grid::region_interval<D> const &interval,
                            bool synchronize = true) {
    // inspect_result(*static_cast<E const *>(&e), interval);

    auto ec = to_cuda_expr(*static_cast<E const *>(&e));
    grid::region_interval<D> region_interval(data.region.origin,
                                             data.region.dims);

    iter_type len = 1;
    for (iter_type i = 0; i < D; ++i) {
      len *= interval.dims[i];
    }

    submit_expr(len, ec, data.values, interval, region_interval,
                data.region.boundary_size);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename E, typename assign_type>
  evaluate_expression_trait(equation_ptr_list_type<assign_type, E> *eq_list,
                            len_type len, bool synchronize = true) {
    for (iter_type i = 0; i < len; ++i) {
      auto *data = eq_list[i].left;
      auto interval = expr::iterable_domain(*data);
      auto ec = to_cuda_expr(*eq_list[i].right);
      submit_expr(len, ec, data->values, interval);
    }
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename E>
  evaluate_expression_trait(
      equation_ptr_list_type<RegionalGridCUDA<T, D>, E> *eq_list, len_type len,
      bool synchronize = true) {
    for (iter_type i = 0; i < len; ++i) {
      auto *data = eq_list[i].left;
      auto interval = expr::iterable_domain(*data);
      grid::region_interval<D> region_interval(data->region.origin,
                                               data->region.dims);
      auto ec = to_cuda_expr(*eq_list[i].right);
      submit_expr(len, ec, data->values, interval, region_interval,
                  data->region.boundary_size);
    }
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  evaluate_expression_trait() {}
};

template <typename T, size_t D>
struct expr::accumulate_expression_trait<GridCUDA<T, D>> {
  template <typename E, typename assign_type, typename interval_type>
  auto submit_expr(len_type len, E const &ec, assign_type &&values,
                   interval_type &&interval) {
    E *ecDevPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&ecDevPtr, sizeof ec));
    CHECK_CUDA_ERROR(
        cudaMemcpy(ecDevPtr, &ec, sizeof ec, cudaMemcpyHostToDevice));

    eval_submit<T>{}.acc(len, ecDevPtr, std::forward<assign_type>(values),
                         std::forward<interval_type>(interval));

    CHECK_CUDA_ERROR(cudaFree(ecDevPtr));
  }

  template <typename E, typename assign_type>
  auto submit_expr(len_type len, E const &ec, assign_type &&values,
                   grid::region_interval<D> const &interval,
                   grid::region_interval<D> const &region_interval,
                   len_type boundary_size) {
    E *ecDevPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&ecDevPtr, sizeof(E)));
    CHECK_CUDA_ERROR(
        cudaMemcpy(ecDevPtr, &ec, sizeof(E), cudaMemcpyHostToDevice));

    eval_submit<T>{}.acc(len, ecDevPtr, std::forward<assign_type>(values),
                         interval, region_interval, boundary_size);
    CHECK_CUDA_ERROR(cudaFree(ecDevPtr));
  }

  //! Add the result of the expression into the underlying data member.
  /*!
   * The expression is evaluated and the result is added to the existing
   * values in the data array.
   *
   * \param e Expression that is evaluated.
   * \param data The array of data.
   * \param len The length of the array.
   */
  template <typename E>
  accumulate_expression_trait(OpEvaluable<E> const &e, GridCUDA<T, D> &data,
                              len_type len, bool synchronize = true) {
    auto ec = to_cuda_expr(*static_cast<E const *>(&e));
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    submit_expr(len, ec, data.values, len);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename E>
  accumulate_expression_trait(OpEvaluable<E> const &e, GridCUDA<T, D> &data,
                              grid::region_interval<D> const &interval,
                              bool synchronize = true) {
    auto ec = to_cuda_expr(*static_cast<E const *>(&e));

    iter_type len = 1;
    for (iter_type i = 0; i < D; ++i) {
      len *= interval.intervals[i][1] - interval.intervals[i][0];
    }
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    submit_expr(len, ec, data.values, interval);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  template <typename E>
  accumulate_expression_trait(OpEvaluable<E> const &e,
                              RegionalGridCUDA<T, D> &data,
                              grid::region_interval<D> const &interval,
                              bool synchronize = true) {
    // inspect_result(*static_cast<E const *>(&e), interval);

    auto ec = to_cuda_expr(*static_cast<E const *>(&e));
    grid::region_interval<D> region_interval(data.region.origin,
                                             data.region.dims);

    iter_type len = 1;
    for (iter_type i = 0; i < D; ++i) {
      len *= interval.dims[i];
    }

    submit_expr(len, ec, data.values, interval, region_interval,
                data.region.boundary_size);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  accumulate_expression_trait() {}
};

template <typename T, size_t D>
struct expr::result_sum_trait<GridCUDA<T, D>> {
  template <typename E, typename assign_type, typename interval_type>
  void submit_expr(len_type len, E const &ec, assign_type &&output,
                   interval_type &&interval, bool synchronize = true) {
    E *ecDevPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&ecDevPtr, sizeof(E)));
    CHECK_CUDA_ERROR(
        cudaMemcpy(ecDevPtr, &ec, sizeof(E), cudaMemcpyHostToDevice));

    sum_submit<T>{}.eval(len, ecDevPtr, std::forward<assign_type>(output),
                         std::forward<interval_type>(interval));
    CHECK_CUDA_ERROR(cudaFree(ecDevPtr));
  }

  template <typename E, typename assign_type>
  void result(OpEvaluable<E> const &e, assign_type &&output,
              grid::region_interval<D> const &interval,
              bool synchronize = true) {
    // inspect_result(*static_cast<E const *>(&e), interval);

    auto ec = to_cuda_expr(*static_cast<E const *>(&e));
    iter_type len = 1;
    for (iter_type i = 0; i < D; ++i) {
      len *= interval.dims[i];
    }
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;

    T *devPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&devPtr, sizeof(T) * numBlocks));
    submit_expr(len, ec, devPtr, interval);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    CHECK_CUDA_ERROR(cudaMemcpy(result, devPtr, sizeof(T) * numBlocks,
                                cudaMemcpyDeviceToHost));
    T sum = T{};
    for (iter_type i = 0; i < numBlocks; ++i) {
      sum += result[i];
    }
    delete[] result;
    CHECK_CUDA_ERROR(cudaFree(devPtr));
    std::forward<assign_type>(output) = sum;
  }

  template <typename E, typename assign_type>
  void result(OpEvaluable<E> const &e, assign_type &&output, len_type len,
              bool synchronize = true) {
    auto ec = to_cuda_expr(*static_cast<E const *>(&e));
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    T *devPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&devPtr, sizeof(T) * numBlocks));
    submit_expr(len, ec, devPtr, len);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    CHECK_CUDA_ERROR(cudaMemcpy(result, devPtr, sizeof(T) * numBlocks,
                                cudaMemcpyDeviceToHost));
    T sum = T{};
    for (iter_type i = 0; i < numBlocks; ++i) {
      sum += result[i];
    }
    delete[] result;
    CHECK_CUDA_ERROR(cudaFree(devPtr));
    std::forward<assign_type>(output) = sum;
  }

  template <typename... Es, typename assign_type, size_t... Is>
  void result(OpAdd<Es...> const &e, assign_type &&output, len_type len,
              std::index_sequence<Is...>, bool synchronize = true) {
    auto ec = to_cuda_expr((expr::get<Is>(e) + ...));
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    T *devPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&devPtr, sizeof(T) * numBlocks));
    submit_expr(len, ec, devPtr, len);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    CHECK_CUDA_ERROR(cudaMemcpy(result, devPtr, sizeof(T) * numBlocks,
                                cudaMemcpyDeviceToHost));
    T sum = T{};
    for (iter_type i = 0; i < numBlocks; ++i) {
      sum += result[i];
    }
    delete[] result;
    CHECK_CUDA_ERROR(cudaFree(devPtr));
    std::forward<assign_type>(output) = sum;
  }

  template <typename... Es, typename assign_type, size_t... Is>
  void result(OpAdd<Es...> const &e, assign_type &&output,
              grid::region_interval<D> const &interval,
              std::index_sequence<Is...>, bool synchronize = true) {
    // inspect_result(e, interval);

    auto ec = to_cuda_expr((expr::get<Is>(e) + ...));
    iter_type len = 1;
    for (iter_type i = 0; i < D; ++i) {
      len *= interval.dims[i];
    }
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    T *devPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&devPtr, sizeof(T) * numBlocks));
    submit_expr(len, ec, devPtr, interval);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    CHECK_CUDA_ERROR(cudaMemcpy(result, devPtr, sizeof(T) * numBlocks,
                                cudaMemcpyDeviceToHost));
    T sum = T{};
    for (iter_type i = 0; i < numBlocks; ++i) {
      sum += result[i];
    }
    delete[] result;
    CHECK_CUDA_ERROR(cudaFree(devPtr));
    std::forward<assign_type>(output) = sum;
  }
};

template <typename T, size_t D>
struct expr::result_sum_trait<GridCUDA<any_vector_t<T, D>, D>> {
  template <typename E, typename assign_type, typename interval_type>
  auto submit_expr(len_type len, E const &ec, assign_type &&output,
                   interval_type &&interval, bool synchronize = true) {
    E *ecDevPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&ecDevPtr, sizeof(E)));
    CHECK_CUDA_ERROR(
        cudaMemcpy(ecDevPtr, &ec, sizeof(E), cudaMemcpyHostToDevice));

    sum_submit<any_vector_t<T, D>>{}.eval(
        len, ecDevPtr, std::forward<assign_type>(output),
        std::forward<interval_type>(interval));
    CHECK_CUDA_ERROR(cudaFree(ecDevPtr));
  }

  template <typename E, typename assign_type>
  void result(OpEvaluable<E> const &e, assign_type &&output,
              grid::region_interval<D> const &interval,
              bool synchronize = true) {
    // inspect_result(*static_cast<E const *>(&e), interval);

    auto ec = to_cuda_expr(*static_cast<E const *>(&e));
    iter_type len = 1;
    for (iter_type i = 0; i < D; ++i) {
      len *= interval.dims[i];
    }
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;

    T *devPtrArr[D];
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr[i], sizeof(T) * numBlocks));
    }
    submit_expr(len, ec, devPtrArr, interval);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    any_vector_t<T, D> sum{};
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMemcpy(result, devPtrArr[i], sizeof(T) * numBlocks,
                                  cudaMemcpyDeviceToHost));
      sum[i] = T{};
      for (iter_type n = 0; n < numBlocks; ++n) {
        sum[i] += result[n];
      }
    }
    delete[] result;
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaFree(devPtrArr[i]));
    }
    std::forward<assign_type>(output) = sum;
  }

  template <typename E, typename assign_type>
  void result(OpEvaluable<E> const &e, assign_type &&output, len_type len,
              bool synchronize = true) {
    auto ec = to_cuda_expr(*static_cast<E const *>(&e));
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    T *devPtrArr[D];
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr[i], sizeof(T) * numBlocks));
    }
    CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr, sizeof(T) * numBlocks));
    submit_expr(len, ec, devPtrArr, len);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    any_vector_t<T, D> sum{};
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMemcpy(result, devPtrArr[i], sizeof(T) * numBlocks,
                                  cudaMemcpyDeviceToHost));
      sum[i] = T{};
      for (iter_type n = 0; n < numBlocks; ++n) {
        sum[i] += result[n];
      }
    }
    delete[] result;
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaFree(devPtrArr[i]));
    }
    std::forward<assign_type>(output) = sum;
  }

  template <typename... Es, typename assign_type, size_t... Is>
  void result(OpAdd<Es...> const &e, assign_type &&output, len_type len,
              std::index_sequence<Is...>, bool synchronize = true) {
    auto ec = to_cuda_expr((expr::get<Is>(e) + ...));
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    T *devPtrArr[D];
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr[i], sizeof(T) * numBlocks));
    }
    CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr, sizeof(T) * numBlocks));
    submit_expr(len, ec, devPtrArr, len);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    any_vector_t<T, D> sum{};
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMemcpy(result, devPtrArr[i], sizeof(T) * numBlocks,
                                  cudaMemcpyDeviceToHost));
      sum[i] = T{};
      for (iter_type n = 0; n < numBlocks; ++n) {
        sum[i] += result[n];
      }
    }
    delete[] result;
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaFree(devPtrArr[i]));
    }
    std::forward<assign_type>(output) = sum;
  }

  template <typename... Es, typename assign_type, size_t... Is>
  void result(OpAdd<Es...> const &e, assign_type &&output,
              grid::region_interval<D> const &interval,
              std::index_sequence<Is...>, bool synchronize = true) {
    // inspect_result(e, interval);

    auto ec = to_cuda_expr((expr::get<Is>(e) + ...));
    iter_type len = 1;
    for (iter_type i = 0; i < D; ++i) {
      len *= interval.dims[i];
    }
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    T *devPtrArr[D];
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr[i], sizeof(T) * numBlocks));
    }
    CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr, sizeof(T) * numBlocks));
    submit_expr(len, ec, devPtrArr, interval);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    any_vector_t<T, D> sum{};
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMemcpy(result, devPtrArr[i], sizeof(T) * numBlocks,
                                  cudaMemcpyDeviceToHost));
      sum[i] = T{};
      for (iter_type n = 0; n < numBlocks; ++n) {
        sum[i] += result[n];
      }
    }
    delete[] result;
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaFree(devPtrArr[i]));
    }
    std::forward<assign_type>(output) = sum;
  }
};

template <typename T, size_t D>
struct expr::result_sum_only_trait<GridCUDA<any_vector_t<T, D>, D>> {
  template <typename E, typename assign_type, typename interval_type>
  void submit_expr(len_type len, E const &ec, assign_type &&output,
                   interval_type &&interval, bool synchronize = true) {
    E *ecDevPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&ecDevPtr, sizeof(E)));
    CHECK_CUDA_ERROR(
        cudaMemcpy(ecDevPtr, &ec, sizeof(E), cudaMemcpyHostToDevice));

    sum_submit<any_vector_t<T, D>>{}.eval(
        len, ecDevPtr, std::forward<assign_type>(output),
        std::forward<interval_type>(interval));
    CHECK_CUDA_ERROR(cudaFree(ecDevPtr));
  }
  template <typename... Es, typename assign_type, size_t... Is>
  void result(OpAdd<Es...> const &e, assign_type &&output, len_type len,
              std::index_sequence<Is...>, bool synchronize = true) {
    auto ec = to_cuda_expr((expr::get<Is>(e) + ...));
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    T *devPtrArr[D];
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr[i], sizeof(T) * numBlocks));
    }
    CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr, sizeof(T) * numBlocks));
    submit_expr(len, ec, devPtrArr, len);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    any_vector_t<T, D> sum{};
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMemcpy(result, devPtrArr[i], sizeof(T) * numBlocks,
                                  cudaMemcpyDeviceToHost));
      sum[i] = T{};
      for (iter_type n = 0; n < numBlocks; ++n) {
        sum[i] += result[n];
      }
    }
    delete[] result;
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaFree(devPtrArr[i]));
    }
    std::forward<assign_type>(output) = sum;
  }

  template <typename... Es, typename assign_type, size_t... Is>
  void result(OpAdd<Es...> const &e, assign_type &&output,
              grid::region_interval<D> const &interval,
              std::index_sequence<Is...>, bool synchronize = true) {
    // inspect_result(e, interval);

    auto ec = to_cuda_expr((expr::get<Is>(e) + ...));
    iter_type len = 1;
    for (iter_type i = 0; i < D; ++i) {
      len *= interval.dims[i];
    }
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    T *devPtrArr[D];
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr[i], sizeof(T) * numBlocks));
    }
    submit_expr(len, ec, devPtrArr, interval);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    any_vector_t<T, D> sum{};
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMemcpy(result, devPtrArr[i], sizeof(T) * numBlocks,
                                  cudaMemcpyDeviceToHost));
      sum[i] = T{};
      for (iter_type n = 0; n < numBlocks; ++n) {
        sum[i] += result[n];
      }
    }
    delete[] result;
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaFree(devPtrArr[i]));
    }
    std::forward<assign_type>(output) = sum;
  }
};

template <typename T, size_t D>
struct expr::result_sum_only_trait<GridCUDA<T, D>> {
  template <typename E, typename assign_type, typename interval_type>
  void submit_expr(len_type len, E const &ec, assign_type &&output,
                   interval_type &&interval, bool synchronize = true) {
    E *ecDevPtr;
    CHECK_CUDA_ERROR(cudaMalloc(&ecDevPtr, sizeof(E)));
    CHECK_CUDA_ERROR(
        cudaMemcpy(ecDevPtr, &ec, sizeof(E), cudaMemcpyHostToDevice));

    sum_submit<T>{}.eval(len, ecDevPtr, std::forward<assign_type>(output),
                         std::forward<interval_type>(interval));
    CHECK_CUDA_ERROR(cudaFree(ecDevPtr));
  }
  template <typename... Es, typename assign_type, size_t... Is>
  void result(OpAdd<Es...> const &e, assign_type &&output, len_type len,
              std::index_sequence<Is...>, bool synchronize = true) {
    auto ec = to_cuda_expr((expr::get<Is>(e) + ...));
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    T *devPtrArr[D];
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr[i], sizeof(T) * numBlocks));
    }
    submit_expr(len, ec, devPtrArr, len);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    any_vector_t<T, D> sum{};
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMemcpy(result, devPtrArr[i], sizeof(T) * numBlocks,
                                  cudaMemcpyDeviceToHost));
      sum[i] = T{};
      for (iter_type n = 0; n < numBlocks; ++n) {
        sum[i] += result[n];
      }
    }
    delete[] result;
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaFree(devPtrArr[i]));
    }
    std::forward<assign_type>(output) = sum;
  }

  template <typename... Es, typename assign_type, size_t... Is>
  void result(OpAdd<Es...> const &e, assign_type &&output,
              grid::region_interval<D> const &interval,
              std::index_sequence<Is...>, bool synchronize = true) {
    // inspect_result((expr::get<Is>(e) + ...), interval);

    auto ec = to_cuda_expr((expr::get<Is>(e) + ...));
    iter_type len = 1;
    for (iter_type i = 0; i < D; ++i) {
      len *= interval.dims[i];
    }
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    T *devPtrArr[D];
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr[i], sizeof(T) * numBlocks));
    }
    CHECK_CUDA_ERROR(cudaMalloc(&devPtrArr, sizeof(T) * numBlocks));
    submit_expr(len, ec, devPtrArr, interval);
    if (synchronize) CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    T *result = new T[numBlocks];
    any_vector_t<T, D> sum{};
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMemcpy(result, devPtrArr[i], sizeof(T) * numBlocks,
                                  cudaMemcpyDeviceToHost));
      sum[i] = T{};
      for (iter_type n = 0; n < numBlocks; ++n) {
        sum[i] += result[n];
      }
    }
    delete[] result;
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaFree(devPtrArr[i]));
    }
    std::forward<assign_type>(output) = sum;
  }
};

template <typename T, size_t D>
struct expr::result_sum_only_trait<BoundaryGridCUDA<T, D>>
    : expr::result_sum_only_trait<GridCUDA<T, D>> {};
template <typename T, size_t D>
struct expr::result_sum_only_trait<RegionalGridCUDA<T, D>>
    : expr::result_sum_only_trait<GridCUDA<T, D>> {};
template <typename T, size_t D>
struct expr::result_sum_trait<BoundaryGridCUDA<T, D>>
    : expr::result_sum_trait<GridCUDA<T, D>> {};
template <typename T, size_t D>
struct expr::result_sum_trait<RegionalGridCUDA<T, D>>
    : expr::result_sum_trait<GridCUDA<T, D>> {};

template <typename T, size_t D>
struct expr::evaluate_expression_trait<BoundaryGridCUDA<T, D>>
    : expr::evaluate_expression_trait<GridCUDA<T, D>> {
  using evaluate_expression_trait<GridCUDA<T, D>>::evaluate_expression_trait;
};
template <typename T, size_t D>
struct expr::evaluate_expression_trait<RegionalGridCUDA<T, D>>
    : expr::evaluate_expression_trait<GridCUDA<T, D>> {
  using evaluate_expression_trait<GridCUDA<T, D>>::evaluate_expression_trait;
};

template <typename T, size_t D>
struct expr::accumulate_expression_trait<BoundaryGridCUDA<T, D>>
    : expr::accumulate_expression_trait<GridCUDA<T, D>> {
  using accumulate_expression_trait<
      GridCUDA<T, D>>::accumulate_expression_trait;
};
template <typename T, size_t D>
struct expr::accumulate_expression_trait<RegionalGridCUDA<T, D>>
    : expr::accumulate_expression_trait<GridCUDA<T, D>> {
  using accumulate_expression_trait<
      GridCUDA<T, D>>::accumulate_expression_trait;
};

template <typename T, size_t D>
struct expr::result_only_trait<GridCUDA<T, D>> {
  template <typename... Es, typename assign_type, size_t... Is>
  result_only_trait(OpAdd<Es...> const &e, assign_type &&data,
                    grid::region_interval<D> const &interval,
                    std::index_sequence<Is...>) {
    evaluate_expression_trait<GridCUDA<T, D>>(
        (expr::get<Is>(e) + ...), std::forward<assign_type>(data), interval);
  }

  template <typename... Es, typename assign_type, size_t... Is>
  result_only_trait(OpAdd<Es...> const &e, assign_type &&data, len_type len,
                    std::index_sequence<Is...>) {
    evaluate_expression_trait<GridCUDA<T, D>>(
        (expr::get<Is>(e) + ...), std::forward<assign_type>(data), len);
  }
};

#endif
