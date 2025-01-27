
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
 * MODULE:  lib
 * PURPOSE: Adds additional necessary features for applying Fourier transforms.
 *
 * ***************************************************************************
 */

#pragma once

#include "spslib.h"

#ifdef USING_CUDA
#include <cuda_runtime.h>
#endif

namespace grid {

#ifdef USING_CUDA

template <typename T>
__global__ void scale_arr_on_device(T* device_data_y, len_type len, double r) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < len) {
    device_data_y[idx] /= r;
  }
}

template <typename T>
void scale_cuda(T* device_data_y, len_type len) {
  // Calculate the number of blocks and threads per block
  int threadsPerBlock = 256;
  int blocksPerGrid = (len + threadsPerBlock - 1) / threadsPerBlock;

  scale_arr_on_device<<<blocksPerGrid, threadsPerBlock>>>(device_data_y, len,
                                                          1.0 / len);

  // Wait for GPU to finish before accessing on host
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

#endif

}  // namespace grid
