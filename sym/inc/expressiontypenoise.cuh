
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
 * PURPOSE: Defines random noise, which represents, for instance, thermal
 * fluctuations in a field.
 *
 * ***************************************************************************
 */

#pragma once

#include <curand_kernel.h>

#include <ctime>

#include "expressionproperties.cuh"
#include "expressiontypenoise.h"

__global__ void generate_kernel(curandState *state, int len, int seed) {
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  if (id < len) {
    curand_init(seed, id, 0, &state[id]);
  }
}

namespace expr {

template <typename T, size_t D>
struct random_state<GridCUDA<T, D>> {
  len_type len;
  curandState *states;

  random_state() : len{0}, states{nullptr} {}

  random_state(iter_type len) : len{len}, states{nullptr} {
    cudaMalloc(&states, len * sizeof(curandState));
    auto numBlocks = (len + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned long long seed = static_cast<unsigned long long>(time(NULL));

    generate_kernel CUDA_KERNEL(numBlocks, BLOCK_SIZE)(states, len, seed);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  }

  random_state(random_state &&other) : states() { swap(*this, other); }
  random_state(random_state const &other) : random_state(other.len) {
    cudaMemcpy(states, other.states, len * sizeof(curandState),
               cudaMemcpyDeviceToDevice);
  }

  friend void swap(random_state &a, random_state &b) {
    using std::swap;
    swap(a.len, b.len);
    swap(a.states, b.states);
  }

  ~random_state() { CHECK_CUDA_ERROR(cudaFree(states)); }
};

}  // namespace expr
