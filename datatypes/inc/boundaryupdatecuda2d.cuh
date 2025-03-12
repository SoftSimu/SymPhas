
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

/* 2 dimensional boundaries
 *
 */

template <Side side0, Side side1>
struct update_periodic_cuda {
  template <typename T>
  __device__ static void update(T* grid, int N, int M);
  template <typename T>
  __device__ static void update(T* grid0, T* grid1, int N, int M);
};

// Kernel to update left ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::LEFT, Side::LEFT>::update(T* grid,
                                                                     int N,
                                                                     int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= BOUNDARY_DEPTH && j < M - BOUNDARY_DEPTH) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid[j * N + i] = grid[j * N + (N - 2 * BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update right ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::RIGHT, Side::RIGHT>::update(T* grid,
                                                                       int N,
                                                                       int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= BOUNDARY_DEPTH && j < M - BOUNDARY_DEPTH) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid[j * N + (N - BOUNDARY_DEPTH + i)] =
          grid[j * N + (BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update top ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::TOP, Side::TOP>::update(T* grid,
                                                                   int N,
                                                                   int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= BOUNDARY_DEPTH && i < N - BOUNDARY_DEPTH) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      grid[j * N + i] = grid[(M - 2 * BOUNDARY_DEPTH + j) * N + i];
    }
  }
}

// Kernel to update bottom ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::BOTTOM, Side::BOTTOM>::update(
    T* grid, int N, int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= BOUNDARY_DEPTH && i < N - BOUNDARY_DEPTH) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      grid[(M - BOUNDARY_DEPTH + j) * N + i] =
          grid[(BOUNDARY_DEPTH + j) * N + i];
    }
  }
}

// Kernel to update left ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::TOP, Side::LEFT>::update(T* grid,
                                                                    int N,
                                                                    int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= M - BOUNDARY_DEPTH && j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid[j * N + i] = grid[j * N + (N - 2 * BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update left ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::BOTTOM, Side::LEFT>::update(T* grid,
                                                                       int N,
                                                                       int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= 0 && j < BOUNDARY_DEPTH) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid[j * N + i] = grid[j * N + (N - 2 * BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update right ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::TOP, Side::RIGHT>::update(T* grid,
                                                                     int N,
                                                                     int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= M - BOUNDARY_DEPTH && j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid[j * N + (N - BOUNDARY_DEPTH + i)] =
          grid[j * N + (BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update right ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::BOTTOM, Side::RIGHT>::update(T* grid,
                                                                        int N,
                                                                        int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= 0 && j < BOUNDARY_DEPTH) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid[j * N + (N - BOUNDARY_DEPTH + i)] =
          grid[j * N + (BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update left ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::LEFT, Side::LEFT>::update(T* grid0,
                                                                     T* grid1,
                                                                     int N,
                                                                     int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= BOUNDARY_DEPTH && j < M - BOUNDARY_DEPTH) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid0[j * N + i] = grid0[j * N + (N - 2 * BOUNDARY_DEPTH + i)];
      grid1[j * N + i] = grid1[j * N + (N - 2 * BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update right ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::RIGHT, Side::RIGHT>::update(T* grid0,
                                                                       T* grid1,
                                                                       int N,
                                                                       int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= BOUNDARY_DEPTH && j < M - BOUNDARY_DEPTH) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid0[j * N + (N - BOUNDARY_DEPTH + i)] =
          grid0[j * N + (BOUNDARY_DEPTH + i)];
      grid1[j * N + (N - BOUNDARY_DEPTH + i)] =
          grid1[j * N + (BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update top ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::TOP, Side::TOP>::update(T* grid0,
                                                                   T* grid1,
                                                                   int N,
                                                                   int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= BOUNDARY_DEPTH && i < N - BOUNDARY_DEPTH) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      grid0[j * N + i] = grid0[(M - 2 * BOUNDARY_DEPTH + j) * N + i];
      grid1[j * N + i] = grid1[(M - 2 * BOUNDARY_DEPTH + j) * N + i];
    }
  }
}

// Kernel to update bottom ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::BOTTOM, Side::BOTTOM>::update(
    T* grid0, T* grid1, int N, int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= BOUNDARY_DEPTH && i < N - BOUNDARY_DEPTH) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      grid0[(M - BOUNDARY_DEPTH + j) * N + i] =
          grid0[(BOUNDARY_DEPTH + j) * N + i];
      grid1[(M - BOUNDARY_DEPTH + j) * N + i] =
          grid1[(BOUNDARY_DEPTH + j) * N + i];
    }
  }
}

// Kernel to update left ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::TOP, Side::LEFT>::update(T* grid0,
                                                                    T* grid1,
                                                                    int N,
                                                                    int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= M - BOUNDARY_DEPTH && j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid0[j * N + i] = grid0[j * N + (N - 2 * BOUNDARY_DEPTH + i)];
      grid1[j * N + i] = grid1[j * N + (N - 2 * BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update left ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::BOTTOM, Side::LEFT>::update(T* grid0,
                                                                       T* grid1,
                                                                       int N,
                                                                       int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= 0 && j < BOUNDARY_DEPTH) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid0[j * N + i] = grid0[j * N + (N - 2 * BOUNDARY_DEPTH + i)];
      grid1[j * N + i] = grid1[j * N + (N - 2 * BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update right ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::TOP, Side::RIGHT>::update(T* grid0,
                                                                     T* grid1,
                                                                     int N,
                                                                     int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= M - BOUNDARY_DEPTH && j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid0[j * N + (N - BOUNDARY_DEPTH + i)] =
          grid0[j * N + (BOUNDARY_DEPTH + i)];
      grid1[j * N + (N - BOUNDARY_DEPTH + i)] =
          grid1[j * N + (BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update right ghost cells
template <>
template <typename T>
__device__ void update_periodic_cuda<Side::BOTTOM, Side::RIGHT>::update(
    T* grid0, T* grid1, int N, int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= 0 && j < BOUNDARY_DEPTH) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid0[j * N + (N - BOUNDARY_DEPTH + i)] =
          grid0[j * N + (BOUNDARY_DEPTH + i)];
      grid1[j * N + (N - BOUNDARY_DEPTH + i)] =
          grid1[j * N + (BOUNDARY_DEPTH + i)];
    }
  }
}

// Kernel to update left ghost cells
template <Side side0, Side side1, typename T>
__global__ void call_update_periodic_cuda(T* grid, int N, int M) {
  update_periodic_cuda<side0, side1>::update(grid, N, M);
}
// Kernel to update left ghost cells
template <Side side0, Side side1, typename T>
__global__ void call_update_periodic_cuda_vec(T* grid0, T* grid1, int N,
                                              int M) {
  update_periodic_cuda<side0, side1>::update(grid0, grid1, N, M);
}

namespace grid {

// Kernel to update left ghost cells
template <Side side0, Side side1, typename T>
struct update_boundary_2d {
  static void update(T* grid, int N, int M) {
    int nn = std::max(M, N);
    int numBlocks = (nn + BLOCK_SIZE - 1) / BLOCK_SIZE;
    call_update_periodic_cuda<side0, side1> CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        grid, N, M);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }
};

// Kernel to update left ghost cells
template <Side side0, Side side1, typename T>
struct update_boundary_2d<side0, side1, any_vector_t<T, 2>> {
  static void update(T* (&grid)[2], int N, int M) {
    int nn = std::max(M, N);
    int numBlocks = (nn + BLOCK_SIZE - 1) / BLOCK_SIZE;
    call_update_periodic_cuda_vec<side0, side1> CUDA_KERNEL(
        numBlocks, BLOCK_SIZE)(grid[0], grid[1], N, M);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }
};

}  // namespace grid

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 1>::
operator()(const grid::Boundary<T, 1>*, GridCUDA<T, 2>& grid) {
  grid::update_boundary_2d<Side::LEFT, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::TOP, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::BOTTOM, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::RIGHT,
                                        1>::operator()(const grid::Boundary<T,
                                                                            1>*,
                                                       GridCUDA<T, 2>& grid) {
  grid::update_boundary_2d<Side::RIGHT, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::TOP, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::BOTTOM, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::TOP, 1>::
operator()(const grid::Boundary<T, 1>*, GridCUDA<T, 2>& grid) {
  grid::update_boundary_2d<Side::TOP, Side::TOP, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BOTTOM,
                                        1>::operator()(const grid::Boundary<T,
                                                                            1>*,
                                                       GridCUDA<T, 2>& grid) {
  grid::update_boundary_2d<Side::BOTTOM, Side::BOTTOM, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::LEFT,
                                        1>::operator()(const grid::Boundary<T,
                                                                            1>*,
                                                       GridCUDA<T, 2>& grid) {
  grid::update_boundary_2d<Side::LEFT, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::TOP, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::BOTTOM, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::RIGHT,
                                        1>::operator()(const grid::Boundary<T,
                                                                            1>*,
                                                       GridCUDA<T, 2>& grid) {
  grid::update_boundary_2d<Side::RIGHT, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::TOP, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::BOTTOM, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::TOP, 1>::
operator()(const grid::Boundary<T, 1>*, GridCUDA<T, 2>& grid) {
  grid::update_boundary_2d<Side::TOP, Side::TOP, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::TOP, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::TOP, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::BOTTOM,
                                        1>::operator()(const grid::Boundary<T,
                                                                            1>*,
                                                       GridCUDA<T, 2>& grid) {
  grid::update_boundary_2d<Side::BOTTOM, Side::BOTTOM, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::BOTTOM, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
  grid::update_boundary_2d<Side::BOTTOM, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1]);
}

// Regional Boundaries

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 1>::
operator()(const grid::Boundary<T, 1>*, RegionalGridCUDA<T, 2>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT>{},
                           grid);
  regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::TOP>{},
                           grid);
  regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::BOTTOM>{},
                           grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC, Side::RIGHT,
    1>::operator()(const grid::Boundary<T, 1>*, RegionalGridCUDA<T, 2>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT>{},
                           grid);
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::TOP>{},
                           grid);
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::BOTTOM>{},
                           grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::TOP, 1>::
operator()(const grid::Boundary<T, 1>*, RegionalGridCUDA<T, 2>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::TOP>{},
                           grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC, Side::BOTTOM,
    1>::operator()(const grid::Boundary<T, 1>*, RegionalGridCUDA<T, 2>& grid) {
  regional_update_boundary(
      symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC0, Side::LEFT,
    1>::operator()(const grid::Boundary<T, 1>*, RegionalGridCUDA<T, 2>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC0, Side::RIGHT,
    1>::operator()(const grid::Boundary<T, 1>*, RegionalGridCUDA<T, 2>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::TOP, 1>::
operator()(const grid::Boundary<T, 1>*, RegionalGridCUDA<T, 2>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::BOTTOM>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC0, Side::BOTTOM,
    1>::operator()(const grid::Boundary<T, 1>*, RegionalGridCUDA<T, 2>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::BOTTOM>{}, grid);
}

// *********************************************************************
/* DEFAULT BOUNDARY ALGORITHMS
 */

template <typename T>
__global__ void copyLeftBoundaryFromGrid(const T* grid, T* boundaryDevice,
                                         int N, int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      boundaryDevice[j * BOUNDARY_DEPTH + i] = grid[j * N + i];
    }
  }
}

template <typename T>
__global__ void copyRightBoundaryFromGrid(const T* grid, T* boundaryDevice,
                                          int N, int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      boundaryDevice[j * BOUNDARY_DEPTH + i] =
          grid[j * N + (N - BOUNDARY_DEPTH + i)];
    }
  }
}

template <typename T>
__global__ void copyTopBoundaryFromGrid(const T* grid, T* boundaryDevice, int N,
                                        int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      boundaryDevice[i * BOUNDARY_DEPTH + j] = grid[j * N + i];
    }
  }
}

template <typename T>
__global__ void copyBottomBoundaryFromGrid(const T* grid, T* boundaryDevice,
                                           int N, int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      boundaryDevice[i * BOUNDARY_DEPTH + j] =
          grid[(M - BOUNDARY_DEPTH + j) * N + i];
    }
  }
}

template <typename T>
__global__ void copyLeftBoundaryToGrid(const T* boundaryDevice, T* grid, int N,
                                       int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid[j * N + i] = boundaryDevice[j * BOUNDARY_DEPTH + i];
    }
  }
}

template <typename T>
__global__ void copyRightBoundaryToGrid(const T* boundaryDevice, T* grid, int N,
                                        int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid[j * N + (N - BOUNDARY_DEPTH + i)] =
          boundaryDevice[j * BOUNDARY_DEPTH + i];
    }
  }
}

template <typename T>
__global__ void copyTopBoundaryToGrid(const T* boundaryDevice, T* grid, int N,
                                      int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      grid[j * N + i] = boundaryDevice[i * BOUNDARY_DEPTH + j];
    }
  }
}

template <typename T>
__global__ void copyBottomBoundaryToGrid(const T* boundaryDevice, T* grid,
                                         int N, int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      grid[(M - BOUNDARY_DEPTH + j) * N + i] =
          boundaryDevice[i * BOUNDARY_DEPTH + j];
    }
  }
}

template <typename T>
__global__ void copyLeftBoundaryFromGrid(T* grid0, T* grid1,
                                         any_vector_t<T, 2>* boundaryDevice,
                                         int N, int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      boundaryDevice[j * BOUNDARY_DEPTH + i][0] = grid0[j * N + i];
      boundaryDevice[j * BOUNDARY_DEPTH + i][1] = grid1[j * N + i];
    }
  }
}

template <typename T>
__global__ void copyRightBoundaryFromGrid(T* grid0, T* grid1,
                                          any_vector_t<T, 2>* boundaryDevice,
                                          int N, int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      boundaryDevice[j * BOUNDARY_DEPTH + i][0] =
          grid0[j * N + (N - BOUNDARY_DEPTH + i)];
      boundaryDevice[j * BOUNDARY_DEPTH + i][1] =
          grid1[j * N + (N - BOUNDARY_DEPTH + i)];
    }
  }
}

template <typename T>
__global__ void copyTopBoundaryFromGrid(T* grid0, T* grid1,
                                        any_vector_t<T, 2>* boundaryDevice,
                                        int N, int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      boundaryDevice[i * BOUNDARY_DEPTH + j][0] = grid0[j * N + i];
      boundaryDevice[i * BOUNDARY_DEPTH + j][1] = grid1[j * N + i];
    }
  }
}

template <typename T>
__global__ void copyBottomBoundaryFromGrid(T* grid0, T* grid1,
                                           any_vector_t<T, 2>* boundaryDevice,
                                           int N, int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      boundaryDevice[i * BOUNDARY_DEPTH + j][0] =
          grid0[(M - BOUNDARY_DEPTH + j) * N + i];
      boundaryDevice[i * BOUNDARY_DEPTH + j][1] =
          grid1[(M - BOUNDARY_DEPTH + j) * N + i];
    }
  }
}

template <typename T>
__global__ void copyLeftBoundaryToGrid(const any_vector_t<T, 2>* boundaryDevice,
                                       T* grid0, T* grid1, int N, int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid0[j * N + i] = boundaryDevice[j * BOUNDARY_DEPTH + i][0];
      grid1[j * N + i] = boundaryDevice[j * BOUNDARY_DEPTH + i][1];
    }
  }
}

template <typename T>
__global__ void copyRightBoundaryToGrid(
    const any_vector_t<T, 2>* boundaryDevice, T* grid0, T* grid1, int N,
    int M) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j < M) {
    for (int i = 0; i < BOUNDARY_DEPTH; ++i) {
      grid0[j * N + (N - BOUNDARY_DEPTH + i)] =
          boundaryDevice[j * BOUNDARY_DEPTH + i][0];
      grid1[j * N + (N - BOUNDARY_DEPTH + i)] =
          boundaryDevice[j * BOUNDARY_DEPTH + i][1];
    }
  }
}

template <typename T>
__global__ void copyTopBoundaryToGrid(const any_vector_t<T, 2>* boundaryDevice,
                                      T* grid0, T* grid1, int N, int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      grid0[j * N + i] = boundaryDevice[i * BOUNDARY_DEPTH + j][0];
      grid1[j * N + i] = boundaryDevice[i * BOUNDARY_DEPTH + j][1];
    }
  }
}

template <typename T>
__global__ void copyBottomBoundaryToGrid(
    const any_vector_t<T, 2>* boundaryDevice, T* grid0, T* grid1, int N,
    int M) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    for (int j = 0; j < BOUNDARY_DEPTH; ++j) {
      grid0[(M - BOUNDARY_DEPTH + j) * N + i] =
          boundaryDevice[i * BOUNDARY_DEPTH + j][0];
      grid1[(M - BOUNDARY_DEPTH + j) * N + i] =
          boundaryDevice[i * BOUNDARY_DEPTH + j][1];
    }
  }
}

template <Side side>
struct do_boundary_copy_to_grid {
  template <typename T>
  void operator()(const T* boundaryDevice, T* grid, int N, int M);
  template <typename T>
  void operator()(const any_vector_t<T, 2>* boundaryDevice, T* (&grid)[2],
                  int N, int M);
};

template <>
template <typename T>
void do_boundary_copy_to_grid<Side::LEFT>::operator()(const T* boundaryDevice,
                                                      T* grid, int N, int M) {
  int numBlocks = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyLeftBoundaryToGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(boundaryDevice,
                                                            grid, N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_to_grid<Side::LEFT>::operator()(
    const any_vector_t<T, 2>* boundaryDevice, T* (&grid)[2], int N, int M) {
  int numBlocks = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyLeftBoundaryToGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      boundaryDevice, grid[0], grid[1], N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_to_grid<Side::RIGHT>::operator()(const T* boundaryDevice,
                                                       T* grid, int N, int M) {
  int numBlocks = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyRightBoundaryToGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(boundaryDevice,
                                                             grid, N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_to_grid<Side::RIGHT>::operator()(
    const any_vector_t<T, 2>* boundaryDevice, T* (&grid)[2], int N, int M) {
  int numBlocks = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyRightBoundaryToGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      boundaryDevice, grid[0], grid[1], N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_to_grid<Side::TOP>::operator()(const T* boundaryDevice,
                                                     T* grid, int N, int M) {
  int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyTopBoundaryToGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(boundaryDevice, grid,
                                                           N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_to_grid<Side::TOP>::operator()(
    const any_vector_t<T, 2>* boundaryDevice, T* (&grid)[2], int N, int M) {
  int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyTopBoundaryToGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      boundaryDevice, grid[0], grid[1], N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_to_grid<Side::BOTTOM>::operator()(const T* boundaryDevice,
                                                        T* grid, int N, int M) {
  int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyBottomBoundaryToGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(boundaryDevice,
                                                              grid, N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_to_grid<Side::BOTTOM>::operator()(
    const any_vector_t<T, 2>* boundaryDevice, T* (&grid)[2], int N, int M) {
  int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyBottomBoundaryToGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      boundaryDevice, grid[0], grid[1], N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <Side side>
struct do_boundary_copy_from_grid {
  template <typename T>
  void operator()(const T* grid, T* boundaryDevice, int N, int M);
  template <typename T>
  void operator()(T* const (&grid)[2], any_vector_t<T, 2>* boundaryDevice,
                  int N, int M);
};

template <>
template <typename T>
void do_boundary_copy_from_grid<Side::LEFT>::operator()(const T* grid,
                                                        T* boundaryDevice,
                                                        int N, int M) {
  int numBlocks = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyLeftBoundaryFromGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      grid, boundaryDevice, N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_from_grid<Side::LEFT>::operator()(
    T* const (&grid)[2], any_vector_t<T, 2>* boundaryDevice, int N, int M) {
  int numBlocks = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyLeftBoundaryFromGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      grid[0], grid[1], boundaryDevice, N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_from_grid<Side::RIGHT>::operator()(const T* grid,
                                                         T* boundaryDevice,
                                                         int N, int M) {
  int numBlocks = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyRightBoundaryFromGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      grid, boundaryDevice, N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_from_grid<Side::RIGHT>::operator()(
    T* const (&grid)[2], any_vector_t<T, 2>* boundaryDevice, int N, int M) {
  int numBlocks = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyRightBoundaryFromGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      grid[0], grid[1], boundaryDevice, N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_from_grid<Side::TOP>::operator()(const T* grid,
                                                       T* boundaryDevice, int N,
                                                       int M) {
  int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyTopBoundaryFromGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      grid, boundaryDevice, N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_from_grid<Side::TOP>::operator()(
    T* const (&grid)[2], any_vector_t<T, 2>* boundaryDevice, int N, int M) {
  int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyTopBoundaryFromGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      grid[0], grid[1], boundaryDevice, N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_from_grid<Side::BOTTOM>::operator()(const T* grid,
                                                          T* boundaryDevice,
                                                          int N, int M) {
  int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyBottomBoundaryFromGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      grid, boundaryDevice, N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

template <>
template <typename T>
void do_boundary_copy_from_grid<Side::BOTTOM>::operator()(
    T* const (&grid)[2], any_vector_t<T, 2>* boundaryDevice, int N, int M) {
  int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
  copyBottomBoundaryFromGrid CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
      grid[0], grid[1], boundaryDevice, N, M);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

// 2 dimension

template <Side side, typename T>
void update_default_boundary(const grid::Boundary<T, 1>* b,
                             GridCUDA<T, 2>& grid, double time) {
  auto* bd =
      static_cast<grid::BoundaryApplied<T, 1, BoundaryType::DEFAULT> const*>(b);

  const double v0 = bd->v[0];
  const double v1 = bd->v[1];
  double h[2];

  // go backwards or forwards in iteration depending on the interval
  int fx = (v0 < v1) ? 1 : -1;
  double delta_h = bd->h * fx;

  const double x0 = (side == Side::LEFT || side == Side::RIGHT) ? 0 : v0;
  const double x1 = (side == Side::LEFT || side == Side::RIGHT) ? 0 : v1;
  const double y0 = (side == Side::TOP || side == Side::BOTTOM) ? 0 : v0;
  const double y1 = (side == Side::TOP || side == Side::BOTTOM) ? 0 : v1;

  h[0] = (side == Side::LEFT || side == Side::RIGHT) ? 0 : delta_h;
  h[1] = (side == Side::TOP || side == Side::BOTTOM) ? 0 : delta_h;

  iter_type L =
      (side == Side::LEFT || side == Side::RIGHT) ? grid.dims[0] : grid.dims[1];
  iter_type M =
      (side == Side::TOP || side == Side::BOTTOM) ? grid.dims[0] : grid.dims[1];

  T* boundaryHost = new T[BOUNDARY_DEPTH * M]{};

  T* boundaryDevice;
  CHECK_CUDA_ERROR(cudaMalloc(&boundaryDevice, BOUNDARY_DEPTH * M * sizeof(T)));

  do_boundary_copy_from_grid<side>{}(grid.values, boundaryDevice, grid.dims[0],
                                     grid.dims[1]);

  CHECK_CUDA_ERROR(cudaMemcpy(boundaryHost, boundaryDevice,
                              BOUNDARY_DEPTH * M * sizeof(T),
                              cudaMemcpyDeviceToHost));

  for (iter_type i = 0; i < BOUNDARY_DEPTH; ++i) {
    for (iter_type j = 0; j < BOUNDARY_DEPTH; ++j) {
      bd->update(boundaryHost[j * BOUNDARY_DEPTH + i], x0, y0, time);
    }
    for (iter_type j = BOUNDARY_DEPTH; j < M - BOUNDARY_DEPTH; ++j) {
      bd->update(boundaryHost[j * BOUNDARY_DEPTH + i],
                 x0 + (j - BOUNDARY_DEPTH) * h[0],
                 y0 + (j - BOUNDARY_DEPTH) * h[1], time);
    }
    for (iter_type j = M - BOUNDARY_DEPTH; j < M; ++j) {
      bd->update(boundaryHost[j * BOUNDARY_DEPTH + i], x1, y1, time);
    }
  }

  CHECK_CUDA_ERROR(cudaMemcpy(boundaryDevice, boundaryHost,
                              BOUNDARY_DEPTH * M * sizeof(T),
                              cudaMemcpyHostToDevice));

  do_boundary_copy_to_grid<side>{}(boundaryDevice, grid.values, grid.dims[0],
                                   grid.dims[1]);

  CHECK_CUDA_ERROR(cudaFree(boundaryDevice));
  delete[] boundaryHost;
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 1>::
operator()(const grid::Boundary<T, 1>* b, GridCUDA<T, 2>& grid, double time) {
  update_default_boundary<Side::LEFT>(b, grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 1>::
operator()(const grid::Boundary<T, 1>* b, GridCUDA<T, 2>& grid, double time) {
  update_default_boundary<Side::RIGHT>(b, grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::TOP, 1>::
operator()(const grid::Boundary<T, 1>* b, GridCUDA<T, 2>& grid, double time) {
  update_default_boundary<Side::TOP>(b, grid, time);
}

template <>
template <typename T>
void symphas::internal::
    update_boundary<BoundaryType::DEFAULT, Side::BOTTOM, 1>::operator()(
        const grid::Boundary<T, 1>* b, GridCUDA<T, 2>& grid, double time) {
  update_default_boundary<Side::BOTTOM>(b, grid, time);
}

// 2 dimension

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::LEFT, 1>::
operator()(const grid::Boundary<T, 1>* b, RegionalGridCUDA<T, 2>& grid,
           double time) {
  regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT>{}, b,
                           grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::RIGHT, 1>::
operator()(const grid::Boundary<T, 1>* b, RegionalGridCUDA<T, 2>& grid,
           double time) {
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT>{},
                           b, grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::DEFAULT, Side::TOP, 1>::
operator()(const grid::Boundary<T, 1>* b, RegionalGridCUDA<T, 2>& grid,
           double time) {
  regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::TOP>{}, b,
                           grid, time);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::DEFAULT, Side::BOTTOM,
    1>::operator()(const grid::Boundary<T, 1>* b, RegionalGridCUDA<T, 2>& grid,
                   double time) {
  regional_update_boundary(
      symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM>{}, b, grid, time);
}
