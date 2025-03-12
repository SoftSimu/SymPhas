
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
/* an implementation of all BOUNDARY ITERATION ALGORITHMS
 */

template <Side side0, Side side1, Side side2>
struct update_periodic_cuda_3d {
  template <typename T>
  __device__ static void update(T* grid, int N, int M, int P) {}
  template <typename T>
  __device__ static void update(T* grid0, T* grid1, T* grid2, int N, int M,
                                int P) {}
};

////////////////////////////////////////////////////////////////////////////////////////
// Important disclaimer
//
// the 0,0,0 point starts at LEFT, BOTTOM, FRONT
// x index increases from left to right
// y index increases bottom to top
// z index increases front to back
//
// This is meant to be just like a cartesian grid
//
////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
__device__ void general_periodic_iteration(T* grid, len_type (&dims)[2],
                                           iter_type (&stride)[3],
                                           iter_type (&current_offset)[3],
                                           iter_type (&periodic_offset)[3]) {
  int n = blockIdx.x * blockDim.x + threadIdx.x;
  int i = n % dims[0];
  int j = n / dims[0];
  if (i < dims[0] && j < dims[1]) {
    for (int k = 0; k < BOUNDARY_DEPTH; ++k) {
      iter_type current_index = (current_offset[2] + k) * stride[2] +
                                (current_offset[1] + j) * stride[1] +
                                (current_offset[0] + i) * stride[0];
      iter_type offset =
          (current_offset[2] + periodic_offset[2] + k) * stride[2] +
          (current_offset[1] + periodic_offset[1] + j) * stride[1] +
          (current_offset[0] + periodic_offset[0] + i) * stride[0];
      grid[current_index] = grid[offset];
    }
  }
}

template <>
template <typename T>
__device__ void
update_periodic_cuda_3d<Side::TOP, Side::TOP, Side::TOP>::update(T* grid, int N,
                                                                 int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {N - BOUNDARY_DEPTH * 2, P - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH,
                                 M - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {0, 0, -(M - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void
update_periodic_cuda_3d<Side::TOP, Side::TOP, Side::FRONT>::update(T* grid,
                                                                   int N, int M,
                                                                   int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {N - BOUNDARY_DEPTH * 2, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {BOUNDARY_DEPTH, 0, M - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {0, P - 2 * BOUNDARY_DEPTH,
                                  -(M - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::TOP, Side::LEFT,
                                        Side::FRONT>::update(T* grid, int N,
                                                             int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {0, 0, M - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {N - 2 * BOUNDARY_DEPTH,
                                  P - 2 * BOUNDARY_DEPTH,
                                  -(M - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::TOP, Side::RIGHT,
                                        Side::FRONT>::update(T* grid, int N,
                                                             int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {N - BOUNDARY_DEPTH, 0, M - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {-(N - 2 * BOUNDARY_DEPTH),
                                  P - 2 * BOUNDARY_DEPTH,
                                  -(M - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void
update_periodic_cuda_3d<Side::TOP, Side::LEFT, Side::BACK>::update(T* grid,
                                                                   int N, int M,
                                                                   int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {0, P - BOUNDARY_DEPTH, M - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {N - 2 * BOUNDARY_DEPTH,
                                  -(P - 2 * BOUNDARY_DEPTH),
                                  -(M - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::TOP, Side::RIGHT,
                                        Side::BACK>::update(T* grid, int N,
                                                            int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {N - BOUNDARY_DEPTH, P - BOUNDARY_DEPTH,
                                 M - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {-(N - 2 * BOUNDARY_DEPTH),
                                  -(P - 2 * BOUNDARY_DEPTH),
                                  -(M - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void
update_periodic_cuda_3d<Side::TOP, Side::TOP, Side::BACK>::update(T* grid,
                                                                  int N, int M,
                                                                  int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {N - BOUNDARY_DEPTH * 2, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {BOUNDARY_DEPTH, P - BOUNDARY_DEPTH,
                                 M - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {0, -(P - 2 * BOUNDARY_DEPTH),
                                  -(M - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void
update_periodic_cuda_3d<Side::TOP, Side::TOP, Side::LEFT>::update(T* grid,
                                                                  int N, int M,
                                                                  int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, P - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {0, BOUNDARY_DEPTH, M - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {N - 2 * BOUNDARY_DEPTH, 0,
                                  -(M - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void
update_periodic_cuda_3d<Side::TOP, Side::TOP, Side::RIGHT>::update(T* grid,
                                                                   int N, int M,
                                                                   int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, P - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {N - BOUNDARY_DEPTH, BOUNDARY_DEPTH,
                                 M - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {-(N - 2 * BOUNDARY_DEPTH), 0,
                                  -(M - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BOTTOM, Side::BOTTOM,
                                        Side::BOTTOM>::update(T* grid, int N,
                                                              int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {N - BOUNDARY_DEPTH * 2, P - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH, 0};
  iter_type periodic_offset[3] = {0, 0, M - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BOTTOM, Side::BOTTOM,
                                        Side::FRONT>::update(T* grid, int N,
                                                             int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {N - BOUNDARY_DEPTH * 2, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {BOUNDARY_DEPTH, 0, 0};
  iter_type periodic_offset[3] = {0, P - 2 * BOUNDARY_DEPTH,
                                  M - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BOTTOM, Side::LEFT,
                                        Side::FRONT>::update(T* grid, int N,
                                                             int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {0, 0, 0};
  iter_type periodic_offset[3] = {
      N - 2 * BOUNDARY_DEPTH, P - 2 * BOUNDARY_DEPTH, M - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BOTTOM, Side::RIGHT,
                                        Side::FRONT>::update(T* grid, int N,
                                                             int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {N - BOUNDARY_DEPTH, 0, 0};
  iter_type periodic_offset[3] = {-(N - 2 * BOUNDARY_DEPTH),
                                  P - 2 * BOUNDARY_DEPTH,
                                  M - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BOTTOM, Side::LEFT,
                                        Side::BACK>::update(T* grid, int N,
                                                            int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {0, P - BOUNDARY_DEPTH, 0};
  iter_type periodic_offset[3] = {N - 2 * BOUNDARY_DEPTH,
                                  -(P - 2 * BOUNDARY_DEPTH),
                                  M - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BOTTOM, Side::RIGHT,
                                        Side::BACK>::update(T* grid, int N,
                                                            int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {N - BOUNDARY_DEPTH, P - BOUNDARY_DEPTH, 0};
  iter_type periodic_offset[3] = {-(N - 2 * BOUNDARY_DEPTH),
                                  -(P - 2 * BOUNDARY_DEPTH),
                                  M - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BOTTOM, Side::BOTTOM,
                                        Side::BACK>::update(T* grid, int N,
                                                            int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {N - BOUNDARY_DEPTH * 2, BOUNDARY_DEPTH};
  iter_type current_offset[3] = {BOUNDARY_DEPTH, P - BOUNDARY_DEPTH, 0};
  iter_type periodic_offset[3] = {0, -(P - 2 * BOUNDARY_DEPTH),
                                  M - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BOTTOM, Side::BOTTOM,
                                        Side::LEFT>::update(T* grid, int N,
                                                            int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, P - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {0, BOUNDARY_DEPTH, 0};
  iter_type periodic_offset[3] = {N - 2 * BOUNDARY_DEPTH, 0,
                                  M - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BOTTOM, Side::BOTTOM,
                                        Side::RIGHT>::update(T* grid, int N,
                                                             int M, int P) {
  iter_type stride[3] = {1, N * M, N};
  iter_type dims[2] = {BOUNDARY_DEPTH, P - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {N - BOUNDARY_DEPTH, BOUNDARY_DEPTH, 0};
  iter_type periodic_offset[3] = {-(N - 2 * BOUNDARY_DEPTH), 0,
                                  M - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::LEFT, Side::LEFT,
                                        Side::LEFT>::update(T* grid, int N,
                                                            int M, int P) {
  iter_type stride[3] = {N, N * M, 1};
  iter_type dims[2] = {M - BOUNDARY_DEPTH * 2, P - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH, 0};
  iter_type periodic_offset[3] = {0, 0, N - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::RIGHT, Side::RIGHT,
                                        Side::RIGHT>::update(T* grid, int N,
                                                             int M, int P) {
  iter_type stride[3] = {N, N * M, 1};
  iter_type dims[2] = {M - BOUNDARY_DEPTH * 2, P - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH,
                                 N - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {0, 0, -(N - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::FRONT, Side::FRONT,
                                        Side::FRONT>::update(T* grid, int N,
                                                             int M, int P) {
  iter_type stride[3] = {1, N, N * M};
  iter_type dims[2] = {N - BOUNDARY_DEPTH * 2, M - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH, 0};
  iter_type periodic_offset[3] = {0, 0, P - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::FRONT, Side::FRONT,
                                        Side::LEFT>::update(T* grid, int N,
                                                            int M, int P) {
  iter_type stride[3] = {1, N, N * M};
  iter_type dims[2] = {BOUNDARY_DEPTH, M - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {0, BOUNDARY_DEPTH, 0};
  iter_type periodic_offset[3] = {N - 2 * BOUNDARY_DEPTH, 0,
                                  P - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::FRONT, Side::FRONT,
                                        Side::RIGHT>::update(T* grid, int N,
                                                             int M, int P) {
  iter_type stride[3] = {1, N, N * M};
  iter_type dims[2] = {BOUNDARY_DEPTH, M - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {N - BOUNDARY_DEPTH, BOUNDARY_DEPTH, 0};
  iter_type periodic_offset[3] = {-(N - 2 * BOUNDARY_DEPTH), 0,
                                  P - 2 * BOUNDARY_DEPTH};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BACK, Side::BACK,
                                        Side::BACK>::update(T* grid, int N,
                                                            int M, int P) {
  iter_type stride[3] = {1, N, N * M};
  iter_type dims[2] = {N - BOUNDARY_DEPTH * 2, M - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {BOUNDARY_DEPTH, BOUNDARY_DEPTH,
                                 P - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {0, 0, -(P - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BACK, Side::BACK,
                                        Side::LEFT>::update(T* grid, int N,
                                                            int M, int P) {
  iter_type stride[3] = {1, N, N * M};
  iter_type dims[2] = {BOUNDARY_DEPTH, M - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {0, BOUNDARY_DEPTH, P - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {N - 2 * BOUNDARY_DEPTH, 0,
                                  -(P - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

template <>
template <typename T>
__device__ void update_periodic_cuda_3d<Side::BACK, Side::BACK,
                                        Side::RIGHT>::update(T* grid, int N,
                                                             int M, int P) {
  iter_type stride[3] = {1, N, N * M};
  iter_type dims[2] = {BOUNDARY_DEPTH, M - BOUNDARY_DEPTH * 2};
  iter_type current_offset[3] = {N - BOUNDARY_DEPTH, BOUNDARY_DEPTH,
                                 P - BOUNDARY_DEPTH};
  iter_type periodic_offset[3] = {-(N - 2 * BOUNDARY_DEPTH), 0,
                                  -(P - 2 * BOUNDARY_DEPTH)};
  general_periodic_iteration(grid, dims, stride, current_offset,
                             periodic_offset);
}

// Kernel to update left ghost cells
template <Side side0, Side side1, Side side2, typename T>
__global__ void call_update_periodic_cuda_3d(T* grid, int N, int M, int P) {
  update_periodic_cuda_3d<side0, side1, side2>::update(grid, N, M, P);
}
// Kernel to update left ghost cells
template <Side side0, Side side1, Side side2, typename T>
__global__ void call_update_periodic_cuda_3d_vec(T* grid0, T* grid1, T* grid2,
                                                 int N, int M, int P) {
  update_periodic_cuda_3d<side0, side1, side2>::update(grid0, grid1, grid2, N,
                                                       M, P);
}

template <Side side0, Side side1, Side side2>
struct boundary_surface_area;

template <>
struct boundary_surface_area<Side::LEFT, Side::LEFT, Side::LEFT> {
  int operator()(int N, int M, int P) {
    return (M - BOUNDARY_DEPTH * 2) * (P - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::RIGHT, Side::RIGHT, Side::RIGHT> {
  int operator()(int N, int M, int P) {
    return (M - BOUNDARY_DEPTH * 2) * (P - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::TOP, Side::TOP, Side::TOP> {
  int operator()(int N, int M, int P) {
    return (N - BOUNDARY_DEPTH * 2) * (P - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::TOP, Side::TOP, Side::LEFT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * (P - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::TOP, Side::TOP, Side::RIGHT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * (P - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::TOP, Side::TOP, Side::FRONT> {
  int operator()(int N, int M, int P) {
    return (N - BOUNDARY_DEPTH * 2) * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::TOP, Side::LEFT, Side::FRONT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::TOP, Side::RIGHT, Side::FRONT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::TOP, Side::TOP, Side::BACK> {
  int operator()(int N, int M, int P) {
    return (N - BOUNDARY_DEPTH * 2) * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::TOP, Side::LEFT, Side::BACK> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::TOP, Side::RIGHT, Side::BACK> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::BOTTOM, Side::BOTTOM, Side::BOTTOM> {
  int operator()(int N, int M, int P) {
    return (N - BOUNDARY_DEPTH * 2) * (P - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::BOTTOM, Side::BOTTOM, Side::LEFT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * (P - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::BOTTOM, Side::BOTTOM, Side::RIGHT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * (P - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::BOTTOM, Side::BOTTOM, Side::FRONT> {
  int operator()(int N, int M, int P) {
    return (N - BOUNDARY_DEPTH * 2) * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::BOTTOM, Side::LEFT, Side::FRONT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::BOTTOM, Side::RIGHT, Side::FRONT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::BOTTOM, Side::BOTTOM, Side::BACK> {
  int operator()(int N, int M, int P) {
    return (N - BOUNDARY_DEPTH * 2) * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::BOTTOM, Side::LEFT, Side::BACK> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::BOTTOM, Side::RIGHT, Side::BACK> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * BOUNDARY_DEPTH;
  }
};

template <>
struct boundary_surface_area<Side::FRONT, Side::FRONT, Side::FRONT> {
  int operator()(int N, int M, int P) {
    return (N - BOUNDARY_DEPTH * 2) * (M - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::FRONT, Side::FRONT, Side::LEFT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * (M - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::FRONT, Side::FRONT, Side::RIGHT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * (M - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::BACK, Side::BACK, Side::BACK> {
  int operator()(int N, int M, int P) {
    return (N - BOUNDARY_DEPTH * 2) * (M - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::BACK, Side::BACK, Side::LEFT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * (M - BOUNDARY_DEPTH * 2);
  }
};

template <>
struct boundary_surface_area<Side::BACK, Side::BACK, Side::RIGHT> {
  int operator()(int N, int M, int P) {
    return BOUNDARY_DEPTH * (M - BOUNDARY_DEPTH * 2);
  }
};

namespace grid {

// Kernel to update left ghost cells
template <Side side0, Side side1, Side side2, typename T>
struct update_boundary_3d {
  static void update(T* grid, int N, int M, int P) {
    int surfaceArea = boundary_surface_area<side0, side1, side2>{}(N, M, P);
    int numBlocks = (surfaceArea + BLOCK_SIZE - 1) / BLOCK_SIZE;
    call_update_periodic_cuda_3d<side0, side1, side2> CUDA_KERNEL(
        numBlocks, BLOCK_SIZE)(grid, N, M, P);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }
};

// Kernel to update left ghost cells
template <Side side0, Side side1, Side side2, typename T>
struct update_boundary_3d<side0, side1, side2, any_vector_t<T, 3>> {
  static void update(T* (&grid)[3], int N, int M, int P) {
    int surfaceArea = boundary_surface_area<side0, side1, side2>{}(N, M, P);
    int numBlocks = (surfaceArea + BLOCK_SIZE - 1) / BLOCK_SIZE;
    call_update_periodic_cuda_3d_vec<side0, side1, side2> CUDA_KERNEL(
        numBlocks, BLOCK_SIZE)(grid[0], grid[1], grid[2], N, M, P);
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
  }
};

}  // namespace grid

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::FRONT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::FRONT, Side::FRONT, Side::FRONT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::FRONT, Side::FRONT, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::FRONT, Side::FRONT, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BACK, 2>::
operator()(const grid::Boundary<T, 2>*, GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::BACK, Side::BACK, Side::BACK, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::BACK, Side::BACK, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::BACK, Side::BACK, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 2>::
operator()(const grid::Boundary<T, 2>*, GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::LEFT, Side::LEFT, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::RIGHT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::RIGHT, Side::RIGHT, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::TOP, 2>::
operator()(const grid::Boundary<T, 2>*, GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::TOP, Side::TOP, Side::TOP, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::TOP, Side::TOP, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::TOP, Side::TOP, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::TOP, Side::TOP, Side::FRONT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::TOP, Side::LEFT, Side::FRONT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::TOP, Side::RIGHT, Side::FRONT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::TOP, Side::TOP, Side::BACK, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::TOP, Side::LEFT, Side::BACK, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::TOP, Side::RIGHT, Side::BACK, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BOTTOM,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::BOTTOM, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::FRONT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::BOTTOM, Side::LEFT, Side::FRONT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::BOTTOM, Side::RIGHT, Side::FRONT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::BACK, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::BOTTOM, Side::LEFT, Side::BACK, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  grid::update_boundary_3d<Side::BOTTOM, Side::RIGHT, Side::BACK, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::FRONT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::FRONT, Side::FRONT, Side::FRONT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::TOP, Side::TOP, Side::LEFT, T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::LEFT,
  // T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::BACK,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::BACK, Side::BACK, Side::BACK, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::TOP, Side::TOP, Side::LEFT, T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::LEFT,
  // T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::LEFT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::LEFT, Side::LEFT, Side::LEFT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::TOP, Side::TOP, Side::LEFT, T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::LEFT,
  // T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::RIGHT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::RIGHT, Side::RIGHT, Side::RIGHT, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::TOP, Side::TOP, Side::LEFT, T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::LEFT,
  // T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::TOP, 2>::
operator()(const grid::Boundary<T, 2>*, GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::TOP, Side::TOP, Side::TOP, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::TOP, Side::TOP, Side::LEFT, T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::LEFT,
  // T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::BOTTOM,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::BOTTOM, T>::update(
      grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::TOP, Side::TOP, Side::LEFT, T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
  // grid::update_boundary_3d<Side::BOTTOM, Side::BOTTOM, Side::LEFT,
  // T>::update(
  //     grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::FRONT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //          (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_FRONT_3A(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2]);

  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //              (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
  //          (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_FRONT_LEFT_ALL(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2]);

  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //              (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
  //          (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_FRONT_RIGHT_ALL(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::BACK,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //          (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_BACK_3A(
  //     { grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);

  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //              (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
  //          (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_BACK_RIGHT_ALL(
  //     { grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);

  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //              (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
  //          (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_BACK_LEFT_ALL(
  //     { grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::LEFT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  // ITER_GRID3_LEFT_3A(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::RIGHT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  // ITER_GRID3_RIGHT_3A(
  //     { grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::TOP,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  // BOUNDARY_DEPTH); ITER_GRID3_TOP_3A(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);

  // offset = grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  // BOUNDARY_DEPTH) +
  //          (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_LEFT_TOP_ALL(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);

  // offset = grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  // BOUNDARY_DEPTH) -
  //          (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_RIGHT_TOP_ALL(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::BOTTOM,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  // BOUNDARY_DEPTH); ITER_GRID3_BOTTOM_3A(
  //     { grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);

  // offset = -grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  // BOUNDARY_DEPTH) +
  //          (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_LEFT_BOTTOM_ALL(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);

  // offset = -grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  // BOUNDARY_DEPTH) -
  //          (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_RIGHT_BOTTOM_ALL(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::FRONT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //          (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_FRONT_3AA(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2]);

  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //              (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
  //          grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  //          BOUNDARY_DEPTH);
  // ITER_GRID3_FRONT_TOP_ALL(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2]);

  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //              (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
  //          grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  //          BOUNDARY_DEPTH);
  // ITER_GRID3_FRONT_BOTTOM_ALL(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::BACK,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //          (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  // ITER_GRID3_BACK_3AA(
  //     { grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);

  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //              (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
  //          grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  //          BOUNDARY_DEPTH);
  // ITER_GRID3_BACK_TOP_ALL(
  //     { grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);

  // offset = grid.dims[0] * grid.dims[1], grid.dims[2] *
  //              (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
  //          grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  //          BOUNDARY_DEPTH);
  // ITER_GRID3_BACK_BOTTOM_ALL(
  //     { grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::LEFT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  // ITER_GRID3_LEFT_3AA(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::RIGHT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  // ITER_GRID3_RIGHT_3AA(
  //     { grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::TOP,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  // BOUNDARY_DEPTH); ITER_GRID3_TOP_3AA(
  //     { grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::BOTTOM,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  // iter_type offset;
  // offset = grid.dims[0] * (grid.dims[1], grid.dims[2] - BOUNDARY_DEPTH -
  // BOUNDARY_DEPTH); ITER_GRID3_BOTTOM_3AA(
  //     { grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0], grid.dims[1],
  //     grid.dims[2], grid.dims[2]);
}
