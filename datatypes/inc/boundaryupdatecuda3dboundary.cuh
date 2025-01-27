
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

#include "boundaryupdate.h"

#ifdef USING_CUDA

#include <cuda_runtime.h>

// *********************************************************************
/* an implementation of all BOUNDARY ITERATION ALGORITHMS
 */

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::FRONT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * grid.dims[1] *
           (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                   grid.dims[1]);

  // edges

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_TOP({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                       grid.dims[1]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_BOTTOM({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                          grid.dims[1]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BACK, 2>::
operator()(const grid::Boundary<T, 2>*, GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * grid.dims[1] *
           (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                  grid.dims[1], grid.dims[2]);

  // edges

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_TOP({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                      grid.dims[1], grid.dims[2]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_BOTTOM({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                         grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 2>::
operator()(const grid::Boundary<T, 2>*, GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  ITER_GRID3_LEFT({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                  grid.dims[1], grid.dims[2]);

  // edges

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_LEFT({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                       grid.dims[1], grid.dims[2]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_LEFT({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                        grid.dims[1]);

  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_LEFT_TOP({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                      grid.dims[1], grid.dims[2]);

  offset = -grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_LEFT_BOTTOM({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                         grid.dims[1], grid.dims[2]);

  // corners

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_LEFT_TOP_FRONT({ grid[INDEX] = grid[INDEX + offset]; },
                            grid.dims[0], grid.dims[1]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_LEFT_BOTTOM_FRONT({ grid[INDEX] = grid[INDEX + offset]; },
                               grid.dims[0], grid.dims[1]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_TOP_LEFT({ grid[INDEX] = grid[INDEX - offset]; },
                           grid.dims[0], grid.dims[1], grid.dims[2]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_BOTTOM_LEFT({ grid[INDEX] = grid[INDEX - offset]; },
                              grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::RIGHT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  ITER_GRID3_RIGHT({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                   grid.dims[1], grid.dims[2]);

  // edges

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_RIGHT({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                         grid.dims[1]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_RIGHT({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                        grid.dims[1], grid.dims[2]);

  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_RIGHT_TOP({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                       grid.dims[1], grid.dims[2]);

  offset = -grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_RIGHT_BOTTOM({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                          grid.dims[1], grid.dims[2]);

  // corners

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_TOP_RIGHT({ grid[INDEX] = grid[INDEX + offset]; },
                             grid.dims[0], grid.dims[1]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_BOTTOM_RIGHT({ grid[INDEX] = grid[INDEX + offset]; },
                                grid.dims[0], grid.dims[1]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_RIGHT_TOP_BACK({ grid[INDEX] = grid[INDEX - offset]; },
                            grid.dims[0], grid.dims[1], grid.dims[2]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_RIGHT_BOTTOM_BACK({ grid[INDEX] = grid[INDEX - offset]; },
                               grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::TOP, 2>::
operator()(const grid::Boundary<T, 2>*, GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_TOP({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                 grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BOTTOM,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BOTTOM({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                    grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::FRONT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * grid.dims[1] *
           (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_ALL({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                       grid.dims[1]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::BACK,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * grid.dims[1] *
           (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_ALL({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                      grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::LEFT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  ITER_GRID3_LEFT_ALL({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                      grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::RIGHT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  ITER_GRID3_RIGHT_ALL({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                       grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::TOP, 2>::
operator()(const grid::Boundary<T, 2>*, GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_TOP_ALL({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                     grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::BOTTOM,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BOTTOM_ALL({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                        grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::FRONT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * grid.dims[1] *
           (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_3A({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                      grid.dims[1]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_LEFT_ALL({ grid[INDEX] = grid[INDEX + offset]; },
                            grid.dims[0], grid.dims[1]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_RIGHT_ALL({ grid[INDEX] = grid[INDEX + offset]; },
                             grid.dims[0], grid.dims[1]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::BACK,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * grid.dims[1] *
           (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_3A({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                     grid.dims[1], grid.dims[2]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_RIGHT_ALL({ grid[INDEX] = grid[INDEX - offset]; },
                            grid.dims[0], grid.dims[1], grid.dims[2]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_LEFT_ALL({ grid[INDEX] = grid[INDEX - offset]; },
                           grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::LEFT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  ITER_GRID3_LEFT_3A({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                     grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::RIGHT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  ITER_GRID3_RIGHT_3A({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                      grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::TOP,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_TOP_3A({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                    grid.dims[1], grid.dims[2]);

  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_LEFT_TOP_ALL({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                          grid.dims[1], grid.dims[2]);

  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_RIGHT_TOP_ALL({ grid[INDEX] = grid[INDEX + offset]; },
                           grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XZ, Side::BOTTOM,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BOTTOM_3A({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                       grid.dims[1], grid.dims[2]);

  offset = -grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_LEFT_BOTTOM_ALL({ grid[INDEX] = grid[INDEX + offset]; },
                             grid.dims[0], grid.dims[1], grid.dims[2]);

  offset = -grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           (grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_RIGHT_BOTTOM_ALL({ grid[INDEX] = grid[INDEX + offset]; },
                              grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::FRONT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * grid.dims[1] *
           (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_3AA({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                       grid.dims[1]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_TOP_ALL({ grid[INDEX] = grid[INDEX + offset]; },
                           grid.dims[0], grid.dims[1]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_FRONT_BOTTOM_ALL({ grid[INDEX] = grid[INDEX + offset]; },
                              grid.dims[0], grid.dims[1]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::BACK,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * grid.dims[1] *
           (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_3AA({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                      grid.dims[1], grid.dims[2]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) -
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_TOP_ALL({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                          grid.dims[1], grid.dims[2]);

  offset = grid.dims[0] * grid.dims[1] *
               (grid.dims[2] - BOUNDARY_DEPTH - BOUNDARY_DEPTH) +
           grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BACK_BOTTOM_ALL({ grid[INDEX] = grid[INDEX - offset]; },
                             grid.dims[0], grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::LEFT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  ITER_GRID3_LEFT_3AA({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                      grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::RIGHT,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] - BOUNDARY_DEPTH - BOUNDARY_DEPTH;
  ITER_GRID3_RIGHT_3AA({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                       grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::TOP,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_TOP_3AA({ grid[INDEX] = grid[INDEX + offset]; }, grid.dims[0],
                     grid.dims[1], grid.dims[2]);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC3XY, Side::BOTTOM,
                                        2>::operator()(const grid::Boundary<T,
                                                                            2>*,
                                                       GridCUDA<T, 3>& grid) {
  iter_type offset;
  offset = grid.dims[0] * (grid.dims[1] - BOUNDARY_DEPTH - BOUNDARY_DEPTH);
  ITER_GRID3_BOTTOM_3AA({ grid[INDEX] = grid[INDEX - offset]; }, grid.dims[0],
                        grid.dims[1], grid.dims[2]);
}

#endif
