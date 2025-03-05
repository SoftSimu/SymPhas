
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
 * PURPOSE: Defines the boundary class used for the finite difference
 * implementation of solvers, i.e. when central differencing is used for
 * all derivative approximations.
 *
 * ***************************************************************************
 */

#pragma once

#include "boundary.h"
#include "gridcuda.cuh"

template <typename T>
template <typename vector_type, typename T0>
void grid::BoundaryApplied<T, 0, BoundaryType::DEFAULT>::update(
    multi_value_cuda<1, T0> val, axis_coord_t, axis_coord_t,
    double time) const {
  any_vector_t<T0, 1> vector = val;
  this->update(vector, 0, 0, time);
  val = vector;
}

template <typename T>
void grid::BoundaryApplied<T, 0, BoundaryType::DEFAULT>::update(
    carry_value_cuda<T> val, axis_coord_t, axis_coord_t, double time) const {
  if (val.is_valid()) {
    T value;
    CHECK_CUDA_ERROR(
        cudaMemcpy(&value, val.value, sizeof(T), cudaMemcpyDeviceToHost));
    update(value, 0, 0, time);
    CHECK_CUDA_ERROR(
        cudaMemcpy(val.value, &value, sizeof(T), cudaMemcpyHostToDevice));
  }
}

template <typename T>
template <typename vector_type, typename T0>
void grid::BoundaryApplied<T, 1, BoundaryType::DEFAULT>::update(
    multi_value_cuda<2, T0> val, axis_coord_t x, axis_coord_t,
    double time) const {
  any_vector_t<T0, 2> vector = val;
  this->update(vector, x, 0, time);
  val = vector;
}

template <typename T>
void grid::BoundaryApplied<T, 1, BoundaryType::DEFAULT>::update(
    carry_value_cuda<T> val, axis_coord_t x, axis_coord_t, double time) const {
  if (val.is_valid()) {
    T value;
    CHECK_CUDA_ERROR(
        cudaMemcpy(&value, val.value, sizeof(T), cudaMemcpyDeviceToHost));
    update(value, x, 0, time);
    CHECK_CUDA_ERROR(
        cudaMemcpy(val.value, &value, sizeof(T), cudaMemcpyHostToDevice));
  }
}

template <typename T>
template <typename vector_type, typename T0>
void grid::BoundaryApplied<T, 2, BoundaryType::DEFAULT>::update(
    multi_value_cuda<3, T0> val, axis_coord_t x, axis_coord_t y,
    double time) const {
  any_vector_t<T0, 2> vector = val;
  this->update(vector, x, y, time);
  val = vector;
}

template <typename T>
void grid::BoundaryApplied<T, 2, BoundaryType::DEFAULT>::update(
    carry_value_cuda<T> val, axis_coord_t x, axis_coord_t y,
    double time) const {
  if (val.is_valid()) {
    T value;
    CHECK_CUDA_ERROR(
        cudaMemcpy(&value, val.value, sizeof(T), cudaMemcpyDeviceToHost));
    update(value, x, y, time);
    CHECK_CUDA_ERROR(
        cudaMemcpy(val.value, &value, sizeof(T), cudaMemcpyHostToDevice));
  }
}
