
/* ***************************************************************************
 * This file is part of the SymPhas library, a framework for implementing
 * solvers for phase-field problems with compile-time symbolic algebra.
 *
 * Copyright (c) 2018-2020 by Steven A. Silber and Mikko Karttunen
 *
 * SymPhas is free software, which can be redistributed or modified under
 * the terms of the GNU Lesser General Public License (LGPL) as published
 * by the Free Software Foundation; LGPL version 3, or later versions at
 * your choice
 *
 * SymPhas is distributed with the faith that it will be helpful and
 * practical but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details
 *
 * ***************************************************************************
 *
 * MODULE:  sol
 * PURPOSE:
 *
 * The stencil object defines the three finite difference functions, and
 * uses templates to specialize each one by dimension and order, for the
 * gradient, laplacian, bilaplacian and gradlaplacian.
 *
 * A specialized template for the stencil object in turn specializes each
 * function to take advantage of generic programming.
 *
 * A generic stencil can also be specialized for any stencils which are not
 * included in the specialization of first to fourth order derivatives, and
 * allows higher order derivatives to be programmatically or manually computed
 * by specializing the dimension, accuracy and derivative order.
 *
 * ***************************************************************************
 */

#pragma once

#include "gridcuda.cuh"
#include "stencil.h"

template <typename Sp>
template <Axis ax, size_t O, typename T, size_t D>
__host__ __device__ auto
Stencil<Sp>::applied_generalized_directional_derivative(
    RegionalGridCUDA<T, D> const &grid, iter_type n) const {
  auto v = grid[n];
  if (v.is_valid()) {
    return cast().template apply_directional<ax, O>(v.value, grid.region.dims);
  } else {
    return *v.value;
  }
}

template <typename Sp>
template <size_t... Os, typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_generalized_mixed_derivative(
    RegionalGridCUDA<T, D> const &grid, iter_type n) const {
  auto v = grid[n];
  if (v.is_valid()) {
    return cast().template apply_mixed<Os...>(v.value, grid.region.dims);
  } else {
    return *v.value;
  }
}

template <typename Sp>
template <Axis ax, size_t O, typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_generalized_derivative(
    RegionalGridCUDA<T, D> const &grid, iter_type n) const {
  auto v = grid[n];
  if (v.is_valid()) {
    len_type stride[D];
    grid::get_stride<ax>(stride, grid.region.dims);
    return cast().template apply<O>(v.value, stride);
  } else {
    return *v.value;
  }
}

template <typename Sp>
template <typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_laplacian(
    RegionalGridCUDA<T, D> const &grid, iter_type n) const {
  auto v = grid[n];
  if (v.is_valid()) {
    return laplacian(v.value, grid.region.stride);
  } else {
    return *v.value;
  }
}

template <typename Sp>
template <typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_bilaplacian(
    RegionalGridCUDA<T, D> const &grid, iter_type n) const {
  auto v = grid[n];
  if (v.is_valid()) {
    return bilaplacian(v.value, grid.region.stride);
  } else {
    return *v.value;
  }
}

template <typename Sp>
template <Axis ax, typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_gradlaplacian(
    RegionalGridCUDA<T, D> const &grid, iter_type n) const {
  auto v = grid[n];
  if (v.is_valid()) {
    len_type stride[D];
    grid::get_stride<ax>(stride, grid.region.dims);
    return gradlaplacian(v.value, stride);
  } else {
    return *v.value;
  }
}

template <typename Sp>
template <Axis ax, typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_gradient(
    RegionalGridCUDA<T, D> const &grid, iter_type n) const {
  auto v = grid[n];
  if (v.is_valid()) {
    len_type stride[D];
    grid::get_stride<ax>(stride, grid.region.dims);
    return gradient(v.value, stride);
  } else {
    return *v.value;
  }
}

template <typename Sp>
template <Axis ax, size_t O, typename T, size_t D>
__host__ __device__ auto
Stencil<Sp>::applied_generalized_directional_derivative(
    RegionGridDataCUDA<T, D> const &grid, iter_type n) const {
  return cast().template apply_directional<ax, O>(&grid[n], grid.dims);
}

template <typename Sp>
template <size_t... Os, typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_generalized_mixed_derivative(
    RegionGridDataCUDA<T, D> const &grid, iter_type n) const {
  return cast().template apply_mixed<Os...>(&grid[n], grid.dims);
}

template <typename Sp>
template <Axis ax, size_t O, typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_generalized_derivative(
    RegionGridDataCUDA<T, D> const &grid, iter_type n) const {
  return cast().template apply<O>(&grid[n], grid.stride);
}

template <typename Sp>
template <typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_laplacian(
    RegionGridDataCUDA<T, D> const &grid, iter_type n) const {
  return laplacian(&grid[n], grid.stride);
}

template <typename Sp>
template <typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_bilaplacian(
    RegionGridDataCUDA<T, D> const &grid, iter_type n) const {
  return bilaplacian(&grid[n], grid.stride);
}

template <typename Sp>
template <Axis ax, typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_gradlaplacian(
    RegionGridDataCUDA<T, D> const &grid, iter_type n) const {
  len_type stride[D];
  grid::get_stride<ax>(stride, grid.dims);
  return gradlaplacian(&grid[n], stride);
}

template <typename Sp>
template <Axis ax, typename T, size_t D>
__host__ __device__ auto Stencil<Sp>::applied_gradient(
    RegionGridDataCUDA<T, D> const &grid, iter_type n) const {
  len_type stride[D];
  grid::get_stride<ax>(stride, grid.dims);
  return gradient(&grid[n], stride);
}