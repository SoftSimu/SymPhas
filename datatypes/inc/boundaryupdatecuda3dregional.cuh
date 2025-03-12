
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

// Regional Boundaries

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC, Side::FRONT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(
      symphas::lib::side_list<Side::FRONT, Side::FRONT, Side::FRONT>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::LEFT, Side::TOP, Side::FRONT>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::RIGHT, Side::TOP, Side::FRONT>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::LEFT, Side::BOTTOM, Side::FRONT>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::RIGHT, Side::BOTTOM, Side::FRONT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::BACK, 2>::
operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(
      symphas::lib::side_list<Side::BACK, Side::BACK, Side::BACK>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::LEFT, Side::TOP, Side::BACK>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::RIGHT, Side::TOP, Side::BACK>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::LEFT, Side::BOTTOM, Side::BACK>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::RIGHT, Side::BOTTOM, Side::BACK>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::LEFT, 2>::
operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(
      symphas::lib::side_list<Side::LEFT, Side::LEFT, Side::LEFT>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::LEFT, Side::LEFT, Side::TOP>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::LEFT, Side::LEFT, Side::BOTTOM>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::LEFT, Side::LEFT, Side::BACK>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::LEFT, Side::LEFT, Side::FRONT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC, Side::RIGHT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(
      symphas::lib::side_list<Side::RIGHT, Side::RIGHT, Side::RIGHT>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::RIGHT, Side::RIGHT, Side::TOP>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::RIGHT, Side::RIGHT, Side::BOTTOM>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::RIGHT, Side::RIGHT, Side::BACK>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::RIGHT, Side::RIGHT, Side::FRONT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC, Side::TOP, 2>::
operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(
      symphas::lib::side_list<Side::TOP, Side::TOP, Side::TOP>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::TOP, Side::TOP, Side::FRONT>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::TOP, Side::TOP, Side::BACK>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC, Side::BOTTOM,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(
      symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM, Side::BOTTOM>{},
      grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM, Side::FRONT>{}, grid);
  regional_update_boundary(
      symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM, Side::BACK>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3YZ, Side::FRONT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::FRONT, Side::FRONT>{},
                           grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3YZ, Side::BACK,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::BACK, Side::BACK>{},
                           grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3YZ, Side::LEFT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3YZ, Side::RIGHT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3YZ, Side::TOP,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3YZ, Side::BOTTOM,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XZ, Side::FRONT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::FRONT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XZ, Side::BACK,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::BACK>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XZ, Side::LEFT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::LEFT, Side::LEFT>{},
                           grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XZ, Side::RIGHT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT, Side::RIGHT>{},
                           grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XZ, Side::TOP,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XZ, Side::BOTTOM,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XY, Side::FRONT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XY, Side::BACK,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XY, Side::LEFT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XY, Side::RIGHT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XY, Side::TOP,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::TOP, Side::TOP>{},
                           grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC3XY, Side::BOTTOM,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(
      symphas::lib::side_list<Side::BOTTOM, Side::BOTTOM>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC0, Side::FRONT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::FRONT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC0, Side::BACK,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::BACK>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC0, Side::LEFT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::LEFT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC0, Side::RIGHT,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::RIGHT>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<BoundaryType::PERIODIC0, Side::TOP, 2>::
operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::TOP>{}, grid);
}

template <>
template <typename T>
void symphas::internal::update_boundary<
    BoundaryType::PERIODIC0, Side::BOTTOM,
    2>::operator()(const grid::Boundary<T, 2>*, RegionalGridCUDA<T, 3>& grid) {
  regional_update_boundary(symphas::lib::side_list<Side::BOTTOM>{}, grid);
}
