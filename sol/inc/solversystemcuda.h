
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
 * MODULE:  sol
 * PURPOSE: Manages a group of provisional systems, used by the phase field
 * model in order to manage the provisional variables.
 *
 * ***************************************************************************
 */

#pragma once

#include "solversystem.h"

#ifdef USING_CUDA

//! The default phase field system.
/*!
 * The representation of a phase field system, storing the values of the
 * order parameter that is defined in a phase field problem.
 *
 * Unless explicitly specified, this phase field system type will be
 * used by the solver. It does not manage boundaries or other data.
 *
 * \tparam T The order parameter type.
 * \tparam D The order parameter dimension.
 */
template <typename T, size_t D>
using SolverSystemCUDA = SystemCUDA<T, D>;

template <typename T, size_t D>
struct SolverSystemFDCUDA : BoundarySystemCUDA<T, D> {
  using BoundarySystemCUDA<T, D>::dims;

  BoundaryGridCUDA<T, D> dframe;  // the working grid for the solver
  SolverSystemFDCUDA(symphas::init_data_type const& tdata,
                     symphas::interval_data_type const& vdata,
                     symphas::b_data_type const& bdata, size_t id = 0)
      : BoundarySystemCUDA<T, D>(tdata, vdata, bdata, id), dframe{dims} {}
  SolverSystemFDCUDA() : BoundarySystemCUDA<T, D>(), dframe{0} {}
};

template <typename T, size_t D>
struct SolverSystemFDwSDCUDA : RegionalSystemCUDA<T, D> {
  using RegionalSystemCUDA<T, D>::dims;

  RegionalGridCUDA<T, D> dframe;  // the working grid for the solver
  SolverSystemFDwSDCUDA(symphas::init_data_type const& tdata,
                        symphas::interval_data_type const& vdata,
                        symphas::b_data_type const& bdata, size_t id = 0)
      : RegionalSystemCUDA<T, D>(tdata, vdata, bdata, id), dframe{dims} {}
  SolverSystemFDwSDCUDA() : RegionalSystemCUDA<T, D>(), dframe{0} {}

  inline void update(iter_type index, double time) {
    RegionalSystemCUDA<T, D>::update(index, time);
    dframe.adjust(RegionalSystemCUDA<T, D>::region);
  }
};

DEFINE_BASE_DATA_INHERITED((typename T, size_t D), (SolverSystemFDCUDA<T, D>),
                           (BoundaryGridCUDA<T, D>))

DEFINE_BASE_DATA_INHERITED((typename T, size_t D),
                           (SolverSystemFDwSDCUDA<T, D>),
                           (RegionalGridCUDA<T, D>))

#endif
