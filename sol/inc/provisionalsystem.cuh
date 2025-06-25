
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
 * PURPOSE: Defines a provisional system. A provisional system stores values
 * of the provisional variables which can constitute part of the model
 * definition. Provisional variables can act as virtual variables; they
 * are computed at each step and can be used in the equations of motion, but
 * are not physical variables like the order parameters.
 *
 * ***************************************************************************
 */

#pragma once

#include "boundarysystem.cuh"

#ifdef USING_CUDA
template <typename T, size_t D>
struct ProvisionalSystemFDCUDA : BoundarySystemCUDA<T, D> {
  //! Create a provisional system.
  /*!
   * Create a provisional system using the given intervals and boundary data.
   */
  ProvisionalSystemFDCUDA(const symphas::interval_data_type vdata,
                          const symphas::b_data_type bdata)
      : BoundarySystemCUDA<T, D>(provisional_init, vdata, bdata) {}
};
#endif
