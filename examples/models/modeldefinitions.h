
/* ***************************************************************************
 * This file is part of the SymPhas package, containing a framework for
 * implementing solvers for phase-field problems with compile-time symbolic
 * algebra.
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
 * This file is responsible for including all the phase-field model examples.
 *
 * ***************************************************************************
 */

#pragma once

#include "modelmacros.h"

#define dpsi dop(1)
#define psi op(1)
#define drho dop(2)
#define rho op(2)
#define j rho
#define dj drho


#ifndef BASIC_MODELS

#ifdef MODEL_SET_1
//! Model A with noise.
MODEL(MAWN, (SCALAR),
      EVOLUTION(dpsi = lap(psi) + (c(1) - c(2) * psi * psi) * psi +
                       _nW(SCALAR)))
LINK_WITH_NAME(MAWN, MODELA_WN)

//! Model A.
MODEL(MA, (SCALAR),
      EVOLUTION(dpsi = lap(psi) + (c(1) - c(2) * psi * psi) * psi))
LINK_WITH_NAME(MA, MODELA)

//! Model B.
MODEL(MB, (SCALAR),
      EVOLUTION(dpsi = -bilap(psi) - lap((c(1) - c(2) * psi * psi) * psi)))
LINK_WITH_NAME(MB, MODELB)

//! Model B.
MODEL(MBWN, (SCALAR),
      EVOLUTION(dpsi = -c(3) * bilap(psi) -
                       lap((c(1) - c(2) * psi * psi) * psi + _nW(SCALAR))))
LINK_WITH_NAME(MBWN, MODELB_WN)

//! Model C.
MODEL(MCWN, (SCALARS(2)),
      EVOLUTION(dpsi = -bilap(psi) -
                       lap((c(1) - c(2) * psi * psi + _nW(SCALAR)) * psi +
                           c(5) * rho * rho),
                drho = lap(rho) +
                       (c(3) - c(4) * rho * rho + _nW(SCALAR)) * rho +
                       2_n * c(5) * psi * rho))
LINK_WITH_NAME(MCWN, MODELC)
DEFINE_MODEL_FIELD_NAMES(MCWN, ("psi", "m"))

//! Model C.
MODEL(MC, (SCALARS(2)),
      EVOLUTION(dpsi = -bilap(psi) -
                       lap((c(1) - c(2) * psi * psi) * psi + c(5) * rho * rho),
                drho = lap(rho) + (c(3) - c(4) * rho * rho) * rho +
                       2_n * c(5) * psi * rho))
LINK_WITH_NAME(MC, MODELC)
DEFINE_MODEL_FIELD_NAMES(MC, ("psi", "m"))

// Model H
MODEL(MH, (SCALAR, VECTOR),
      EVOLUTION_PREAMBLE(
          (auto f = lap(psi) + (c(1) - c(2) * psi * psi) * psi;),
          dpsi = -lap(f) - c(3) * grad(psi) * j + lap(_nW(SCALAR, -2)),
          dj = lap(j) + c(3) * grad(psi) * f + lap(_nW(VECTOR, 2))))
LINK_WITH_NAME(MH, MODELH)
DEFINE_MODEL_FIELD_NAMES(MH, ("m", "j"))

#endif

#ifdef MODEL_SET_2

//! Example of coupling through the free energy using an iterated sum.
MODEL(MH_FE, (SCALAR, VECTOR),
      FREE_ENERGY((EQUATION_OF(1)(-lap(-DF(1)) - grad(op(1)) * DF(2)),
                   EQUATION_OF(2)(lap(DF(2)) + grad(op(1)) * -DF(1))),
                  INT(LANDAU_FE(op(1)) + _2 * pow<2>(op(2)))))
LINK_WITH_NAME(MH_FE, MODELH_FE)
DEFINE_MODEL_FIELD_NAMES(MH_FE, ("m", "j"))

//! Model A defined by the free energy, where the free energy uses
//! a general index for any number of fields that are defined.
MODEL(MA_FE, (SCALAR),
      FREE_ENERGY((NONCONSERVED), INT(SUM(ii)(LANDAU_FE(op_ii, c(1), c(2))))))
LINK_WITH_NAME(MA_FE, MODELA_FE)

//! Example of coupling through the free energy using an iterated sum.

//! Model B defined by the free energy.
MODEL(MB_FE, (SCALAR), FREE_ENERGY((CONSERVED), INT(SUM(ii)(LANDAU_FE(op_ii)))))
LINK_WITH_NAME(MB_FE, MODELB_FE)

//! Model C defined by the free energy.
MODEL(MC_FE, (SCALAR, SCALAR),
      FREE_ENERGY((NONCONSERVED, CONSERVED),
                  INT(LANDAU_FE(psi) + LANDAU_FE(rho) + psi * psi * rho)))
LINK_WITH_NAME(MC_FE, MODELC_FE)

#endif

#ifdef MODEL_SET_3

//! The Gray-Scott phase field model.
MODEL(GRAYSCOTT, (SCALAR, SCALAR),
      EVOLUTION_PREAMBLE((auto prr = psi * rho * rho;),
                         dpsi = c(1) * lap(psi) - prr + c(3) * (one - psi),
                         drho = c(2) * lap(rho) + prr - (c(3) + c(4)) * rho))
LINK_WITH_NAME(GRAYSCOTT, GRAYSCOTT)

// Turing Model
MODEL(Turing, (SCALAR, SCALAR),
      EVOLUTION(dpsi = c(1) * lap(psi) +
                       c(2) * (psi + c(3) * rho - psi * rho * rho -
                               c(4) * psi * rho),
                drho = lap(rho) + c(2) * (c(5) * rho + c(6) * psi +
                                          psi * rho * rho + c(4) * psi * rho)))
LINK_WITH_NAME(Turing, TURING)

#endif

#ifdef MODEL_SET_4

MODEL(COUPLING4, (SCALARS(4)),
      FREE_ENERGY((ALL_NONCONSERVED(ii)),
                  INT(SUM(ii)(LANDAU_FE(op_ii)) +
                      _4 * SUM(ii, jj != ii)(op_ii * op_jj * op_jj) + op(1))))
LINK_WITH_NAME(COUPLING4, MODEL_COUPLING4)
DEFINE_MODEL_FIELD_NAMES(COUPLING4, ("A", "B", "C", "D"))

//! Example of provisional variables, Model B. Note the <= sign when
//! creating an expression for setting the provisional variable.
MODEL(MBB, (SCALAR),
      PROVISIONAL_DEF((SCALAR),
                      var(1) <= -(c(1) - lit(4.) * c(2) * psi * psi) * psi)
          EVOLUTION(dpsi = -bilap(psi) + lap(var(1))))
LINK_WITH_NAME(MBB, MODELBB)

#endif

#endif

#undef j
#undef dj
#undef rho
#undef psi
#undef dpsi
#undef drho

// *******************************************************************************
