
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
 * This file is responsible for including all the model definitions and
 * forward declaring the function that will be used to run the model
 * simulation workflow. This is enabled and specified through the CMake
 * configuration.
 *
 * ***************************************************************************
 */

#pragma once

// By default a basic set of models is compiled. Define one or more of
//   USE_EXTENDED_MODELS  - the H&H family (Models A-F, FE variants, etc.)
//   USE_PFC_MODELS       - the PFC family (PFC_C, PFC_NC, PFC_NC4)
//   USE_GRID_HH          - minimal H&H subset (MA, MB, MC, MH, MF) used for
//                          the validation grid; cheaper to compile than the
//                          full USE_EXTENDED_MODELS set.
// to opt into the extended model definitions. Either may be used alone or
// together. When neither is set the basic in-line definitions below are used.
#if !defined(BASIC_MODELS) && !defined(USE_EXTENDED_MODELS) && \
    !defined(USE_PFC_MODELS) && !defined(USE_GRID_HH)
#define BASIC_MODELS
#endif

#ifdef BASIC_MODELS

#include "modelmacros.h"

#define dpsi dop(1)
#define psi op(1)
#define drho dop(2)
#define rho op(2)

MODEL(NOCHANGE, (SCALAR), EVOLUTION(dpsi = 0_n))
LINK_WITH_NAME(NOCHANGE, NOCHANGE)

MODEL(MA, (SCALAR),
      EVOLUTION(dpsi = lap(psi) + (c(0) - 4_n * c(1) * psi * psi) * psi))
LINK_WITH_NAME(MA, MODELA)

//! Model B (Cahn-Hilliard).
MODEL(MB, (SCALAR),
      EVOLUTION(dpsi = -bilap(psi) - lap((c(1) - c(2) * psi * psi) * psi)))
LINK_WITH_NAME(MB, MODELB)

//! Model C (two-field: conserved + non-conserved).
MODEL(MC, (SCALARS(2)),
      EVOLUTION(dpsi = -bilap(psi) -
                       lap((c(1) - c(2) * psi * psi) * psi + c(5) * rho * rho),
                drho = lap(rho) + (c(3) - c(4) * rho * rho) * rho +
                       2_n * c(5) * psi * rho))
LINK_WITH_NAME(MC, MODELC)
DEFINE_MODEL_FIELD_NAMES(MC, ("psi", "m"))

#undef dpsi
#undef psi
#undef drho
#undef rho

#else

#ifdef USE_GRID_HH
// Minimal HH subset for the validation grid: MA, MB, MC, MH, MF only.
#include "modelmacros.h"
#define dpsi dop(1)
#define psi op(1)
#define dj dop(2)
#define j op(2)

MODEL(MA, (SCALAR),
      EVOLUTION(dpsi = lap(psi) + (c(0) - 4_n * c(1) * psi * psi) * psi))
LINK_WITH_NAME(MA, MODELA)

MODEL(MB, (SCALAR),
      EVOLUTION(dpsi = -bilap(psi) - lap((c(1) - c(2) * psi * psi) * psi)))
LINK_WITH_NAME(MB, MODELB)

MODEL(MC, (SCALARS(2)),
      EVOLUTION(dpsi = -bilap(psi) -
                       lap((c(1) - c(2) * psi * psi) * psi + c(5) * j * j),
                dj    = lap(j) + (c(3) - c(4) * j * j) * j +
                        2_n * c(5) * psi * j))
LINK_WITH_NAME(MC, MODELC)
DEFINE_MODEL_FIELD_NAMES(MC, ("psi", "m"))

#ifndef USE_GRID_HH_SCALAR_ONLY
MODEL(MH, (SCALAR, VECTOR),
      EVOLUTION_PREAMBLE((auto f = lap(psi) + (c(1) - c(2) * psi * psi) * psi;),
                         dpsi = -lap(f) - c(3) * grad * (psi * j),
                         dj = lap(j) - c(3) * psi * grad(f)))
LINK_WITH_NAME(MH, MODELH)
DEFINE_MODEL_FIELD_NAMES(MH, ("m", "j"))

MODEL(MF, (SCALAR, VECTOR),
      EVOLUTION(
            dpsi = c(4) * lap(psi) + (c(1) - c(2) * psi * psi) * psi
                   - c(3) * div(j),
            dj = lap(j) + c(3) * grad(psi)
      ))
LINK_WITH_NAME(MF, MODELF)
DEFINE_MODEL_FIELD_NAMES(MF, ("n", "g"))
#endif

#undef dpsi
#undef psi
#undef dj
#undef j
#endif

#ifdef USE_EXTENDED_MODELS
// #include "advancedmodeldefs.h"
#include "modeldefinitions.h"
// #include "modelacmms.h"
#endif

#ifdef USE_PFC_MODELS
#include "pfcdefs.h"
#endif

#endif
