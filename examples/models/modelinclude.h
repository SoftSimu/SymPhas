
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

// By default a basic set of models is compiled. Define USE_EXTENDED_MODELS
// (optionally together with MODEL_SET_n and USE_PFC_MODELS) to opt into the
// extended model definitions in modeldefinitions.h / pfcdefs.h.
#if !defined(BASIC_MODELS) && !defined(USE_EXTENDED_MODELS)
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

// #include "advancedmodeldefs.h"
#include "modeldefinitions.h"
#ifdef USE_PFC_MODELS
#include "pfcdefs.h"
#endif
// #include "modelacmms.h"
#endif
