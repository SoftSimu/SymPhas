
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

#ifdef BASIC_MODELS

#include "modelmacros.h"

#define dpsi dop(1)
#define psi op(1)

MODEL(NOCHANGE, (SCALAR), EVOLUTION(dpsi = 0_n))
LINK_WITH_NAME(NOCHANGE, NOCHANGE)

MODEL(MA, (SCALAR),
      EVOLUTION(dpsi = lap(psi) + (c(0) - 4_n * c(1) * psi * psi) * psi))
LINK_WITH_NAME(MA, MODELA)

MODEL(CONV, (SCALAR), EVOLUTION(dpsi = -smoothing(psi)))
LINK_WITH_NAME(CONV, CONVOLUTION)

#else

// #include "advancedmodeldefs.h"
#include "modeldefinitions.h"
// #include "pfcdefs.h"
// #include "modelacmms.h"
#endif
