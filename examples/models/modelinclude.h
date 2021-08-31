
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


namespace simulate
{
	template<typename M>
	int simulate(double const *coeff, size_t num_coeff);
}


#ifdef DEBUG

#include "modelmacros.h"

#define dpsi dop(1)
#define psi op(1)

MODEL(NOCHANGE, (SCALAR),
	MODEL_DEF(dpsi = lit(0))
)
LINK_WITH_NAME(NOCHANGE, NOCHANGE)

MODEL(MA, (SCALAR),
	MODEL_DEF(
		dpsi = lap(psi) + (c1 - lit(4.) * c2 * psi * psi) * psi)
)
LINK_WITH_NAME(MA, MODELA)

MODEL(CONV, (SCALAR),
	MODEL_DEF(
		dpsi = -smoothing(psi))
)
LINK_WITH_NAME(CONV, CONVOLUTION)


#else

#include "pfcdefs.h"
#include "modeldefs.h"
#endif


