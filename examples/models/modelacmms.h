
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

#define eta op(1)
#define deta dop(1)

#define kappa c(1)
#define A1 c(2)
#define B1 c(3)
#define A2 c(4)
#define B2 c(5)
#define C2 c(6)

MODEL(AC_MMS, (SCALAR),
	EVOLUTION_PREAMBLE(
		(
			auto th = tanh(sqrt(2) * (A1 * t * sin(B1 * x) + A2 * sin(B2 * x + C2 * t) - y + 0.25) / (2 * sqrt(kappa)));
			auto a0 = (A1 * B1 * t * cos(B1 * x) + A2 * B2 * cos(B2 * x + C2 * t));
			auto S = -(th * th - one) * 
				(0.5 * sqrt(kappa) * (a0 * a0 + one) * tanh(sqrt(2) * (A1 * t * sin(B1 * x) + A2 * sin(B2 * x + C2 * t) - y + 0.25) / (2 * sqrt(kappa))) 
					- 0.5 * sqrt(kappa) * tanh(sqrt(2) * (A1 * t * sin(B1 * x) + A2 * sin(B2 * x + C2 * t) - y + 0.25) / (2 * sqrt(kappa))) 
					+ 0.25 * sqrt(2) * kappa * (A1 * B1 * B1 * t * sin(B1 * x) + A2 * B2 * B2 * sin(B2 * x + C2 * t)) 
					+ 0.25 * sqrt(2) * (A1 * sin(B1 * x) + A2 * C2 * cos(B2 * x + C2 * t))
					) / sqrt(kappa);
		),
		deta = -(4 * eta * (eta - 1) * (eta - 0.5) - kappa * lap(eta)) + S
	)
	RESTRICT_DIMENSIONS(2)
)
LINK_WITH_NAME(AC_MMS, AC_MMS)

#define end_time lit(8.)
INITIAL_CONDITION_EQUATION(
	MMS_AC_INIT, (2, 3), 
	-0.5 * tanh(sqrt(2) * (-A2 * sin(B2 * x) + y - 0.25) / (2 * sqrt(kappa))) + 0.5)
INITIAL_CONDITION_EQUATION(
	MMS_AC_FINAL, (2, 3), 
	-0.5 * tanh(sqrt(2) * (-A1 * end_time * sin(B1 * x) - A2 * sin(B2 * x + C2 * end_time) + y - 0.25) / (2 * sqrt(kappa))) + 0.5)

#undef kappa
#undef A1
#undef B1
#undef A2
#undef B2
#undef C2
#undef SQRT2
