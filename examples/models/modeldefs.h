
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

//! Model A.
MODEL(MA, (SCALAR),
	MODEL_DEF(
		dpsi = lap(psi) + (c1 - c2 * psi * psi) * psi)
)
LINK_WITH_NAME(MA, MODELA)

//! The Gray-Scott phase field model.
MODEL(GRAYSCOTT, (SCALAR, SCALAR),
	MODEL_DEF_PREAMBLE((auto prr = psi * rho * rho;),
		dpsi = c1 * lap(psi) - prr + c3 * (one - psi),
		drho = c2 * lap(rho) + prr - (c3 + c4) * rho)
)
LINK_WITH_NAME(GRAYSCOTT, GRAYSCOTT)

//! Model B.
MODEL(MB, (SCALAR),
	MODEL_DEF(
		dpsi = -bilap(psi) - lap((c1 - c2 * psi * psi) * psi))
)
LINK_WITH_NAME(MB, MODELB)

#ifndef BASIC_MODELS

//! Model C.
MODEL(MC, (SCALAR, SCALAR),
	MODEL_DEF_PREAMBLE(
		(
			auto psi3 = psi * psi * psi;
			auto rho3 = rho * rho * rho;
		),
		dpsi = lap(psi) + c1 * psi - c2 * psi3 + lit(2.) * c5 * psi * rho,
		drho = -bilap(rho) - lap(c3 * rho - c4 * rho3 + c5 * psi * psi))
)
LINK_WITH_NAME(MC, MODELC)
DEFINE_MODEL_FIELD_NAMES(MC, ("psi", "m"))

/*

//! Model C with provisional variables.
MODEL(MCC, (SCALAR, SCALAR),
	PROVISIONAL_DEF((SCALAR, SCALAR),
		var(1) = -(c1 - lit(4.) * c2 * psi * psi) * psi - lit(2.) * c5 * psi * rho,
		var(2) = -(c3 - lit(4.) * c4 * rho * rho) * rho - c5 * psi * psi)
	MODEL_DEF(
		dpsi = lap(psi) - var(1),
		drho = -bilap(rho) + lap(var(2)))
)
LINK_WITH_NAME(MCC, MODELCC)

//! Model B with provisional variables.
MODEL(MBB, (SCALAR),
	PROVISIONAL_DEF((SCALAR),
		var(1) = -(c1 - lit(4.) * c2 * psi * psi) * psi)
	MODEL_DEF(
		dpsi = -bilap(psi) + lap(var(1)))
)
LINK_WITH_NAME(MBB, MODELBB)

// Model D
MODEL(ME, (COMPLEX, SCALAR),
	PROVISIONAL_DEF((COMPLEX, SCALAR), 
		var(1) = -(c1 - lit(4.) * c2 * psi * conj(psi)) * psi,
		var(2) = -(c3 - lit(4.) * c4 * rho * rho) * rho)
	MODEL_DEF(
		dpsi = lit(2.) * (lap(psi) - var(1)) - c5 * Ii * psi * (-c4 * lap(rho) + var(2)),
		drho = -c6 * (c4 * bilap(rho) - lap(var(2))) - lit(2.) * c5 * Im(conj(psi) * (-lap(psi) + var(1))))
)
LINK_WITH_NAME(ME, MODELE)

// Model H
MODEL(MH, (SCALAR, VECTOR),
	PROVISIONAL_DEF((SCALAR, VECTOR), 
		var(1) = -(c1 - lit(4.) * c2 * psi * psi) * psi,
		var(2) = -(c3 - lit(4.) * c4 * rho * rho) * rho)
	MODEL_DEF(
		dpsi = -bilap(psi) + lap(var(1)) - lit(2.) * c5 * grad(psi) * (-lap(rho) + var(2)),
		drho = -bilap(rho) + lap(var(2)) + c5 * grad(psi) * (-lap(psi) + var(1)))
)
LINK_WITH_NAME(MH, MODELH)


// Model F
MODEL(MF, (COMPLEX, SCALAR),
	PROVISIONAL_DEF((COMPLEX, SCALAR), 
		var(1) = -(c1 - lit(4.) * c2 * psi * conj(psi)) * psi - lit(2.) * c5 * psi * rho,
		var(2) = -(c3 - lit(4.) * c4 * rho * rho) * rho - c5 * Re(psi * conj(psi)))
	MODEL_DEF(
		dpsi = lit(2.) * (lap(psi) - var(1)) - c6 * Ii * psi * (-lap(rho) + var(2)),
		drho = -c7 * (bilap(rho) - lap(var(2))) - lit(2.) * c6 * Im(conj(psi) * (-lap(psi) + var(1))))
)
LINK_WITH_NAME(MF, MODELF)

// Model G
MODEL(MG, (VECTOR, VECTOR),
	PROVISIONAL_DEF((VECTOR, VECTOR), 
		var(1) = -(c1 - lit(4.) * c2 * psi * psi) * psi,
		var(2) = -(c3 - lit(4.) * c4 * rho * rho) * rho)
	MODEL_DEF(
		dpsi = lap(psi) - var(1) + psi % (lap(rho) - var(2)),
		drho = -bilap(rho) + lap(var(2)) + psi % (lap(psi) - var(1)) + rho % (lap(rho) - var(2)))
)
LINK_WITH_NAME(MG, MODELG)

// Model J
MODEL(MJ, (VECTOR),
	PROVISIONAL_DEF((VECTOR), 
		var(1) = -(c1 - lit(4.) * c2 * psi * psi) * psi)
	MODEL_DEF(
		dpsi = -bilap(psi) + lap(var(1)) + psi % (lap(psi) - var(1)))
)
LINK_WITH_NAME(MJ, MODELJ)
*/


// Turing Model
MODEL(Turing, (SCALAR, SCALAR),
	MODEL_DEF(
		dpsi = c1 * lap(psi) + c2 * (psi + c3 * rho - psi * rho * rho - c4 * psi * rho),
		drho = lap(rho) + c2 * (c5 * rho + c6 * psi + psi * rho * rho + c4 * psi * rho))
)
LINK_WITH_NAME(Turing, TURING)

#endif

#undef rho
#undef psi
#undef dpsi
#undef drho

// *******************************************************************************


