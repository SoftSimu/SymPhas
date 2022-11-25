
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

#ifndef BASIC_MODELS

//! Model A.
MODEL(MA, (SCALAR),
	EVOLUTION(
		dpsi = lap(psi) + (c1 - c2 * psi * psi) * psi)
)
LINK_WITH_NAME(MA, MODELA)


//! Model B.
MODEL(MB, (SCALAR),
	EVOLUTION(
		dpsi = -bilap(psi) - lap((c1 - c2 * psi * psi) * psi))
)
LINK_WITH_NAME(MB, MODELB)

//! Model C.
MODEL(MC, (SCALARS(2)),
	EVOLUTION(
		dpsi = -bilap(psi) - lap((c1 - c2 * psi * psi) * psi + c5 * rho * rho), 
		drho = lap(rho) + (c3 - c4 * rho * rho) * rho + integer(2) * c5 * psi * rho)
)
LINK_WITH_NAME(MC, MODELC)
DEFINE_MODEL_FIELD_NAMES(MC, ("psi", "m"))

#define j rho
#define dj drho

// Model H
MODEL(MH, (SCALAR, VECTOR),
	EVOLUTION_PREAMBLE(
		(auto f = lap(psi) + (c1 - c2 * psi * psi) * psi; ),
		dpsi = -lap(f) - c3 * grad(psi) * j,
		dj = lap(j) + c4 * grad(psi) * f)
)
LINK_WITH_NAME(MH, MODELH)
DEFINE_MODEL_FIELD_NAMES(MH, ("m", "j"))


#undef j
#undef dj


//! Model A defined by the free energy.
MODEL(MA_FE, (SCALAR),
	FREE_ENERGY((NONCONSERVED),
		SUM(ii)(LANDAU_FE(op_ii)))
)
LINK_WITH_NAME(MA_FE, MODELA_FE)


//! Model B defined by the free energy.
MODEL(MB_FE, (SCALAR),
	FREE_ENERGY((CONSERVED),
		SUM(ii)(LANDAU_FE(op_ii)))
)
LINK_WITH_NAME(MB_FE, MODELB_FE)

//! Model B defined by the free energy.
MODEL(MC_FE, (SCALAR, SCALAR),
	FREE_ENERGY((NONCONSERVED, CONSERVED),
		SUM(ii)(LANDAU_FE(op_ii)) + psi * psi * rho)
)
LINK_WITH_NAME(MC_FE, MODELC_FE)


//! The Gray-Scott phase field model.
MODEL(GRAYSCOTT, (SCALAR, SCALAR),
	EVOLUTION_PREAMBLE((auto prr = psi * rho * rho;),
		dpsi = c1 * lap(psi) - prr + c3 * (one - psi),
		drho = c2 * lap(rho) + prr - (c3 + c4) * rho)
)
LINK_WITH_NAME(GRAYSCOTT, GRAYSCOTT)

// Turing Model
MODEL(Turing, (SCALAR, SCALAR),
	EVOLUTION(
		dpsi = c1 * lap(psi) + c2 * (psi + c3 * rho - psi * rho * rho - c4 * psi * rho),
		drho = lap(rho) + c2 * (c5 * rho + c6 * psi + psi * rho * rho + c4 * psi * rho))
)
LINK_WITH_NAME(Turing, TURING)

/*
//! Example of provisional variables, Model B.
MODEL(MBB, (SCALAR),
	PROVISIONAL_DEF((SCALAR),
		var(1) = -(c1 - lit(4.) * c2 * psi * psi) * psi)
	EVOLUTION(
		dpsi = -bilap(psi) + lap(var(1)))
)
LINK_WITH_NAME(MBB, MODELBB)
*/


#endif

#undef rho
#undef psi
#undef dpsi
#undef drho

// *******************************************************************************


