
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
 * This file is responsible for defining advanced phase-field models.
 *
 * ***************************************************************************
 */

#pragma once

#include "modelmacros.h"

#define lambda c1
#define lambda2 (lambda * lambda)
#define mu_n (c2 / 10.0)
#define R (c3 * LENGTH(X) / 10)
#define R2 (R * R)
#define kappa c4


#define dpsi dop(1)
#define psi op(1)
#define drho dop(2)
#define rho op(2)

//! Cellular migration
//MODEL(CELL_MIGRATION, (SCALARS(4)),
//	FREE_ENERGY((EQUATION_OF_EACH(ii)(-_2 * DF_ii - (op_ii * grad(op_ii) * SUM(jj != ii)(op_jj * op_jj)) * grad(op_ii))),
//		SUM(ii)(CELLULAR_FE(op_ii, lambda)) + mu_n / (Pi * R2) * pow<2>(Pi * R2 - SUM(ii)(op_ii)) 
//		+ integer(30) * kappa / lambda2 * SUM(ii, jj != ii)(op_ii * op_ii * op_jj * op_jj))
//)
//LINK_WITH_NAME(CELL_MIGRATION, CELL_MODEL)
//DEFINE_MODEL_FIELD_NAMES_FORMAT(CELL_MIGRATION, "\\psi_{%d}")



//MODEL(MBE, (SCALAR),
//	EVOLUTION(dpsi = -c1 * bilap(psi) - grad * (one - grad(psi) * grad(psi)) * grad(psi))
//)
//LINK_WITH_NAME(MBE, MBE_MODEL)

MODEL(FLOCK, (VECTOR, SCALAR),
	EVOLUTION(
		dpsi = c1 * psi - c2 * psi * psi * psi - grad(c6 * (rho - _2) + c7 * POW(2)(rho - _2)) + c3 * grad(grad * psi) + c4 * lap(psi) + c5 * POW(2)(psi * grad) * psi - (psi * grad) * psi,
		drho = -grad * (psi * rho))
)
LINK_WITH_NAME(FLOCK, FLOCK_MODEL)


#undef lambda
#undef lambda2
#undef mu_n
#undef R2
#undef R

