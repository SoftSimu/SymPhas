
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

#define cell_params param_matrix(1)
#define gamma_n_ cell_params(1)

#define lambda (c(1))
#define R (c(3)) //(c(3) * AREA / (Pi * NUM_FIELDS))
#define R2 (c(3) * c(3)) //(c(3) * AREA / (Pi * NUM_FIELDS))
#define mu (c(2))
#define kappa c(4)
#define xi c(5)
#define gamma (c(6))
#define gammaC (c(7))
#define vel c(8)
#define tau c(9)

#define dpsi dop(1)
#define psi op(1)
#define drho dop(2)
#define rho op(2)

PARAMETERIZED_TYPE(NORMAL_CELL, SCALAR, (gamma))
PARAMETERIZED_TYPE(CANCER_CELL, SCALAR, (gammaC))

MODEL(CELL_MIGRATION_NO_MOTILITY, (MANY(CANCER_CELL, CONFIGURATION), MANY(NORMAL_CELL, CONFIGURATION)),
	FREE_ENERGY(
		(EQUATION_OF(ii)(
			-_2 * DF_(ii) -
			(60_n * kappa / (xi * lambda * lambda) * INT(op_ii * grad(op_ii) * SUM(jj != ii)(op_jj * op_jj))) * grad(op_ii)
			)),
		SUM(ii)(gamma_n_[ii] * INT(CELLULAR_FE(op_ii, lambda)) + mu / (Pi * R2) * pow<2>(Pi * R2 - INT(op_ii * op_ii)))
		+ INT(SUM(ii, jj != ii)(30_n * kappa / (lambda * lambda) * op_ii * op_ii * op_jj * op_jj)))
)
LINK_WITH_NAME(CELL_MIGRATION_NO_MOTILITY, CELL_MODEL_NO_MOT)
DEFINE_MODEL_FIELD_NAMES_FORMAT(CELL_MIGRATION_NO_MOTILITY, "\\phi_{%d}")

#undef lambda
#undef kappa
#undef R
#undef R2
#undef mu
#undef xi
#undef gamma
#undef gammaC
#undef vel
#undef tau

