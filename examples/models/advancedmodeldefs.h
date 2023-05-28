
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

#define cell_params C(6)


#define gamma_n_ cell_params(6)

#define lambda_ cell_params(1)
#define lambda2_(I) (lambda_(I) * lambda_(I))

#define R2 (c3 * c3) //(c3 * AREA / (Pi * NUM_FIELDS))

//#define R2 (c3 * c3)
//#define R (sqrt(R2))

#define lambda (c1)
#define R (c3) //(c3 * AREA / (Pi * NUM_FIELDS))
#define R2 (c3 * c3) //(c3 * AREA / (Pi * NUM_FIELDS))
#define mu (c2)
#define kappa c4
#define xi c5
#define gamma (c6)
#define gammaC (c7)
#define vel c8
#define tau c9

#define dpsi dop(1)
#define psi op(1)
#define drho dop(2)
#define rho op(2)

PARAMETERIZED_TYPE(NORMAL_CELL, SCALAR, (lambda, mu, R, kappa, xi, gamma))
PARAMETERIZED_TYPE(CANCER_CELL, SCALAR, (lambda, mu, R, kappa, xi, gammaC))

MODEL(CELL_MIGRATION, (MANY(CANCER_CELL, 1), MANY(NORMAL_CELL, CONFIGURATION)),
	FREE_ENERGY(
		(EQUATION_OF(ii)(
			-_2 * DF_(ii) - 
			(vel * e(2 * Pi * ARRAY(ii)(Poisson(1. / tau)(t, ii))) 
				+ val<60> * kappa / (xi * lambda_[ii]) * integral(op_ii * grad(op_ii) * SUM(jj != ii)(op_jj * op_jj))) * grad(op_ii)
			)),		
		SUM(ii)(gamma_n_[ii] * integral(CELLULAR_FE(op_ii, lambda_[ii])) + mu / (Pi * R2) * pow<2>(Pi * R2 - integral(op_ii * op_ii)))
		+ integral(SUM(ii, jj != ii)(integer(30) * kappa / lambda_[ii] * op_ii * op_ii * op_jj * op_jj)))
)
LINK_WITH_NAME(CELL_MIGRATION, CELL_MODEL)
DEFINE_MODEL_FIELD_NAMES_FORMAT(CELL_MIGRATION, "\\phi_{%d}")

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

//
//MODEL(CELL_MIGRATION, (MANY(CANCER_CELL, 1), MANY(NORMAL_CELL, CONFIGURATION)),
//	FREE_ENERGY((EQUATION_OF(ii)(-_2 * DF_(ii) - val<60> * kappa / (xi * lambda2_(ii)) * integral(op_ii * grad(op_ii) * SUM(jj != ii)(op_jj * op_jj)) * grad(op_ii))),
//		SUM(ii)(gamma_n_(ii) * integral(CELLULAR_FE(op_ii, lambda_(ii))) + mu / (Pi * R2) * pow<2>(Pi * R2 - integral(op_ii * op_ii)))
//		+ integral(SUM(ii, jj != ii)(integer(30) * kappa / lambda2_(ii) * op_ii * op_ii * op_jj * op_jj)))
//)
//LINK_WITH_NAME(CELL_MIGRATION, CELL_MODEL)
//DEFINE_MODEL_FIELD_NAMES_FORMAT(CELL_MIGRATION, "\\phi_{%d}")

//MODEL(MODEL_B_MANY, (SCALARS(RUNTIME)),
//	FREE_ENERGY((EQUATION_OF(ii)(lap * DF_(ii))),
//		SUM(ii)(integral(LANDAU_FE(op_ii))))
//)
//LINK_WITH_NAME(MODEL_B_MANY, MODELB_MANY)

//MODEL(FLOCK, (VECTOR, SCALAR),
//	EVOLUTION(
//		dpsi = c1 * psi - c2 * psi * psi * psi - grad(c6 * (rho - STATS.mean(rho)) + c7 * POW(2)(rho - STATS.mean(rho))) + c3 * grad(grad * psi) + c4 * lap(psi) + c5 * POW(2)(psi * grad) * psi - (psi * grad) * psi + _nW(VECTOR),
//		drho = -grad * (psi * rho))
//)
//LINK_WITH_NAME(FLOCK, FLOCK_MODEL)

//MODEL(FLOCK_B, (VECTOR, SCALAR),
//	EVOLUTION_PREAMBLE(
//		(
//			auto lam = c1;
//			auto sig2 = c2 * c2;
//			auto sig2_0 = c3 * c3;
//			auto nu = frac<1, 4> / (lam * (1 - exp(-2 * sig2_0)) + val<4> / Pi * rho * (frac<14, 15> + frac<2, 3> * exp(-2 * sig2)));
//			auto gg = val<8> / Pi * nu * (frac<16, 15> + val<2> * exp(-2 * sig2) - exp(-sig2 / 2));
//			auto pp = val<8> / Pi * nu * (frac<4, 15> + val<2> * exp(-2 * sig2) - exp(-sig2 / 2));
//			auto mu = val<4> / Pi * rho * (exp(-sig2/2) - frac<2, 3>) - lam * (1 - exp(-sig2_0 / 2));
//			auto xi = val<64> / (Pi * Pi) * nu * (exp(-sig2 / 2) - frac<2, 5>) * (frac<1, 3> + exp(-2 * sig2));
//		),
//		dpsi = -gg * (psi * grad) * psi - _2 * grad * (rho - pp * psi * psi) + (mu - xi * psi * psi) * psi + nu * lap(psi) - pp * (grad * psi) * psi,
//		drho = -grad * (psi))
//)
//LINK_WITH_NAME(FLOCK_B, BERTIN_MODEL)

//#define w op(1)
//#define dw dop(1)
//
//MODEL(FLOCK_B2, (VECTOR, SCALAR),
//	EVOLUTION_PREAMBLE(
//		(
//			auto lam = c1;
//			auto sig2 = c2 * c2;
//			auto sig2_0 = c3 * c3;
//
//			auto v0 = c4;
//			auto v02 = v0 * v0;
//			auto d0 = one;
//
//			auto nu = v02 * frac<1, 4> / (lam * (1 - exp(-2 * sig2_0)) + frac<16, 3> / Pi * d0 * v0 * rho * (frac<7, 5> + exp(-2 * sig2)));
//			auto dnu = diff(nu, rho);
//			auto gg = val<16> / Pi * nu * (d0 / v0) * (frac<16, 15> + 2 * exp(-2 * sig2) - exp(-sig2 / 2));
//			auto pp = val<16> / Pi * nu * (d0 / v0) * (frac<4, 15> + val<2> * exp(-2 * sig2) - exp(-sig2 / 2));
//			auto mu = val<8> / Pi * d0 * v0 * rho * (exp(-sig2/2) - frac<2, 3>) - lam * (1 - exp(-sig2_0 / 2));
//			auto xi = val<256> / (Pi * Pi) * (d0 * d0 / v02) * nu * (exp(-sig2 / 2) - frac<2, 5>) * (frac<1, 3> + exp(-2 * sig2));
//			auto M = _2(grad(w) + T(grad(w)));
//		),
//		dw = -gg * (w * grad) * w - v02 / 2 * grad(rho) + pp / 2 * grad(w * w) + (mu - xi * w * w) * w + nu * lap(w) - pp(grad * w) * w + 2 * dnu * M * grad(rho) - dnu * (grad * w) * grad(rho),
//		drho = -grad * (w))
//)
//LINK_WITH_NAME(FLOCK_B2, BERTIN2_MODEL)

#undef lambda
#undef lambda2
#undef mu_n
#undef R2
#undef R

