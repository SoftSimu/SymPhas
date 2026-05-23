
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

#include "modelmacros.h"

#define PoissonSolver(E) expr::poisson_solver(E)

// TODO: Solvers (especially JFNK) need the ability to "block" provisional
// variables from appearing in evolution equations. Provisional variables are
// evaluated once per timestep in update(), but implicit solvers need to
// re-evaluate all field-dependent quantities at trial states during Newton
// iterations. Models using provisionals should either inline the provisional
// expressions (as done here) or the solver interface should be extended so
// that form_expr receives the provisional expressions alongside the evolution
// equations, allowing re-evaluation during the solve.

// MagneticPFC2013: Magnetic Phase-Field Crystal (Faghihi et al., PRE 88, 2013)
//
// The demagnetizing field B_ind = curl(A_z) where A_z = PoissonSolver(curl(M))
// is inlined directly into the magnetization evolution equation rather than
// using provisional variables. This is mathematically equivalent to the
// original formulation but ensures the demagnetizing field is a pure functional
// of M, which is essential for implicit solvers (JFNK) that must re-evaluate
// the RHS at arbitrary trial states.
//
// B_ind_x =  dAz/dy = grady(PoissonSolver(curl(M)))
// B_ind_y = -dAz/dx = -gradx(PoissonSolver(curl(M)))

MODEL(MagneticPFC2013, (SCALAR, VECTOR),
      EVOLUTION(
            dop(1) = lap(c(1) * op(1) + c(2) * op(1) +
                        c(2) * 2_n * lap(op(1)) + c(2) * bilap(op(1)) -
                        c(3) * power(op(1), 2) + c(4) * power(op(1), 3) -
                        c(8) * op(1) * dot(op(2), op(2))) +
                  lap(c(10) * dot(op(2), grad(op(1))) * div(op(2))) +
                  lap(c(10) * dot(op(2), grad(dot(grad(op(1)), op(2))))),
            dop(2) = c(6) * c(6) * lap(op(2)) - c(7) * op(2) +
                  c(8) * power(op(1), 2) * op(2) -
                  c(9) * op(2) * dot(op(2), op(2)) +
                  c(10) * grad(op(1)) * dot(op(2), grad(op(1))) +
                  grady(PoissonSolver(curl(op(2)))) * e_x -
                  gradx(PoissonSolver(curl(op(2)))) * e_y
      )
)
LINK_WITH_NAME(MagneticPFC2013, MAGNETICPFC2013)

MODEL(FMPFCLinearField, (SCALAR, VECTOR),
      PROVISIONAL_DEF((SCALAR, VECTOR, SCALAR), 
        var(1) <= PoissonSolver(curl(op(2))),
        var(2) <= grady(var(1)) * e_x - gradx(var(1)) * e_y,
        var(3) <= c(11) * asin(sin(2 * pi_n * t / c(12)))
      )
      EVOLUTION(
            dop(1) = lap(c(1) * op(1) + c(2) * op(1) +
                        c(2) * 2_n * lap(op(1)) + c(2) * bilap(op(1)) -
                        c(3) * power(op(1), 2) + c(4) * power(op(1), 3) -
                        c(8) * op(1) * dot(op(2), op(2))) +
                  lap(c(10) * dot(op(2), grad(op(1))) * div(op(2))) +
                  lap(c(10) * dot(op(2), grad(dot(grad(op(1)), op(2))))),
            dop(2) = c(6) * c(6) * lap(op(2)) - c(7) * op(2) +
                  c(8) * power(op(1), 2) * op(2) -
                  c(9) * op(2) * dot(op(2), op(2)) +
                  c(10) * grad(op(1)) * dot(op(2), grad(op(1))) + var(2) + var(3) * e_y
      )
)
LINK_WITH_NAME(FMPFCLinearField, FMPFCLINEARFIELD)

// AnisotropicFMPFC: MagneticPFC2013 plus higher-order magnetostriction.
//
// Adds free-energy contributions
//   F += -omega * alpha3 * (m . grad n)^4 / 4
//   F += -omega * alpha5 * (m . grad n)^6 / 6
// and a uniform applied magnetic field (c(13), c(14)).  In the EOM the
// functional derivatives become the c(11) and c(12) terms below; each
// power is split into two pieces by the product rule
//   div( m * (m.grad n)^k ) = (m.grad n)^k div(m)
//                           + k * (m.grad n)^(k-1) * m.grad(m.grad n)
//
// Coefficient slots:
//   c(1)..c(10) match MagneticPFC2013 (DeltaB, Bs, t, v, unused, W0,
//   omega*r_c, omega*beta, omega*gamma, omega*alpha).
//   c(11) = cubic magnetostriction (alpha3).
//   c(12) = quintic magnetostriction (alpha5).
//   c(13), c(14) = uniform applied field (B_app_x, B_app_y).
MODEL(AnisotropicFMPFC, (SCALAR, VECTOR),
      EVOLUTION(
            dop(1) = lap(c(1) * op(1) + c(2) * op(1) +
                        c(2) * 2_n * lap(op(1)) + c(2) * bilap(op(1)) -
                        c(3) * power(op(1), 2) + c(4) * power(op(1), 3) -
                        c(8) * op(1) * dot(op(2), op(2))) +
                  lap(c(10) * dot(op(2), grad(op(1))) * div(op(2))) +
                  lap(c(10) * dot(op(2), grad(dot(grad(op(1)), op(2))))) +
                  lap(c(11) * power(dot(op(2), grad(op(1))), 3) * div(op(2))) +
                  lap(c(11) * 3_n * power(dot(op(2), grad(op(1))), 2) * dot(op(2), grad(dot(grad(op(1)), op(2))))) +
                  lap(c(12) * power(dot(op(2), grad(op(1))), 5) * div(op(2))) +
                  lap(c(12) * 5_n * power(dot(op(2), grad(op(1))), 4) * dot(op(2), grad(dot(grad(op(1)), op(2))))),
            dop(2) = c(6) * c(6) * lap(op(2)) - c(7) * op(2) +
                  c(8) * power(op(1), 2) * op(2) -
                  c(9) * op(2) * dot(op(2), op(2)) +
                  c(10) * grad(op(1)) * dot(op(2), grad(op(1))) +
                  c(11) * grad(op(1)) * power(dot(op(2), grad(op(1))), 3) +
                  c(12) * grad(op(1)) * power(dot(op(2), grad(op(1))), 5) +
                  grady(PoissonSolver(curl(op(2)))) * e_x -
                  gradx(PoissonSolver(curl(op(2)))) * e_y +
                  c(13) * e_x + c(14) * e_y
      )
)
LINK_WITH_NAME(AnisotropicFMPFC, ANISOTROPICFMPFC)

#if defined(ANISO_T1) || defined(ANISO_T2) || defined(ANISO_T3) || defined(ANISO_T4) || defined(ANISO_T5) || defined(ANISO_T6) || \
    defined(ANISO_T3a) || defined(ANISO_T3b) || defined(ANISO_T3c) || defined(ANISO_T3d) || defined(ANISO_T3e) || \
    defined(ANISO_T7) || defined(ANISO_T8) || \
    defined(ANISO_L6) || defined(ANISO_L7) || defined(ANISO_L8)
// Bisection-ladder test fixtures used by tests/testanisocompile.  These
// stay guarded because they intentionally exercise many overlapping
// model templates that drive cc1plus to several GB just by instantiation.
#include "aniso_fmpfc_test_models.h"
#endif

// #include "advancedmodeldefs.h"
// #include "modeldefinitions.h"
// #include "pfcdefs.h"
// #include "modelacmms.h"
#endif
