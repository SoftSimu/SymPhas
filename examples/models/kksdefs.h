
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
 * This file includes the Kim-Kim-Suzuki (KKS) two-phase coupled
 * phase-field model.
 *
 * ***************************************************************************
 */

#pragma once

#include "modelmacros.h"

//! Kim-Kim-Suzuki phase-field model with equal parabolic phase free energies.
//!   Field 1 (op(1)) = phi, non-conserved phase indicator (~1 in solid, ~0 in liquid).
//!   Field 2 (op(2)) = c, conserved composition.
//! With equal parabolic curvatures A_alpha = A_beta = A, the equal-mu constraint
//! and lever rule give closed-form
//!     c_alpha(phi, c) = c - (1 - h(phi)) * (c_eq_beta - c_eq_alpha),
//!     c_beta (phi, c) = c +     h(phi)  * (c_eq_beta - c_eq_alpha),
//!     mu             = 2 A * (c - h(phi)*c_eq_alpha - (1-h(phi))*c_eq_beta),
//!     f_alpha - f_beta = 0   (cancels for equal A).
//! Equations of motion (Karma-Rappel sign convention):
//!     d_t phi = (1/tau) * [ eps^2 * lap(phi) - W * g'(phi)
//!                          - h'(phi) * mu * (c_eq_beta - c_eq_alpha) ]
//!     d_t c   = M_c * lap(mu)
//! with h(phi)  = phi^3 (10 - 15 phi + 6 phi^2),
//!      h'(phi) = 30 phi^2 (1 - phi)^2,
//!      g(phi)  = phi^2 (1-phi)^2,
//!      g'(phi) = 2 phi (1 - phi) (1 - 2 phi).
//! Coefficients (config "coefficients.cN"):
//!   c(1) = c_eq_alpha, c(2) = c_eq_beta, c(3) = A,
//!   c(4) = W, c(5) = epsilon^2, c(6) = 1/tau, c(7) = M_c.
MODEL(KKS, (SCALARS(2)),
      EVOLUTION_PREAMBLE(
            (auto h_phi   = pow<3>(op(1)) * (10_n - 15_n * op(1) + 6_n * pow<2>(op(1)));
             auto h_prime = 30_n * pow<2>(op(1)) * pow<2>(one - op(1));
             auto g_prime = 2_n * op(1) * (one - op(1)) * (one - 2_n * op(1));
             auto mu      = 2_n * c(3) * (op(2) - h_phi * c(1) - (one - h_phi) * c(2));
            ),
            dop(1) = c(6) * (c(5) * lap(op(1)) - c(4) * g_prime
                             - h_prime * mu * (c(2) - c(1))),
            dop(2) = c(7) * lap(mu)
      ))
LINK_WITH_NAME(KKS, KKS)
DEFINE_MODEL_FIELD_NAMES(KKS, ("phi", "c"))
