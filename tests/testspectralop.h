#pragma once
//
// Regression test for the spectral solver's linear-operator construction.
//
// Goal: isolate whether the apparent sign/term issue in the printed
// "spectral scheme" for Model B (Cahn-Hilliard, dpsi = -bilap(psi) - lap(psi)
//  + lap(psi^3)) is:
//   (1) purely a display/printing artifact, or
//   (2) a real symbolic-algebra regression where l_op actually evaluates to
//       the wrong numeric value.
//
// We construct the Model-B RHS symbolically, run it through the SP solver's
// get_l_op / form_A_op / form_B_op pipeline, print the symbolic form, AND
// numerically evaluate l_op at hand-picked k-grid indices and compare against
// the closed-form Cahn-Hilliard symbol L(k) = -|k|^4 + |k|^2.
//
// Build: -DUSE_TESTS=ON -DCHOSEN_TESTS=testspectralop -DUSE_FFTW3=ON

void testspectralop();
