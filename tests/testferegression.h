#pragma once
//
// FE-pattern regression test.
//
// Exercises the operator chains that free-energy-derived models rely on:
//   - lap of polynomial (Model A FE: lap(c1*psi + c4*psi^3))
//   - bilap of scalar (Model B FE)
//   - lap of scalar * scalar (Model C FE coupling: lap(psi^2 * rho))
//   - grad of scalar * scalar (Model H FE interface: grad(psi*rho))
//   - scalar * grad(lap(scalar)) (Model H residual: psi*grad(lap(psi)))
//   - grad of sum-of-products with a derivative inside (catches the
//     combine_mixed_derivatives O1==1 regression where the cross
//     derivative gains a spurious factor-of-2)
//
// Each block prints its POST-`apply_operators` form via expr::printe and
// asserts at runtime that the substring of forbidden patterns (e.g.
// "2d2/dx2" doubled coefficient, or duplicate terms like "*X*X*Y + X*X*Y")
// does not appear.
//
// Build:  -DUSE_TESTS=ON -DCHOSEN_TESTS=testferegression

void testferegression();
