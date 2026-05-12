#pragma once
//
// Symbolic-algebra runtime regression tests for the AnisotropicFMPFC
// investigation.  Each function returns the number of failures it
// detected; the master entry point sums them and exits nonzero if any
// triggered.
//
// Covered:
//   - Bug 1: lap(a + b) must evaluate identically to lap(a) + lap(b)
//   - Bug 2: OpAdd containing div(...) -> correct divergence values
//   - Bug 3: pow<N>(expr) evaluates to expr^N for N up to 5, with both
//            scalar and deep-subtree inner expressions
//   - Sanity:  laplacian of sin(k*x) -> -k^2 sin(k*x) at interior points
//   - Sanity:  product rule d/dx(f*g) == df*g + f*dg
//
// Build:  add testanisofmpfc.cpp to USE_TESTS; the runner is the
// existing tests/ harness from tests.cpp.in.

#include "conf.h"
#include "expressiontypeincludes.h"
#include "modelspecialized.h"
#include "solverinclude.h"

// Public entry point.  Aggregates per-bug subtest counts and prints a
// summary line.  Returns nothing; sets a static counter that the
// running process inspects via testanisofmpfc_failures().
void testanisofmpfc();

// Total number of subtest assertions that failed during the most recent
// testanisofmpfc() call.  Zero means the symbolic algebra is preserved.
int testanisofmpfc_failures();
