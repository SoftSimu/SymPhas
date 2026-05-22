#pragma once
//
// Compile-only regression tests for the AnisotropicFMPFC investigation.
//
// Each model defined in examples/models/aniso_fmpfc_test_models.h
// exercises a specific class of expression-template interactions.
// Including the header forces the templates to be instantiated; any
// regression in the symbolic-algebra rules manifests as a build
// failure that the CMake test harness picks up before runtime.
//
// Build:  -DCHOSEN_TESTS=testanisocompile  via the existing USE_TESTS
//         mechanism.  Each model variant gates on its own ANISO_*
//         define, set here so they're all active when this test
//         translation unit is compiled.

#include "modelspecialized.h"
#include "solverinclude.h"

// Force instantiation of every bisection-ladder model.  Surrounding
// code paths in aniso_fmpfc_test_models.h are gated on these.
#define ANISO_T1
#define ANISO_T2
#define ANISO_T3
#define ANISO_T4
#define ANISO_T5
#define ANISO_T6
#define ANISO_T3a
#define ANISO_T3b
#define ANISO_T3c
#define ANISO_T3d
#define ANISO_T3e
#define ANISO_T7
#define ANISO_T8

#include "modelmacros.h"

// Reach the model-definition macros' expectation of a PoissonSolver
// shorthand (set up the same way modelinclude.h does).
#define PoissonSolver(E) expr::poisson_solver(E)

// Bug-4 ladder model gates: lighter ones (L6) on by default; heavier
// (L7, L8) opt-in via SYMPHAS_TEST_BUG4_HEAVY because L8 alone takes
// ~5 GB to compile.
#define ANISO_L6
#ifdef SYMPHAS_TEST_BUG4_HEAVY
#define ANISO_L7
#define ANISO_L8
#endif

// AnisotropicFMPFC is now an unconditional model defined in
// modelinclude.h alongside MagneticPFC2013, so no opt-in is required
// to compile its template.  SYMPHAS_TEST_ANISO_FMPFC_ORIGINAL remains
// recognized as a no-op for back-compat with build scripts.

// Pull in the bisection-ladder model definitions (ANISO_T*/ANISO_L*)
// and the production AnisotropicFMPFC via modelinclude.h.
#include "modelinclude.h"
#include "aniso_fmpfc_test_models.h"

void testanisocompile();

