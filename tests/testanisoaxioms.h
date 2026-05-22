#pragma once
//
// Compile-time axiom checks for the SymPhas symbolic-algebra core.
//
// Each block below builds a minimal expression matching one closure or
// shape axiom (see also: review notes -- "C1, S1..S9 axiom list") and
// asserts at compile time that `expr::eval_type<E>::rank` reports the
// algebraically-correct value.  A build failure here pinpoints exactly
// which axiom is currently violated by the framework.
//
// Each axiom appears with a one-line comment describing what shape the
// math requires.  Failing static_asserts are the regression contract.
//
// Why a separate file (vs. testanisocompile.cpp):
//   testanisocompile only asks "does the model TEMPLATE instantiate".
//   This file asks the stricter question "does eval_type report the
//   shape the algebra denotes".  C1 (rank query terminates) is implicit
//   in the assertions compiling at all; S1..S9 are explicit.
//
// Build:  -DUSE_TESTS=ON -DCHOSEN_TESTS=testanisoaxioms
//
// Known-failing axioms (Bug 4 / Bug ?-equiv) are guarded behind
// SYMPHAS_TEST_AXIOMS_STRICT so a default build still completes; flip
// the define on to surface them as compile failures.

#include "conf.h"
#include "expressiontypeincludes.h"
#include "modelspecialized.h"
#include "solverinclude.h"

// Public entry: trivial, just exists so the test harness can call it.
// The actual checking is at compile time inside testanisoaxioms.cpp.
void testanisoaxioms();
