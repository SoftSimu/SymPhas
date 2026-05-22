#pragma once
//
// Regression test for the printed equation form of each registered
// model.  Constructs the model and captures the "given equation"
// string emitted by MakeEquation::make_equations / printe, then
// compares against an expected canonical form.
//
// Why this is useful as a regression contract:
//   The printed form is a faithful read-out of the expression tree the
//   framework actually built.  If any symbolic rule (algebra closure,
//   eval_type rank, derivative distribution) regresses, the printed
//   string changes -- a fast, deterministic signal that doesn't
//   require running the integrator.
//
// Build:  -DCHOSEN_TESTS=testmodelequations
//
// Returns 0 on success, otherwise nonzero count of failing models.

int testmodelequations();
