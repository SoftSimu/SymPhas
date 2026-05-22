#pragma once
//
// Focused regression tests for the autogen stencil dictionary builder
// (sym/inc/expressionstencils.h).  Each subtest:
//   1. Calls expr::get_central_space_stencil<O, OA, D>() (or the
//      mixed-derivative variant) for a specific (O, OA, D) tuple.
//   2. Wraps the result in a GeneratedStencilApply.
//   3. Applies it to a known polynomial sample p(x, y) = x^a * y^b and
//      compares the output to the analytic derivative at the center.
//
// The test asserts numeric correctness, which catches both:
//   (a) dictionary entries that fail to reduce (`S_symbol` residual,
//       which evaluates to expr::symbols::Symbol{} and propagates a
//       runtime "Symbol" through computation), and
//   (b) wrong coefficient magnitudes from under-determined moment
//       systems.
//
// Build:  -DUSE_TESTS=ON -DCHOSEN_TESTS=testautogenstencil
//         -DAVAILABLE_STENCILS_AUTOGENERATION=ON
//
// Counts failures and reports them; entry point returns int.

int testautogenstencil();
