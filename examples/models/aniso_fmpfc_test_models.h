// Test fixtures for the AnisotropicFMPFC investigation.
// See review/spectral_sanity/_anisofmpfc_plan.md for context.
// Each ladder model is guarded by its own ANISO_T*/ANISO_L* define so
// the testanisocompile harness can opt in to one or more at a time.
//
// The production AnisotropicFMPFC definition now lives in
// examples/models/modelinclude.h alongside MagneticPFC2013; the
// ENABLE_ANISO_FMPFC_* alternative formulations below are kept here
// for investigation purposes only.
#pragma once

#ifdef ENABLE_ANISO_FMPFC_REWRITTEN
// Provisional CSE + merged-lap rewrite.  Now compiles for the merged lap
// thanks to the OpAdd distribution fix in initialize_derivative_order;
// dop(2) still fails because PROVISIONAL_DEF + inline PoissonSolver chain
// hits a separate eval_type bug (T7/T8 territory).
MODEL(AnisotropicFMPFC, (SCALAR, VECTOR),
      PROVISIONAL_DEF((SCALAR, SCALAR, SCALAR),
            var(1) <= dot(op(2), grad(op(1))),
            var(2) <= dot(op(2), grad(dot(grad(op(1)), op(2)))),
            var(3) <= div(op(2))
      )
      EVOLUTION(
            dop(1) = lap(
                        c(1) * op(1) + c(2) * op(1)
                      + c(2) * 2_n * lap(op(1)) + c(2) * bilap(op(1))
                      - c(3) * power(op(1), 2) + c(4) * power(op(1), 3)
                      - c(8) * op(1) * dot(op(2), op(2))
                      + c(10) * var(1) * var(3) + c(10) * var(2)
                      + c(11) * power(var(1), 3) * var(3)
                      + c(11) * 3_n * power(var(1), 2) * var(2)
                      + c(12) * power(var(1), 5) * var(3)
                      + c(12) * 5_n * power(var(1), 4) * var(2)
                  ),
            dop(2) = c(6) * c(6) * lap(op(2)) - c(7) * op(2)
                  + c(8) * power(op(1), 2) * op(2)
                  - c(9) * op(2) * dot(op(2), op(2))
                  + c(10) * grad(op(1)) * var(1)
                  + grady(PoissonSolver(curl(op(2)))) * e(x)
                  - gradx(PoissonSolver(curl(op(2)))) * e(y)
      )
)
LINK_WITH_NAME(AnisotropicFMPFC, ANISOTROPICFMPFC)
#endif

#ifdef ENABLE_ANISO_FMPFC_PROV2
// Same as ENABLE_ANISO_FMPFC_PROV, but div(op(2)) is NOT captured as a
// provisional; it stays inlined inside each product so the surrounding
// multiplication keeps the type system from triggering the lap(div(...))
// auto-deduction recursion (Bug 1 in the progress doc).
MODEL(AnisotropicFMPFC, (SCALAR, VECTOR),
      PROVISIONAL_DEF((SCALAR, VECTOR, SCALAR),
            var(1) <= PoissonSolver(curl(op(2))),
            var(2) <= grady(var(1)) * e(x) - gradx(var(1)) * e(y),
            var(3) <= dot(op(2), grad(op(1)))
      )
      EVOLUTION(
            dop(1) = lap(
                        c(1) * op(1) + c(2) * op(1)
                      + c(2) * 2_n * lap(op(1)) + c(2) * bilap(op(1))
                      - c(3) * power(op(1), 2) + c(4) * power(op(1), 3)
                      - c(8) * op(1) * dot(op(2), op(2))
                      + c(10) * var(3) * div(op(2))
                      + c(10) * dot(op(2), grad(var(3)))
                      + c(11) * power(var(3), 3) * div(op(2))
                      + c(11) * 3_n * power(var(3), 2) * dot(op(2), grad(var(3)))
                      + c(12) * power(var(3), 5) * div(op(2))
                      + c(12) * 5_n * power(var(3), 4) * dot(op(2), grad(var(3)))
                  ),
            dop(2) = c(6) * c(6) * lap(op(2)) - c(7) * op(2)
                  + c(8) * power(op(1), 2) * op(2)
                  - c(9) * op(2) * dot(op(2), op(2))
                  + c(10) * grad(op(1)) * var(3)
                  + c(11) * grad(op(1)) * power(var(3), 3)
                  + c(12) * grad(op(1)) * power(var(3), 5)
                  + var(2)
                  + c(13) * e(x) + c(14) * e(y)
      )
)
LINK_WITH_NAME(AnisotropicFMPFC, ANISOTROPICFMPFC)
#endif
// provisional block (var(1)=A_z, var(2)=B_ind), and the deep dot/grad

// ---------- Bisection ladder ----------

#ifdef ANISO_T1
MODEL(AnisoT1, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(c(1) * op(1)),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT1, ANISOT1)
#endif

#ifdef ANISO_T2
MODEL(AnisoT2, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(c(1) * op(1) + c(2) * lap(op(1))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT2, ANISOT2)
#endif

#ifdef ANISO_T3
MODEL(AnisoT3, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(c(1) * op(1) + c(10) * dot(op(2), grad(op(1))) * div(op(2))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT3, ANISOT3)
#endif

#ifdef ANISO_T4
MODEL(AnisoT4, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(c(1) * op(1) + c(10) * dot(op(2), grad(dot(grad(op(1)), op(2))))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT4, ANISOT4)
#endif

#ifdef ANISO_T5
MODEL(AnisoT5, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(c(10) * dot(op(2), grad(op(1))) * div(op(2))
                 + c(10) * dot(op(2), grad(dot(grad(op(1)), op(2))))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT5, ANISOT5)
#endif

#ifdef ANISO_T6
MODEL(AnisoT6, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(c(1) * op(1))
             + lap(c(10) * dot(op(2), grad(op(1))) * div(op(2)))
             + lap(c(10) * dot(op(2), grad(dot(grad(op(1)), op(2))))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT6, ANISOT6)
#endif

#ifdef ANISO_T3a
MODEL(AnisoT3a, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(c(10) * dot(op(2), grad(op(1))) * div(op(2))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT3a, ANISOT3A)
#endif

#ifdef ANISO_T3b
MODEL(AnisoT3b, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(c(1) * op(1) + c(10) * div(op(2))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT3b, ANISOT3B)
#endif

#ifdef ANISO_T3c
MODEL(AnisoT3c, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(c(1) * op(1) + c(10) * dot(op(2), grad(op(1)))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT3c, ANISOT3C)
#endif

#ifdef ANISO_T3d
MODEL(AnisoT3d, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = c(1) * op(1) + c(10) * dot(op(2), grad(op(1))) * div(op(2)),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT3d, ANISOT3D)
#endif

#ifdef ANISO_T3e
MODEL(AnisoT3e, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(c(1) * op(1) + c(10) * div(op(2)) * dot(op(2), grad(op(1)))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT3e, ANISOT3E)
#endif

#ifdef ANISO_T7
// lap(c(10) * div(op(2))) alone -- exposes a separate pre-existing bug
// with lap(div(...)) auto-deduction.  Not addressed by the OpAdd
// distribution fix.
MODEL(AnisoT7, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(c(10) * div(op(2))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT7, ANISOT7)
#endif

#ifdef ANISO_T8
MODEL(AnisoT8, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(div(op(2))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoT8, ANISOT8)
#endif

// Bug-4 ladder models.  Each progressively adds one piece of the
// user's failing AnisotropicFMPFC term so a compile failure pinpoints
// the layer whose simplification rule is broken.  All models reuse
// the same primary fields (psi, M) as the bisection ladder above.

#ifdef ANISO_L6
// L6: lap(power(dot(M, grad psi), 3)).  Adds an outer Laplacian over
// a scalar power of a dot product -- the simplest shape where a
// per-component decomposition tensor V interacts with the OpPow node.
MODEL(AnisoL6, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(power(dot(op(2), grad(op(1))), 3)),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoL6, ANISOL6)
#endif

#ifdef ANISO_L7
// L7: lap(power(dot(M, grad psi), 3) * div(M)).  Adds the * div(M)
// factor -- the exact shape the user model has.
MODEL(AnisoL7, (SCALAR, VECTOR), EVOLUTION(
      dop(1) = lap(power(dot(op(2), grad(op(1))), 3) * div(op(2))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoL7, ANISOL7)
#endif

#ifdef ANISO_L8
// L8: lap of the full user-model dop(1) inner term at coefficient 3.
// (The user model also has 5 elsewhere; this is the first one that
// fails in the original definition.)
MODEL(AnisoL8, (SCALAR, VECTOR), EVOLUTION(
      dop(1) =
          lap(c(11) * power(dot(op(2), grad(op(1))), 3) * div(op(2))) +
          lap(c(11) * 3_n * power(dot(op(2), grad(op(1))), 2)
              * dot(op(2), grad(dot(grad(op(1)), op(2))))),
      dop(2) = c(6) * c(6) * lap(op(2))))
LINK_WITH_NAME(AnisoL8, ANISOL8)
#endif
